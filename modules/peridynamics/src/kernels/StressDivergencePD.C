/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergencePD.h"
#include "MooseMesh.h"
#include "Material.h"
#include "Assembly.h"

template<>
InputParameters validParams<StressDivergencePD>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the variable this kernel acts on (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements", "The displacement variables");
  params.addCoupledVar("temp", "The temperature");
  params.addCoupledVar("bond_status", "Auxiliary variable for bond status");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

StressDivergencePD::StressDivergencePD(const InputParameters & parameters)
  :Kernel(parameters),
  _bond_force(getMaterialProperty<Real>("bond_force")),
  _bond_dfdU(getMaterialProperty<Real>("bond_dfdU")),
  _bond_dfdT(getMaterialProperty<Real>("bond_dfdT")),
  _component(getParam<unsigned int>("component")),
  _ndisp(coupledComponents("displacements")),
  _temp_coupled(isCoupled("temp")),
  _temp_var(_temp_coupled ? coupled("temp") : 0),
  _bond_status(coupledValue("bond_status")),
  _orientation(NULL)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var.push_back(coupled("displacements", i));
}

void
StressDivergencePD::initialSetup()
{
  _orientation = &_subproblem.assembly(_tid).getFE(FEType(), 1)->get_dxyzdxi();
}

void
StressDivergencePD::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  mooseAssert(re.size() == 2, "Truss elements must have two nodes");
  _local_re.resize(re.size());
  _local_re.zero();

  RealGradient ori((*_orientation)[0]);
  ori /= ori.size();
  VectorValue<Real> force_local = _bond_force[0] * ori * _bond_status[0];
  _local_re(0) = - force_local(_component);
  _local_re(1) = - _local_re(0);

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _save_in.size(); i++)
      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
  }
}

void
StressDivergencePD::computeStiffness(DenseVector<Real> & stiff_elem)
{
  RealGradient ori((*_orientation)[0]);
  Real dist = 2.0 * ori.size(); //ori.size() only gives half of the actual distance between two nodes
  ori /= ori.size();

  //the effect of truss ori change has been accounted for
  stiff_elem(0) = ori(0) * ori(0) * _bond_dfdU[0] + _bond_force[0] * (1.0 - ori(0) * ori(0)) / dist;
  stiff_elem(1) = ori(1) * ori(1) * _bond_dfdU[0] + _bond_force[0] * (1.0 - ori(1) * ori(1)) / dist;
  stiff_elem(2) = ori(2) * ori(2) * _bond_dfdU[0] + _bond_force[0] * (1.0 - ori(2) * ori(2)) / dist;
}

void
StressDivergencePD::computeOffDiagStiffness(DenseMatrix<Real> & off_stiff_elem)
{
  RealGradient ori((*_orientation)[0]);
  Real dist = 2.0 * ori.size(); //ori.size() only gives half of the actual distance between two nodes
  ori /= ori.size();

  //the effect of truss ori change has been accounted for
  off_stiff_elem(0, 0) = ori(0) * ori(0) * _bond_dfdU[0] + _bond_force[0] * (1.0 - ori(0) * ori(0)) / dist;
  off_stiff_elem(0, 1) = ori(0) * ori(1) * _bond_dfdU[0] - _bond_force[0] * ori(0) * ori(1) / dist;
  off_stiff_elem(0, 2) = ori(0) * ori(2) * _bond_dfdU[0] - _bond_force[0] * ori(0) * ori(2) / dist;
  off_stiff_elem(0, 3) = ori(0) * _bond_dfdT[0];
  off_stiff_elem(1, 0) = ori(1) * ori(0) * _bond_dfdU[0] - _bond_force[0] * ori(1) * ori(0) / dist;
  off_stiff_elem(1, 1) = ori(1) * ori(1) * _bond_dfdU[0] + _bond_force[0] * (1.0 - ori(1) * ori(1)) / dist;
  off_stiff_elem(1, 2) = ori(1) * ori(2) * _bond_dfdU[0] - _bond_force[0] * ori(1) * ori(2) / dist;
  off_stiff_elem(1, 3) = ori(1) * _bond_dfdT[0];
  off_stiff_elem(2, 0) = ori(2) * ori(0) * _bond_dfdU[0] - _bond_force[0] * ori(2) * ori(0) / dist;
  off_stiff_elem(2, 1) = ori(2) * ori(1) * _bond_dfdU[0] - _bond_force[0] * ori(2) * ori(1) / dist;
  off_stiff_elem(2, 2) = ori(2) * ori(2) * _bond_dfdU[0] + _bond_force[0] * (1.0 - ori(2) * ori(2)) / dist;
  off_stiff_elem(2, 3) = ori(2) * _bond_dfdT[0];
}

void
StressDivergencePD::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  DenseVector<Real> stiff_elem(3);
  computeStiffness(stiff_elem);

  for (unsigned int i = 0; i < _test.size(); ++i)
    for (unsigned int j = 0; j < _phi.size(); ++j)
       _local_ke(i, j) += (i == j ? 1 : -1) * stiff_elem(_component) * _bond_status[0];

  ke += _local_ke;

  if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Real> diag(rows);
    for (unsigned int i = 0; i < rows; ++i)
      diag(i) = _local_ke(i,i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
    _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
  }
}

void
StressDivergencePD::computeOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _var.number())
    computeJacobian();
  else
  {
    unsigned int coupled_component = 0;
    bool active = false;
    for (unsigned int i = 0; i < _ndisp; ++i)
      if (jvar == _disp_var[i])
      {
        coupled_component = i;
        active = true;
      }

    if (_temp_coupled && jvar == _temp_var)
    {
      coupled_component = 3;
      active = true;
    }

    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);
    if (active)
    {
      DenseMatrix<Real> off_stiff_elem(3, 4);
      computeOffDiagStiffness(off_stiff_elem);
      for (unsigned int i = 0; i < _test.size(); ++i)
        for (unsigned int j = 0; j < _phi.size(); ++j)
          ke(i, j) += (i == j ? 1 : -1) * off_stiff_elem(_component, coupled_component) * _bond_status[0];
    }
  }
}
