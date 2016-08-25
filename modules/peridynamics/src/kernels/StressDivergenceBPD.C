/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergenceBPD.h"
#include "MooseMesh.h"
#include "Material.h"
#include "Assembly.h"

template<>
InputParameters validParams<StressDivergenceBPD>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("BPD Stress divergence kernel");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the variable this kernel acts on (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements", "Coupled variable for the displacements");
  params.addCoupledVar("temp", "Coupled variable for the temperature");
  params.addCoupledVar("strain_zz", "Coupled variable for the strain_zz");
  params.addCoupledVar("bond_status", "Auxiliary variable for failure status of each bond");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

StressDivergenceBPD::StressDivergenceBPD(const InputParameters & parameters) :
  Kernel(parameters),
  _bond_force(getMaterialProperty<Real>("bond_force")),
  _bond_dfdU(getMaterialProperty<Real>("bond_dfdU")),
  _bond_dfdE(getMaterialProperty<Real>("bond_dfdE")),
  _bond_dfdT(getMaterialProperty<Real>("bond_dfdT")),
  _aux(_fe_problem.getAuxiliarySystem()),
  _component(getParam<unsigned int>("component")),
  _ndisp(coupledComponents("displacements")),
  _temp_coupled(isCoupled("temp")),
  _temp_var(_temp_coupled ? coupled("temp") : 0),
  _strain_zz_coupled(isCoupled("strain_zz")),
  _strain_zz_var(_strain_zz_coupled ? coupled("strain_zz") : 0),
  _bond_status(coupledValue("bond_status")),
  _bond_status_var(getVar("bond_status", 0)),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh)),
  _orientation(NULL)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var.push_back(coupled("displacements", i));
}

void
StressDivergenceBPD::initialSetup()
{
  _orientation = &_assembly.getFE(FEType(), 1)->get_dxyzdxi();
}

void
StressDivergenceBPD::computeResidual()
{
  NumericVector<Number> & sln = _aux.solution();
  dof_id_type bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  unsigned int bond_status = sln(bs_dof);

//std::cout<< _bond_status[0] <<" : " << bond_status << std::endl;

  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  mooseAssert(re.size() == 2, "Truss element has only two nodes");
  _local_re.resize(re.size());
  _local_re.zero();

  RealGradient ori = (*_orientation)[0];
  ori /= ori.size();
  _local_re(0) = - _bond_force[0] * ori(_component) * bond_status;
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
StressDivergenceBPD::computeJacobian()
{
  NumericVector<Number> & sln = _aux.solution();
  long int bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  unsigned int bond_status = sln(bs_dof);

  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  RealGradient ori = (*_orientation)[0];
  double current_length = 2.0 * ori.size(); // ori.size() only gives half current length
  ori /= ori.size();

  double stiff_elem = ori(_component) * ori(_component) * _bond_dfdU[0] + _bond_force[0] * (1.0 - ori(_component) * ori(_component)) / current_length;

  for (unsigned int i = 0; i < _test.size(); ++i)
    for (unsigned int j = 0; j < _phi.size(); ++j)
       _local_ke(i, j) += (i == j ? 1 : -1) * stiff_elem * bond_status;

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
StressDivergenceBPD::computeOffDiagJacobian(unsigned int jvar)
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

    if (_strain_zz_coupled && jvar == _strain_zz_var)
    {
      coupled_component = 4;
      active = true;
    }

    NumericVector<Number> & sln = _aux.solution();
    long int bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
    unsigned int bond_status = sln(bs_dof);

    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);
    if (active)
    {
      RealGradient ori = (*_orientation)[0];
      double current_length = 2.0 * ori.size(); // ori.size() only gives half current length
      ori /= ori.size();

      double off_stiff_elem;
      if (coupled_component == 3)
      {
        off_stiff_elem = ori(_component) * _bond_dfdT[0];
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            ke(i, j) += (i == 1 ? 1 : -1) * off_stiff_elem * bond_status;
      }
      else if (coupled_component == 4)
      {
        off_stiff_elem = ori(_component) * _bond_dfdE[0];
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            ke(i, j) += (i == 1 ? 1 : -1) * off_stiff_elem * bond_status;
      }
      else
      {
        off_stiff_elem = ori(_component) * ori(coupled_component) * _bond_dfdU[0] - _bond_force[0] * ori(_component) * ori(coupled_component) / current_length;
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            ke(i, j) += (i == j ? 1 : -1) * off_stiff_elem * bond_status;
      }
    }
  }
}
