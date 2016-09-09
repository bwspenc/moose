/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergenceBPD.h"
#include "Assembly.h"

template<>
InputParameters validParams<StressDivergenceBPD>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("BPD Stress divergence kernel");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the variable this kernel acts on (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements", "Variable for the displacements");
  params.addCoupledVar("temp", "Variable for the temperature");
  params.addRequiredParam<NonlinearVariableName>("bond_status", "Auxiliary variable for the bond failure status");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

StressDivergenceBPD::StressDivergenceBPD(const InputParameters & parameters) :
  Kernel(parameters),
  _bond_force_ij(getMaterialProperty<Real>("bond_force_ij")),
  _bond_dfdU_ij(getMaterialProperty<Real>("bond_dfdU_ij")),
  _bond_dfdT_ij(getMaterialProperty<Real>("bond_dfdT_ij")),
  _aux(_fe_problem.getAuxiliarySystem()),
  _aux_sln(*_aux.currentSolution()),
  _component(getParam<unsigned int>("component")),
  _ndisp(coupledComponents("displacements")),
  _temp_coupled(isCoupled("temp")),
  _temp_var(_temp_coupled ? coupled("temp") : 0),
  _bond_status_var(&_fe_problem.getVariable(_tid, getParam<NonlinearVariableName>("bond_status"))),
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
  // bond status for bond ij
  dof_id_type bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  Number bond_status_ij = _aux_sln(bs_dof);

  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  mooseAssert(re.size() == 2, "Truss element has only two nodes");
  _local_re.resize(re.size());
  _local_re.zero();

  RealGradient ori = (*_orientation)[0];
  ori /= ori.size();
  _local_re(0) = - _bond_force_ij[0] * ori(_component) * bond_status_ij;
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
  dof_id_type bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  Number bond_status_ij = _aux_sln(bs_dof);

  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  RealGradient ori = (*_orientation)[0];
  double current_length = 2.0 * ori.size(); // ori.size() only gives half current length
  ori /= ori.size();

  double diag = ori(_component) * ori(_component) * _bond_dfdU_ij[0] + _bond_force_ij[0] * (1.0 - ori(_component) * ori(_component)) / current_length;

  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < _phi.size(); _j++)
       _local_ke(_i, _j) += (_i == _j ? 1 : -1) * diag * bond_status_ij;

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
    {
      if (jvar == _disp_var[i])
      {
        coupled_component = i;
        active = true;
      }
    }
    if (_temp_coupled && jvar == _temp_var)
    {
      coupled_component = 3;
      active = true;
    }

    dof_id_type bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
    Number bond_status_ij = _aux_sln(bs_dof);

    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);
    if (active)
    {
      RealGradient ori = (*_orientation)[0];
      double current_length = 2.0 * ori.size(); // ori.size() only gives half current length
      ori /= ori.size();

      if (coupled_component == 3)
      {
        for (_i = 0; _i < _test.size(); _i++)
          for (_j = 0; _j < _phi.size(); _j++)
            ke(_i, _j) += (_i == 1 ? 1 : -1) * ori(_component) * _bond_dfdT_ij[0] * bond_status_ij;
      }
      else
      {
        double off_diag = ori(_component) * ori(coupled_component) * _bond_dfdU_ij[0] - _bond_force_ij[0] * ori(_component) * ori(coupled_component) / current_length;
        for (_i = 0; _i < _test.size(); _i++)
          for (_j = 0; _j < _phi.size(); _j++)
            ke(_i, _j) += (_i == _j ? 1 : -1) * off_diag * bond_status_ij;
      }
    }
  }
}
