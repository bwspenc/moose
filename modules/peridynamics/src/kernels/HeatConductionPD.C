/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "HeatConductionPD.h"
#include "MooseMesh.h"
#include "Material.h"
#include "Assembly.h"

template<>
InputParameters validParams<HeatConductionPD>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<NonlinearVariableName>("bond_status", "Auxiliary variable for the bond failure status");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

HeatConductionPD::HeatConductionPD(const InputParameters & parameters) :
  Kernel(parameters),
  _bond_response(getMaterialProperty<Real>("bond_response")),
  _bond_drdT(getMaterialProperty<Real>("bond_drdT")),
  _aux(_fe_problem.getAuxiliarySystem()),
  _aux_sln(*_aux.currentSolution()),
  _bond_status_var(&_fe_problem.getVariable(_tid, getParam<NonlinearVariableName>("bond_status"))),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh))
{
}

HeatConductionPD::~HeatConductionPD()
{
}

void
HeatConductionPD::computeResidual()
{
  dof_id_type bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  Number bond_status = _aux_sln(bs_dof);

  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

  _local_re(0) = - _bond_response[0] * bond_status;
  _local_re(1) = - _local_re(0);

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _save_in.size(); ++i)
      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
  }
}

void
HeatConductionPD::computeJacobian()
{
  dof_id_type bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  Number bond_status = _aux_sln(bs_dof);

  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  for (unsigned int i = 0; i < _test.size(); ++i)
    for (unsigned int j = 0; j < _phi.size(); ++j)
      _local_ke(i, j) += (i == j ? 1 : -1) * _bond_drdT[0] * bond_status;

  ke += _local_ke;

  if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Number> diag(rows);
    for (unsigned int i = 0; i < rows; ++i)
      diag(i) = _local_ke(i, i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
      _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
  }
}
