/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                         Peridynamics                         */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "FailureIndex.h"
#include "PeridynamicMesh.h"

template<>
InputParameters validParams<FailureIndex>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addCoupledVar("intact_bonds", "Auxiliary variable for number of intact bonds at each node");
  params.addCoupledVar("bond_status", "Auxiliary variable for bond status");
  params.addCoupledVar("bond_critical_strain", "Auxiliary variable for bond critical strain");
  return params;
}

FailureIndex :: FailureIndex(const InputParameters & parameters) :
  ElementUserObject(parameters),
  _aux(_fe_problem.getAuxiliarySystem()),
  _intact_bonds_var(getVar("intact_bonds",0)),
  _bond_mechanic_strain(getMaterialProperty<Real>("bond_mechanic_strain")),
  _bond_critical_strain(coupledValue("bond_critical_strain")),
  _bond_status(coupledValue("bond_status"))
{
}

FailureIndex::~FailureIndex()
{
}

void
FailureIndex::initialize()
{
  std::vector<std::string> zero_vars;
  zero_vars.push_back("intact_bonds");
  _aux.zeroVariables(zero_vars);
}

void
FailureIndex::execute()
{
  long int ib_dof0 = _current_elem->get_node(0)->dof_number(_aux.number(), _intact_bonds_var->number(), 0);
  long int ib_dof1 = _current_elem->get_node(1)->dof_number(_aux.number(), _intact_bonds_var->number(), 0);

  if (std::abs(_bond_status[0] - 1.0) < 0.01 && std::abs(_bond_mechanic_strain[0]) < _bond_critical_strain[0])
  {
    _aux.solution().add(ib_dof0, 1.0);
    _aux.solution().add(ib_dof1, 1.0);
  }
  else
  {
    _aux.solution().add(ib_dof0, 0.0);
    _aux.solution().add(ib_dof1, 0.0);
  }
}

void
FailureIndex::threadJoin(const UserObject & u )
{
}

void
FailureIndex::finalize()
{
  _aux.solution().close();
}

Real
FailureIndex::computeFailureIndex(unsigned int node_id) const
{
  NumericVector<Number> & sln = _aux.solution();
  long int ib_dof = _mesh.nodePtr(node_id)->dof_number(_aux.number(), _intact_bonds_var->number(), 0);
  // get the total bonds number for node node_id
  PeridynamicMesh & pdmesh = dynamic_cast<PeridynamicMesh &>(_mesh);
  unsigned int tb = pdmesh.n_neighbors(node_id);

  return 1.0 - sln(ib_dof)/tb;
}

