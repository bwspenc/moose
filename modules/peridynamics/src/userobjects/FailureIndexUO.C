/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                         Peridynamics                         */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "FailureIndexUO.h"
#include "PeridynamicMesh.h"

template<>
InputParameters validParams<FailureIndexUO>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addCoupledVar("intact_bonds", "Auxiliary variable for number of intact bonds at each node");
  params.addCoupledVar("bond_status", "Auxiliary variable for bond status");
  return params;
}

FailureIndexUO::FailureIndexUO(const InputParameters & parameters) :
  ElementUserObject(parameters),
  _aux(_fe_problem.getAuxiliarySystem()),
  _intact_bonds_var(getVar("intact_bonds", 0)),
  _bond_status_var(getVar("bond_status", 0))
{
}

void
FailureIndexUO::initialize()
{
  std::vector<std::string> zero_vars;
  zero_vars.push_back("intact_bonds");
  _aux.zeroVariables(zero_vars);
}

void
FailureIndexUO::execute()
{
  NumericVector<Number> & sln = _aux.solution();
  long int bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  unsigned int bond_status = sln(bs_dof);

  long int ib_dof0 = _current_elem->get_node(0)->dof_number(_aux.number(), _intact_bonds_var->number(), 0);
  long int ib_dof1 = _current_elem->get_node(1)->dof_number(_aux.number(), _intact_bonds_var->number(), 0);

  if (std::abs(bond_status - 1.0) < 0.01)
  {
    sln.add(ib_dof0, 1.0);
    sln.add(ib_dof1, 1.0);
  }
  else
  {
    sln.add(ib_dof0, 0.0);
    sln.add(ib_dof1, 0.0);
  }
}

void
FailureIndexUO::threadJoin(const UserObject & uo)
{
}

void
FailureIndexUO::finalize()
{
  _aux.solution().close();
}

Real
FailureIndexUO::computeFailureIndex(unsigned int node_id) const
{
  NumericVector<Number> & sln = _aux.solution();
  long int ib_dof = _mesh.nodePtr(node_id)->dof_number(_aux.number(), _intact_bonds_var->number(), 0);
  // get the total bonds number for node node_id
  PeridynamicMesh & pdmesh = dynamic_cast<PeridynamicMesh &>(_mesh);
  unsigned int tb = pdmesh.n_neighbors(node_id);

  return 1.0 - sln(ib_dof) / tb;
}
