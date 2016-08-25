/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                         Peridynamics                         */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "BondStatusUO.h"
#include "PeridynamicMesh.h"

template<>
InputParameters validParams<BondStatusUO>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addCoupledVar("bond_status", "Auxiliary variable for failure status of each bond");
  return params;
}

BondStatusUO::BondStatusUO(const InputParameters & parameters) :
  ElementUserObject(parameters),
  _aux(_fe_problem.getAuxiliarySystem()),
  _bond_elastic_strain(getMaterialProperty<Real>("bond_elastic_strain")),
  _bond_critical_strain(getMaterialProperty<Real>("bond_critical_strain")),
  _bond_status_var(getVar("bond_status", 0))
{
}

void
BondStatusUO::initialize()
{
}

void
BondStatusUO::execute()
{
  NumericVector<Number> & sln = _aux.solution();
  long int bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  unsigned int bond_status = sln(bs_dof);

//std::cout << bond_status << _bond_elastic_strain[0] << " : " << _bond_critical_strain[0] << std::endl;

  if (std::abs(bond_status - 1.0) < 0.01 && _bond_elastic_strain[0] < _bond_critical_strain[0])
    _aux.solution().set(bs_dof, 1.0);
  else
    _aux.solution().set(bs_dof, 0.0);
}

void
BondStatusUO::threadJoin(const UserObject & uo)
{
}

void
BondStatusUO::finalize()
{
  _aux.solution().close();
}
