/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "BondStatusAux.h"

template<>
InputParameters validParams<BondStatusAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addCoupledVar("bond_status", "Auxiliary variable for bond status");
  params.addCoupledVar("bond_critical_strain", "Auxiliary variable for bond critical strain");
  return params;
}

BondStatusAux::BondStatusAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _bond_mechanic_strain(getMaterialProperty<Real>("bond_mechanic_strain")),
  _bond_critical_strain(coupledValue("bond_critical_strain")),
  _bond_status(coupledValue("bond_status"))
{
}

Real
BondStatusAux::computeValue()
{
  if (std::abs(_bond_status[0] - 1.0) < 0.01 && std::abs(_bond_mechanic_strain[0]) < _bond_critical_strain[0])
    return 1.0;
  else
    return 0.0;
}
