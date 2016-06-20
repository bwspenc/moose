/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "BondContactAux.h"

template<>
InputParameters validParams<BondContactAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addCoupledVar("bond_critical_strain", "Auxiliary variable for bond critical strain");
  return params;
}

BondContactAux::BondContactAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _bond_mechanic_strain(getMaterialProperty<Real>("bond_mechanic_strain")),
  _bond_critical_strain(coupledValue("bond_critical_strain"))
{
}

Real
BondContactAux::computeValue()
{
  if (_bond_mechanic_strain[0] < _bond_critical_strain[0] / 10.0)
    return 1.0;
  else
    return 0.0;
}
