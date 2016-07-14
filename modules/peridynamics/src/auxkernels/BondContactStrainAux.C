/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "BondContactStrainAux.h"

template<>
InputParameters validParams<BondContactStrainAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addCoupledVar("bond_critical_strain", "Auxiliary variable for bond critical strain");
  return params;
}

BondContactStrainAux::BondContactStrainAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _bond_critical_strain(coupledValue("bond_critical_strain"))
{
}

Real
BondContactStrainAux::computeValue()
{
// Critical strain for contact
  return _bond_critical_strain[0] / 5.0;
}
