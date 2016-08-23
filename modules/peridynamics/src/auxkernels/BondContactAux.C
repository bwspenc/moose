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
  params.addCoupledVar("bond_status", "Auxiliary variable for bond status");
  params.addCoupledVar("bond_critical_strain", "Auxiliary variable for bond critical strain");
  return params;
}

BondContactAux::BondContactAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _bond_elastic_strain(getMaterialProperty<Real>("bond_elastic_strain")),
  _bond_critical_strain(coupledValue("bond_critical_strain")),
  _bond_status(coupledValue("bond_status"))
{
}

Real
BondContactAux::computeValue()
{
//  if (_bond_status[0] < 0.01 && _bond_elastic_strain[0] < _bond_critical_strain[0] / 5.0)
  if (_bond_elastic_strain[0] < _bond_critical_strain[0] / 5.0)
    return 1.0;
  else
    return 0.0;
}
