/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "VThermalPDMaterial.h"

template<>
InputParameters validParams<VThermalPDMaterial>()
{
  InputParameters params = validParams<ThermalPDMaterial>();
  return params;
}

VThermalPDMaterial::VThermalPDMaterial(const InputParameters & parameters) :
  ThermalPDMaterial(parameters)
{
}

Real
VThermalPDMaterial::computeBondModulus()
{
  double val = _pddim * _kappa * (1.0 / _nvsum_i + 1.0 / _nvsum_j) / _origin_length;

  return val;
}
