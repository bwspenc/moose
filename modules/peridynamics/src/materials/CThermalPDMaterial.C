/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CThermalPDMaterial.h"

template<>
InputParameters validParams<CThermalPDMaterial>()
{
  InputParameters params = validParams<ThermalPDMaterial>();
  return params;
}

CThermalPDMaterial::CThermalPDMaterial(const InputParameters & parameters) :
  ThermalPDMaterial(parameters)
{
}

double
CThermalPDMaterial::computeBondModulus()
{
  double Kij = (6.0 * _kappa / (3.14159265358 * std::pow(_horizon_i, _pddim + 1)) + 6.0 * _kappa / (3.14159265358 * std::pow(_horizon_j, _pddim + 1))) / 2.0;

  return Kij;
}
