/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CElasticBPDMaterial.h"

template<>
InputParameters validParams<CElasticBPDMaterial>()
{
  InputParameters params = validParams<ElasticBPDMaterial>();
  return params;
}

CElasticBPDMaterial::CElasticBPDMaterial(const InputParameters & parameters) :
  ElasticBPDMaterial(parameters)
{
}

Real
CElasticBPDMaterial::computeBondModulus()
{
  double Cij = (6.0 * _pddim * _bulk_modulus / (3.14159265358 * std::pow(_horizon_i, _pddim + 1)) + 6.0 * _pddim * _bulk_modulus / (3.14159265358 * std::pow(_horizon_j, _pddim + 1))) / 2.0;

  return Cij;
}
