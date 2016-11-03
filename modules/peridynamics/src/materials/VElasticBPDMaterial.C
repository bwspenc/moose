/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "VElasticBPDMaterial.h"

template<>
InputParameters validParams<VElasticBPDMaterial>()
{
  InputParameters params = validParams<ElasticBPDMaterial>();
  return params;
}

VElasticBPDMaterial::VElasticBPDMaterial(const InputParameters & parameters) :
  ElasticBPDMaterial(parameters)
{
}

Real
VElasticBPDMaterial::computeBondModulus()
{
  double val = _pddim * _pddim * _bulk_modulus * (1.0 / _nvsum_i + 1.0 / _nvsum_j) / _origin_length;

  return val;
}
