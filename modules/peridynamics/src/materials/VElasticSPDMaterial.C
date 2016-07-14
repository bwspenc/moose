/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "VElasticSPDMaterial.h"

template<>
InputParameters validParams<VElasticSPDMaterial>()
{
  InputParameters params = validParams<ElasticSPDMaterial>();
  return params;
}

VElasticSPDMaterial::VElasticSPDMaterial(const InputParameters & parameters) :
  ElasticSPDMaterial(parameters)
{
}

Real
VElasticSPDMaterial::computeBondModulus()
{
  _a = 0.5 * (_bulk_modulus  -  (8.0 - _pddim) / 3.0 * _shear_modulus);

  // _b = _bij * _horizon_i + _bji * _horizon_j
  _b = _pddim * _pddim * (_bulk_modulus / 2.0 - _a) * (1.0 / _nvsum_i + 1.0 / _nvsum_j);

  // _d_i = _di * _horizon_i = _pddim / _nvsum_i
  _d_i = _pddim / _nvsum_i;
  _d_j = _pddim / _nvsum_j;

  return 0;
}
