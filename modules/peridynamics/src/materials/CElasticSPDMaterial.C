/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CElasticSPDMaterial.h"

template<>
InputParameters validParams<CElasticSPDMaterial>()
{
  InputParameters params = validParams<ElasticSPDMaterial>();
  return params;
}

CElasticSPDMaterial::CElasticSPDMaterial(const InputParameters & parameters) :
  ElasticSPDMaterial(parameters)
{
}

Real
CElasticSPDMaterial::computeBondModulus()
{
  _a = 0.5 * (_bulk_modulus  -  (8.0 - _pddim) / 3.0 * _shear_modulus);

  // _b = 2 * _b * _horizon_(i/j) * _origin_length //_origin_length will be cancelled out in parent material model
  _b = _origin_length * (3.0 * _pddim + 6.0) * _shear_modulus / 3.1415926 / std::pow(_horizon_i, _pddim + 1.0);

  // _d_i = _di * _horizon_(i/j)
  _d_i = (_pddim / 4.0 + 1.5) / 3.1415926 / std::pow(_horizon_i, _pddim);
  _d_j = _d_i;

  return 0;
}
