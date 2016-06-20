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
  double Cij;
  double h0 = _pdmesh.horizon(_current_elem->get_node(0)->id());
  double h1 = _pdmesh.horizon(_current_elem->get_node(1)->id());

  Cij = (6.0 * _pddim * _bulk_modulus / (3.14159265358 * std::pow(h0, _pddim + 1)) + 6.0 * _pddim * _bulk_modulus / (3.14159265358 * std::pow(h1, _pddim + 1))) / 2.0;

  return Cij;
}
