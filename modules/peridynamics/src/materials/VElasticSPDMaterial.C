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
  dof_id_type node_i = _current_elem->get_node(0)->id();
  dof_id_type node_j = _current_elem->get_node(1)->id();

  std::vector<dof_id_type> i_neighbors = _pdmesh.neighbors(node_i);
  std::vector<dof_id_type> j_neighbors = _pdmesh.neighbors(node_j);
  unsigned int i_nneighbor = _pdmesh.n_neighbors(node_i);
  unsigned int j_nneighbor = _pdmesh.n_neighbors(node_j);

  double val1 = 0, val2 = 0;
  for (unsigned int k = 0; k < i_nneighbor; ++k)
    val1 += _pdmesh.volume(i_neighbors[k]);

  for (unsigned int k = 0; k < j_nneighbor; ++k)
    val2 += _pdmesh.volume(j_neighbors[k]);

  double val = (1.0 / val1 + 1.0 / val2) / _origin_length;

  return _pddim * _pddim * _bulk_modulus * val;
}
