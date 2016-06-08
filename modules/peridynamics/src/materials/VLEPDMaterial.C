/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "VLEPDMaterial.h"

template<>
InputParameters validParams<VLEPDMaterial>()
{
  InputParameters params = validParams<MechanicPDMaterial>();
  return params;
}

VLEPDMaterial::VLEPDMaterial(const InputParameters & parameters)
  :MechanicPDMaterial(parameters)
{
}

double
VLEPDMaterial::computeBondModulus()
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

  double Cij;
  if (_pddim == 2)
  {
    if (_plane_strain)
      Cij = 4.0 * (_youngs_modulus / 2.0 / (1.0 + _poissons_ratio) / (1.0 - 2.0 * _poissons_ratio)) * val;// plane strain
    else
      Cij = 4.0 * (_youngs_modulus / 2.0 / (1.0 - _poissons_ratio)) * val;// plane stress
  }
  else if (_pddim == 3)
    Cij = 9.0 * (_youngs_modulus / 3.0 / (1.0 - 2.0 * _poissons_ratio)) * val;

  return Cij;
}

void
VLEPDMaterial::computeQpStrain()
{
  _bond_total_strain[_qp] = _current_length / _origin_length - 1.0;
  // bond temperature is taken as the average of two end nodes
  _bond_mechanic_strain[_qp] = _bond_total_strain[_qp] - _alpha * ((_temp[0] + _temp[1]) / 2.0 - _temp_ref);
  _bond_elastic_strain[_qp] = _bond_mechanic_strain[_qp];
}

void
VLEPDMaterial::computeQpForce()
{
  double Cij = computeBondModulus();
  // bond_force and bond_dfdU, bond_dfdT
  _bond_force[_qp] = Cij * _bond_mechanic_strain[_qp] * _nv_j * _nv_i;
  _bond_dfdU[_qp] =  Cij / _origin_length * _nv_j * _nv_i;
  _bond_dfdT[_qp] = - Cij * _alpha / 2.0 * _nv_j * _nv_i;
}

