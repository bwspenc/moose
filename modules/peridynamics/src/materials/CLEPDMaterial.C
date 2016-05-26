/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CLEPDMaterial.h"

template<>
InputParameters validParams<CLEPDMaterial>()
{
  InputParameters params = validParams<MechanicPDMaterial>();
  return params;
}

CLEPDMaterial::CLEPDMaterial(const InputParameters & parameters)
  :MechanicPDMaterial(parameters)
{
}

void
CLEPDMaterial::computeQpStrain()
{
  _bond_total_strain[_qp] = _current_length / _origin_length - 1.0;
  // bond temperature is taken as the average of two end nodes
  _bond_mechanic_strain[_qp] = _bond_total_strain[_qp] - _alpha * ((_temp[0] + _temp[1]) / 2.0 - _temp_ref);
  _bond_elastic_strain[_qp] = _bond_mechanic_strain[_qp];
}

void
CLEPDMaterial::computeQpForce()
{
  double bulk_modulus;
  // bond_force and bond_dfdU, bond_dfdT
  if (_pddim == 2)
  {
    if (_plane_strain)
      bulk_modulus = _youngs_modulus / 2.0 / (1.0 + _poissons_ratio) / (1.0 - 2.0 * _poissons_ratio);// plane strain
    else
      bulk_modulus = _youngs_modulus / 2.0 / (1.0 - _poissons_ratio);// plane stress

    _bond_force[_qp] = 12.0 * bulk_modulus / (3.14159265358 * std::pow(_horizon, 3)) * _bond_mechanic_strain[_qp] * _nv_i * _nv_j;
    _bond_dfdU[_qp] = 12.0 * bulk_modulus / (3.14159265358 * std::pow(_horizon, 3)) / _origin_length * _nv_i * _nv_j;
    _bond_dfdT[_qp] = - 12.0 * bulk_modulus / (3.14159265358 * std::pow(_horizon, 3)) * _alpha / 2.0 * _nv_i * _nv_j;
  }
  else if (_pddim == 3)
  {
    bulk_modulus = _youngs_modulus / 3.0 / (1.0 - 2.0 * _poissons_ratio);
    _bond_force[_qp] = 18.0 * bulk_modulus / (3.14159265358 * std::pow(_horizon, 4)) * _bond_mechanic_strain[_qp] * _nv_i * _nv_j;
    _bond_dfdU[_qp] = 18.0 * bulk_modulus / (3.14159265358 * std::pow(_horizon, 4)) / _origin_length * _nv_i * _nv_j;
    _bond_dfdT[_qp] = - 18.0 * bulk_modulus / (3.14159265358 * std::pow(_horizon, 4)) * _alpha / 2.0 * _nv_i * _nv_j;
  }
}

