/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ElasticBPDMaterial.h"

template<>
InputParameters validParams<ElasticBPDMaterial>()
{
  InputParameters params = validParams<MechanicPDMaterial>();
  return params;
}

ElasticBPDMaterial::ElasticBPDMaterial(const InputParameters & parameters) :
  MechanicPDMaterial(parameters)
{
}

void
ElasticBPDMaterial::computeQpStrain()
{
  _bond_total_strain[_qp] = _current_length / _origin_length - 1.0;
  // bond temperature is taken as the average of two end nodes
  _bond_mechanic_strain[_qp] = _bond_total_strain[_qp] - _alpha * ((_temp[0] + _temp[1]) / 2.0 - _temp_ref);
  _bond_elastic_strain[_qp] = _bond_mechanic_strain[_qp];
}

void
ElasticBPDMaterial::computeQpForce()
{
  double Cij = computeBondModulus();
  // bond_force, bond_dfdU and bond_dfdT
//  _bond_force[_qp] = Cij * (_bond_mechanic_strain[_qp] - _bond_contact[0] * _bond_contact_strain[0]) * _nv_i * _nv_j * _bond_sign;
  _bond_force[_qp] = Cij * _bond_mechanic_strain[_qp] * _nv_i * _nv_j * _bond_sign;
  _bond_dfdU[_qp] = Cij / _origin_length * _nv_i * _nv_j * _bond_sign;
  _bond_dfdT[_qp] = - Cij * _alpha / 2.0 * _nv_i * _nv_j * _bond_sign;
}
