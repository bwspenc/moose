/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ElasticSPDMaterial.h"

template<>
InputParameters validParams<ElasticSPDMaterial>()
{
  InputParameters params = validParams<MechanicPDMaterial>();
  return params;
}

ElasticSPDMaterial::ElasticSPDMaterial(const InputParameters & parameters) :
  MechanicPDMaterial(parameters),
  _bond_force_ij(declareProperty<Real>("bond_force_ij")),
  _bond_force_i(declareProperty<Real>("bond_force_i")),
  _bond_force_j(declareProperty<Real>("bond_force_j")),
  _bond_dfdU_ij(declareProperty<Real>("bond_dfdU_ij")),
  _bond_dfdU_i(declareProperty<Real>("bond_dfdU_i")),
  _bond_dfdU_j(declareProperty<Real>("bond_dfdU_j")),
  _bond_dfdT_ij(declareProperty<Real>("bond_dfdT_ij"))
{
}

void
ElasticSPDMaterial::computeQpForce()
{
  // assign bond constants values
  double zero = computeBondModulus();

  // need to use the temperature at the two end nodes
  _bond_force_ij[_qp] = (- 2.0 * _pddim * _alpha * _a * ( _d_i * (_temp[0] - _temp_ref) + _d_j * (_temp[1] - _temp_ref)) + 2.0 * _b * _bond_mechanic_strain[_qp] ) / _origin_length * _nv_i * _nv_j;
  _bond_force_i[_qp] = 2.0 * _a * _d_i * _d_i * _bond_total_strain[_qp] * _nv_i * _nv_j;
  _bond_force_j[_qp] = 2.0 * _a * _d_j * _d_j * _bond_total_strain[_qp] * _nv_i * _nv_j;

  _bond_dfdU_ij[_qp] = 2.0 * _b / _origin_length / _origin_length * _nv_i * _nv_j;
  _bond_dfdU_i[_qp] = 2.0 * _a * _d_i * _d_i / _origin_length * _nv_i * _nv_j;
  _bond_dfdU_j[_qp] = 2.0 * _a * _d_j * _d_j / _origin_length * _nv_i * _nv_j;

  _bond_dfdT_ij[_qp] = - 2.0 * _b / _origin_length * _alpha / 2.0 * _nv_i * _nv_j;
}
