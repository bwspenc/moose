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
  _bond_force_i_j(declareProperty<Real>("bond_force_i_j")),
  _bond_dfdU_ij(declareProperty<Real>("bond_dfdU_ij")),
  _bond_dfdU_i_j(declareProperty<Real>("bond_dfdU_i_j")),
  _bond_dfdE_ij(declareProperty<Real>("bond_dfdE_ij")),
  _bond_dfdE_i_j(declareProperty<Real>("bond_dfdE_i_j")),
  _bond_dfdT_ij(declareProperty<Real>("bond_dfdT_ij")),
  _bond_dfdT_i_j(declareProperty<Real>("bond_dfdT_i_j"))
{
}

void
ElasticSPDMaterial::computeQpForce()
{
  // assign bond constants values
  double zero = computeBondModulus();

  _bond_force_ij[_qp] = 2.0 * _b * (_bond_elastic_strain[_qp] + _poissons_ratio * (_strain_zz_i + _strain_zz_j) / 2.0) / _origin_length * _nv_i * _nv_j;
  _bond_force_i_j[0] = 2.0 * _a * _d_i * _d_i * (_bond_elastic_strain[_qp] + _alpha * ((_temp_i + _temp_j) / 2.0 - _temp_ref) - _alpha * (_temp_i - _temp_ref) + _poissons_ratio * _strain_zz_i) * _nv_i * _nv_j;
  _bond_force_i_j[1] = 2.0 * _a * _d_j * _d_j * (_bond_elastic_strain[_qp] + _alpha * ((_temp_i + _temp_j) / 2.0 - _temp_ref) - _alpha * (_temp_j - _temp_ref) + _poissons_ratio * _strain_zz_j) * _nv_i * _nv_j;

  _bond_dfdU_ij[_qp] = 2.0 * _b / _origin_length / _origin_length * _nv_i * _nv_j;
  _bond_dfdU_i_j[0] = 2.0 * _a * _d_i * _d_i / _origin_length * _nv_i * _nv_j;
  _bond_dfdU_i_j[1] = 2.0 * _a * _d_j * _d_j / _origin_length * _nv_i * _nv_j;

  _bond_dfdE_ij[_qp] = 2.0 * _b * _poissons_ratio / 2.0 / _origin_length * _nv_i * _nv_j;
  _bond_dfdE_i_j[0] = 2.0 * _a * _d_i * _d_i * _poissons_ratio * _nv_i * _nv_j;
  _bond_dfdE_i_j[1] = 2.0 * _a * _d_j * _d_j * _poissons_ratio * _nv_i * _nv_j;

  _bond_dfdT_ij[_qp] = - 2.0 * _b * (1.0 + _poissons_ratio) * _alpha / 2.0 / _origin_length * _nv_i * _nv_j;
  _bond_dfdT_i_j[0] = - 2.0 * _a * _d_i * _d_i * (1.0 + _poissons_ratio) * _alpha * _nv_i * _nv_j;
  _bond_dfdT_i_j[1] = - 2.0 * _a * _d_j * _d_j * (1.0 + _poissons_ratio) * _alpha * _nv_i * _nv_j;
}
