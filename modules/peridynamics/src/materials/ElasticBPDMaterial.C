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
ElasticBPDMaterial::computeQpForce()
{
  double Cij = computeBondModulus();

  // residuals
  if (_strain_zz_coupled)
    _bond_force_ij[_qp] = Cij * (_bond_elastic_strain[_qp] + _poissons_ratio * (_strain_zz[0] - _alpha * ((_temp_i + _temp_j) / 2.0 - _temp_ref))) * _nv_i * _nv_j;
  else
    _bond_force_ij[_qp] = Cij * _bond_elastic_strain[_qp] * _nv_i * _nv_j;

  _bond_force_i_j[0] = 0;
  _bond_force_i_j[1] = 0;

  // derivatives of residuals
  _bond_dfdU_ij[_qp] = Cij / _origin_length * _nv_i * _nv_j;
  _bond_dfdU_i_j[0] = 0;
  _bond_dfdU_i_j[1] = 0;

  _bond_dfdE_ij[_qp] = Cij * _poissons_ratio * _nv_i * _nv_j;
  _bond_dfdE_i_j[0] = 0;
  _bond_dfdE_i_j[1] = 0;

  if (_strain_zz_coupled)
    _bond_dfdT_ij[_qp] = - Cij * (1.0 + _poissons_ratio) * _alpha / 2.0 * _nv_i * _nv_j;
  else
    _bond_dfdT_ij[_qp] = - Cij * _alpha / 2.0 * _nv_i * _nv_j;

  _bond_dfdT_i_j[0] = 0;
  _bond_dfdT_i_j[1] = 0;
}
