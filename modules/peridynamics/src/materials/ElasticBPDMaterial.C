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
  MechanicPDMaterial(parameters),
  _bond_force(declareProperty<Real>("bond_force")),
  _bond_dfdU(declareProperty<Real>("bond_dfdU")),
  _bond_dfdE(declareProperty<Real>("bond_dfdE")),
  _bond_dfdT(declareProperty<Real>("bond_dfdT"))
{
}

void
ElasticBPDMaterial::computeQpForce()
{
  double Cij = computeBondModulus();

  _bond_force[_qp] = Cij * (_bond_elastic_strain[_qp] + _poissons_ratio * (_strain_zz_i + _strain_zz_j) / 2.0) * _nv_i * _nv_j;

  _bond_dfdU[_qp] = Cij / _origin_length * _nv_i * _nv_j;
  _bond_dfdE[_qp] = Cij * _poissons_ratio / 2.0 * _nv_i * _nv_j;
  _bond_dfdT[_qp] = - Cij * (1.0 + _poissons_ratio) * _alpha / 2.0 * _nv_i * _nv_j;
}
