/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
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
  // bond_force and bond_dfdU, bond_dfdT
  if (_pddim == 2)
  {
    double lambda = computeLambda2D(_poissons_ratio);
    _bond_force[_qp] = (_youngs_modulus * std::exp(-_origin_length / _horizon) / _mesh_spacing / lambda) * _bond_mechanic_strain[_qp] * _nv_i * _nv_j;
    _bond_dfdU[_qp] = (_youngs_modulus * std::exp(-_origin_length / _horizon) / _mesh_spacing / lambda) / _origin_length * _nv_i * _nv_j;
    _bond_dfdT[_qp] = -(_youngs_modulus * std::exp(-_origin_length / _horizon) / _mesh_spacing / lambda) * _alpha / 2.0 * _nv_i * _nv_j;
  }
  else if (_pddim == 3)
  {
    double lambda = computeLambda3D(_poissons_ratio);
    _bond_force[_qp] = (_youngs_modulus * std::exp(-_origin_length / _horizon) / _mesh_spacing / lambda) * _bond_mechanic_strain[_qp] * _nv_i * _nv_j;
    _bond_dfdU[_qp] = (_youngs_modulus * std::exp(-_origin_length / _horizon) / _mesh_spacing / lambda) / _origin_length * _nv_i * _nv_j;
    _bond_dfdT[_qp] = -(_youngs_modulus * std::exp(-_origin_length / _horizon) / _mesh_spacing / lambda) * _alpha / 2.0 * _nv_i * _nv_j;
  }
}

double
VLEPDMaterial::computeLambda2D(double poissons_ratio)
{
  double e = 0.0001, lambda = 0, ftemp, s;
  for (unsigned int i = -3; i < 4; ++i)
    for (unsigned int j = 1; j < 4; ++j)
    {
      ftemp = std::sqrt(i * i + j * j);
      if (ftemp <= 3.0001)
      {
        s = (std::sqrt(std::pow(j * (1 + e), 2) + std::pow(i * (1 - poissons_ratio * e), 2)) - std::sqrt(i * i + j * j)) / std::sqrt(i * i + j * j);
        lambda += j * (std::exp(-ftemp / 3) * j / ftemp * s);
      }
    }
  lambda /= e;

  return lambda;
}

double
VLEPDMaterial::computeLambda3D(double poissons_ratio)
{
  double e = 0.0001, lambda = 0, ftemp, s;
  for (unsigned int i = -3; i < 4; ++i)
    for (unsigned int j = -3; j < 4; ++j)
      for (unsigned int k = 1; k < 4; ++k)
      {
        ftemp = std::sqrt(i * i + j * j + k * k);
        if (ftemp <= 3.0001)
        {
          s = (std::sqrt(std::pow(k * (1.0 + e), 2) + std::pow(j * (1.0 - poissons_ratio * e), 2) + std::pow(i * (1.0 - poissons_ratio * e), 2)) - std::sqrt(i * i + j * j + k * k)) / std::sqrt(i * i + j * j + k * k);
          lambda += k * (std::exp(-ftemp / 3.0) * k / ftemp * s);
        }
      }
  lambda /= e;

  return lambda;
}
