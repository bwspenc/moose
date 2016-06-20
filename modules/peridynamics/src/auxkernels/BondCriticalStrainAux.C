/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "BondCriticalStrainAux.h"
#include "PeridynamicMesh.h"

template<>
InputParameters validParams<BondCriticalStrainAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<std::string>("plane_strain", "Plane strain problem");
  params.addRequiredParam<Real>("critical_energy_release_rate", "critical strain energy release rate");
  params.addRequiredParam<Real>("Youngs_modulus", "Youngs modulus");
  params.addRequiredParam<Real>("Poissons_ratio", "Poisson's ratio");
  return params;
}

BondCriticalStrainAux::BondCriticalStrainAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _Gc(getParam<Real>("critical_energy_release_rate")),
  _E(getParam<Real>("Youngs_modulus")),
  _mu(getParam<Real>("Poissons_ratio")),
  _pddim(dynamic_cast<PeridynamicMesh &>(_mesh).dim())
{
  if (isParamValid("plane_strain"))
    _kappa = _E / 2.0 / (1.0 + _mu) / (1.0 - 2.0 * _mu);
  else
    _kappa = _E / _pddim / (1.0 - (_pddim - 1.0) * _mu);

  setRandomResetFrequency(EXEC_INITIAL);
}

Real
BondCriticalStrainAux::computeValue()
{
  // Generate randomized critical stretch by Box-Muller method: randomized_critical_strain = small_random_number + _input_critical_strain
  double val = std::sqrt(2.0 * _Gc / _kappa / std::pow(_pddim, 2) / _current_elem_volume);
//  double val = std::sqrt(8.0 * 3.1415926 * _Gc / 27.0 / _E / 0.17697842);
//  double val = 0.000446185; // 3x
//  double val = 0.00034595; // 5x
  return  (std::sqrt(- 2.0 * std::log(getRandomReal())) * std::cos(2.0 * 3.14159265358 * getRandomReal()) * 0.05 + 1.0) * val;
}
