/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "BondCriticalStrainPDAux.h"

template<>
InputParameters validParams<BondCriticalStrainPDAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("critical_strain",1.0, "critical_strain");
  params.addParam<Real>("standard_deviation",0.0, "standard_deviation");
  return params;
}

BondCriticalStrainPDAux::BondCriticalStrainPDAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _critical_strain(isParamValid("critical_strain") ? getParam<Real>("critical_strain") : 1.0),
  _standard_deviation(isParamValid("standard_deviation") ? getParam<Real>("standard_deviation") : 0.0)
{
  setRandomResetFrequency(EXEC_INITIAL);
}

Real
BondCriticalStrainPDAux::computeValue()
{
// Generate randomized critical stretch by Box-Muller method
  return  std::sqrt(-2.0 * std::log(getRandomReal())) * std::cos(2.0 * 3.14159265358 * getRandomReal()) * _standard_deviation + _critical_strain;

}
