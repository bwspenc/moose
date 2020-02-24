//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VoltageLimitedFluxBC.h"
#include "Function.h"

registerMooseObject("MooseApp", VoltageLimitedFluxBC);

defineLegacyParams(VoltageLimitedFluxBC);

InputParameters
VoltageLimitedFluxBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<PostprocessorName>("flux_postprocessor", "PP that computes the flux");
  params.addRequiredParam<PostprocessorName>("voltage_postprocessor", "PP that computes the voltage");
  params.addRequiredParam<Real>("flux_limit", "Value of flux at which to turn on this BC");
  params.addRequiredParam<Real>("switching_voltage", "Value of voltage at which the flux switches from a low to a high value");
  params.addRequiredParam<Real>("flux_high", "Imposed flux after the voltage drops below switching_voltage");
  params.addRequiredParam<Real>("flux_low", "Imposed flux before the voltage drops below switching_voltage");
  params.addClassDescription("Imposes the integrated boundary condition "
                             "$\\frac{\\partial u}{\\partial n}=h(t,\\vec{x})$, "
                             "where $h$ is a (possibly) time and space-dependent MOOSE Function.");
  return params;
}

VoltageLimitedFluxBC::VoltageLimitedFluxBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
  _flux_postprocessor(getPostprocessorValueOld("flux_postprocessor")),
  _voltage_postprocessor(getPostprocessorValueOld("voltage_postprocessor")),
  _flux_limit(getParam<Real>("flux_limit")),
  _switching_voltage(getParam<Real>("switching_voltage")),
  _flux_high(getParam<Real>("flux_high")),
  _flux_low(getParam<Real>("flux_low")),
  _exceeded_flux_limit(false),
  _below_switching_voltage(false)
{
}

Real
VoltageLimitedFluxBC::computeQpResidual()
{
  if (_below_switching_voltage)
    return -_test[_i][_qp]*_flux_high;
  else
    return -_test[_i][_qp]*_flux_low;
}

void
VoltageLimitedFluxBC::timestepSetup()
{
  if (_exceeded_flux_limit)
  {
    if (_below_switching_voltage == false)
      if (_voltage_postprocessor < _switching_voltage)
        _below_switching_voltage = true;
  }
  else
  {
    if (_flux_postprocessor > _flux_limit)
      _exceeded_flux_limit = true;
  }
}

bool
VoltageLimitedFluxBC::shouldApply()
{
  return (_exceeded_flux_limit);
}
