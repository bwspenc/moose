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
  params.addRequiredParam<PostprocessorName>("switch_postprocessor", "PP that switches BCs");
  params.addRequiredParam<Real>("flux_high_voltage", "Imposed flux in the high voltage regime");
  params.addRequiredParam<Real>("flux_low_voltage", "Imposed flux in the low voltage regime");
  params.addClassDescription("Imposes the integrated boundary condition "
                             "$\\frac{\\partial u}{\\partial n}=h(t,\\vec{x})$, "
                             "where $h$ is a (possibly) time and space-dependent MOOSE Function.");
  return params;
}

VoltageLimitedFluxBC::VoltageLimitedFluxBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
  _flux_high_voltage(getParam<Real>("flux_high_voltage")),
  _flux_low_voltage(getParam<Real>("flux_low_voltage")),
  _switch_postprocessor(getPostprocessorValue("switch_postprocessor"))
{
}

Real
VoltageLimitedFluxBC::computeQpResidual()
{
  if (_switch_postprocessor == 1)
    return -_test[_i][_qp]*_flux_high_voltage;
  else if (_switch_postprocessor == 3)
    return -_test[_i][_qp]*_flux_low_voltage;
  else
    mooseError("Shouldn't get here");
}

bool
VoltageLimitedFluxBC::shouldApply()
{
  return (_switch_postprocessor == 1 || _switch_postprocessor == 3);
}
