//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "USCVoltageLimitedFluxBC.h"
#include "Function.h"

registerMooseObject("MooseApp", USCVoltageLimitedFluxBC);

defineLegacyParams(USCVoltageLimitedFluxBC);

InputParameters
USCVoltageLimitedFluxBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<PostprocessorName>("switch_postprocessor", "PP that switches BCs");
  params.addRequiredParam<FunctionName>("flux", "Imposed flux");
  params.addClassDescription("Imposes the integrated boundary condition "
                             "$\\frac{\\partial u}{\\partial n}=h(t,\\vec{x})$, "
                             "where $h$ is a (possibly) time and space-dependent MOOSE Function.");
  return params;
}

USCVoltageLimitedFluxBC::USCVoltageLimitedFluxBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
  _flux(getFunction("flux")),
  _switch_postprocessor(getPostprocessorValue("switch_postprocessor"))
{
}

Real
USCVoltageLimitedFluxBC::computeQpResidual()
{
  if (_switch_postprocessor == 0)
    return -_test[_i][_qp]*_flux.value(_t, _q_point[_qp]);
  else if (_switch_postprocessor == 2)
    return -_test[_i][_qp]*_flux.value(_t, _q_point[_qp]);
  else
    mooseError("Shouldn't get here");
}

bool
USCVoltageLimitedFluxBC::shouldApply()
{
  return (_switch_postprocessor == 0 || _switch_postprocessor == 2);
}
