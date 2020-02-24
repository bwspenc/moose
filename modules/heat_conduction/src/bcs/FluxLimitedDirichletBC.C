//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FluxLimitedDirichletBC.h"
#include "Function.h"

registerMooseObject("MooseApp", FluxLimitedDirichletBC);

defineLegacyParams(FluxLimitedDirichletBC);

InputParameters
FluxLimitedDirichletBC::validParams()
{
  InputParameters params = DirichletBCBase::validParams();
  params.addRequiredParam<FunctionName>("function", "The forcing function.");
  params.addRequiredParam<PostprocessorName>("flux_postprocessor", "PP that computes the flux");
  params.addRequiredParam<Real>("flux_limit", "Value of flux at which to shut off this BC");
  params.addClassDescription(
      "Imposes the essential boundary condition $u=g(t,\\vec{x})$, where $g$ "
      "is a (possibly) time and space-dependent MOOSE Function.");
  return params;
}

FluxLimitedDirichletBC::FluxLimitedDirichletBC(const InputParameters & parameters)
  : DirichletBCBase(parameters),
  _func(getFunction("function")),
  _flux_postprocessor(getPostprocessorValueOld("flux_postprocessor")),
  _flux_limit(getParam<Real>("flux_limit")),
  _exceeded_flux_limit(false)
{
}

Real
FluxLimitedDirichletBC::computeQpValue()
{
  return _func.value(_t, *_current_node);
}

void
FluxLimitedDirichletBC::timestepSetup()
{
  if (_exceeded_flux_limit == false)
  {
    if (_flux_postprocessor > _flux_limit)
      _exceeded_flux_limit = true;
  }
}

bool
FluxLimitedDirichletBC::shouldApply()
{
  return (_exceeded_flux_limit == false);
}
