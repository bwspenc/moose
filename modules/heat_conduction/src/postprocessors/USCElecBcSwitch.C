//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "USCElecBcSwitch.h"
#include "Function.h"

registerMooseObject("HeatConductionApp", USCElecBcSwitch);

defineLegacyParams(USCElecBcSwitch);

InputParameters
USCElecBcSwitch::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<PostprocessorName>("flux_postprocessor", "PP that provides the flux");
  params.addRequiredParam<PostprocessorName>("voltage_postprocessor", "PP that provides the voltage");
  params.addRequiredParam<FunctionName>("flux_switch", "Value of the flux at which the voltage BC is switched to a flux BC");
  params.addRequiredParam<Real>("voltage_switch", "Value of the voltage at which the flux BC is switched to a voltage Dirichlet BC");
  params.addClassDescription("This postprocessor computes the stage of the ANL experiment that used resistive heating on a stack of 7 pellets");
  return params;
}

USCElecBcSwitch::USCElecBcSwitch(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
  _flux_postprocessor(getPostprocessorValue("flux_postprocessor")),
  _voltage_postprocessor(getPostprocessorValue("voltage_postprocessor")),
  _voltage_switch(getParam<Real>("voltage_switch")),
  _flux_switch(getFunction("flux_switch")),
  _stage(0)
{
}

void
USCElecBcSwitch::timestepSetup()
{
  if (_stage == 0)
  {
    if (_voltage_postprocessor >= _voltage_switch)
      _stage = 1;
  }
  else if (_stage == 1)
  {
    if (_flux_postprocessor >= _flux_switch.value(_t, Point()))
      _stage = 2;
  }
}

PostprocessorValue
USCElecBcSwitch::getValue()
{
  return _stage;
}

