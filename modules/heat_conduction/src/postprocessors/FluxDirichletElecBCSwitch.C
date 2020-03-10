//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FluxDirichletElecBCSwitch.h"

registerMooseObject("HeatConductionApp", FluxDirichletElecBCSwitch);

defineLegacyParams(FluxDirichletElecBCSwitch);

InputParameters
FluxDirichletElecBCSwitch::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<PostprocessorName>("flux_postprocessor", "PP that provides the flux");
  params.addRequiredParam<PostprocessorName>("voltage_postprocessor", "PP that provides the voltage");
  params.addRequiredParam<Real>("flux_switch_high_voltage", "Value of the flux at which the initial high voltage BC is switched to a flux BC");
  params.addRequiredParam<Real>("voltage_switch_low_voltage", "Value of the voltage at which the flux BC is switched to a low voltage Dirichlet BC");
  params.addRequiredParam<Real>("flux_switch_low_voltage", "Value of the flux at which the low voltage BC is switched to a high flux BC");
  params.addClassDescription("This postprocessor computes the stage of the ANL experiment that used resistive heating on a stack of 7 pellets");
  return params;
}

FluxDirichletElecBCSwitch::FluxDirichletElecBCSwitch(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
  _flux_postprocessor(getPostprocessorValue("flux_postprocessor")),
  _voltage_postprocessor(getPostprocessorValue("voltage_postprocessor")),
  _flux_switch_high_voltage(getParam<Real>("flux_switch_high_voltage")),
  _voltage_switch_low_voltage(getParam<Real>("voltage_switch_low_voltage")),
  _flux_switch_low_voltage(getParam<Real>("flux_switch_low_voltage")),
  _stage(0)
{
}

void
FluxDirichletElecBCSwitch::timestepSetup()
{
  if (_stage == 0)
  {
    if (_flux_postprocessor >= _flux_switch_high_voltage)
      _stage = 1;
  }
  else if (_stage == 1)
  {
    if (_voltage_postprocessor <= _voltage_switch_low_voltage)
      _stage = 2;
  }
  else if (_stage == 2)
  {
    if (_flux_postprocessor >= _flux_switch_low_voltage)
      _stage = 3;
  }
}

PostprocessorValue
FluxDirichletElecBCSwitch::getValue()
{
  return _stage;
}

