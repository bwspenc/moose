//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

class FluxDirichletElecBCSwitch;

template <>
InputParameters validParams<FluxDirichletElecBCSwitch>();

/**
 * This postprocessor computes the stage of the ANL experiment
 * that used resistive heating on a stack of 7 pellets
 */
class FluxDirichletElecBCSwitch : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  FluxDirichletElecBCSwitch(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual void timestepSetup() override;
  virtual PostprocessorValue getValue() override;

protected:
  /// The postprocessor that provides the flux
  const PostprocessorValue & _flux_postprocessor;
  /// The postprocessor that provides the voltage
  const PostprocessorValue & _voltage_postprocessor;
  /// Value of the flux at which the initial high voltage BC is switched to a flux BC
  const Real _flux_switch_high_voltage;
  /// Value of the voltage at which the flux BC is switched to a low voltage Dirichlet BC
  const Real _voltage_switch_low_voltage;
  /// Value of the flux at which low voltage BC is switched to a high flux BC
  const Real _flux_switch_low_voltage;

  /* Stage of the experiment:
   * 0: Initial high voltage, dirichlet BC
   * 1: Initial low current, flux BC
   * 2: Ligh voltage, dirichlet BC
   * 3: High current, flux BC
   */
  unsigned int _stage;
};
