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

class USCElecBcSwitch;

template <>
InputParameters validParams<USCElecBcSwitch>();

/**
 * This postprocessor computes the stage of the ANL experiment
 * that used resistive heating on a stack of 7 pellets
 */
class USCElecBcSwitch : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  USCElecBcSwitch(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual void timestepSetup() override;
  virtual PostprocessorValue getValue() override;

protected:
  /// The postprocessor that provides the flux
  const PostprocessorValue & _flux_postprocessor;
  /// The postprocessor that provides the voltage
  const PostprocessorValue & _voltage_postprocessor;
  /// Value of the voltage at which the flux BC is switched to a Dirichlet BC
  const Real _voltage_switch;
  /// Value of the flux at which the voltage BC is switched back to a flux BC
  const Function & _flux_switch;

  /* Stage of the experiment:
   * 0: Initial high voltage, dirichlet BC
   * 1: Initial low current, flux BC
   * 2: Ligh voltage, dirichlet BC
   * 3: High current, flux BC
   */
  unsigned int _stage;
};
