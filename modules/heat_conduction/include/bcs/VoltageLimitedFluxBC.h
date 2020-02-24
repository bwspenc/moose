//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBC.h"

// Forward Declarations
class VoltageLimitedFluxBC;
class Function;

template <>
InputParameters validParams<VoltageLimitedFluxBC>();

/**
 * BC that prescribes a flux
 */
class VoltageLimitedFluxBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  VoltageLimitedFluxBC(const InputParameters & parameters);

  virtual void timestepSetup() override;
  virtual bool shouldApply() override;

protected:
  virtual Real computeQpResidual() override;

  /// The postprocessor that provides the flux
  const PostprocessorValue & _flux_postprocessor;
  /// The postprocessor that provides the voltage
  const PostprocessorValue & _voltage_postprocessor;
  /// Limiting value of the flux above which this BC turns on
  const Real & _flux_limit;
  /// Voltage at which the flux switches to a higher value
  const Real & _switching_voltage;
  /// Imposed flux after the voltage drops below switching_voltage
  const Real & _flux_high;
  /// Imposed flux before the voltage drops below switching_voltage
  const Real & _flux_low;
  /// Whether the limiting value of the flux has been exceeded
  bool _exceeded_flux_limit;
  /// Whether we're below the voltage at which the flux changes
  bool _below_switching_voltage;
};
