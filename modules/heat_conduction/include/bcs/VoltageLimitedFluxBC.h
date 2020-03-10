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

  virtual bool shouldApply() override;

protected:
  virtual Real computeQpResidual() override;

  /// Imposed flux in the high voltage regime
  const Real & _flux_high_voltage;
  /// Imposed flux in the low voltage regime
  const Real & _flux_low_voltage;
  /// The postprocessor that switches between stages
  const PostprocessorValue & _switch_postprocessor;
};
