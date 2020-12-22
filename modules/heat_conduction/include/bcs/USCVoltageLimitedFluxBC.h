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
class USCVoltageLimitedFluxBC;
class Function;

template <>
InputParameters validParams<USCVoltageLimitedFluxBC>();

/**
 * BC that prescribes a flux
 */
class USCVoltageLimitedFluxBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  USCVoltageLimitedFluxBC(const InputParameters & parameters);

  virtual bool shouldApply() override;

protected:
  virtual Real computeQpResidual() override;

  /// Imposed flux
  const Function & _flux;
  /// The postprocessor that switches between stages
  const PostprocessorValue & _switch_postprocessor;
};
