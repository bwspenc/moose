//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DirichletBCBase.h"

// Forward Declarations
class FluxLimitedDirichletBC;
class Function;

template <>
InputParameters validParams<FluxLimitedDirichletBC>();

/**
 * Defines a boundary condition that forces the value to be a user specified
 * function at the boundary.
 */
class FluxLimitedDirichletBC : public DirichletBCBase
{
public:
  static InputParameters validParams();

  FluxLimitedDirichletBC(const InputParameters & parameters);

  virtual void timestepSetup() override;
  virtual bool shouldApply() override;

protected:
  virtual Real computeQpValue() override;

  /// The function being used for evaluation
  const Function & _func;
  /// The postprocessor that provides the flux
  const PostprocessorValue & _flux_postprocessor;
  /// Limiting value of the flux above which this BC shuts off
  const Real & _flux_limit;
  /// Whether the limiting value of the flux has been exceeded
  bool _exceeded_flux_limit;
};
