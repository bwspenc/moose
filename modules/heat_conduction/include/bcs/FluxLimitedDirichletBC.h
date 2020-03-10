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

  virtual bool shouldApply() override;

protected:
  virtual Real computeQpValue() override;

  /// The high voltage function
  const Function & _func_high;
  /// The low voltage function
  const Function & _func_low;
  /// The postprocessor that switches between stages
  const PostprocessorValue & _switch_postprocessor;
};
