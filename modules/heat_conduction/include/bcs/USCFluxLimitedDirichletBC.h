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
class USCFluxLimitedDirichletBC;
class Function;

template <>
InputParameters validParams<USCFluxLimitedDirichletBC>();

/**
 * Defines a boundary condition that forces the value to be a user specified
 * function at the boundary.
 */
class USCFluxLimitedDirichletBC : public DirichletBCBase
{
public:
  static InputParameters validParams();

  USCFluxLimitedDirichletBC(const InputParameters & parameters);

  virtual bool shouldApply() override;

protected:
  virtual Real computeQpValue() override;

  /// The voltage function
  const Function & _func;
  /// The postprocessor that switches between stages
  const PostprocessorValue & _switch_postprocessor;
};
