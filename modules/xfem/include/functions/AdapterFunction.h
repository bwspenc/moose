//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "Function.h"
#include "FunctionInterface.h"

/**
 * Adapter function that...
 */
class AdapterFunction : public Function,
  protected FunctionInterface
{
public:
  static InputParameters validParams();

  AdapterFunction(const InputParameters & parameters);

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

protected:
  const Function * _x;
  const Function * _y;
  const Function * _z;
  const Function * _t;
  const Function & _function;
};
