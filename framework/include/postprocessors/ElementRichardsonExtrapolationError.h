//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralVariablePostprocessor.h"

class Function;

class ElementRichardsonExtrapolationError : public ElementIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  ElementRichardsonExtrapolationError(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const Function & _function_1;
  const Function & _function_2;
  const Function & _function_n;
};
