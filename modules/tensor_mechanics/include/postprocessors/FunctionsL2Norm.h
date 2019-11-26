//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef FUNCTIONSL2NORM_H
#define FUNCTIONSL2NORM_H

#include "ElementIntegralPostprocessor.h"

class FunctionsL2Norm;
class Function;

template <>
InputParameters validParams<FunctionsL2Norm>();

class FunctionsL2Norm : public ElementIntegralPostprocessor
{
public:
  FunctionsL2Norm(const InputParameters & parameters);

  virtual Real getValue() override;

protected:
  virtual Real computeQpIntegral() override;

  const bool _has_func_0;
  const Function * _func_0;

  const bool _has_func_1;
  const Function * _func_1;

  const bool _has_func_2;
  const Function * _func_2;
};

#endif // FUNCTIONSL2NORM
