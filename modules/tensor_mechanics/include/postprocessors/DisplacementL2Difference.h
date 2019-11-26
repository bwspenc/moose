//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DISPLACEMENTL2DIFFERENCE_H
#define DISPLACEMENTL2DIFFERENCE_H

#include "ElementIntegralPostprocessor.h"

class DisplacementL2Difference;
class Function;

template <>
InputParameters validParams<DisplacementL2Difference>();

class DisplacementL2Difference : public ElementIntegralPostprocessor
{
public:
  DisplacementL2Difference(const InputParameters & parameters);

  virtual Real getValue() override;

protected:
  virtual Real computeQpIntegral() override;

  const unsigned int _ndisp;
  std::vector<const VariableValue *> _disp;

  const bool _has_func_0;
  const Function * _func_0;

  const bool _has_func_1;
  const Function * _func_1;

  const bool _has_func_2;
  const Function * _func_2;
};

#endif // DISPLACEMENTL2DIFFERENCE
