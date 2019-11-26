//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef RANKTWOCOMPONENTL2DIFFERENCE_H
#define RANKTWOCOMPONENTL2DIFFERENCE_H

#include "ElementIntegralPostprocessor.h"
#include "RankTwoTensor.h"

class RankTwoComponentL2Difference;
class Function;

template <>
InputParameters validParams<RankTwoComponentL2Difference>();

class RankTwoComponentL2Difference : public ElementIntegralPostprocessor
{
public:
  RankTwoComponentL2Difference(const InputParameters & parameters);

  virtual Real getValue() override;

protected:
  virtual Real computeQpIntegral() override;

  /// The function being used for evaluation
  const Function & _func;

  const MaterialProperty<RankTwoTensor> & _tensor;

  const unsigned int _i;
  const unsigned int _j;
};

#endif // RANKTWOCOMPONENTL2DIFFERENCE
