//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RankTwoComponentL2Difference.h"
#include "Function.h"

registerMooseObject("TensorMechanicsApp", RankTwoComponentL2Difference);

template <>
InputParameters
validParams<RankTwoComponentL2Difference>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredParam<FunctionName>(
      "function", "The known function for specific component of rank_two_tensor");
  params.addRequiredParam<MaterialPropertyName>(
      "rank_two_tensor", "the name of the rank two tensor we are going to use");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i",
      "index_i >= 0 & index_i <= 2",
      "The index i of ij for the tensor to output (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j",
      "index_j >= 0 & index_j <= 2",
      "The index j of ij for the tensor to output (0, 1, 2)");

  return params;
}

RankTwoComponentL2Difference::RankTwoComponentL2Difference(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _func(getFunction("function")),
    _tensor(getMaterialProperty<RankTwoTensor>("rank_two_tensor")),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j"))
{
}

Real
RankTwoComponentL2Difference::getValue()
{
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}

Real
RankTwoComponentL2Difference::computeQpIntegral()
{
  Real diff = _tensor[_qp](_i, _j) - _func.value(_t, _q_point[_qp]);

  return diff * diff;
}
