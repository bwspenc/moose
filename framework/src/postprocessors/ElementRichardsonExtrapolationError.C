//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementRichardsonExtrapolationError.h"
#include "Function.h"

registerMooseObject("MooseApp", ElementRichardsonExtrapolationError);

InputParameters
ElementRichardsonExtrapolationError::validParams()
{
  InputParameters params = ElementIntegralVariablePostprocessor::validParams();
  params.addRequiredParam<FunctionName>("function_1", "The solution at the most refined mesh");
  params.addRequiredParam<FunctionName>("function_1", "The solution at the second-most refined mesh");
  params.addRequiredParam<FunctionName>("function_n", "The solution at the n-th refined mesh");
  params.addClassDescription(
      "Computes L2 error between a field variable and an analytical function");
  return params;
}

ElementRichardsonExtrapolationError::ElementRichardsonExtrapolationError(const InputParameters & parameters)
  : ElementIntegralVariablePostprocessor(parameters),
  _function_1(getFunction("function_1")),
  _function_2(getFunction("function_2")),
  _function_n(getFunction("function_n"))
{
}

Real
ElementRichardsonExtrapolationError::computeQpIntegral()
{
  //TODO: This is where you plug in the Richardson extraploation function
  Real err = _function_n.value(_t, _q_point[_qp]) - _function_2.value(_t, _q_point[_qp]) - _function_1.value(_t, _q_point[_qp]);
  return err;
}
