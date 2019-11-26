//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DisplacementL2Difference.h"
#include "Function.h"

registerMooseObject("TensorMechanicsApp", DisplacementL2Difference);

template <>
InputParameters
validParams<DisplacementL2Difference>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addParam<FunctionName>("function_0", "The known function for displacement component 0.");
  params.addParam<FunctionName>("function_1", "The known function for displacement component 1.");
  params.addParam<FunctionName>("function_2", "The known function for displacement component 2.");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");

  return params;
}

DisplacementL2Difference::DisplacementL2Difference(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp(_ndisp),
    _has_func_0(isParamValid("function_0")),
    _func_0(_has_func_0 ? &getFunction("function_0") : NULL),
    _has_func_1(isParamValid("function_1")),
    _func_1(_has_func_1 ? &getFunction("function_1") : NULL),
    _has_func_2(isParamValid("function_2")),
    _func_2(_has_func_2 ? &getFunction("function_2") : NULL)
{
  if (_has_func_0 + _has_func_1 + _has_func_2 == 0)
    mooseError("Must provide at least one displacement function for integral error check!");

  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp[i] = &coupledValue("displacements", i);
}

Real
DisplacementL2Difference::getValue()
{
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}

Real
DisplacementL2Difference::computeQpIntegral()
{
  Real diff = 0;
  if (_has_func_0)
  {
    diff += ((*_disp[0])[_qp] - _func_0->value(_t, _q_point[_qp])) *
            ((*_disp[0])[_qp] - _func_0->value(_t, _q_point[_qp]));
  }
  if (_has_func_1)
  {
    diff += ((*_disp[1])[_qp] - _func_1->value(_t, _q_point[_qp])) *
            ((*_disp[1])[_qp] - _func_1->value(_t, _q_point[_qp]));
  }
  if (_has_func_2)
  {
    diff += ((*_disp[2])[_qp] - _func_2->value(_t, _q_point[_qp])) *
            ((*_disp[2])[_qp] - _func_2->value(_t, _q_point[_qp]));
  }

  return diff;
}
