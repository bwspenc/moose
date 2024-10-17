//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AdapterFunction.h"
#include "FunctionInterface.h"

registerMooseObject("XFEMApp", AdapterFunction);

InputParameters
AdapterFunction::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Function that calls another function using coorinates and time from other functions");
  params.addParam<FunctionName>("x", "The function that provides the x coordinate");
  params.addParam<FunctionName>("y", "The function that provides the y coordinate");
  params.addParam<FunctionName>("z", "The function that provides the z coordinate");
  params.addParam<FunctionName>("t", "The function that provides the time");
  params.addRequiredParam<FunctionName>("function", "The function is evaluated using overwritten values of x, y, z, and t");

  return params;
}

AdapterFunction::AdapterFunction(const InputParameters & parameters)
  : Function(parameters),
  FunctionInterface(this),
    _x(isParamSetByUser("x") ? & getFunction("x") : nullptr),
    _y(isParamSetByUser("y") ? & getFunction("y") : nullptr),
    _z(isParamSetByUser("z") ? & getFunction("z") : nullptr),
    _t(isParamSetByUser("t") ? & getFunction("t") : nullptr),
    _function(getFunction("function"))
{
}

Real
AdapterFunction::value(Real t, const Point & p) const
{
  auto p_local = p;
  auto t_local = t;
  if (_x)
    p_local(0) = _x->value(t, p);
  if (_y)
    p_local(1) = _y->value(t, p);
  if (_z)
    p_local(2) = _z->value(t, p);
  if (_t)
    t_local = _t->value(t, p);

  return _function.value(t_local, p_local);
}
