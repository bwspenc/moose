//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MaterialRealAux.h"

registerMooseObject("MooseApp", MaterialRealAux);

defineLegacyParams(MaterialRealAux);

InputParameters
MaterialRealAux::validParams()
{
  InputParameters params = NodalPatchRecovery::validParams();
  params.addClassDescription("Outputs element volume-averaged material properties");
  params.addRequiredParam<MaterialPropertyName>("property", "Name of the Real property to output");
  return params;
}

MaterialRealAux::MaterialRealAux(const InputParameters & parameters)
  : NodalPatchRecovery(parameters),
  _prop(getMaterialProperty<Real>("property"))
{
}

Real
MaterialRealAux::computeValue()
{
  return _prop[_qp];
}
