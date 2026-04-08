//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OldValueCoupledNodalAux.h"

registerMooseObject("SolidMechanicsApp", OldValueCoupledNodalAux);

InputParameters
OldValueCoupledNodalAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Skeleton nodal auxiliary kernel that computes an auxiliary value "
                             "from its old value and a coupled variable's current and old values.");
  params.addRequiredCoupledVar("coupled_var", "The variable whose current and old values are used.");
  return params;
}

OldValueCoupledNodalAux::OldValueCoupledNodalAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _u_old(uOld()),
    _coupled_value(coupledValue("coupled_var")),
    _coupled_value_old(coupledValueOld("coupled_var"))
{
  if (!isNodal())
    paramError("variable", "OldValueCoupledNodalAux must be used with a nodal auxiliary variable.");
}

Real
OldValueCoupledNodalAux::computeValue()
{
  // Placeholder update rule for customization.
  return _u_old[_qp] + _coupled_value[_qp] - _coupled_value_old[_qp];
}
