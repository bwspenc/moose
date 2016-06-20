/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "NodalVolSumAux.h"

template<>
InputParameters validParams<NodalVolSumAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<UserObjectName>("pd_nodal", "The name of the PDNodal user object");
  return params;
}

NodalVolSumAux::NodalVolSumAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _pd_nodal(&getUserObject<PDNodalUO>("pd_nodal"))
{
}

Real
NodalVolSumAux::computeValue()
{
  if (!isNodal())
    mooseError("must run on a nodal variable");

  return _pd_nodal->computeVolSum(_current_node->id());
}
