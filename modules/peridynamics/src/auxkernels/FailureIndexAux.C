/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "FailureIndexAux.h"

template<>
InputParameters validParams<FailureIndexAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<UserObjectName>("failure_index_uo", "The name of the FailureIndexUO UserObject");
  return params;
}

FailureIndexAux::FailureIndexAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _failure_index_uo(&getUserObject<FailureIndexUO>("failure_index_uo"))
{
}

Real
FailureIndexAux::computeValue()
{
  if (!isNodal())
    mooseError("must run on a nodal variable");

  return _failure_index_uo->computeFailureIndex(_current_node->id());
}
