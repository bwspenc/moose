/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "NodalDilatationAux.h"

template<>
InputParameters validParams<NodalDilatationAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<UserObjectName>("pd_nodal", "The name of the PDNodal user object");
  return params;
}

NodalDilatationAux::NodalDilatationAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _pd_nodal(&getUserObject<PDNodalUO>("pd_nodal"))
{
}

Real
NodalDilatationAux::computeValue()
{
  if (!isNodal())
    mooseError("must run on a nodal variable");

  return _pd_nodal->computeDilatation(_current_node->id());
}
