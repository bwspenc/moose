/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NODALVOLSUMAUX_H
#define NODALVOLSUMAUX_H

#include "AuxKernel.h"
#include "PDNodalUO.h"

class NodalVolSumAux;

template<>
InputParameters validParams<NodalVolSumAux>();

class NodalVolSumAux : public AuxKernel
{
public:
  NodalVolSumAux(const InputParameters & parameters);
  virtual ~NodalVolSumAux() {}

protected:
  virtual Real computeValue();

  const PDNodalUO * const _pd_nodal;
};

#endif //NODALVOLSUMAUX_H
