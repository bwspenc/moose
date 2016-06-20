/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NODALDILATATIONAUX_H
#define NODALDILATATIONAUX_H

#include "AuxKernel.h"
#include "PDNodalUO.h"

class NodalDilatationAux;

template<>
InputParameters validParams<NodalDilatationAux>();

class NodalDilatationAux : public AuxKernel
{
public:
  NodalDilatationAux(const InputParameters & parameters);
  virtual ~NodalDilatationAux() {}

protected:
  virtual Real computeValue();

  const PDNodalUO * const _pd_nodal;
};

#endif //NODALDILATATIONAUX_H
