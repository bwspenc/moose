/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FAILUREINDEXAUX_H
#define FAILUREINDEXAUX_H

#include "AuxKernel.h"
#include "FailureIndexUO.h"

class FailureIndexAux;

template<>
InputParameters validParams<FailureIndexAux>();

class FailureIndexAux : public AuxKernel
{
public:
  FailureIndexAux(const InputParameters & parameters);
  virtual ~FailureIndexAux() {}

protected:
  const FailureIndexUO * const _failure_index;
  virtual Real computeValue();
};

#endif //FAILUREINDEXAUX_H
