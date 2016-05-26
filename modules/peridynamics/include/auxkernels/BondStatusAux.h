/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                        Peridynamics                          */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef BONDSTATUSAUX_H
#define BONDSTATUSAUX_H

#include "AuxKernel.h"

class BondStatusAux;

template<>
InputParameters validParams<BondStatusAux>();

class BondStatusAux : public AuxKernel
{
public:
  BondStatusAux(const InputParameters & parameters);
  virtual ~BondStatusAux() {}

protected:
  virtual Real computeValue();

  const MaterialProperty<Real> & _bond_mechanic_strain;

  const VariableValue & _bond_critical_strain;
  const VariableValue & _bond_status;
};

#endif //BONDSTATUSAUX_H
