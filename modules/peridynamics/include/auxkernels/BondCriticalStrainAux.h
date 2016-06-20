/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef BONDCRITICALSTRAINAUX_H
#define BONDCRITICALSTRAINAUX_H

#include "AuxKernel.h"

class BondCriticalStrainAux;

template<>
InputParameters validParams<BondCriticalStrainAux>();

class BondCriticalStrainAux : public AuxKernel
{
public:
  BondCriticalStrainAux(const InputParameters & parameters);
  virtual ~BondCriticalStrainAux() {}

protected:
  virtual Real computeValue();

  const double _Gc;
  const double _E;
  const double _mu;

  const int _pddim;

  double _kappa;
};

#endif //BONDCRITICALSTRAINAUX_H
