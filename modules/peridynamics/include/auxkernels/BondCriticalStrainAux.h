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
  const Real _critical_strain;
  const Real _standard_deviation;

  virtual Real computeValue();
};

#endif //BONDCRITICALSTRAINAUX_H
