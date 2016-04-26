/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#ifndef BONDCRITICALSTRAINPDAUX_H
#define BONDCRITICALSTRAINPDAUX_H

#include "AuxKernel.h"

//Forward Declarations
class BondCriticalStrainPDAux;

template<>
InputParameters validParams<BondCriticalStrainPDAux>();

class BondCriticalStrainPDAux : public AuxKernel
{
public:

  BondCriticalStrainPDAux(const InputParameters & parameters);

  virtual ~BondCriticalStrainPDAux() {}

protected:

  const Real _critical_strain;
  const Real _standard_deviation;

  virtual Real computeValue();

};

#endif //BONDCRITICALSTRAINPDAUX_H
