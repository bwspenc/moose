/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef BONDCONTACTSTRAINAUX_H
#define BONDCONTACTSTRAINAUX_H

#include "AuxKernel.h"

class BondContactStrainAux;

template<>
InputParameters validParams<BondContactStrainAux>();

class BondContactStrainAux : public AuxKernel
{
public:
  BondContactStrainAux(const InputParameters & parameters);
  virtual ~BondContactStrainAux() {}

protected:
  const VariableValue & _bond_critical_strain;

  virtual Real computeValue();
};

#endif //BONDCONTACTSTRAINAUX_H
