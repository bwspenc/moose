/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                        Peridynamics                          */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef BONDCONTACTAUX_H
#define BONDCONTACTAUX_H

#include "AuxKernel.h"

class BondContactAux;

template<>
InputParameters validParams<BondContactAux>();

class BondContactAux :
public AuxKernel
{
public:
  BondContactAux(const InputParameters & parameters);
  virtual ~BondContactAux() {}

protected:
  virtual Real computeValue();

  const MaterialProperty<Real> & _bond_mechanic_strain;

  const VariableValue & _bond_critical_strain;
//  const VariableValue & _bond_contact_strain;
};

#endif //BONDCONTACTAUX_H
