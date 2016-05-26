/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FAILUREINDEX_H
#define FAILUREINDEX_H

#include "ElementUserObject.h"

class FailureIndex;

template<>
InputParameters validParams<FailureIndex>();

class FailureIndex :
  public ElementUserObject
{
public:
  FailureIndex(const InputParameters & parameters);

  virtual ~FailureIndex();

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & u );
  virtual void finalize();
  virtual Real computeFailureIndex(unsigned int nodeid) const;

protected:

  AuxiliarySystem & _aux;

  MooseVariable * _intact_bonds_var;
  MooseVariable * _total_bonds_var;

  const MaterialProperty<Real> & _bond_mechanic_strain;

  const VariableValue & _bond_critical_strain;
  const VariableValue & _bond_status;
};

#endif // FAILUREINDEX_H
