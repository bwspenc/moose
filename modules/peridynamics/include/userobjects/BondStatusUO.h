/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef BONDSTATUSUO_H
#define BONDSTATUSUO_H

#include "ElementUserObject.h"

class BondStatusUO;

template<>
InputParameters validParams<BondStatusUO>();

class BondStatusUO : public ElementUserObject
{
public:
  BondStatusUO(const InputParameters & parameters);
  virtual ~BondStatusUO() {}

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & uo);
  virtual void finalize();

protected:
  AuxiliarySystem & _aux;

  const MaterialProperty<Real> & _bond_elastic_strain;
  const MaterialProperty<Real> & _bond_critical_strain;

  MooseVariable * _bond_status_var;
};

#endif // BONDSTATUSUO_H
