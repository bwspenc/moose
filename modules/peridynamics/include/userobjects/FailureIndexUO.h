/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FAILUREINDEXUO_H
#define FAILUREINDEXUO_H

#include "ElementUserObject.h"

class FailureIndexUO;

template<>
InputParameters validParams<FailureIndexUO>();

class FailureIndexUO : public ElementUserObject
{
public:
  FailureIndexUO(const InputParameters & parameters);
  virtual ~FailureIndexUO() {}

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & uo);
  virtual void finalize();
  virtual Real computeFailureIndex(unsigned int node_id) const;

protected:
  AuxiliarySystem & _aux;

  MooseVariable * _intact_bonds_var;

  MooseVariable * _bond_status_var;
};

#endif // FAILUREINDEXUO_H
