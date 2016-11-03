/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PERIDYNAMICSACTION_H
#define PERIDYNAMICSACTION_H

#include "Action.h"

class PeridynamicsAction;

template<>
InputParameters validParams<PeridynamicsAction>();

class PeridynamicsAction : public Action
{
public:
  PeridynamicsAction(const InputParameters & params);

  virtual void act();
  virtual void addkernel(const std::string & name, InputParameters & params);
};

#endif //PERIDYNAMICSACTION_H
