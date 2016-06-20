/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DILATATIONUO_H
#define DILATATIONUO_H

#include "ElementUserObject.h"
#include "PeridynamicMesh.h"

class DilatationUO;

template<>
InputParameters validParams<DilatationUO>();

class DilatationUO : public ElementUserObject
{
public:
  DilatationUO(const InputParameters & parameters);
  virtual ~DilatationUO() {}

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & uo);
  virtual void finalize();

protected:
  AuxiliarySystem & _aux;

  PeridynamicMesh & _pdmesh;
  const unsigned int _pddim;

  const Real _temp_ref;
  const VariableValue & _temp;

  MooseVariable * _nodal_dila_var;

  std::vector<MooseVariable *> _disp_var;

  const VariableValue & _bond_status;
  const VariableValue & _bond_contact;

  double _alpha;
};

#endif // DILATATIONUO_H
