/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PDNODALUO_H
#define PDNODALUO_H

#include "ElementUserObject.h"
#include "PeridynamicMesh.h"

class PDNodalUO;

template<>
InputParameters validParams<PDNodalUO>();

class PDNodalUO : public ElementUserObject
{
public:
  PDNodalUO(const InputParameters & parameters);
  virtual ~PDNodalUO() {}

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & uo);
  virtual void finalize();
  virtual Real computeVolSum(unsigned int node_id) const;
  virtual Real computeDilatation(unsigned int node_id) const;

protected:
  AuxiliarySystem & _aux;

  PeridynamicMesh & _pdmesh;
  const unsigned int _pddim;

  const Real _temp_ref;
  const VariableValue & _temp;

  MooseVariable * _nodal_volsum_var;
  MooseVariable * _nodal_dila_var;

  std::vector<MooseVariable *> _disp_var;

  const VariableValue & _bond_status;
  const VariableValue & _bond_contact;

  double _alpha;
};

#endif // PDNODALUO_H
