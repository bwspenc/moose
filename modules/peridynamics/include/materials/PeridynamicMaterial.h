/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PERIDYNAMICMATERIAL_H
#define PERIDYNAMICMATERIAL_H

#include "Material.h"
#include "Assembly.h"
#include "PeridynamicMesh.h"

class PeridynamicMaterial;

template<>
InputParameters validParams<PeridynamicMaterial>();

class PeridynamicMaterial : public Material
{
public:
  PeridynamicMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() = 0;
  virtual void computeProperties();
  virtual void computeQpProperties() = 0;

  virtual Real computeBondModulus() = 0;

  virtual Real computeVolSum(dof_id_type node_id);

  AuxiliarySystem & _aux;

  const VariableValue & _bond_status;
  const VariableValue & _bond_contact;
  const VariableValue & _bond_contact_strain;

  PeridynamicMesh & _pdmesh;

  const unsigned int _pddim;

  double _nv_i;
  double _nv_j;
  double _nvsum_i;
  double _nvsum_j;

  double _origin_length;
  double _bond_sign;
};

#endif //PERIDYNAMICMATERIAL_H
