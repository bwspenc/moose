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
#include "PeridynamicMesh.h"
#include "NonlinearSystem.h"
#include "RankTwoTensor.h"

class PeridynamicMaterial;

template<>
InputParameters validParams<PeridynamicMaterial>();

class PeridynamicMaterial : public Material
{
public:
  PeridynamicMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties(){};
  virtual void computeProperties();
  virtual void computeQpStrain() = 0;
  virtual void computeQpForce() = 0;
  virtual void computeNodalTemp() = 0;
  virtual void computeElasticStrainTensor(){};
  virtual void computeStressTensor(){};

  virtual Real computeBondModulus(){return 0;}
  virtual Real computeBondCurrentLength(){return 0;}

  NonlinearSystem & _nsys;

  const VariableValue & _bond_status;
  const VariableValue & _bond_contact;
  const VariableValue & _bond_contact_strain;

  PeridynamicMesh & _pdmesh;

  const unsigned int _pddim;

  double _horizon_i;
  double _horizon_j;
  double _nv_i;
  double _nv_j;
  double _nvsum_i;
  double _nvsum_j;
  double _temp_i;
  double _temp_j;

  double _origin_length;
  double _current_length;
  double _bond_sign;
};

#endif //PERIDYNAMICMATERIAL_H
