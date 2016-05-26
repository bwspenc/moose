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

  PeridynamicMesh & _pdmesh;

  const double _mesh_spacing;
  const unsigned int _pddim;

  const double _horizon;

  double _nv_i;
  double _nv_j;

  double _origin_length;
};

#endif //PERIDYNAMICMATERIAL_H
