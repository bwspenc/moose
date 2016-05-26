/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PeridynamicMaterial.h"
#include "NonlinearSystem.h"

// libmesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<PeridynamicMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("horizon", "The horizon size");
  return params;
}

PeridynamicMaterial::PeridynamicMaterial(const InputParameters & parameters)
  :Material(parameters),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh)),
  _mesh_spacing(_pdmesh.mesh_spacing()),
  _pddim(_pdmesh.dim()),
  _horizon(getParam<Real>("horizon"))
{
}

void
PeridynamicMaterial::initQpStatefulProperties()
{
}

void
PeridynamicMaterial::computeProperties()
{
  // fetch the volume for the two end nodes
  _nv_i = _pdmesh.volume(_current_elem->get_node(0)->id());
  _nv_j = _pdmesh.volume(_current_elem->get_node(1)->id());
  // original length of a truss element
  _origin_length = _current_elem->volume();

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    computeQpProperties();
  }
}

