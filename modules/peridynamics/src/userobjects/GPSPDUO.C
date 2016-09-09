/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                         Peridynamics                         */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "GPSPDUO.h"

template<>
InputParameters validParams<GPSPDUO>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("Generalized Plane Strain UserObject to provide Residual and diagonal Jacobian entry");
  return params;
}

GPSPDUO::GPSPDUO(const InputParameters & parameters) :
  ElementUserObject(parameters),
  _Cijkl(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
  _stress(getMaterialProperty<RankTwoTensor>("stress")),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh))
{
}

void
GPSPDUO::initialize()
{
  _residual = 0;
  _jacobian = 0;
}

void
GPSPDUO::execute()
{
  // nodal area for node i and j
  double nv_i = _pdmesh.volume(_current_elem->get_node(0)->id());
  double nv_j = _pdmesh.volume(_current_elem->get_node(1)->id());

  // number of neighbors for node i and j
  unsigned int nn_i = _pdmesh.neighbors(_current_elem->get_node(0)->id()).size();
  unsigned int nn_j = _pdmesh.neighbors(_current_elem->get_node(1)->id()).size();

  // residual
  _residual += _stress[0](2, 2) * nv_i / nn_i + _stress[1](2, 2) * nv_j / nn_j;

  // diagonal jacobian
  _jacobian += _Cijkl[0](2, 2, 2, 2) * nv_i / nn_i + _Cijkl[0](2, 2, 2, 2) * nv_j / nn_j;
}

Real
GPSPDUO::returnResidual() const
{
  return _residual;
}

Real
GPSPDUO::returnJacobian() const
{
  return _jacobian;
}
