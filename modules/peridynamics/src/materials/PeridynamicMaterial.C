/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PeridynamicMaterial.h"

// libmesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<PeridynamicMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addCoupledVar("bond_status", "Auxiliary variable for bond failure status");
  params.addCoupledVar("bond_contact", "Auxiliary variable for bond contact status");
  params.addCoupledVar("bond_contact_strain", "Auxiliary variable for bond contact strain");
  return params;
}

PeridynamicMaterial::PeridynamicMaterial(const InputParameters & parameters) :
  Material(parameters),
  _nsys(_fe_problem.getNonlinearSystem()),
  _bond_status(coupledValue("bond_status")),
  _bond_contact(coupledValue("bond_contact")),
  _bond_contact_strain(coupledValue("bond_contact_strain")),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh)),
  _pddim(_pdmesh.dim())
{
}

void
PeridynamicMaterial::computeProperties()
{
  // the horizon size of two end nodes
  _horizon_i = _pdmesh.horizon(_current_elem->get_node(0)->id());
  _horizon_j = _pdmesh.horizon(_current_elem->get_node(1)->id());

  // the volume of two end nodes
  _nv_i = _pdmesh.volume(_current_elem->get_node(0)->id());
  _nv_j = _pdmesh.volume(_current_elem->get_node(1)->id());

  // the volume sum from all neighbors of two end nodes
  _nvsum_i = _pdmesh.volumesum(_current_elem->get_node(0)->id());
  _nvsum_j = _pdmesh.volumesum(_current_elem->get_node(1)->id());

  // the temperature of two end nodes
  computeNodalTemp();

  // the nodal strain tensor
  computeElasticStrainTensor();

  // the nodal stress tensor
  computeStressTensor();

  // original length of a truss element
  _origin_length = (*_current_elem->get_node(1) - *_current_elem->get_node(0)).size();

  // current length of a truss element
  _current_length = computeBondCurrentLength();

  // check whether the bond is broken and/or in contact
  if (std::abs(_bond_status[0] - 1.0) < 0.01)
    _bond_sign = 1.0;
  else
    if (std::abs(_bond_contact[0] - 1.0) < 0.01)
      _bond_sign = 1.0;
    else
      _bond_sign = 0.0;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    computeQpStrain();
    computeQpForce();
  }
}
