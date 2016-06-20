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
  _aux(_fe_problem.getAuxiliarySystem()),
  _bond_status(coupledValue("bond_status")),
  _bond_contact(coupledValue("bond_contact")),
  _bond_contact_strain(coupledValue("bond_contact_strain")),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh)),
  _pddim(_pdmesh.dim())
{
}

void
PeridynamicMaterial::initQpStatefulProperties()
{
}

void
PeridynamicMaterial::computeProperties()
{
  // the volume for the two end nodes
  _nv_i = _pdmesh.volume(_current_elem->get_node(0)->id());
  _nv_j = _pdmesh.volume(_current_elem->get_node(1)->id());
  // the volume sum of all neighbors for the two end nodes
  _nvsum_i = computeVolSum(_current_elem->get_node(0)->id());
  _nvsum_j = computeVolSum(_current_elem->get_node(1)->id());

  // original length of a truss element
  RealGradient dxyz;
  for (unsigned int i = 0; i < _pddim; ++i)
    dxyz(i) = (*_current_elem->get_node(1))(i) - (*_current_elem->get_node(0))(i);

  _origin_length = dxyz.norm();

  // check whether the bond is broken and/or in contact
  if (std::abs(_bond_status[0] - 1.0) < 0.01)
    _bond_sign = 1.0;
  else
    _bond_sign = 0.0;
//    if (std::abs(_bond_contact[0] - 1.0) < 0.01)
//      _bond_sign = 1.0;
//    else
//      _bond_sign = 0.0;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    computeQpProperties();
}

Real
PeridynamicMaterial::computeVolSum(dof_id_type node_id)
{
  double val = 0;
  unsigned int nneighbors = _pdmesh.n_neighbors(node_id);
  std::vector<dof_id_type> neighbors = _pdmesh.neighbors(node_id);

  for (unsigned int k = 0; k < nneighbors; ++k)
    val += _pdmesh.volume(neighbors[k]);

  return val;
}
