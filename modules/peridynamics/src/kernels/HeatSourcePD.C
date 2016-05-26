/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "HeatSourcePD.h"
#include "Function.h"
#include "Assembly.h"

template<>
InputParameters validParams<HeatSourcePD>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<Real>("power_density", 0, "Volumetric heat source density");
  params.addParam<FunctionName>("power_density_function", "", "Function describing the volumetric heat source density");
  return params;
}

HeatSourcePD::HeatSourcePD(const InputParameters & parameters)
  :Kernel(parameters),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh)),
  _power_density(getParam<Real>("power_density")),
  _power_density_function(getParam<FunctionName>("power_density_function") != "" ? &getFunction("power_density_function") : NULL)
{
}

void
HeatSourcePD::computeResidual()
{
//get the volume and total_bonds for each node
  double nv0 = _pdmesh.volume(_current_elem->get_node(0)->id());
  double nv1 = _pdmesh.volume(_current_elem->get_node(1)->id());

  unsigned int tb0 = _pdmesh.n_neighbors(_current_elem->get_node(0)->id());
  unsigned int tb1 = _pdmesh.n_neighbors(_current_elem->get_node(1)->id());

  if(_power_density_function)
  {
    Point p;
    _power_density = _power_density_function->value(_t, p);
  }

  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  mooseAssert(re.size() == 2, "Truss elements must have two nodes");
  _local_re.resize(re.size());
  _local_re.zero();

  _local_re(0) = -_power_density * nv0 / tb0;
  _local_re(1) = -_power_density * nv1 / tb1;

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _save_in.size(); i++)
      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
  }
}
