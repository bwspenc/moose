/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contens are licensed under LGPL V2.1            */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef PERIDYNAMICMESH_H
#define PERIDYNAMICMESH_H

#include "MooseMesh.h"

class PeridynamicMesh;

template<>
InputParameters validParams<PeridynamicMesh>();

struct pd_node
{
  Point coord;
  double volume;
  unsigned int bond_num;
};

class PeridynamicMesh : public MooseMesh
{
public:
  PeridynamicMesh(const InputParameters & parameters);
  virtual ~PeridynamicMesh();

  virtual void buildMesh(){};

  virtual void find_neighbor();

  virtual std::vector<unsigned int> neighbors(dof_id_type node_id);
  virtual double volume(dof_id_type node_id);
  virtual unsigned int n_neighbors(dof_id_type node_id);
  virtual int dim();
  virtual double mesh_spacing();
  virtual unsigned int n_nodes();
  virtual unsigned int n_bonds();

protected:
  int _pddim;
  double _horizon;

  struct pd_node * _node;

  double _mesh_spacing;
  unsigned int _total_nodes;
  unsigned int _total_bonds;

  std::vector<std::vector<unsigned int> > _neighbors;
};

#endif /* PERIDYNAMICMESH_H */
