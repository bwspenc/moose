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
  /*
   * function for neighbor search with given horizon
   */
  virtual void find_neighbor();

  /*
   * return neighbor nodes indices for node node_id
   */
  virtual std::vector<dof_id_type> neighbors(dof_id_type node_id);

  /*
   * return coordinates for node node_id
   */
  virtual Point coord(dof_id_type node_id);

  /*
   * return nodal volume for node node_id
   */
  virtual double volume(dof_id_type node_id);

  /*
   * return neighbor number for node node_id
   */
  virtual unsigned int n_neighbors(dof_id_type node_id);

  /*
   * return PD dimension
   */
  virtual int dim();

  /*
   * return mesh_spacing
   */
  virtual double mesh_spacing();

  /*
   * return total number of nodes
   */
  virtual unsigned int n_nodes();

  /*
   * return total number of bonds
   */
  virtual unsigned int n_bonds();

protected:
  int _pddim;
  double _horizon;

  struct pd_node * _node;

  double _mesh_spacing;
  unsigned int _total_nodes;
  unsigned int _total_bonds;

  std::vector<std::vector<dof_id_type> > _neighbors;
};

#endif // PERIDYNAMICMESH_H
