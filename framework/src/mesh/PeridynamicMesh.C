/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                        Peridynamics                          */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "PeridynamicMesh.h"

// libmesh includes
//#include "libmesh/nanoflann.hpp"

template<>
InputParameters validParams<PeridynamicMesh>()
{
  InputParameters params = validParams<MooseMesh>();
  params.addRequiredParam<Real>("horizon", "The horizon size");
  return params;
}

PeridynamicMesh::PeridynamicMesh(const InputParameters & parameters) :
  MooseMesh(parameters),
  _horizon(getParam<Real>("horizon"))
{
}

PeridynamicMesh::~PeridynamicMesh()
{
}

void
PeridynamicMesh::find_neighbor()
{
// TODO: use nanoflann class for efficient neighbor search
//  // convert the data structure
//  std::vector<Point> cloud(_total_nodes);
//  for (unsigned int i = 0; i < _total_nodes; ++i)
//    cloud[i] = _node[i].coord;

//  // construct a kd-tree index
//  typedef nanoflann::KDTreeSingleIndexAdaptor<
//          nanoflann::L2_Simple_Adaptor<double, Point>,
//          Point,
//          /* dim */ 3> kd_tree;
//  kd_tree index(/* dim */ 3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
//  index.build_index();

//  // search distance with tolerance, nanoflann use distance rather than radius
//  const double search_dist = 1.0001 * _horizon * _horizon;
//  std::vector<std::pair<int, double> > matches;
//  nanoflann::SearchParams params;

//  for (unsigned int i = 0; i < _total_nodes; ++i)
//  {
//    unsigned int nMatches = index.radiusSearch(&_node[i].coord, search_dist, matches, params);
//    _node[i].bond_num = nMatches;
//    for (unsigned int j = 0; j < nMatches; ++j)
//      _neighbors[i].push_back(matches[j].first);
//  }

  double dis = 0;
  for (unsigned int i = 0; i < _total_nodes; ++i)
  {
    unsigned int num = 0;
    for (unsigned int j = 0; j < _total_nodes; ++j)
    {
      dis = std::sqrt(std::pow((_node[i].coord)(0) - (_node[j].coord)(0), 2) + std::pow((_node[i].coord)(1) - (_node[j].coord)(1), 2) + std::pow((_node[i].coord)(2) - (_node[j].coord)(2), 2));
      if (dis <= 1.0001 * _horizon && j != i)
      {
        _neighbors[i].push_back(j);
        ++num;
      }
    }
    _node[i].bond_num = num;
  }
}

std::vector<dof_id_type>
PeridynamicMesh::neighbors(dof_id_type node_id)
{
  return _neighbors[node_id];
}

Point
PeridynamicMesh::coord(dof_id_type node_id)
{
  return _node[node_id].coord;
}

double
PeridynamicMesh::volume(dof_id_type node_id)
{
  return _node[node_id].volume;
}

unsigned int
PeridynamicMesh::n_neighbors(dof_id_type node_id)
{
  return _neighbors[node_id].size();
}

int
PeridynamicMesh::dim()
{
  return _pddim;
}

double
PeridynamicMesh::mesh_spacing()
{
  return _mesh_spacing;
}

unsigned int
PeridynamicMesh::n_nodes()
{
  return _total_nodes;
}

unsigned int
PeridynamicMesh::n_bonds()
{
  return _total_bonds;
}

