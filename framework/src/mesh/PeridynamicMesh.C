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
  params.addParam<Real>("horizon_size", "Horizon size");
  params.addParam<unsigned int>("horizon_number", "The ratio number of horizon radius to the effective mesh spacing");
  return params;
}

PeridynamicMesh::PeridynamicMesh(const InputParameters & parameters) :
  MooseMesh(parameters),
  _horizon_size(isParamValid("horizon_size") ? getParam<Real>("horizon_size") : 0),
  _horizon_number(isParamValid("horizon_number") ? getParam<unsigned int>("horizon_number") : 0)
{
  if (isParamValid("horizon_size") && isParamValid("horizon_number"))
    mooseError("You can specify only one option: horizon_size or horizon_number in the mesh block!");

  if (!isParamValid("horizon_size") && !isParamValid("horizon_number"))
    mooseError("You must specify one option: horizon_size or horizon_number in the mesh block!");
}

PeridynamicMesh::~PeridynamicMesh()
{
}

double
PeridynamicMesh::computeHorizon(double spacing)
{
  if (isParamValid("horizon_number"))
    return _horizon_number * spacing;
  else
    return _horizon_size;
}

void
PeridynamicMesh::findNodeNeighbor()
{
// TODO: use nanoflann class for efficient neighbor search
//  // convert the data structure
//  std::vector<Point> cloud(_total_nodes);
//  for (unsigned int i = 0; i < _total_nodes; ++i)
//    cloud[i] = _pdnode[i].coord;

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
//    unsigned int nMatches = index.radiusSearch(&_pdnode[i].coord, search_dist, matches, params);
//    _pdnode[i].bond_num = nMatches;
//    for (unsigned int j = 0; j < nMatches; ++j)
//      _node_neighbors[i].push_back(matches[j].first);
//  }

  for (unsigned int i = 0; i < _total_nodes; ++i)
  {
    double dis = 0;
    for (unsigned int j = 0; j < _total_nodes; ++j)
    {
      dis = std::sqrt(std::pow((_pdnode[i].coord)(0) - (_pdnode[j].coord)(0), 2) + std::pow((_pdnode[i].coord)(1) - (_pdnode[j].coord)(1), 2) + std::pow((_pdnode[i].coord)(2) - (_pdnode[j].coord)(2), 2));
      if (dis <= 1.0001 * _pdnode[i].horizon && j != i)
      {
        // check whether j was already considered as a neighbor of i, if not, add j to i's neighborlist
        if (std::find(_node_neighbors[i].begin(), _node_neighbors[i].end(), j) == _node_neighbors[i].end())
        {
          _node_neighbors[i].push_back(j);
          _pdnode[i].volumesum += _pdnode[j].volume;
        }
        // check whether i was also considered as a neighbor of j, if not, add i to j's neighborlist
        if (std::find(_node_neighbors[j].begin(), _node_neighbors[j].end(), i) == _node_neighbors[j].end())
        {
          _node_neighbors[j].push_back(i);
          _pdnode[j].volumesum += _pdnode[i].volume;
        }
      }
    }
  }
}

std::vector<dof_id_type>
PeridynamicMesh::neighbors(dof_id_type node_id)
{
  return _node_neighbors[node_id];
}

std::vector<dof_id_type>
PeridynamicMesh::bonds(dof_id_type node_id)
{
  return _node_bonds[node_id];
}

Point
PeridynamicMesh::coord(dof_id_type node_id)
{
  return _pdnode[node_id].coord;
}

double
PeridynamicMesh::volume(dof_id_type node_id)
{
  return _pdnode[node_id].volume;
}

double
PeridynamicMesh::volumesum(dof_id_type node_id)
{
  return _pdnode[node_id].volumesum;
}

unsigned int
PeridynamicMesh::n_neighbors(dof_id_type node_id)
{
  return _node_neighbors[node_id].size();
}

int
PeridynamicMesh::dim()
{
  return _pddim;
}

double
PeridynamicMesh::mesh_spacing(dof_id_type node_id)
{
  return _pdnode[node_id].mesh_spacing;
}

double
PeridynamicMesh::horizon(dof_id_type node_id)
{
  return _pdnode[node_id].horizon;
}

unsigned int
PeridynamicMesh::total_nodes()
{
  return _total_nodes;
}

unsigned int
PeridynamicMesh::total_bonds()
{
  return _total_bonds;
}
