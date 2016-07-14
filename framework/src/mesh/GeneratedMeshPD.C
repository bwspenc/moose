/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "GeneratedMeshPD.h"

// libMesh includes
#include "libmesh/edge_edge2.h"

template<>
InputParameters validParams<GeneratedMeshPD>()
{
  InputParameters params = validParams<PeridynamicMesh>();
  params.addRequiredParam<int>("pddim", "PD mesh dimension");
  params.addRequiredParam<int>("shape", "1: Rectangular; 0: Circular");// the cross section shape
  params.addRangeCheckedParam<int>("nx", 1, "nx > 0", "Number of elements in the X/R direction");
  params.addParam<Real>("xmin", 0.0, "Lower X Coordinate of the generated mesh");
  params.addParam<Real>("ymin", 0.0, "Lower Y Coordinate of the generated mesh");
  params.addParam<Real>("zmin", 0.0, "Lower Z Coordinate of the generated mesh");
  params.addParam<Real>("xmax", 1.0, "Upper X Coordinate of the generated mesh");
  params.addParam<Real>("ymax", 1.0, "Upper Y Coordinate of the generated mesh");
  params.addParam<Real>("zmax", 1.0, "Upper Z Coordinate of the generated mesh");
  params.addParam<Real>("R", 1.0, "Radius of the circular domain if applicable");
  return params;
}

GeneratedMeshPD::GeneratedMeshPD(const InputParameters & parameters) :
  PeridynamicMesh(parameters),
  _nx(getParam<int>("nx")),
  _shape(getParam<int>("shape")),
  _xmin(getParam<Real>("xmin")),
  _ymin(getParam<Real>("ymin")),
  _zmin(getParam<Real>("zmin")),
  _xmax(getParam<Real>("xmax")),
  _ymax(getParam<Real>("ymax")),
  _zmax(getParam<Real>("zmax")),
  _R(getParam<Real>("R"))
{
  _pddim = getParam<int>("pddim");
}

GeneratedMeshPD::~GeneratedMeshPD()
{
}

MooseMesh &
GeneratedMeshPD::clone() const
{
  return *(new GeneratedMeshPD(*this));
}

void
GeneratedMeshPD::buildMesh()
{
  UnstructuredMesh & mesh = dynamic_cast<UnstructuredMesh &> (getMesh());
  mesh.clear();
  mesh.set_mesh_dimension(1);
  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  _total_bonds = 0;
  _total_nodes = 0;

  if (_pddim == 2 && _shape == 1)// 2D rectangular domain
    build2DRectangular(mesh, boundary_info);
  else if (_pddim == 2 && _shape == 0)// 2D circular domain
    build2DCircular(mesh, boundary_info);
  else if (_pddim == 3 && _shape == 1)// 3D rectangular domain
    build3DRectangular(mesh, boundary_info);
  else if (_pddim == 3 && _shape == 0)// 3D cylindrical domain
    build3DCylindrical(mesh, boundary_info);
  else
   mooseError("The domain cross section shape can only be 0: circular and 1: rectangular!");

//-----------------------
//  double val = 0;
//  for (unsigned int i = 0; i < _total_nodes; ++i)
//  {
//    if ((_node[i].coord)(1) < -0.0001)
//    {   
//      double voli = 0;
//      for(unsigned int k = 0; k < _neighbors[i].size(); ++k)
//        voli += _node[_neighbors[i][k]].volume;

//      for(unsigned int j = 0; j < _neighbors[i].size(); ++j)
//      {   
//        if ((_node[_neighbors[i][j]].coord)(1) > 0.0001)
//        {   
//          double volj = 0;
//          for(unsigned int k = 0; k < _neighbors[_neighbors[i][j]].size(); ++k)
//            volj += _node[_neighbors[_neighbors[i][j]][k]].volume;
//          double dx = (_node[i].coord)(0) - (_node[_neighbors[i][j]].coord)(0);
//          double dy = (_node[i].coord)(1) - (_node[_neighbors[i][j]].coord)(1);
//          double cosij = dy / std::sqrt(dx * dx + dy * dy);
//          val += _node[i].volume * _node[_neighbors[i][j]].volume * (1 / voli + 1 / volj) * cosij * cosij;
//        }   
//      }   
//    }   
//  }
//  double s = std::sqrt(2 * 0.004 * 8.2 / (239836.911 * 4 * val));
//  std::cout << s << std::endl;
//-----------------------

  std::cout << "Total Node Number: " << _total_nodes << std::endl;
  std::cout << "Total Bond Number: " << _total_bonds << std::endl;

  // prepare for use
  mesh.prepare_for_use (/*skip_renumber =*/ true);
}

void
GeneratedMeshPD::build2DRectangular(UnstructuredMesh & mesh, BoundaryInfo & boundary_info)
{
  double X, Y;
  unsigned int k = 0;
  double spacing = (_xmax - _xmin) / _nx;
  double horizon = PeridynamicMesh::computeHorizon(spacing);
  unsigned int ny = static_cast<int> ((_ymax - _ymin) / spacing);
  _total_nodes = _nx * ny;
  _node = (struct pd_node*) malloc(_total_nodes * sizeof(struct pd_node));
  _neighbors.resize(_total_nodes);
  mesh.reserve_nodes(_total_nodes);

  // define nodal coordinates
  for (unsigned int j = 0; j < ny; ++j)
    for (unsigned int i = 0; i < _nx; ++i)
    {
      X = _xmin + i * spacing + spacing / 2;
      Y = _ymin + j * spacing + spacing / 2;
      mesh.add_point(Point(X, Y, 0.0), k);
      _node[k].coord = Point(X, Y, 0.0);
      _node[k].mesh_spacing = spacing;
      _node[k].horizon = horizon;
      _node[k].volume = spacing * spacing;
      _node[k].volumesum = 0.0;
      ++k;
    }

  // search node neighbor
  PeridynamicMesh::find_neighbor();

  // generate PD mesh
  for (unsigned int i = 0; i < _total_nodes; ++i)
    _total_bonds += _neighbors[i].size();

  _total_bonds /= 2;
  mesh.reserve_elem(_total_bonds);
  for (unsigned int i = 0; i < _total_nodes; ++i)
    for (unsigned int j = 0; j < _neighbors[i].size(); ++j)
      if (_neighbors[i][j] > i)
      {
        Elem* elem = mesh.add_elem(new Edge2);
        elem->set_node(0) = mesh.node_ptr(i);
        elem->set_node(1) = mesh.node_ptr(_neighbors[i][j]);
        elem->subdomain_id() = 0;
      }

  // define boundary nodeset
  for (unsigned int i = 0; i < _total_nodes; ++i)
  {
    X = (_node[i].coord)(0);
    Y = (_node[i].coord)(1);
    if (X < _xmin + spacing)
      boundary_info.add_node(mesh.node_ptr(i), 0);
    if (X > _xmax - spacing)
      boundary_info.add_node(mesh.node_ptr(i), 1);
    if (Y < _ymin + spacing)
      boundary_info.add_node(mesh.node_ptr(i), 2);
    if (Y > _ymax - spacing)
      boundary_info.add_node(mesh.node_ptr(i), 3);
    if (std::abs((_xmax - _xmin) * (_ymin - Y) - (_xmin - X) * (_ymax - _ymin)) / std::sqrt(std::pow(_xmax - _xmin, 2) + std::pow(_ymax - _ymin, 2)) < 0.1 * spacing)
      boundary_info.add_node(mesh.node_ptr(i), 4);
  }
  boundary_info.nodeset_name(0) = "Left";
  boundary_info.nodeset_name(1) = "Right";
  boundary_info.nodeset_name(2) = "Bottom";
  boundary_info.nodeset_name(3) = "Top";
  boundary_info.nodeset_name(4) = "UpDiag";
}

void
GeneratedMeshPD::build2DCircular(UnstructuredMesh & mesh, BoundaryInfo & boundary_info)
{
  double X, Y;
  unsigned int k = 0;
  double spacing = 2 * _R / (2 * _nx - 1);
  double horizon = PeridynamicMesh::computeHorizon(spacing);
  // calculate the total node number by cutting out a circular domain from square patch
  for (unsigned int j = 0; j < 2 * _nx - 1; ++j)
    for (unsigned int i = 0; i < 2 * _nx - 1; ++i)
    {
      X = - _R + i * spacing + spacing / 2;
      Y = - _R + j * spacing + spacing / 2;
      if (X * X + Y * Y < _R * _R)
        ++_total_nodes;
    }
  _node = (struct pd_node*) malloc(_total_nodes * sizeof(struct pd_node));
  _neighbors.resize(_total_nodes);
  mesh.reserve_nodes(_total_nodes);

  // define nodal coordinates
  for (unsigned int j = 0; j < 2 * _nx - 1; ++j)
    for (unsigned int i = 0; i < 2 * _nx - 1; ++i)
    {
      X = - _R + i * spacing + spacing / 2;
      Y = - _R + j * spacing + spacing / 2;
      if (X * X + Y * Y < _R * _R)
      {
        mesh.add_point(Point(X, Y, 0.0), k);
        _node[k].coord = Point(X, Y, 0.0);
        _node[k].mesh_spacing = spacing;
        _node[k].horizon = horizon;
        _node[k].volume = spacing * spacing;
        _node[k].volumesum = 0.0;
        ++k;
      }
    }

  // search node neighbor
  PeridynamicMesh::find_neighbor();

  // generate PD mesh
  for(unsigned int i = 0; i < _total_nodes; ++i)
    _total_bonds += _neighbors[i].size();

  _total_bonds /= 2;
  mesh.reserve_elem(_total_bonds);
  for (unsigned int i = 0; i < _total_nodes; ++i)
    for (unsigned int j = 0; j < _neighbors[i].size(); ++j)
      if (_neighbors[i][j] > i)
      {
        Elem* elem = mesh.add_elem(new Edge2);
        elem->set_node(0) = mesh.node_ptr(i);
        elem->set_node(1) = mesh.node_ptr(_neighbors[i][j]);
        elem->subdomain_id() = 0;
      }

  // define boundary nodeset
  for (unsigned int i = 0; i < _total_nodes; ++i)
  {
    X = (_node[i].coord)(0);
    Y = (_node[i].coord)(1);
    double dis = std::sqrt(X * X + Y * Y);
    if (dis > _R - spacing)
      boundary_info.add_node(mesh.node_ptr(i), 0);
    if (dis < 0.001 * spacing)
      boundary_info.add_node(mesh.node_ptr(i), 1);
    if (std::abs(Y) < 0.001 * spacing && dis > _R - spacing && X < 0.0)
      boundary_info.add_node(mesh.node_ptr(i), 2);
    if (std::abs(Y) < 0.001 * spacing && dis > _R - spacing && X > 0.0)
      boundary_info.add_node(mesh.node_ptr(i), 3);
    if (std::abs(X) < 0.001 * spacing && dis > _R - spacing && Y < 0.0)
      boundary_info.add_node(mesh.node_ptr(i), 4);
    if (std::abs(X) < 0.001 * spacing && dis > _R - spacing && Y > 0.0)
      boundary_info.add_node(mesh.node_ptr(i), 5);
  }
  boundary_info.nodeset_name(0) = "Periphery";
  boundary_info.nodeset_name(1) = "CenterPoint";
  boundary_info.nodeset_name(2) = "LeftPoint";
  boundary_info.nodeset_name(3) = "RightPoint";
  boundary_info.nodeset_name(4) = "BottomPoint";
  boundary_info.nodeset_name(5) = "TopPoint";
}

void
GeneratedMeshPD::build3DRectangular(UnstructuredMesh & mesh, BoundaryInfo & boundary_info)
{
  double X, Y, Z;
  unsigned int k = 0;
  double spacing = (_xmax - _xmin) / _nx;
  double horizon = PeridynamicMesh::computeHorizon(spacing);
  unsigned int ny = static_cast<int> ((_ymax - _ymin) / spacing);
  unsigned int nz = static_cast<int> ((_zmax - _zmin) / spacing);
  _total_nodes = _nx * ny * nz;
  _node = (struct pd_node*) malloc(_total_nodes * sizeof(struct pd_node));
  _neighbors.resize(_total_nodes);
  mesh.reserve_nodes(_total_nodes);

  // define nodal coordinates
  for (unsigned int n = 0; n < nz; ++n)
    for (unsigned int j = 0; j < ny; ++j)
      for (unsigned int i = 0; i < _nx; ++i)
      {
        X = _xmin + i * spacing + spacing / 2;
        Y = _ymin + j * spacing + spacing / 2;
        Z = _zmin + n * spacing + spacing / 2;
        mesh.add_point(Point(X, Y, Z), k);
        _node[k].coord = Point(X, Y, Z);
        _node[k].mesh_spacing = spacing;
        _node[k].horizon = horizon;
        _node[k].volume = spacing * spacing * spacing;
        _node[k].volumesum = 0.0;
        ++k;
      }

  // search node neighbor
  PeridynamicMesh::find_neighbor();

  // generate PD mesh
  for (unsigned int i = 0; i < _total_nodes; ++i)
    _total_bonds += _neighbors[i].size();

  _total_bonds /= 2;
  mesh.reserve_elem(_total_bonds);
  for (unsigned int i = 0; i < _total_nodes; ++i)
    for (unsigned int j = 0; j < _neighbors[i].size(); ++j)
      if (_neighbors[i][j] > i)
      {
        Elem* elem = mesh.add_elem (new Edge2);
        elem->set_node(0) = mesh.node_ptr(i);
        elem->set_node(1) = mesh.node_ptr(_neighbors[i][j]);
        elem->subdomain_id() = 0;
      }

  // define boundary nodeset
  for (unsigned int i = 0; i < _total_nodes; ++i)
  {
    X = (_node[i].coord)(0);
    Y = (_node[i].coord)(1);
    Z = (_node[i].coord)(2);
    if (Y < _ymin + spacing)
      boundary_info.add_node(mesh.node_ptr(i), 0);
    if (Y > _ymax - spacing)
      boundary_info.add_node(mesh.node_ptr(i), 1);
    if (Z < _zmin + spacing)
      boundary_info.add_node(mesh.node_ptr(i), 2);
    if (Z > _zmax - spacing)
      boundary_info.add_node(mesh.node_ptr(i), 3);
    if (X < _xmin + spacing)
      boundary_info.add_node(mesh.node_ptr(i), 4);
    if (X > _xmax - spacing)
      boundary_info.add_node(mesh.node_ptr(i), 5);
    if (std::sqrt(((std::pow(_xmin - X, 2) + std::pow(_ymin - Y, 2) + std::pow(_zmin - Z, 2)) * (std::pow(_xmax - _xmin, 2) + std::pow(_ymax - _ymin, 2) + std::pow(_zmax - _zmin, 2)) - std::pow((_xmin - X) * (_xmax - _xmin) + (_ymin - Y) * (_ymax - _ymin) + (_zmin - Z) * (_zmax - _zmin), 2)) / (std::pow(_xmax - _xmin, 2) + std::pow(_ymax - _ymin, 2) + std::pow(_zmax - _zmin, 2))) < 0.1 * spacing)
      boundary_info.add_node(mesh.node_ptr(i), 6);
  }
  boundary_info.nodeset_name(0) = "Left";
  boundary_info.nodeset_name(1) = "Right";
  boundary_info.nodeset_name(2) = "Bottom";
  boundary_info.nodeset_name(3) = "Top";
  boundary_info.nodeset_name(4) = "Back";
  boundary_info.nodeset_name(5) = "Front";
  boundary_info.nodeset_name(6) = "UpDiag";
}

void
GeneratedMeshPD::build3DCylindrical(UnstructuredMesh & mesh, BoundaryInfo & boundary_info)
{
  double X, Y, Z;
  unsigned int k = 0;
  double spacing = 2.0 * _R / (2 * _nx - 1);
  double horizon = PeridynamicMesh::computeHorizon(spacing);
  unsigned int nz = static_cast<int> ((_zmax - _zmin) / spacing);

  // calculate the total node number by cutting out a cylindrical domain from rectangular domain
  for (unsigned int n = 0; n < nz; ++n)
    for (unsigned int j = 0; j < 2 * _nx - 1; ++j)
      for (unsigned int i = 0; i < 2 * _nx - 1; ++i)
      {
        X = - _R + i * spacing + spacing / 2;
        Y = - _R + j * spacing + spacing / 2;
        if (X * X + Y * Y < _R * _R)
          ++_total_nodes;
      }
  _node = (struct pd_node*) malloc(_total_nodes * sizeof(struct pd_node));
  _neighbors.resize(_total_nodes);
  mesh.reserve_nodes(_total_nodes);

  // define nodal coordinates
  for (unsigned int n = 0; n < nz; ++n)
    for (unsigned int j = 0; j < 2 * _nx - 1; ++j)
      for (unsigned int i = 0; i < 2 * _nx - 1; ++i)
      {
        X = - _R + i * spacing + spacing / 2;
        Y = - _R + j * spacing + spacing / 2;
        Z = _zmin + n * spacing + spacing / 2;
        if (X * X + Y * Y < _R * _R)
        {
          mesh.add_point(Point(X, Y, Z), k);
          _node[k].coord = Point(X, Y, Z);
          _node[k].mesh_spacing = spacing;
          _node[k].horizon = horizon;
          _node[k].volume = spacing * spacing * spacing;
          _node[k].volumesum = 0.0;
          ++k;
        }
      }

  // earch node neighbor
  PeridynamicMesh::find_neighbor();

  // generate PD mesh
  for (unsigned int i = 0; i < _total_nodes; ++i)
    _total_bonds += _neighbors[i].size();

  _total_bonds /= 2;
  mesh.reserve_elem(_total_bonds);
  for (unsigned int i = 0; i < _total_nodes; ++i)
    for (unsigned int j = 0; j < _neighbors[i].size(); ++j)
      if (_neighbors[i][j] > i)
      {
        Elem* elem = mesh.add_elem(new Edge2);
        elem->set_node(0) = mesh.node_ptr(i);
        elem->set_node(1) = mesh.node_ptr(_neighbors[i][j]);
        elem->subdomain_id() = 0;
      }

  // define boundary nodeset
  for (unsigned int i = 0; i < _total_nodes; ++i)
  {
    X = (_node[i].coord)(0);
    Y = (_node[i].coord)(1);
    Z = (_node[i].coord)(2);
    double dis = std::sqrt(X * X + Y * Y);
    if (dis > _R - spacing)
      boundary_info.add_node(mesh.node_ptr(i), 0);
    if (dis < 0.001 * spacing)
      boundary_info.add_node(mesh.node_ptr(i), 1);
    if (dis < 0.001 * spacing && std::abs(Z - (_zmax + _zmin) / 2) < spacing)
      boundary_info.add_node(mesh.node_ptr(i), 2);
    if (std::abs(Y) < 0.001 * spacing && dis > _R - spacing && X < 0.0)
      boundary_info.add_node(mesh.node_ptr(i), 3);
    if (std::abs(Y) < 0.001 * spacing && dis > _R - spacing && X > 0.0)
      boundary_info.add_node(mesh.node_ptr(i), 4);
    if (std::abs(X) < 0.001 * spacing && dis > _R - spacing && Y < 0.0)
      boundary_info.add_node(mesh.node_ptr(i), 5);
    if (std::abs(X) < 0.001 * spacing && dis > _R - spacing && Y > 0.0)
      boundary_info.add_node(mesh.node_ptr(i), 6);
    if (Z < _zmin + spacing)
      boundary_info.add_node(mesh.node_ptr(i), 7);
    if (Z > _zmax - spacing)
      boundary_info.add_node(mesh.node_ptr(i), 8);
  }
  boundary_info.nodeset_name(0) = "Periphery";
  boundary_info.nodeset_name(1) = "CenterLine";
  boundary_info.nodeset_name(2) = "Center";
  boundary_info.nodeset_name(3) = "LeftLine";
  boundary_info.nodeset_name(4) = "RightLine";
  boundary_info.nodeset_name(5) = "BackLine";
  boundary_info.nodeset_name(6) = "FrontLine";
  boundary_info.nodeset_name(7) = "Bottom";
  boundary_info.nodeset_name(8) = "Top";
}
