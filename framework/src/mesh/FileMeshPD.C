/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                         Peridynamics                         */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "FileMeshPD.h"
#include "MooseMesh.h"

// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/edge_edge2.h"

template<>
InputParameters validParams<FileMeshPD>()
{
  InputParameters params = validParams<PeridynamicMesh>();

  params.addRequiredParam<MeshFileName>("file", "Name of the mesh file (must be exodusII file)");
  return params;
}

FileMeshPD::FileMeshPD(const InputParameters & parameters) :
  PeridynamicMesh(parameters)
{
}

FileMeshPD::~FileMeshPD()
{
}

MooseMesh &
FileMeshPD::clone() const
{
  return *(new FileMeshPD(*this));
}

void
FileMeshPD::buildMesh()
{
  // read the temporary mesh from Exodus file
  std::string _file_name = getParam<MeshFileName>("file");
  MooseUtils::checkFileReadable(_file_name);
  MeshBase * fe_mesh = new SerialMesh(_communicator);
  ExodusII_IO * _exodusII_io = new ExodusII_IO(*fe_mesh);
  _exodusII_io->read(_file_name);
  fe_mesh->allow_renumbering(false);
  fe_mesh->prepare_for_use(/*true*/);

  _mesh_spacing = 0;
  _total_nodes = fe_mesh->n_elem();
  _pddim = fe_mesh->mesh_dimension();
  BoundaryInfo & fe_boundary_info = fe_mesh->get_boundary_info();

  // initialize PD mesh
  UnstructuredMesh & pd_mesh = dynamic_cast<UnstructuredMesh &> (getMesh());
  pd_mesh.clear();
  pd_mesh.set_mesh_dimension(_pddim);
  BoundaryInfo & pd_boundary_info = pd_mesh.get_boundary_info();
  pd_mesh.reserve_nodes(_total_nodes);
  _node = (struct pd_node*) malloc(_total_nodes * sizeof(struct pd_node));
  _neighbors.resize(_total_nodes);

  // loop through all fe elements to generate PD nodes structure
  for (MeshBase::element_iterator it = fe_mesh->elements_begin(); it != fe_mesh->elements_end(); ++it)
  {
    Elem *fe_elem = *it;
    _node[fe_elem->id()].coord = fe_elem->centroid();
    _node[fe_elem->id()].volume = fe_elem->volume();
    pd_mesh.add_point(fe_elem->centroid(), fe_elem->id());
  }

  // search node neighbor
  PeridynamicMesh::find_neighbor();

  // generate PD mesh
  _total_bonds = 0;
  for (unsigned int i = 0; i < _total_nodes; ++i)
    _total_bonds += _node[i].bond_num;

  _total_bonds /= 2;
  pd_mesh.reserve_elem(_total_bonds);
  for (unsigned int i = 0; i < _total_nodes; ++i)
    for(unsigned int j = 0; j < _node[i].bond_num; ++j)
      if (_neighbors[i][j] > i)
      {
        Elem* pd_elem = pd_mesh.add_elem(new Edge2);
        pd_elem->set_node(0) = pd_mesh.node_ptr(i);
        pd_elem->set_node(1) = pd_mesh.node_ptr(_neighbors[i][j]);
        pd_elem->subdomain_id() = 0;
      }

  // convert boundary info from fe_boundary_info to pd_boundary_info
  // build element list for user specified nodal boundaries
  std::vector<dof_id_type> elems;
  std::vector<unsigned short int> sides;
  std::vector<boundary_id_type> ids;
  fe_boundary_info.build_side_list_from_node_list();
  fe_boundary_info.build_active_side_list(elems, sides, ids);

  unsigned int n = elems.size();
  // array of boundary elems
  std::vector<BndElement *> _bnd_elems(n);
  // map of set of elem IDs connected to each boundary
  std::map<boundary_id_type, std::set<dof_id_type> > _bnd_elem_ids;
  for (unsigned int i = 0; i < n; ++i)
  {
    _bnd_elems[i] = new BndElement(fe_mesh->elem_ptr(elems[i]), sides[i], ids[i]);
    _bnd_elem_ids[ids[i]].insert(elems[i]);
  }
  // for nodeset only
  std::set<boundary_id_type> fe_node_bid(fe_boundary_info.get_node_boundary_ids());
  std::set<boundary_id_type> fe_edge_bid(fe_boundary_info.get_edge_boundary_ids());
  if (!fe_edge_bid.empty())
    mooseError("FileMeshPD currently only accepts nodesets!");

  for (std::set<boundary_id_type>::iterator bit = fe_node_bid.begin(); bit != fe_node_bid.end(); ++bit)
  {
    pd_boundary_info.nodeset_name(*bit) = fe_boundary_info.get_nodeset_name(*bit);
    for (MeshBase::element_iterator eit = fe_mesh->elements_begin(); eit != fe_mesh->elements_end(); ++eit)
    {
      Elem *fe_elem = *eit;
      std::map<boundary_id_type, std::set<dof_id_type> >::const_iterator it = _bnd_elem_ids.find(*bit);
      if (it != _bnd_elem_ids.end())
        if (it->second.find(fe_elem->id()) != it->second.end())
          pd_boundary_info.add_node(pd_mesh.node_ptr(fe_elem->id()), *bit);
    }
  }

  delete fe_mesh;
  delete _exodusII_io;

  std::cout << "Total Node Number: " << _total_nodes << std::endl;
  std::cout << "Total Bond Number: " << _total_bonds << std::endl;

  // prepare for use
  pd_mesh.prepare_for_use (/*skip_renumber =*/ false);
}
