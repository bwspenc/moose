//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeometricCutUserObject.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/explicit_system.h"
#include "libmesh/equation_systems.h"
#include "Function.h"
#include "libmesh/enum_to_string.h"
#include "XFEMFuncs.h"

/**
 * InterfaceMeshCutUserObjectBase:
 * (1) reads in a mesh describing the interface,
 * (2) uses the mesh to do cutting of 2D/3D elements, and
 * (3) grows the mesh based on interface velocites.
 */

class XFEMMovingInterfaceVelocityBase;

class InterfaceMeshCutUserObjectBase : public GeometricCutUserObject
{
public:
  /**
   * Build the parameter set describing an interface mesh cut user object.
   *
   * @return Input parameters defining the interface mesh cutting configuration.
   */
  static InputParameters validParams();

  InterfaceMeshCutUserObjectBase(const InputParameters & parameters);

  /**
   * Perform initial setup prior to interface mesh cutting.
   */
  virtual void initialSetup() override;

  /**
   * Initialize state for the current execution interval.
   */
  virtual void initialize() override;

  /**
   * Retrieve the ordered crack front points to expose to other XFEM components.
   *
   * @param num_crack_front_points The number of crack front points requested.
   * @return Ordered crack front point coordinates suitable for postprocessing.
   */
  virtual const std::vector<Point>
  getCrackFrontPoints(unsigned int num_crack_front_points) const override;

  /**
   * Provide normals associated with the crack front points.
   *
   * @param num_crack_front_points The number of normals to populate.
   * @return Crack front normals corresponding to the requested point count.
   */
  virtual const std::vector<RealVectorValue>
  getCrackPlaneNormals(unsigned int num_crack_front_points) const override;

  /// Get the cutter mesh pointer
  std::shared_ptr<MeshBase> getCutterMesh() const { return _cutter_mesh; };

  /**
   * Calculate the signed distance for a point relative to the interface surface.
   *
   * @param p Coordinate of the point whose signed distance is requested.
   * @return Signed distance from @p p to the closest interface feature with the correct sign based
   * on the face or edge normal orientation.
   */
  virtual Real calculateSignedDistance(Point p) const = 0;

  /// return the normal at a node in the cutting mesh
  virtual Point nodeNormal(const unsigned int & node_id) = 0;

  /**
   * Calculate the element normal values for all elements in the cutting mesh.
   */
  virtual void calculateNormals() = 0;

  virtual CutSubdomainID getCutSubdomainID(const Node * node) const override;

protected:
  /// The cutter mesh
  std::shared_ptr<MeshBase> _cutter_mesh;
  /// node to element map of cut mesh
  std::unordered_map<dof_id_type, std::vector<dof_id_type>> _node_to_elem_map;
  /// initial nodes location
  std::unordered_map<dof_id_type, Point> _initial_nodes_location;
  /// Pointer to XFEMMovingInterfaceVelocityBase object
  const XFEMMovingInterfaceVelocityBase * _interface_velocity;
  /// Velocity function
  const Function * _func;
  /// The computation domain mesh
  MooseMesh & _mesh;
  /// Pointer to PointLocatorBase object
  std::unique_ptr<PointLocatorBase> _pl;
  /// Exodus for outputing points and values
  std::unique_ptr<ExodusII_IO> _exodus_io;
  /// Local explicit system for exodus output
  ExplicitSystem * _explicit_system;
  /// Local equation system for exodus output
  std::unique_ptr<EquationSystems> _equation_systems;
  /// Output exodus file
  const bool _output_exodus;
  /// The CutSubdomainID for the negative side of the cut
  const CutSubdomainID _negative_id;
  /// The CutSubdomainID for the positive side of the cut
  const CutSubdomainID _positive_id;
  /// Local variable number for displacement x
  unsigned int _var_num_disp_x;
  /// Local variable number for displacement y
  unsigned int _var_num_disp_y;
  /// Local variable number for displacement z
  unsigned int _var_num_disp_z;
};
