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
#include "CrackFrontDefinition.h"

#include <array>

class Function;

/**
 * CrackMeshCut3DUserObject: (1) reads in a mesh describing the crack surface,
 * (2) uses the mesh to do initial cutting of 3D elements, and
 * (3) grows the mesh based on prescribed growth functions.
 */

class CrackMeshCut3DUserObject : public GeometricCutUserObject
{
public:
  static InputParameters validParams();

  CrackMeshCut3DUserObject(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override;

  virtual const std::vector<Point>
  getCrackFrontPoints(unsigned int num_crack_front_points) const override;
  virtual const std::vector<RealVectorValue>
  getCrackPlaneNormals(unsigned int num_crack_front_points) const override;

  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutEdge> & cut_edges,
                                    std::vector<Xfem::CutNode> & cut_nodes) const override;
  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutFace> & cut_faces) const override;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_edges,
                                     std::vector<Xfem::CutEdge> & cut_edges) const override;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_faces,
                                     std::vector<Xfem::CutFace> & cut_faces) const override;

  /**
   * Find all active boundary nodes in the cutter mesh.
   *
   * Boundary nodes that lie outside the structural mesh are marked inactive.
   */
  void findActiveBoundaryNodes();

  /**
   * Get the indices of crack-front points within the active segment.
   *
   * @return Vector of indices where -1 indicates an inactive point and nonnegative values map to
   * entries in the crack front definition.
   */
  std::vector<int> getFrontPointsIndex();

  /**
   * Record the growth size at the active boundary for the mesh cutter.
   *
   * @param growth_size Growth increments to apply at each active boundary node.
   */
  void setSubCriticalGrowthSize(std::vector<Real> & growth_size);

  /**
   * Get the total number of crack front points tracked by the cutter mesh.
   *
   * @return Total number of crack front points available in the mesh cutter.
   */
  unsigned int getNumberOfCrackFrontPoints() const;

protected:
  /// The cutter mesh
  std::unique_ptr<MeshBase> _cut_mesh;

  /// The cutter mesh has triangluar elements only
  const unsigned int _cut_elem_nnode = 3;
  const unsigned int _cut_elem_dim = 2;

  /// The structural mesh
  MooseMesh & _mesh;

  /// The crack front definition
  CrackFrontDefinition * _crack_front_definition;

  /// updated crack front definition
  /// they are in the same order as defined in the input but the number of nodes may increase
  /// its difference from _front is that: _front does not necessarily follow the order of crack front definition
  /// therefore, _crack_front_points is generated from _front with the order of crack front definition
  /// limitation: this approach does not currently support the growth of one crack front into two
  std::vector<dof_id_type> _crack_front_points;

  /// Enum to for crack growth direction
  enum class GrowthDirectionEnum
  {
    MAX_HOOP_STRESS,
    FUNCTION
  };
  /// The direction method for growing mesh at the front
  const GrowthDirectionEnum _growth_dir_method;

  /// Enum to for crack growth rate
  enum class GrowthRateEnum
  {
    FATIGUE,
    FUNCTION
  };
  /// The rate method for growing mesh at the front
  const GrowthRateEnum _growth_rate_method;

  /// The structural mesh must be 3D only
  const unsigned int _elem_dim = 3;

  /// Used to define intersection points
  const Real _const_intersection = 0.01;

  /// Used for cutter mesh refinement and front advancement
  Real _size_control;

  /// Number of steps to grow the mesh
  unsigned int _n_step_growth;

  /// Variables to help control the work flow
  bool _stop;
  bool _grow;

  /// Boundary nodes of the cutter mesh
  std::vector<dof_id_type> _boundary;

  /// Active boundary nodes where growth is allowed
  std::vector<std::vector<dof_id_type>> _active_boundary;

  /// Inactive boundary
  std::vector<unsigned int> _inactive_boundary_pos;

  /// Front nodes that are grown from the crack front definition defined in the input
  /// therefore, they are (1) in the same order as defined in the input and (2) the number of nodes does not change
  std::vector<dof_id_type> _tracked_crack_front_points;

  bool _cfd;

  /// Edges at the boundary
  std::set<Xfem::CutEdge> _boundary_edges;

  /// A map of boundary nodes and their neighbors
  std::map<dof_id_type, std::vector<dof_id_type>> _boundary_map;

  /// Growth direction for active boundaries
  std::vector<std::vector<Point>> _active_direction;

  /// Growth size for the active boundary in a subcritical simulation
  std::vector<Real> _growth_size;

  /// Fatigue life
  std::vector<unsigned long int> _dn;
  std::vector<unsigned long int> _n;

  /// New boundary after growth
  std::vector<std::vector<dof_id_type>> _front;

  /// Indicator that shows if the cutting mesh is modified or not in this calculation step
  bool _is_mesh_modified;

  /// Total number of crack front points in the mesh cutter
  unsigned int _num_crack_front_points;

  /**
   * Check if a line intersects with an element defined by its vertices.
   *
   * @param p1 First endpoint of the line segment.
   * @param p2 Second endpoint of the line segment.
   * @param _vertices Vertices defining the element.
   * @param point Intersection point when an intersection is detected.
   * @return True if the line intersects the element; false otherwise.
   */
  virtual bool intersectWithEdge(const Point & p1,
                                 const Point & p2,
                                 const std::vector<Point> & _vertices,
                                 Point & point) const;

  /**
   * Find the intersection along the positive extension of the vector from @p p1 to @p p2.
   *
   * @param p1 Starting point of the direction vector.
   * @param p2 Ending point defining the direction vector.
   * @param vertices Vertices defining the candidate element.
   * @param point Intersection point when the direction intersects the element.
   * @return True if an intersection is found along the positive direction.
   */
  bool findIntersection(const Point & p1,
                        const Point & p2,
                        const std::vector<Point> & vertices,
                        Point & point) const;

  /**
   * Check if a point lies inside the segment defined by @p p1 and @p p2.
   *
   * @param p1 First endpoint of the segment.
   * @param p2 Second endpoint of the segment.
   * @param p Query point.
   * @return True if @p p lies on the edge; false otherwise.
   */
  bool isInsideEdge(const Point & p1, const Point & p2, const Point & p) const;

  /**
   * Get the relative position of a point measured from @p p1 along the segment to @p p2.
   *
   * @param p1 First endpoint of the segment.
   * @param p2 Second endpoint of the segment.
   * @param p Query point.
   * @return Relative position of @p p measured from @p p1 along the segment.
   */
  Real getRelativePosition(const Point & p1, const Point & p2, const Point & p) const;

  /**
   * Check if a point lies inside a plane defined by a set of vertices.
   *
   * @param _vertices Vertices defining the plane.
   * @param p Query point.
   * @return True if the point is inside the plane; false otherwise.
   */
  bool isInsideCutPlane(const std::vector<Point> & _vertices, const Point & p) const;

  /**
   * Find boundary nodes of the cutter mesh.
   *
   * This simple algorithm is based on checking whether the accumulated angle equals 360 degrees.
   */
  void findBoundaryNodes();

  /**
   * Find boundary edges of the cutter mesh.
   */
  void findBoundaryEdges();

  /**
   * Sort boundary nodes into the correct order along the boundary loop.
   */
  void sortBoundaryNodes();

  /**
   * Find the distance between two boundary nodes.
   *
   * @param node1 Identifier of the first node.
   * @param node2 Identifier of the second node.
   * @return Euclidean distance between the nodes.
   */
  Real findDistance(dof_id_type node1, dof_id_type node2);

  /**
   * Refine the boundary by inserting nodes when the spacing is too large.
   */
  void refineBoundary();

  /**
   * Find the growth direction at each active boundary node.
   */
  void findActiveBoundaryDirection();

  /**
   * Grow the cutter mesh according to the active growth directions and magnitudes.
   */
  void growFront();

  /**
   * Sort the front nodes to maintain a consistent ordering.
   */
  void sortFrontNodes();

  /**
   * Find intersections between the front and the structural mesh.
   */
  void findFrontIntersection();

  /**
   * Refine the cutter mesh near the crack front.
   */
  void refineFront();

  /**
   * Create TRI3 elements between the new front and the old front.
   */
  void triangulation();

  /**
    Join active boundaries and inactive boundaries to be the new boundary
   */
  void joinBoundary();

  /**
    Parsed functions of front growth
   */
  const Function * _func_x;
  const Function * _func_y;
  const Function * _func_z;
  const Function * _func_v;
};
