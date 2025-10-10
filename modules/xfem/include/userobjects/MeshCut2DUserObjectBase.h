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
class MeshCut2DNucleationBase;
/**
 * MeshCut2DUserObjectBase: (1) reads in a mesh describing the crack surface,
 * (2) Fills xfem cut element ojbects.
 * Derived classes modify the class and grow the mesh
 */

class MeshCut2DUserObjectBase : public GeometricCutUserObject
{
public:
  /**
   * Build parameters describing the base 2D mesh cut user object.
   *
   * @return Input parameters defining the base mesh-cut configuration.
   */
  static InputParameters validParams();

  MeshCut2DUserObjectBase(const InputParameters & parameters);

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
   * Retrieve the ordered crack front points managed by the base user object.
   *
   * @param num_crack_front_points Number of crack front points requested by the caller.
   * @return Vector of crack front points suitable for downstream consumers.
   */
  virtual const std::vector<Point>
  getCrackFrontPoints(unsigned int num_crack_front_points) const override;

  /**
   * Retrieve crack-front normals from the XFEM geometric cut user object.
   *
   * CrackFrontDefinition requests normals, so this implementation provides normals for line
   * elements with a tangent direction in the [001] orientation.
   *
   * @param num_crack_front_points Number of crack front normals requested by the caller.
   * @return Vector containing the crack front normals.
   */
  virtual const std::vector<RealVectorValue>
  getCrackPlaneNormals(unsigned int num_crack_front_points) const override;

  MeshBase & getCutterMesh() const;

protected:
  /// The FE solution mesh
  MooseMesh & _mesh;

  /// The xfem cutter mesh
  std::unique_ptr<MeshBase> _cutter_mesh;

  /// 2D UO for nucleating cracks
  const MeshCut2DNucleationBase * _nucleate_uo;

  /// Indicator that shows if the cutting mesh is modified or not in this calculation step
  bool _is_mesh_modified;

  /**
   * This vector of pairs orders crack tips to make the order used in this class the same as those
   * for the fracture integrals vectorpostProcessorsin created by CrackFrontDefinition.
   * The original crack front node ids found in the cutter mesh are put in pair.first and the
   * assosciated current crack front node id that grew from the original crack front node id is in
   * pair.second.  This vector is sorted on pair.first which makes the ordering of this vector the
   * same as that used in the CrackFrontDefinition
   */
  std::vector<std::pair<dof_id_type, dof_id_type>> _original_and_current_front_node_ids;

  /// contains the active node ids and their growth vectors
  std::vector<std::pair<dof_id_type, Point>> _active_front_node_growth_vectors;

  /**
   * Find the growth direction at each active node.
   */
  virtual void findActiveBoundaryGrowth() = 0;

  /**
   * Find the original crack front nodes in the cutter mesh and use to populate
   * _original_and_current_front_node_ids.
   */
  void findOriginalCrackFrontNodes();

  /// grow the cutter mesh
  void growFront();

  /**
   * Calls into MeshCutNucleation UO to add cracks.
   */
  void addNucleatedCracksToMesh();

private:
  /**
   * Remove nucleated cracks that are too close to each other; the lowest map key wins.
   *
   * @param nucleated_elems_map Map from the nucleation user object keyed by element id with the
   * pair of crack-tip nodes.
   * @param nucleationRadius Minimum exclusion distance permitted between nucleated cracks.
   */
  void removeNucleatedCracksTooCloseToEachOther(
      std::map<unsigned int, std::pair<RealVectorValue, RealVectorValue>> & nucleated_elems_map,
      Real nucleationRadius);
  /**
   * Remove nucleated cracks that are too close to a pre-existing crack in the mesh.
   *
   * @param nucleated_elems_map Map from the nucleation user object keyed by element id with the
   * pair of crack-tip nodes.
   * @param nucleationRadius Minimum exclusion distance permitted between nucleated cracks and the
   * existing mesh.
   */
  void removeNucleatedCracksTooCloseToExistingCracks(
      std::map<unsigned int, std::pair<RealVectorValue, RealVectorValue>> & nucleated_elems_map,
      Real nucleationRadius);
};
