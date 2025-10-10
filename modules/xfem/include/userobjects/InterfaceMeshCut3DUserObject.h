//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceMeshCutUserObjectBase.h"

/**
 * Mesh cutter for 3D material interface problems.
 */

class XFEMMovingInterfaceVelocityBase;

class InterfaceMeshCut3DUserObject : public InterfaceMeshCutUserObjectBase
{
public:
  /**
   * Build parameters defining the three-dimensional interface mesh cutter.
   *
   * @return Input parameters describing the 3D interface cutting behavior.
   */
  static InputParameters validParams();

  InterfaceMeshCut3DUserObject(const InputParameters & parameters);

  /**
   * Cut an element using the geometric definition of the 3D interface.
   *
   * @param elem The element being inspected for cuts.
   * @param cut_edges Output vector describing edges intersected by the interface.
   * @param cut_nodes Output vector describing nodes intersected by the interface.
   * @return True when the element is cut by the interface geometry.
   */
  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutEdge> & cut_edges,
                                    std::vector<Xfem::CutNode> & cut_nodes) const override;

  /**
   * Cut an element to identify face intersections with the 3D interface.
   *
   * @param elem Element inspected for face intersections.
   * @param cut_faces Output vector describing faces intersected by the interface.
   * @return True when the element contains faces intersected by the interface.
   */
  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutFace> & cut_faces) const override;

  /**
   * Cut a fragment representation using geometric edge information.
   *
   * @param frag_edges Fragment edge geometry that may be intersected.
   * @param cut_edges Output vector that receives cut edge information.
   * @return True when the fragment is intersected by the interface.
   */
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_edges,
                                     std::vector<Xfem::CutEdge> & cut_edges) const override;

  /**
   * Cut a fragment representation using geometric face information.
   *
   * @param frag_faces Fragment face geometry that may be intersected.
   * @param cut_faces Output vector that receives cut face information.
   * @return True when the fragment contains intersected faces.
   */
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_faces,
                                     std::vector<Xfem::CutFace> & cut_faces) const override;

  /**
   * Compute the signed distance to the interface for a point in the 3D setting.
   *
   * @param p Point for which the signed distance is requested.
   * @return Signed distance relative to the 3D interface surface.
   */
  virtual Real calculateSignedDistance(Point p) const override;

  /**
   * Return the interface normal at a node in the cutting mesh.
   *
   * @param node_id Identifier of the node whose normal is requested.
   * @return Normal vector at the specified node.
   */
  virtual Point nodeNormal(const unsigned int & node_id) override;

  /**
   * Compute the interface normals for all elements of the cutting mesh.
   */
  virtual void calculateNormals() override;

protected:
  /// Map of information defining cutting elements, stored in this order for each element:
  /// pseudo normal, three nodes, and three sides
  std::unordered_map<unsigned int, std::array<Point, 7>> _pseudo_normal;
};
