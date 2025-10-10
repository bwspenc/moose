//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementUserObject.h"

class XFEM;

class MeshCut2DNucleationBase : public ElementUserObject
{
public:
  /**
   * Build parameters describing the base 2D mesh cut nucleation user object.
   *
   * @return Input parameters configuring crack nucleation for XFEM.
   */
  static InputParameters validParams();

  MeshCut2DNucleationBase(const InputParameters & parameters);

  /**
   * Initialize state prior to evaluating crack nucleation criteria for the current execution step.
   */
  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

  /**
   * Provide a map of nucleated cracks for consumption by the mesh-cut user object.
   *
   * @return Map keyed by element id whose values contain the crack-tip node pair that defines the
   * nucleated crack segment on the XFEM cutter mesh.
   */
  std::map<unsigned int, std::pair<RealVectorValue, RealVectorValue>> getNucleatedElemsMap() const
  {
    return _nucleated_elems;
  }

  /**
   * Provide access to the nucleation radius specified in the input file.
   *
   * @return Nucleation exclusion radius configured for the user object.
   */
  Real getNucleationRadius() const { return _nucleation_radius; }

protected:
  /**
   * Determine whether the current element should be cut by a new crack.
   *
   * @param cutterElemNodes Nodes of the line segment that will be used to create the cutter mesh.
   * @return True if the element nucleates a crack.
   */
  virtual bool doesElementCrack(std::pair<RealVectorValue, RealVectorValue> & cutterElemNodes) = 0;

private:
  /// The FE solution mesh
  MooseMesh & _mesh;
  // New cracks can only be nucleated if they are at least this far from another crack
  Real _nucleation_radius;
  /// shared pointer to XFEM
  std::shared_ptr<XFEM> _xfem;
  // Boundaries where cracks can nucleate
  std::vector<BoundaryID> _initiation_boundary_ids;
  // map with key for element id and value is a pair containing the node points for creating a
  // nucleated element on the xfem cutter mesh.
  std::map<unsigned int, std::pair<RealVectorValue, RealVectorValue>> _nucleated_elems;
};
