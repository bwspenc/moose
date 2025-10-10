//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshCut2DUserObjectBase.h"

class CrackFrontDefinition;

/**
 * MeshCut2DFractureUserObject:
 * (1) reads in a mesh describing the crack surface
 * (2) uses the mesh to do initial cutting of 2D elements, and
 * (3) grows the mesh by a fixed growth rate when a fracture-integral-based growth criterion is met.
 */

class MeshCut2DFractureUserObject : public MeshCut2DUserObjectBase
{
public:
  /**
   * Build parameters describing the fracture-driven mesh cut user object.
   *
   * @return Input parameters defining fracture-controlled crack growth.
   */
  static InputParameters validParams();

  MeshCut2DFractureUserObject(const InputParameters & parameters);

  /**
   * Perform initial setup specific to fracture-controlled crack growth.
   */
  virtual void initialSetup() override;

  /**
   * Initialize state prior to applying fracture growth logic for the current step.
   */
  virtual void initialize() override;

protected:
  /**
   * Determine active boundary growth based on fracture integral criteria.
   */
  virtual void findActiveBoundaryGrowth() override;

private:
  /// critical k value for crack growth
  const Real & _k_critical;
  /// amount to grow crack by for each xfem update step
  const Real & _growth_increment;

  const std::string _ring_number_string;

  CrackFrontDefinition * _crack_front_definition;
  /**
   * Compute the squared fracture integral magnitude for each crack tip.
   *
   * @param k1 Fracture integrals from mode-I vector postprocessors.
   * @param k2 Fracture integrals from mode-II vector postprocessors.
   * @return Vector of squared fracture integral magnitudes used for crack growth decisions.
   */
  std::vector<Real> getKSquared(const std::vector<Real> & k1, const std::vector<Real> & k2) const;
};
