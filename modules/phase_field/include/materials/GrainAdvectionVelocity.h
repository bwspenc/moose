/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef GRAINADVECTIONVELOCITY_H
#define GRAINADVECTIONVELOCITY_H

#include "Material.h"
#include "GrainTrackerInterface.h"
#include "GrainForceAndTorqueInterface.h"
#include "DerivativeMaterialInterface.h"

//Forward Declarations
class GrainAdvectionVelocity;

template<>
InputParameters validParams<GrainAdvectionVelocity>();

/**
 * This Material calculates the advection velocity, it's divergence and
 * derivatives acting on a particle/grain
 */
class GrainAdvectionVelocity : public DerivativeMaterialInterface<Material>
{
public:
  GrainAdvectionVelocity(const InputParameters & parameters);

  virtual void initialSetup();
  virtual void residualSetup();

protected:
  virtual void computeQpProperties();
  /// obtain total no. of grains from GrainTracker
  virtual void getTotalNumberOfGrains();

  /// getting userobject for calculating grain centers and volumes
  const GrainTrackerInterface & _grain_tracker;

  /// getting userobject for calculating grain forces and torques
  const GrainForceAndTorqueInterface & _grain_force_torque;
  const std::vector<RealGradient> & _grain_forces;
  const std::vector<RealGradient> & _grain_torques;

private:
  /// constant value corresponding to grain translation
  Real _mt;

  /// constant value corresponding to grain rotation
  Real _mr;

  unsigned int _grain_num;
  unsigned int _op_num;
  std::vector<const VariableValue *> _vals;
  std::vector<const VariableGradient *> _grad_vals;
  VariableName _c_name;
  /// type of force density material
  std::string _base_name;

  /// Material storing advection velocities of grains
  MaterialProperty<std::vector<RealGradient> > & _velocity_advection;

  /// Material storing divergence of advection velocities of grains
  MaterialProperty<std::vector<Real> > & _div_velocity_advection;
};

#endif //GRAINADVECTIONVELOCITY_H
