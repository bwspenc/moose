//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"

class GeneralizedPlaneStrainAction : public Action
{
public:
  static InputParameters validParams();

  GeneralizedPlaneStrainAction(const InputParameters & params);

  void act() override;

protected:
  /// Return the first displacement component in the plane
  unsigned int firstInPlaneDisplacementIndex() const;

  /// Maps the Action's deprecated 'out_of_plane_pressure'/'factor' params onto the kernel's
  /// current 'out_of_plane_pressure_function'/'pressure_factor' params, since the kernel only
  /// carries the current parameter variants. Shared by both the AD and non-AD add_kernel branches.
  void remapDeprecatedPressureParams(InputParameters & params) const;

  std::vector<VariableName> _displacements;

  /// Number of displacement variables
  unsigned int _ndisp;
  const unsigned int _out_of_plane_direction;
  const bool _use_ad;
};
