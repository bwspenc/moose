//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LevelSetBiMaterialBase.h"

/**
 * Compute a Real material property for bi-materials problem (consisting of two
 * different materials) defined by a level set function
 */
template <bool is_ad>
class LevelSetBiMaterialRealTempl : public LevelSetBiMaterialBase
{
public:
  /**
   * Build parameters describing the bi-material real-valued material.
   *
   * @return Input parameters defining the level-set-controlled scalar properties.
   */
  static InputParameters validParams();

  LevelSetBiMaterialRealTempl(const InputParameters & parameters);

protected:
  /**
   * Assign quadrature point properties for the positive level-set region.
   */
  virtual void assignQpPropertiesForLevelSetPositive() override;

  /**
   * Assign quadrature point properties for the negative level-set region.
   */
  virtual void assignQpPropertiesForLevelSetNegative() override;

  /// Real Material properties for the two separate materials in the bi-material system
  std::vector<const GenericMaterialProperty<Real, is_ad> *> _bimaterial_material_prop;

  /// Global Real material property (switch bi-material diffusion coefficient based on level set values)
  GenericMaterialProperty<Real, is_ad> & _material_prop;
};

typedef LevelSetBiMaterialRealTempl<false> LevelSetBiMaterialReal;
typedef LevelSetBiMaterialRealTempl<true> ADLevelSetBiMaterialReal;
