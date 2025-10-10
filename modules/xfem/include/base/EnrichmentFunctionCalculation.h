//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CrackFrontDefinition.h"

/**
 * Perform calculation of enrichment function values and derivatives.
 */
class EnrichmentFunctionCalculation
{
public:
  EnrichmentFunctionCalculation(const CrackFrontDefinition * crack_front_definition);

  /**
   * Calculate the enrichment function values at a spatial point.
   *
   * @param point Location where the enrichment function is evaluated.
   * @param B Vector receiving the evaluated enrichment function values.
   * @return Index of the closest crack front used for the enrichment evaluation.
   */
  virtual unsigned int crackTipEnrichementFunctionAtPoint(const Point & point,
                                                          std::vector<Real> & B);

  /**
   * Calculate the enrichment function derivatives at a spatial point.
   *
   * @param point Location where the enrichment function derivatives are evaluated.
   * @param dB Vector receiving the evaluated enrichment function derivatives.
   * @return Index of the closest crack front used for the enrichment evaluation.
   */
  virtual unsigned int
  crackTipEnrichementFunctionDerivativeAtPoint(const Point & point,
                                               std::vector<RealVectorValue> & dB);

  /**
   * Rotate a vector from the crack-front coordinate system to global coordinates.
   *
   * @param vector Vector expressed in crack-front coordinates.
   * @param rotated_vector Output vector expressed in global coordinates.
   * @param point_index Index of the crack front point defining the local orientation.
   */
  void rotateFromCrackFrontCoordsToGlobal(const RealVectorValue & vector,
                                          RealVectorValue & rotated_vector,
                                          const unsigned int point_index);

private:
  const CrackFrontDefinition & _crack_front_definition;
  Real _r;
  Real _theta;
};
