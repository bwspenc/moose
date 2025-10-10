//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "XFEM.h"

/**
 * Coupled auxiliary value
 */
class XFEMCutPlaneAux : public AuxKernel
{
public:
  /**
   * Create the parameter set describing the XFEM cut plane auxiliary quantity.
   *
   * @return Input parameters required to configure the cut plane auxiliary kernel.
   */
  static InputParameters validParams();

  XFEMCutPlaneAux(const InputParameters & parameters);

  virtual ~XFEMCutPlaneAux() {}

protected:
  /**
   * Evaluate the requested cut plane quantity for the current element.
   *
   * @return The scalar cut-plane value computed at the current quadrature point.
   */
  virtual Real computeValue();

private:
  Xfem::XFEM_CUTPLANE_QUANTITY _quantity;
  std::shared_ptr<XFEM> _xfem;
  unsigned int _plane_id;
};
