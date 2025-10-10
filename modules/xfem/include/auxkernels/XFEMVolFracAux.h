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

class XFEM;

/**
 * Coupled auxiliary value
 */
class XFEMVolFracAux : public AuxKernel
{
public:
  /**
   * Build the parameter set for constructing an XFEM volume fraction auxiliary kernel.
   *
   * @return Input parameters describing the configuration for the XFEM volume fraction kernel.
   */
  static InputParameters validParams();

  XFEMVolFracAux(const InputParameters & parameters);

  virtual ~XFEMVolFracAux() {}

protected:
  /**
   * Compute the auxiliary value reporting the physical volume fraction.
   *
   * @return The computed physical volume fraction for the current quadrature point.
   */
  virtual Real computeValue();

private:
  std::shared_ptr<XFEM> _xfem;
};
