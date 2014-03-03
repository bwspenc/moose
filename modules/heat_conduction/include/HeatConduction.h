/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef HEATCONDUCTIONKERNEL_H
#define HEATCONDUCTIONKERNEL_H

#include "Diffusion.h"
#include "Material.h"

namespace libMesh
{
class XFEM;
}

//Forward Declarations
class HeatConductionKernel;

template<>
InputParameters validParams<HeatConductionKernel>();

/**
 * Note: This class is named HeatConductionKernel instead of HeatConduction
 * to avoid a clash with the HeatConduction namespace.  It is registered
 * as HeatConduction, which means it can be used by that name in the input
 * file.
 */
class HeatConductionKernel : public Diffusion
{
public:

  HeatConductionKernel(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const unsigned _dim;
  const MaterialProperty<Real> & _diffusion_coefficient;
  const MaterialProperty<Real> * const _diffusion_coefficient_dT;
  bool _has_xfem_volfrac;
  VariableValue & _xfem_volfrac;
};

#endif //HEATCONDUCTIONKERNEL_H
