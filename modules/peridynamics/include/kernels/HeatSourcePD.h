/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef HEATSOURCEPD_H
#define HEATSOURCEPD_H

#include "Kernel.h"
#include "PeridynamicMesh.h"

class HeatSourcePD;

template<>
InputParameters validParams<HeatSourcePD>();

class HeatSourcePD : public Kernel
{
public:

  HeatSourcePD(const InputParameters & parameters);

protected:
  virtual void computeResidual();
  virtual Real computeQpResidual() {return 0;}

  PeridynamicMesh & _pdmesh;

  double _power_density;
  Function * _power_density_function;
};

#endif //HEATSOURCEPD_H
