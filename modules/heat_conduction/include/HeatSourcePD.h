/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef HEATSOURCEPD_H
#define HEATSOURCEPD_H

#include "Kernel.h"

//Forward Declaration
class HeatSourcePD;

template<>
InputParameters validParams<HeatSourcePD>();

class HeatSourcePD : public Kernel
{
public:

  HeatSourcePD(const InputParameters & parameters);
  virtual ~HeatSourcePD();

protected:
  virtual void computeResidual();

  virtual Real computeQpResidual() {return 0;}

  Real _power_density;

  Function * _power_density_function;

  const MaterialProperty<Real> & _bond_volume;

};

#endif //HEATSOURCEPD_H
