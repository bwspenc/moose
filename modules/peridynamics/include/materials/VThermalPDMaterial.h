/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef VTHERMALPDMATERIAL_H
#define VTHERMALPDMATERIAL_H

#include "PeridynamicMaterial.h"

class VThermalPDMaterial;
class Function;

template<>
InputParameters validParams<VThermalPDMaterial>();

class VThermalPDMaterial : public PeridynamicMaterial
{
public:
  VThermalPDMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual double computeBondModulus(double temp_avg);
  virtual void computeQpProperties();

  const Real _thermal_conductivity;
  Function * _thermal_conductivity_function;

  MaterialProperty<Real> & _bond_response;
  MaterialProperty<Real> & _bond_drdT;

  MooseVariable * _temp_var;
};

#endif //VTHERMALPDMATERIAL_H
