/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef THERMALPDMATERIAL_H
#define THERMALPDMATERIAL_H

#include "PeridynamicMaterial.h"

class ThermalPDMaterial;
class Function;

template<>
InputParameters validParams<ThermalPDMaterial>();

class ThermalPDMaterial : public PeridynamicMaterial
{
public:
  ThermalPDMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  const Real _thermal_conductivity;
  Function * _thermal_conductivity_function;

  MaterialProperty<Real> & _bond_response;
  MaterialProperty<Real> & _bond_drdT;

  MooseVariable * _temp_var;
};

#endif //THERMALPDMATERIAL_H
