/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef MECHANICPDMATERIAL_H
#define MECHANICPDMATERIAL_H

#include "PeridynamicMaterial.h"

class MechanicPDMaterial;

template<>
InputParameters validParams<MechanicPDMaterial>();

class MechanicPDMaterial : public PeridynamicMaterial
{
public:
  MechanicPDMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  virtual void computeQpStrain() = 0;
  virtual void computeQpForce() = 0;

  MaterialProperty<Real> & _bond_force;
  MaterialProperty<Real> & _bond_dfdU;
  MaterialProperty<Real> & _bond_dfdT;
  MaterialProperty<Real> & _bond_total_strain;
  MaterialProperty<Real> & _bond_mechanic_strain;
  MaterialProperty<Real> & _bond_elastic_strain;

  const bool _plane_strain;
  const Real _youngs_modulus;
  const Real _poissons_ratio;

  const VariableValue & _temp;
  const Real _temp_ref;
  const Real _thermal_expansion_coeff;

  std::vector<MooseVariable *> _disp_var;

  double _alpha;
  double _current_length;
};

#endif //MECHANICPDMATERIAL_H
