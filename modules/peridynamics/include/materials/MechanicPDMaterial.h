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
  virtual void computeQpStrain();
  virtual void computeQpForce() = 0;

  virtual void computeNodalTemp();
  virtual void computeElasticStrainTensor();
  virtual void computeStressTensor();
  virtual Real computeBondCurrentLength();

  MaterialProperty<Real> & _bond_elastic_strain;
  MaterialProperty<RankFourTensor> & _elasticity_tensor;
  MaterialProperty<Real> & _thermal_expansion;
  MaterialProperty<RankTwoTensor> & _shape_tensor;
  MaterialProperty<RankTwoTensor> & _deformation_gradient;
  MaterialProperty<RankTwoTensor> & _elastic_strain;
  MaterialProperty<RankTwoTensor> & _strain;
  MaterialProperty<RankTwoTensor> & _stress;

  const Real _youngs_modulus;
  const Real _poissons_ratio;

  MooseVariable * _strain_zz_var;
  MooseVariable * _temp_var;
  const Real _temp_ref;

  std::vector<MooseVariable *> _disp_var;

  double _alpha;
  double _shear_modulus;
  double _bulk_modulus;

  RankFourTensor _Cijkl;

  double _strain_zz_i;
  double _strain_zz_j;
};

#endif //MECHANICPDMATERIAL_H
