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
  virtual void computeQpStrain();
  virtual void computeQpForce() = 0;

  virtual void computeNodalTemp();
  virtual void computeElasticStrainTensor();
  virtual void computeStressTensor();
  virtual Real computeBondCurrentLength();

  MaterialProperty<Real> & _bond_elastic_strain;
  MaterialProperty<Real> & _bond_critical_strain;
  MaterialProperty<Real> & _bond_critical_strain_old;
  MaterialProperty<RankFourTensor> & _elasticity_tensor;
  MaterialProperty<Real> & _thermal_expansion;
  MaterialProperty<RankTwoTensor> & _shape_tensor;
  MaterialProperty<RankTwoTensor> & _deformation_gradient;
  MaterialProperty<RankTwoTensor> & _elastic_strain;
  MaterialProperty<RankTwoTensor> & _strain;
  MaterialProperty<RankTwoTensor> & _stress;
  MaterialProperty<Real> & _bond_force_ij;
  MaterialProperty<Real> & _bond_force_i_j;
  MaterialProperty<Real> & _bond_dfdU_ij;
  MaterialProperty<Real> & _bond_dfdU_i_j;
  MaterialProperty<Real> & _bond_dfdE_ij;
  MaterialProperty<Real> & _bond_dfdE_i_j;
  MaterialProperty<Real> & _bond_dfdT_ij;
  MaterialProperty<Real> & _bond_dfdT_i_j;

  const Real _youngs_modulus;
  const Real _poissons_ratio;

  bool _strain_zz_coupled;
  VariableValue & _strain_zz;

  MooseVariable * _temp_var;
  MooseVariable * _bond_status_var;

  const Real _temp_ref;

  std::vector<MooseVariable *> _disp_var;

  double _alpha;
  double _shear_modulus;
  double _bulk_modulus;

  RankFourTensor _Cijkl;
};

#endif //MECHANICPDMATERIAL_H
