/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELASTICSPDMATERIAL_H
#define ELASTICSPDMATERIAL_H

#include "MechanicPDMaterial.h"

class ElasticSPDMaterial;

template<>
InputParameters validParams<ElasticSPDMaterial>();

class ElasticSPDMaterial : public MechanicPDMaterial
{
public:
  ElasticSPDMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpForce();

  MaterialProperty<Real> & _bond_force_ij;
  MaterialProperty<Real> & _bond_force_i_j;
  MaterialProperty<Real> & _bond_dfdU_ij;
  MaterialProperty<Real> & _bond_dfdU_i_j;
  MaterialProperty<Real> & _bond_dfdE_ij;
  MaterialProperty<Real> & _bond_dfdE_i_j;
  MaterialProperty<Real> & _bond_dfdT_ij;
  MaterialProperty<Real> & _bond_dfdT_i_j;

  double _a;
  double _b;
  double _d_i;
  double _d_j;
};

#endif //ELASTICSPDMATERIAL_H
