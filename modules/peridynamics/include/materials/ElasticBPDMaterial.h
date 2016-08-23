/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELASTICBPDMATERIAL_H
#define ELASTICBPDMATERIAL_H

#include "MechanicPDMaterial.h"

class ElasticBPDMaterial;

template<>
InputParameters validParams<ElasticBPDMaterial>();

class ElasticBPDMaterial : public MechanicPDMaterial
{
public:
  ElasticBPDMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpForce();

  MaterialProperty<Real> & _bond_force;
  MaterialProperty<Real> & _bond_dfdU;
  MaterialProperty<Real> & _bond_dfdE;
  MaterialProperty<Real> & _bond_dfdT;
};

#endif //ELASTICBPDMATERIAL_H
