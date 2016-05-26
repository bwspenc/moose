/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef VLEPDMATERIAL_H
#define VLEPDMATERIAL_H

#include "MechanicPDMaterial.h"

class VLEPDMaterial;

template<>
InputParameters validParams<VLEPDMaterial>();

class VLEPDMaterial : public MechanicPDMaterial
{
public:
  VLEPDMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpStrain();
  virtual void computeQpForce();

  virtual double computeLambda2D(double poissons_ratio);
  virtual double computeLambda3D(double poissons_ratio);
};

#endif //VLEPDMATERIAL_H
