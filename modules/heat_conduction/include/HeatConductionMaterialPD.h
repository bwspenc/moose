/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef HEATCONDUCTIONMATERIALPD_H
#define HEATCONDUCTIONMATERIALPD_H

#include "Material.h"

//Forward Declarations
class HeatConductionMaterialPD;
class Function;

template<>
InputParameters validParams<HeatConductionMaterialPD>();

class HeatConductionMaterialPD : public Material
{
public:
  HeatConductionMaterialPD(const InputParameters & parameters);
  virtual ~HeatConductionMaterialPD();

protected:

  virtual void computeProperties();

  MooseVariable * _temp_var;
  
  const Real _my_thermal_conductivity;

  MaterialProperty<Real> & _thermal_conductivity;
  Function * _thermal_conductivity_function;

  const int _pddim;
  const Real _mesh_spacing;

  MaterialProperty<Real> & _bond_response;
  MaterialProperty<Real> & _bond_response_dif_temp;
<<<<<<< e48bdee7d4a87c5e6fb10593fa54979f06371c99
  MaterialProperty<Real> & _node_volume;
=======
  MaterialProperty<Real> & _bond_volume;

>>>>>>> completed the heat conduction implementation and the coupled thermomechanical model, wrote a aux kernel to update the bond status at the end of each time step.
};

#endif //HEATCONDUCTIONMATERIALPD_H
