#ifndef TRUSSMATERIAL_H
#define TRUSSMATERIAL_H

#include "Material.h"

//Forward Declarations
class TrussMaterial;

template<>
InputParameters validParams<TrussMaterial>();

/**
 * LinearIsotropic material for use in simple applications that don't need material properties.
 */
class TrussMaterial : public Material
{
public:
  TrussMaterial(const std::string & name,
                InputParameters parameters);

  virtual ~TrussMaterial();

protected:

  VariableValue & _disp_x_nodal_val;
  VariableValue & _disp_y_nodal_val;
  VariableValue & _disp_z_nodal_val;

  virtual void computeProperties();

  MaterialProperty<Real> & _axial_stress;
  MaterialProperty<Real> & _e_over_l;
  Real _youngs_modulus;
  bool _youngs_modulus_coupled;
  VariableValue & _youngs_modulus_var;

  bool _has_temp;
  VariableValue & _temp;
  Real _t_ref;
  Real _alpha;

  unsigned int _dim;
};

#endif //TRUSSMATERIAL_H
