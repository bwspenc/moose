/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STRESSDIVERGENCEPD_H
#define STRESSDIVERGENCEPD_H

#include "Kernel.h"

class StressDivergencePD;

template<>
InputParameters validParams<StressDivergencePD>();

class StressDivergencePD : public Kernel
{
public:
  StressDivergencePD(const InputParameters & parameters);
  virtual ~StressDivergencePD() {}

protected:
  virtual void initialSetup();

  virtual void computeResidual();
  virtual Real computeQpResidual() {return 0;}
  virtual void computeJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);

  virtual void computeStiffness(DenseVector<Real> & stiff_elem);
  virtual void computeOffDiagStiffness(DenseMatrix<Real> & off_stiff_elem);

  const MaterialProperty<Real> & _bond_force;
  const MaterialProperty<Real> & _bond_dfdU;
  const MaterialProperty<Real> & _bond_dfdT;

private:
  const unsigned int _component;

  unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

  const bool _temp_coupled;
  const unsigned int _temp_var;

  const VariableValue & _bond_status;

  const std::vector<RealGradient> * _orientation;

};

#endif //STRESSDIVERGENCEPD_H
