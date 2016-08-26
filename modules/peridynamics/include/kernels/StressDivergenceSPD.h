/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                         Peridynamics                         */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STRESSDIVERGENCESPD_H
#define STRESSDIVERGENCESPD_H

#include "Kernel.h"
#include "PeridynamicMesh.h"

class StressDivergenceSPD;

template<>
InputParameters validParams<StressDivergenceSPD>();

class StressDivergenceSPD : public Kernel
{
public:
  StressDivergenceSPD(const InputParameters & parameters);
  virtual ~StressDivergenceSPD() {}

protected:
  virtual void initialSetup();

  virtual void computeResidual();
  virtual Real computeQpResidual() {return 0;}
  virtual void computeJacobian();
  virtual void computePartialJacobian();
  virtual void computeFullJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);
  virtual void computePartialOffDiagJacobian(unsigned int jvar);
  virtual void computeFullOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _bond_force_ij;
  const MaterialProperty<Real> & _bond_force_i_j;
  const MaterialProperty<Real> & _bond_dfdU_ij;
  const MaterialProperty<Real> & _bond_dfdU_i_j;
  const MaterialProperty<Real> & _bond_dfdE_ij;
  const MaterialProperty<Real> & _bond_dfdE_i_j;
  const MaterialProperty<Real> & _bond_dfdT_ij;
  const MaterialProperty<Real> & _bond_dfdT_i_j;

private:
  AuxiliarySystem & _aux;

  NumericVector<Number> & _aux_sln;

  NonlinearSystem & _nsys;

  const unsigned int _component;

  unsigned int _ndisp;

  std::vector<unsigned int> _disp_var;

  const bool _temp_coupled;
  const unsigned int _temp_var;

  const bool _strain_zz_coupled;
  const unsigned int _strain_zz_var;

  MooseVariable * _bond_status_var;

  PeridynamicMesh & _pdmesh;
  const unsigned int _pddim;

  const std::vector<RealGradient> * _orientation;
};

#endif //STRESSDIVERGENCESPD_H
