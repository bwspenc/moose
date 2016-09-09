/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                        Peridynamics                          */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef GENERALIZEDPLANESTRAINPDDOF_H
#define GENERALIZEDPLANESTRAINPDDOF_H

#include "Kernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "PeridynamicMesh.h"

class GeneralizedPlaneStrainPDDOF;

template<>
InputParameters validParams<GeneralizedPlaneStrainPDDOF>();

class GeneralizedPlaneStrainPDDOF : public Kernel
{
public:
  GeneralizedPlaneStrainPDDOF(const InputParameters & parameters);

protected:
  virtual void computeResidual();
  virtual Real computeQpResidual() {return 0;};
  virtual void computeJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<RankFourTensor> & _Cijkl;
  const MaterialProperty<RankTwoTensor> & _shape_tensor;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _stress;

private:
  AuxiliarySystem & _aux;

  const NumericVector<Number> & _aux_sln;

  unsigned int _ndisp;

  const bool _temp_coupled;

  std::vector<unsigned int> _disp_var;

  const unsigned int _temp_var;

  MooseVariable * _bond_status_var;

  PeridynamicMesh & _pdmesh;
};
#endif //GENERALIZEDPLANESTRAINPDDOF_H
