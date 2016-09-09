/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "GPSPDDiag.h"
#include "Assembly.h"

template<>
InputParameters validParams<GPSPDDiag>()
{
  InputParameters params = validParams<ScalarKernel>();
  params.addClassDescription("Generalized Plane Strain Scalar Kernel");
  params.addRequiredParam<UserObjectName>("gps_uo", "The name of the GPSPDUO UserObject");
  return params;
}

GPSPDDiag::GPSPDDiag(const InputParameters & parameters) :
  ScalarKernel(parameters),
  _gps_uo(getUserObject<GPSPDUO>("gps_uo"))
{
}

void
GPSPDDiag::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  for (_i = 0; _i < re.size(); _i++)
    re(_i) += _gps_uo.returnResidual();
}

/*
  method to provide the diagonal jacobian term for scalar variable using value
  returned from UserObject, off diagonal terms are computed by computeOffDiagJacobianScalar
  in the kernel which needs to couple with the scalar variable
*/
void
GPSPDDiag::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  for (_i = 0; _i < ke.m(); _i++)
    ke(_i, _i) += _gps_uo.returnJacobian();
}
