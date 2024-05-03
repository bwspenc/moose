//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SurfaceTractions.h"

// MOOSE includes
#include "Material.h"

registerMooseObject("SolidMechanicsApp", SurfaceTractions);
//registerMooseObject("SolidMechanicsApp", ADSurfaceTractions); //TODO!!!

InputParameters
SurfaceTractions::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  return params;
}

SurfaceTractions::SurfaceTractions(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),//BWS TODO
    _stress(getMaterialProperty<RankTwoTensor>("stress")),
    _normal(declareProperty<RealVectorValue>("normal")),
    _tractions(declareProperty<RealVectorValue>("tractions")),
    _normal_traction(declareProperty<Real>("normal_traction")),
    _max_shear_traction(declareProperty<Real>("max_shear_traction")),
    _min_shear_traction(declareProperty<Real>("min_shear_traction")),
    _normals(_assembly.normals())
{
}

void
SurfaceTractions::computeQpProperties()
{
  //std::cout<<"Normals: "<<_normals[_qp]<<std::endl;
  _normal[_qp] = _normals[_qp];
  _tractions[_qp] = _stress[_qp] * _normal[_qp];
  _normal_traction[_qp] = _tractions[_qp] * _normal[_qp];
}
