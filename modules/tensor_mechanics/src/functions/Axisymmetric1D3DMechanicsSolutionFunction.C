//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseError.h"
#include "Axisymmetric1D3DMechanicsSolutionFunction.h"
#include "SolutionUserObject.h"

registerMooseObject("MooseApp", Axisymmetric1D3DMechanicsSolutionFunction);

template <>
InputParameters
validParams<Axisymmetric1D3DMechanicsSolutionFunction>()
{
  // Get the Function input parameters
  InputParameters params = validParams<Function>();

  // Add parameters specific to this object
  params.addRequiredParam<UserObjectName>("solution",
                                          "The SolutionUserObject to extract data from.");
  params.addRequiredParam<std::string>(
      "disp_x",
      "The name of the disp_x variable to be extracted from the file");
  params.addRequiredParam<std::string>(
      "strain_zz",
      "The name of the strain_zz variable to be extracted from the file");
  params.addParam<Real>(
      "scale_factor",
      1.0,
      "Scale factor (a) to be applied to the solution (x): ax+b, where b is the 'add_factor'");
  params.addParam<Real>(
      "add_factor",
      0.0,
      "Add this value (b) to the solution (x): ax+b, where a is the 'scale_factor'");

  params.addParam<RealVectorValue>("3d_axis_point1",
                                   RealVectorValue(0, 0, 0),
                                   "Start point for axis of symmetry for the 3d model");
  params.addParam<RealVectorValue>("3d_axis_point2",
                                   RealVectorValue(0, 1, 0),
                                   "End point for axis of symmetry for the 3d model");

  params.addRequiredParam<unsigned int>("component",
                                        "Component of the displacement variable to be computed");

  params.addClassDescription("Function for reading a 1D axisymmetric mechanics solution from file and "
                             "mapping it to a 3D Cartesian mechanics model");

  return params;
}

Axisymmetric1D3DMechanicsSolutionFunction::Axisymmetric1D3DMechanicsSolutionFunction(
    const InputParameters & parameters)
  : Function(parameters),
    _solution_object_ptr(NULL),
    _scale_factor(getParam<Real>("scale_factor")),
    _add_factor(getParam<Real>("add_factor")),
    _3d_axis_point1(getParam<RealVectorValue>("3d_axis_point1")),
    _3d_axis_point2(getParam<RealVectorValue>("3d_axis_point2")),
    _component(getParam<unsigned int>("component")),
    _disp_x_var_name(getParam<std::string>("disp_x")),
    _strain_zz_var_name(getParam<std::string>("strain_zz"))
{
  Point zero;
  Point unit_vec_y;
  unit_vec_y(1) = 1;
  if (_3d_axis_point1 == zero && _3d_axis_point2 == unit_vec_y)
    _default_axes = true;
  else
    _default_axes = false;

  if (_3d_axis_point1.relative_fuzzy_equals(_3d_axis_point2))
    mooseError("3d_axis_point1 and 3d_axis_point2 must be different points");
}

void
Axisymmetric1D3DMechanicsSolutionFunction::initialSetup()
{
  // Get a pointer to the SolutionUserObject. A pointer is used because the UserObject is not
  // available during the
  // construction of the function
  _solution_object_ptr = &getUserObject<SolutionUserObject>("solution");

  _solution_object_var_indices.resize(2);
  _solution_object_var_indices[0] = _solution_object_ptr->getLocalVarIndex(_disp_x_var_name);
  _solution_object_var_indices[1] = _solution_object_ptr->getLocalVarIndex(_strain_zz_var_name);
}

Real
Axisymmetric1D3DMechanicsSolutionFunction::value(Real t, const Point & p) const
{
  Point r_dir_axisym; // direction of the radial coordinate in the axisymmetric model
  Point z_dir_axisym; // direction of the axial coordinate in the axisymmetric model
  Point r_dir_3d;     // direction of the radial coordinate in the 3D model
  Point z_dir_3d;     // direction of the axial coordinate in the 3D model
  bool r_gt_zero = false;

  r_dir_axisym(0) = 1;
  z_dir_axisym(1) = 1;

  Real r;
  Real z;

  if (_default_axes)
  {
    r_dir_3d = p;
    r_dir_3d(1) = 0;
    r = r_dir_3d.norm();
    if (MooseUtils::absoluteFuzzyGreaterThan(r, 0.0))
    {
      r_gt_zero = true;
      r_dir_3d /= r;
    }
    z = p(1);
    z_dir_3d(1) = 1;
  }
  else
  {
    // Find the r, z coordinates of the point in the 3D model relative to the 3D axis
    z_dir_3d = _3d_axis_point2 - _3d_axis_point1;
    z_dir_3d /= z_dir_3d.norm();
    Point v3dp1p(p - _3d_axis_point1);
    z = z_dir_3d * v3dp1p;
    Point axis_proj = _3d_axis_point1 + z * z_dir_3d; // projection of point onto axis
    Point axis_proj_to_p = p - axis_proj;
    r = axis_proj_to_p.norm();
    if (MooseUtils::absoluteFuzzyGreaterThan(r, 0.0))
    {
      r_gt_zero = true;
      r_dir_3d = axis_proj_to_p / r;
    }
  }

  Point axisym_point = r; // point in the axisymmetric model

  // disp_vec_rz contains the displacement vector in RZ coordinates
  Point disp_vec_rz;
  disp_vec_rz(0) = _solution_object_ptr->pointValue(t, axisym_point, _solution_object_var_indices[0]);
  disp_vec_rz(1) = _solution_object_ptr->pointValue(t, axisym_point, _solution_object_var_indices[1]) * z;

  if (!r_gt_zero && !MooseUtils::absoluteFuzzyEqual(disp_vec_rz(0), 0.0))
    mooseError("In Axisymmetric1D3DMechanicsSolutionFunction r=0 and r component of value vector != 0");
  Point disp_vec_3d = disp_vec_rz(0) * r_dir_3d + disp_vec_rz(1) * z_dir_3d;

  Real val = disp_vec_3d(_component);

  return _scale_factor * val + _add_factor;
}
