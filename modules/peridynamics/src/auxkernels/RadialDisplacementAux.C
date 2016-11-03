/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "RadialDisplacementAux.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<RadialDisplacementAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Compute the radial component of the displacement vector for cylindrical or spherical models.");
  params.addRequiredCoupledVar("displacements", "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<RealVectorValue>("axis_point_1", "Start point for line defining cylindrical axis (3D Cartesian models)");
  params.addParam<RealVectorValue>("axis_point_2", "End point for line defining cylindrical axis (3D Cartesian models)");
  params.addParam<RealVectorValue>("origin", "Sphere origin for 3D Cartesian and 2D axisymmetric models, or circle origin for 2D Cartesian models");
  params.set<bool>("use_displaced_mesh") = false;

  return params;
}

RadialDisplacementAux::RadialDisplacementAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _have_axis(false),
  _have_origin(false)
{
  const std::set<SubdomainID> & subdomains = _mesh.meshSubdomains();
  const auto & sbd_begin = *subdomains.begin();
  for (const auto & sbd : subdomains)
  {
    if (sbd == sbd_begin)
      _coord_system = _subproblem.getCoordSystem(sbd);
    else if (_subproblem.getCoordSystem(sbd) != _coord_system)
      mooseError("RadialDisplacementAux requires that all subdomains have the same coordinate type");
  }

  _ndisp = coupledComponents("displacements");
  _disp_vals.resize(_ndisp);

  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_vals[i] = &coupledValue("displacements", i);

//  if (_ndisp != _mesh.dimension())
//    mooseError("The number of displacement variables supplied must match the mesh dimension.");

  if (_coord_system == Moose::COORD_XYZ)
  {
    if (_ndisp == 2)
    {
      if (isParamValid("axis_point_1") || isParamValid("axis_point_2"))
        mooseError("The 'axis_point_1' and 'axis_point_2' parameters are only valid for models " <<
                   "with 3D Cartesian coordinate systems.");
      if (isParamValid("origin"))
      {
        _origin = getParam<RealVectorValue>("origin");
        _have_origin = true;
      }
      else
        mooseError("Must specify 'origin' for models with 2D Cartesian coordinate systems.");
    }
    else if (_ndisp == 3)
    {
      if (isParamValid("axis_point_1") && isParamValid("axis_point_2"))
      {
        _axis_point_1 = getParam<RealVectorValue>("axis_point_1");
        _axis_point_2 = getParam<RealVectorValue>("axis_point_2");
        _have_axis = true;
        const RealVectorValue p1p2(_axis_point_2 - _axis_point_1);
        if (MooseUtils::absoluteFuzzyEqual(p1p2.norm(), 0.0))
          mooseError("Points supplied for 'axis_point_1' and 'axis_point_2' are too close together");
      }
      else if (isParamValid("axis_point_1") || isParamValid("axis_point_2"))
        mooseError("Must specify 'axis_point_1' and 'axis_point_2' together");
      if (isParamValid("origin"))
      {
        if (_have_axis)
          mooseError("Cannot specify both 'origin' and a combination of 'axis_point_1' and 'axis_point_2'");
        _origin = getParam<RealVectorValue>("origin");
        _have_origin = true;
      }
      if (!_have_origin && !_have_axis)
        mooseError("Must specify either 'origin' or a combination of 'axis_point_1' and 'axis_point_2' " <<
                   "for models with 3D Cartesian coordinate systems");
    }
    else
      mooseError("Cannot compute radial displacement for models with 1D Cartesian system");
  }
  else if (_coord_system == Moose::COORD_RZ)
  {
    if (isParamValid("axis_point_1") || isParamValid("axis_point_2"))
      mooseError("The 'axis_point_1' and 'axis_point_2' parameters are only valid for models with " <<
                 "3D Cartesian coordinate systems.");
    if (isParamValid("origin"))
    {
      _origin = getParam<RealVectorValue>("origin");
      _have_origin = true;
    }
  }
  else
  {
    if (isParamValid("axis_point_1") || isParamValid("axis_point_2"))
      mooseError("The 'axis_point_1' and 'axis_point_2' parameters are only valid for models with " <<
                 "3D Cartesian coordinate systems.");
    if (isParamValid("origin"))
      mooseError("The 'origin' parameter is only valid for models with Cartesian coordinate systems.");
  }
}

Real
RadialDisplacementAux::computeValue()
{
  Real rad_disp = 0.0;

  switch (_coord_system)
  {
    case Moose::COORD_XYZ:
    {
      Point current_point (*_current_node);
      if (_ndisp == 2)
      {
        RealVectorValue rad_vec(current_point - _origin);
        Real rad = rad_vec.norm();
        const RealVectorValue disp_vec((*_disp_vals[0])[_qp], (*_disp_vals[1])[_qp]);
        if (rad > 0.0)
        {
          rad_vec /= rad;
          rad_disp = rad_vec * disp_vec;
        }
        else
          rad_disp = disp_vec.norm();
      }
      else if (_ndisp == 3 && _have_origin)
      {
        RealVectorValue rad_vec(current_point - _origin);
        Real rad = rad_vec.norm();
        const RealVectorValue disp_vec((*_disp_vals[0])[_qp], (*_disp_vals[1])[_qp], (*_disp_vals[2])[_qp]);
        if (rad > 0.0)
        {
          rad_vec /= rad;
          rad_disp = rad_vec * disp_vec;
        }
        else
          rad_disp = disp_vec.norm();
      }
      else if (_ndisp == 3 && _have_axis)
      {
        // t is the distance along the axis from point 1 to 2 to the point nearest to the current point.
        const RealVectorValue p1p2(_axis_point_2 - _axis_point_1);
        const RealVectorValue axis_dir = p1p2 / p1p2.norm();
        const RealVectorValue p1pc(current_point - _axis_point_1);
        const Real t = p1pc * axis_dir;

        // The nearest point on the cylindrical axis to current_point is p.
        const RealVectorValue p(_axis_point_1 + t * axis_dir);
        RealVectorValue rad_vec(current_point - p);
        Real rad = rad_vec.norm();
        const RealVectorValue disp_vec((*_disp_vals[0])[_qp], (*_disp_vals[1])[_qp], (*_disp_vals[2])[_qp]);
        if (rad > 0.0)
        {
          rad_vec /= rad;
          rad_disp = rad_vec * disp_vec;
        }
        else
        {
          RealVectorValue rad_disp_vec = disp_vec - (disp_vec * axis_dir) * axis_dir;
          rad_disp = rad_disp_vec.norm();
        }
      }
      else
        mooseError("Cannot compute radial displacement for models with 1D Cartesian system");

      break;
    }

    case Moose::COORD_RZ:
    {
      if (_have_origin)
      {
        Point current_point (*_current_node);
        if (_ndisp == 2)
        {
          RealVectorValue rad_vec(current_point - _origin);
          Real rad = rad_vec.norm();
          const RealVectorValue disp_vec((*_disp_vals[0])[_qp], (*_disp_vals[1])[_qp]);
          if (rad > 0.0)
          {
            rad_vec /= rad;
            rad_disp = rad_vec * disp_vec;
          }
          else
            rad_disp = disp_vec.norm();
        }
      }
      else
        rad_disp = (*_disp_vals[0])[_qp];
      break;
    }

    case Moose::COORD_RSPHERICAL:
    {
      rad_disp = (*_disp_vals[0])[_qp];
      break;
    }

    default:
      mooseError("Unsupported coordinate system");
  }

  return rad_disp;
}
