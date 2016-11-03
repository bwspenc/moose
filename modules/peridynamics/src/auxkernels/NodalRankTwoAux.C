/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "NodalRankTwoAux.h"
#include "RankTwoScalarTools.h"

template<>
InputParameters validParams<NodalRankTwoAux>()
{
  MooseEnum QuantityTypes("Component VonMisesStress Hydrostatic L2norm MaxPrincipal MidPrincipal MinPrincipal VolumetricStrain FirstInvariant SecondInvariant ThirdInvariant AxialStress HoopStress RadialStress TriaxialityStress Direction");
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Access components and scalar quantities of rank two strain/stress tensor");
  params.addRequiredParam<std::vector<NonlinearVariableName> >("displacements", "Variables containing the displacements for axukernel NodalRankTwoAux");
  params.addParam<std::string>("rank_two_tensor", "Which rank two tensor: strain or stress");
  params.addParam<Real>("youngs_modulus", "Young's modulus");
  params.addParam<Real>("poissons_ratio", "Poisson's ratio");
  params.addParam<MooseEnum>("quantity_type", QuantityTypes, "Type of output");
  params.addRequiredRangeCheckedParam<unsigned int>("index_i", "index_i >= 0 & index_i <= 2", "The index i of ij for the tensor to output (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>("index_j", "index_j >= 0 & index_j <= 2", "The index j of ij for the tensor to output (0, 1, 2)");
  params.addParam<Point>("point1", Point(0, 0, 0), "Start point for axis used to calculate some cylinderical material tensor quantities" );
  params.addParam<Point>("point2", Point(0, 1, 0), "End point for axis used to calculate some material tensor quantities");
  params.addParam<Point>("direction", Point(0, 0, 1), "Direction vector");
  return params;
}

NodalRankTwoAux::NodalRankTwoAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _fe_problem(*parameters.get<FEProblem *>("_fe_problem")),
  _nsys(_fe_problem.getNonlinearSystem()),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh)),
  _pddim(_pdmesh.dim()),
  _rank_two_tensor(getParam<std::string>("rank_two_tensor")),
  _youngs_modulus(isParamValid("youngs_modulus") ? getParam<Real>("youngs_modulus") : 1.0),
  _poissons_ratio(isParamValid("poissons_ratio") ? getParam<Real>("poissons_ratio") : 0.0),
  _quantity_type(getParam<MooseEnum>("quantity_type")),
  _i(getParam<unsigned int>("index_i")),
  _j(getParam<unsigned int>("index_j")),
  _point1(parameters.get<Point>("point1")),
  _point2(parameters.get<Point>("point2")),
  _direction(parameters.get<Point>("direction") / parameters.get<Point>("direction").norm())
{
  const std::vector<NonlinearVariableName> & nl_vnames(getParam<std::vector<NonlinearVariableName> >("displacements"));
  for (unsigned int i = 0; i < _pddim; ++i)
    _disp_var.push_back(&_fe_problem.getVariable(_tid, nl_vnames[i]));

  // plane strain case
  if (_rank_two_tensor == "stress" && !isParamValid("youngs_modulus"))
    mooseError("Both Young's modulus and Poisson's ratio must be provided for stress calculation");
  else
  {
    std::vector<Real> iso_const(2);
    iso_const[0] = _youngs_modulus * _poissons_ratio / ((1.0 + _poissons_ratio) * (1.0 - 2.0 * _poissons_ratio));
    iso_const[1] = _youngs_modulus / (2.0 * (1.0 + _poissons_ratio));

    //fill elasticity tensor
    _Cijkl.fillFromInputVector(iso_const, RankFourTensor::symmetric_isotropic);
  }
}

Real
NodalRankTwoAux::computeValue()
{
  Real val = 0.0;
  RankTwoTensor tensor;
  if (_rank_two_tensor == "strain")
    tensor = computeNodalStrain();
  else if (_rank_two_tensor == "stress")
    tensor = computeNodalStress();
  else
    mooseError("NodalRankTwoAux Error: Pass valid rank two tensor type - strain or stress");

  switch (_quantity_type)
  {
    case 0:
      val = RankTwoScalarTools::component(tensor, _i, _j);
      break;
    case 1:
      val = RankTwoScalarTools::vonMisesStress(tensor);
      break;
    case 2:
      val = RankTwoScalarTools::hydrostatic(tensor);
      break;
    case 3:
      val = RankTwoScalarTools::L2norm(tensor);
      break;
    case 4:
      val = RankTwoScalarTools::maxPrinciple(tensor);
      break;
    case 5:
      val = RankTwoScalarTools::midPrinciple(tensor);
      break;
    case 6:
      val = RankTwoScalarTools::minPrinciple(tensor);
      break;
    case 7:
      val = RankTwoScalarTools::volumetricStrain(tensor);
      break;
    case 8:
      val = RankTwoScalarTools::firstInvariant(tensor);
      break;
    case 9:
      val = RankTwoScalarTools::secondInvariant(tensor);
      break;
    case 10:
      val = RankTwoScalarTools::thirdInvariant(tensor);
      break;
    case 11:
      val = RankTwoScalarTools::axialStress(tensor, _point1, _point2, _direction);
      break;
    case 12:
      val = RankTwoScalarTools::hoopStress(tensor, _point1, _point2, *_current_node, _direction);
      break;
    case 13:
      val = RankTwoScalarTools::radialStress(tensor, _point1, _point2, *_current_node, _direction);
      break;
    case 14:
      val = RankTwoScalarTools::triaxialityStress(tensor);
      break;
    case 15:
      val = RankTwoScalarTools::directionValueTensor(tensor, _direction);
      break;
    default:
      mooseError("NodalRankTwoAux Error: Pass valid scalar type - Component VonMisesStress Hydrostatic L2norm MaxPrincipal MidPrincipal MinPrincipal VolumetricStrain FirstInvariant SecondInvariant ThirdInvariant AxialStress HoopStress RadialStress TriaxialityStress Direction");
  }

  return val;
}

RankTwoTensor
NodalRankTwoAux::computeNodalStrain()
{
  std::vector<dof_id_type> neighbors = _pdmesh.neighbors(_current_node->id());
  double horizon_i = _pdmesh.horizon(_current_node->id());

  // calculate the shape tensor and prepare the deformation gradient tensor
  RankTwoTensor shape, dgrad, delta;
  shape.zero();
  dgrad.zero();
  delta.zero();
  delta(0, 0) = delta(1, 1) = delta(2, 2) = 1;
  if (_pddim == 2)
    shape(2, 2) = dgrad(2, 2)  = 1;

  const NumericVector<Number> & sol = *_nsys.currentSolution();
  RealGradient origin_vector(3), current_vector(3);
  origin_vector = 0;
  current_vector = 0;

  for (unsigned int j = 0; j < neighbors.size(); ++j)
  {
    Node * node_j = _mesh.nodePtr(neighbors[j]);
    double vol_j = _pdmesh.volume(neighbors[j]);
    for (unsigned int k = 0; k < _pddim; ++k)
    {
      origin_vector(k) = _pdmesh.coord(neighbors[j])(k) - _pdmesh.coord(_current_node->id())(k);
      current_vector(k) = origin_vector(k) + sol(node_j->dof_number(_nsys.number(), _disp_var[k]->number(), 0)) - sol(_current_node->dof_number(_nsys.number(), _disp_var[k]->number(), 0));
    }
    double origin_length = origin_vector.norm();
    for (unsigned int k = 0; k < _pddim; ++k)
      for (unsigned int l = 0; l < _pddim; ++l)
      {
        shape(k, l) += vol_j * horizon_i / origin_length * origin_vector(k) * origin_vector(l);
        dgrad(k, l) += vol_j * horizon_i / origin_length * current_vector(k) * origin_vector(l);
      }
  }

  //finalize the deformation gradient tensor
  dgrad *= shape.inverse();

  // the green-lagrange strain tensor
  return (dgrad.transpose() * dgrad - delta) / 2.0;
}

RankTwoTensor
NodalRankTwoAux::computeNodalStress()
{
  return _Cijkl * computeNodalStrain();
}
