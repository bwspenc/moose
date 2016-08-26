/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "GeneralizedPlaneStrainPD.h"
#include "Material.h"
#include "Assembly.h"

template<>
InputParameters validParams<GeneralizedPlaneStrainPD>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("BPD Generalized Plane Strain kernel");
  params.addRequiredCoupledVar("displacements", "The coupled variables for displacement (disp_x and disp_y)");
  params.addCoupledVar("temp", "The coupled variable for temperature");
  params.addRequiredParam<NonlinearVariableName>("bond_status", "Auxiliary variable for failure status of each bond");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

GeneralizedPlaneStrainPD::GeneralizedPlaneStrainPD(const InputParameters & parameters) :
  Kernel(parameters),
  _Cijkl(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
  _alpha(getMaterialProperty<Real>("thermal_expansion")),
  _shape_tensor(getMaterialProperty<RankTwoTensor>("shape_tensor")),
  _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
  _stress(getMaterialProperty<RankTwoTensor>("stress")),
  _aux(_fe_problem.getAuxiliarySystem()),
  _aux_sln(*_aux.currentSolution()),
  _ndisp(coupledComponents("displacements")),
  _temp_coupled(isCoupled("temp")),
  _temp_var(_temp_coupled ? coupled("temp") : 0),
  _bond_status_var(&_fe_problem.getVariable(_tid, getParam<NonlinearVariableName>("bond_status"))),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh))
{
  if (_ndisp != 2)
    mooseError("Size of displacements vector can only be 2!");
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var.push_back(coupled("displacements", i));
}

void
GeneralizedPlaneStrainPD::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());

  // nodal area for node i and j
  double nv_i = _pdmesh.volume(_current_elem->get_node(0)->id());
  double nv_j = _pdmesh.volume(_current_elem->get_node(1)->id());

  // number of neighbors for node i and j
  unsigned int nn_i = _pdmesh.neighbors(_current_elem->get_node(0)->id()).size();
  unsigned int nn_j = _pdmesh.neighbors(_current_elem->get_node(1)->id()).size();

  re(0) += - _stress[0](2, 2) * nv_i / nn_i;
  re(1) += - _stress[1](2, 2) * nv_j / nn_j;
}

void
GeneralizedPlaneStrainPD::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());

  // nodal area for node i and j
  double nv_i = _pdmesh.volume(_current_elem->get_node(0)->id());
  double nv_j = _pdmesh.volume(_current_elem->get_node(1)->id());

  // number of neighbors for node i and j
  unsigned int nn_i = _pdmesh.neighbors(_current_elem->get_node(0)->id()).size();
  unsigned int nn_j = _pdmesh.neighbors(_current_elem->get_node(1)->id()).size();

  ke(0, 0) += - _Cijkl[0](0, 0, 0, 0) * nv_i / nn_i;
  ke(1, 1) += - _Cijkl[0](0, 0, 0, 0) * nv_j / nn_j;
}

void
GeneralizedPlaneStrainPD::computeOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _var.number())
    computeJacobian();
  else
  {
    unsigned int coupled_component = 0;
    bool active = false;
    for (unsigned int i = 0; i < _ndisp; ++i)
      if (jvar == _disp_var[i])
      {
        coupled_component = i;
        active = true;
      }

    if (_temp_coupled && jvar == _temp_var)
    {
      coupled_component = 2;
      active = true;
    }

    // bond status for current element
    dof_id_type bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
    Number bond_status = _aux_sln(bs_dof);

    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);
    if (active)
    {
      // nodal area for node i and j
      double nv_i = _pdmesh.volume(_current_elem->get_node(0)->id());
      double nv_j = _pdmesh.volume(_current_elem->get_node(1)->id());

      // horizon radius of node i and j
      double horizon_i = _pdmesh.horizon(_current_elem->get_node(0)->id());
      double horizon_j = _pdmesh.horizon(_current_elem->get_node(1)->id());

      // number of neighbors for node i and j
      unsigned int nn_i = _pdmesh.neighbors(_current_elem->get_node(0)->id()).size();
      unsigned int nn_j = _pdmesh.neighbors(_current_elem->get_node(1)->id()).size();

      RealGradient origin_vector = *_current_elem->get_node(1) - *_current_elem->get_node(0);

      double off_stiff_elem;
      if(coupled_component == 2)
      {
        ke(0, 0) += _alpha[0] * (2.0 * _Cijkl[0](0, 0, 1, 1) + _Cijkl[0](0, 0, 0, 0)) * nv_i / nn_i;
        ke(1, 1) += _alpha[0] * (2.0 * _Cijkl[0](0, 0, 1, 1) + _Cijkl[0](0, 0, 0, 0)) * nv_j / nn_j;
      }
      else if(coupled_component == 1)
      {
        // derivative of strains w.r.t displacements for node i
        double dEidVi = nv_j * horizon_i / origin_vector.size() * ((origin_vector(0) * _shape_tensor[0](0, 0) + origin_vector(1) * _shape_tensor[0](1, 0)) * _deformation_gradient[0](1, 0) + (origin_vector(0) * _shape_tensor[0](0, 1) + origin_vector(1) * _shape_tensor[0](1, 1)) * _deformation_gradient[0](1, 1));
        double dEjdVj = - nv_i * horizon_j / origin_vector.size() * ((origin_vector(0) * _shape_tensor[1](0, 0) + origin_vector(1) * _shape_tensor[1](1, 0)) * _deformation_gradient[1](1, 0) + (origin_vector(0) * _shape_tensor[1](0, 1) + origin_vector(1) * _shape_tensor[1](1, 1)) * _deformation_gradient[1](1, 1));

        ke(0, 0) += _Cijkl[0](0, 0, 1, 1) * dEidVi * nv_i * bond_status;
        ke(0, 1) += - _Cijkl[0](0, 0, 1, 1) * dEidVi * nv_i * bond_status;
        ke(1, 0) += - _Cijkl[0](0, 0, 1, 1) * dEjdVj * nv_j * bond_status;
        ke(1, 1) += _Cijkl[0](0, 0, 1, 1) * dEjdVj * nv_j * bond_status;
      }
      else if(coupled_component == 0)
      {
        // derivative of strains w.r.t displacement u for node i and j
        double dEidUi = nv_j * horizon_i / origin_vector.size() * ((origin_vector(0) * _shape_tensor[0](0, 0) + origin_vector(1) * _shape_tensor[0](1, 0)) * _deformation_gradient[0](0, 0) + (origin_vector(0) * _shape_tensor[0](0, 1) + origin_vector(1) * _shape_tensor[0](1, 1)) * _deformation_gradient[0](0, 1));
        double dEjdUj = - nv_i * horizon_j / origin_vector.size() * ((origin_vector(0) * _shape_tensor[1](0, 0) + origin_vector(1) * _shape_tensor[1](1, 0)) * _deformation_gradient[1](0, 0) + (origin_vector(0) * _shape_tensor[1](0, 1) + origin_vector(1) * _shape_tensor[1](1, 1)) * _deformation_gradient[1](0, 1));

        ke(0, 0) += _Cijkl[0](0, 0, 1, 1) * dEidUi * nv_i * bond_status;
        ke(0, 1) += - _Cijkl[0](0, 0, 1, 1) * dEidUi * nv_i * bond_status;
        ke(1, 0) += - _Cijkl[0](0, 0, 1, 1) * dEjdUj * nv_j * bond_status;
        ke(1, 1) += _Cijkl[0](0, 0, 1, 1) * dEjdUj * nv_j * bond_status;
      }
    }
  }
}
