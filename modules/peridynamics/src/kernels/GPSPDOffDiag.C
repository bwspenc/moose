/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "GPSPDOffDiag.h"
#include "Assembly.h"
#include "NonlinearSystem.h"

template<>
InputParameters validParams<GPSPDOffDiag>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Generalized Plane Strain kernel providing off-diagonal entries for the Jacobian matrix");
  params.addRequiredParam<std::vector<NonlinearVariableName> >("displacements", "Variables for the displacements");
  params.addParam<NonlinearVariableName>("temp", "Variable for the temperature");
  params.addCoupledVar("strain_zz", "Scalar variable for the strain_zz");
  params.addRequiredParam<NonlinearVariableName>("bond_status", "Auxiliary variable for failure status of each bond");
  params.addParam<std::string>("full_jacobian", "whether to use the nonlocal jacobian formulation for the scalar components");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

GPSPDOffDiag::GPSPDOffDiag(const InputParameters & parameters) :
  Kernel(parameters),
  _bond_dfdE_ij(getMaterialProperty<Real>("bond_dfdE_ij")),
  _bond_dfdE_i_j(getMaterialProperty<Real>("bond_dfdE_i_j")),
  _alpha(getMaterialProperty<Real>("thermal_expansion")),
  _Cijkl(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
  _shape(getMaterialProperty<RankTwoTensor>("shape_tensor")),
  _dgrad(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
  _aux(_fe_problem.getAuxiliarySystem()),
  _aux_sln(*_aux.currentSolution()),
  _nsys(_fe_problem.getNonlinearSystem()),
  _temp_var(isParamValid("temp") ? &_fe_problem.getVariable(_tid, getParam<NonlinearVariableName>("temp")) : NULL),
  _strain_zz_var(coupledScalar("strain_zz")),
  _bond_status_var(&_fe_problem.getVariable(_tid, getParam<NonlinearVariableName>("bond_status"))),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh)),
  _pddim(_pdmesh.dim())
{
  const std::vector<NonlinearVariableName> & nl_vnames(getParam<std::vector<NonlinearVariableName> >("displacements"));
  if (nl_vnames.size() != 2)
    mooseError("GeneralizedPlaneStrain only works for two dimensional case!");

  for (unsigned int i = 0; i < _pddim; ++i)
    _disp_var.push_back(&_fe_problem.getVariable(_tid, nl_vnames[i]));
}

void
GPSPDOffDiag::computeOffDiagJacobianScalar(unsigned int jvar)
{
  if (_var.number() == _disp_var[0]->number())
    if (isParamValid("full_jacobian"))
      computeDispFullOffDiagJacobianScalar(0, jvar);
    else
      computeDispPartialOffDiagJacobianScalar(0, jvar);
  else if (_var.number() == _disp_var[1]->number())
    if (isParamValid("full_jacobian"))
      computeDispFullOffDiagJacobianScalar(1, jvar);
    else
      computeDispPartialOffDiagJacobianScalar(1, jvar);
  else if (isParamValid("temp") ? _var.number() == _temp_var->number() : 0)
    computeTempOffDiagJacobianScalar(jvar);
}

void
GPSPDOffDiag::computeDispFullOffDiagJacobianScalar(unsigned int component, unsigned int jvar)
{
  if (jvar == _strain_zz_var)
  {
    DenseMatrix<Number> & ken = _assembly.jacobianBlock(_var.number(), jvar);
    DenseMatrix<Number> & kne = _assembly.jacobianBlock(jvar, _var.number());
    MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);

    const NumericVector<Number> & nsys_sln = *_nsys.currentSolution();
    dof_id_type node_i = _current_elem->get_node(0)->id();
    Node * nd_i = _mesh.nodePtr(node_i);
    dof_id_type node_j = _current_elem->get_node(1)->id();
    Node * nd_j = _mesh.nodePtr(node_j);

    RealGradient ori_ij(3);
    for (unsigned int i = 0; i < _pddim; ++i)
      ori_ij(i) = _pdmesh.coord(node_j)(i) + nsys_sln(nd_j->dof_number(_nsys.number(), _disp_var[i]->number(), 0)) - _pdmesh.coord(node_i)(i) - nsys_sln(nd_i->dof_number(_nsys.number(), _disp_var[i]->number(), 0));
    ori_ij /= ori_ij.size();

    dof_id_type bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
    Number bond_status_ij = _aux_sln(bs_dof);

    // fill in the column corresponding to the scalar variable from bond ij
    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < jv.order(); _j++)
        ken(_i, _j) += (_i == _j ? - 1 : 1) * ori_ij(component) * _bond_dfdE_ij[0] * bond_status_ij;

    _local_ke.resize(ken.m(), ken.n());
    _local_ke.zero();
    std::vector<dof_id_type> dof_ij(2);
    dof_ij[0] = nd_i->dof_number(_sys.number(), _var.number(), 0);
    dof_ij[1] = nd_j->dof_number(_sys.number(), _var.number(), 0);

    // NONLOCAL jacobian contribution
    std::vector<dof_id_type> dof(2), neighbors_i, bonds_i, neighbors_j, bonds_j;

    // calculation of jacobian contribution to node_i's neighbors
    dof[0] = dof_ij[0];
    neighbors_i = _pdmesh.neighbors(node_i);
    bonds_i = _pdmesh.bonds(node_i);
    for (unsigned int k = 0; k < neighbors_i.size(); ++k)
    {
      Node * nd_k = _mesh.nodePtr(neighbors_i[k]);
      dof[1] = nd_k->dof_number(_sys.number(), _var.number(), 0);
      double vol_k = _pdmesh.volume(neighbors_i[k]);

      // obtain bond ik's origin length and current orientation
      double origin_length_ik = 0;
      RealGradient ori_ik(3);
      for (unsigned int j = 0; j < _pddim; ++j)
      {
        origin_length_ik += std::pow(_pdmesh.coord(node_i)(j) - _pdmesh.coord(neighbors_i[k])(j), 2);
        ori_ik(j) = _pdmesh.coord(neighbors_i[k])(j) + nsys_sln(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_i)(j) - nsys_sln(nd_i->dof_number(_nsys.number(), j, 0));
      }
      origin_length_ik = std::sqrt(origin_length_ik);
      ori_ik /= ori_ik.size();

      // bond status for bond ik
      Elem * elem_k = _mesh.elemPtr(bonds_i[k]);
      dof_id_type bs_dof_ik = elem_k->dof_number(_aux.number(), _bond_status_var->number(), 0);
      Number bond_status_ik = _aux_sln(bs_dof_ik);

      _local_ke.zero();
      _local_ke(0, 0) = - ori_ik(component) * _bond_dfdE_i_j[0] / origin_length_ik * vol_k * bond_status_ik * bond_status_ij;
      _local_ke(1, 0) = ori_ik(component) * _bond_dfdE_i_j[0] / origin_length_ik * vol_k * bond_status_ik * bond_status_ij;

     _assembly.cacheJacobianBlock(_local_ke, dof, jv.dofIndices(), _var.scalingFactor());
    }

    // calculation of jacobian contribution to node_j's neighbors
    dof[1] = dof_ij[1];
    neighbors_j = _pdmesh.neighbors(node_j);
    bonds_j = _pdmesh.bonds(node_j);
    for (unsigned int k = 0; k < neighbors_j.size(); ++k)
    {
      Node * nd_k = _mesh.nodePtr(neighbors_j[k]);
      dof[0] = nd_k->dof_number(_sys.number(), _var.number(), 0);
      double vol_k = _pdmesh.volume(neighbors_j[k]);

      // obtain bond ik's origin length and current orientation
      double origin_length_jk = 0;
      RealGradient ori_jk(3);
      for (unsigned int j = 0; j < _pddim; ++j)
      {
        origin_length_jk += std::pow(_pdmesh.coord(node_j)(j) - _pdmesh.coord(neighbors_j[k])(j), 2);
        ori_jk(j) = _pdmesh.coord(neighbors_j[k])(j) + nsys_sln(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_j)(j) - nsys_sln(nd_j->dof_number(_nsys.number(), j, 0));
      }
      origin_length_jk = std::sqrt(origin_length_jk);
      ori_jk /= ori_jk.size();

      // bond status for bond jk
      Elem * elem_k = _mesh.elemPtr(bonds_j[k]);
      dof_id_type bs_dof_jk = elem_k->dof_number(_aux.number(), _bond_status_var->number(), 0);
      Number bond_status_jk = _aux_sln(bs_dof_jk);

      _local_ke.zero();
      _local_ke(0, 0) = ori_jk(component) * _bond_dfdE_i_j[1] / origin_length_jk * vol_k * bond_status_jk * bond_status_ij;
      _local_ke(1, 0) = - ori_jk(component) * _bond_dfdE_i_j[1] / origin_length_jk * vol_k * bond_status_jk * bond_status_ij;

     _assembly.cacheJacobianBlock(_local_ke, dof, jv.dofIndices(), _var.scalingFactor());
    }

    // off-diagonal jacobian entries on the row
    // nodal area for node i and j
    double nv_i = _pdmesh.volume(node_i);
    double nv_j = _pdmesh.volume(node_j);

    // horizon radius of node i and j
    double h_i = _pdmesh.horizon(node_i);
    double h_j = _pdmesh.horizon(node_j);

    RealGradient ov = *_current_elem->get_node(1) - *_current_elem->get_node(0);
    double dEidUi = - nv_j * h_i / ov.size() * (_Cijkl[0](2, 2, 0, 0) * (ov(0) * _shape[0](0, 0) + ov(1) * _shape[0](1, 0)) * _dgrad[0](component, 0) + _Cijkl[0](2, 2, 1, 1) * (ov(0) * _shape[0](0, 1) + ov(1) * _shape[0](1, 1)) * _dgrad[0](component, 1));
    double dEjdUj = nv_i * h_j / ov.size() * (_Cijkl[0](2, 2, 0, 0) * (ov(0) * _shape[1](0, 0) + ov(1) * _shape[1](1, 0)) * _dgrad[1](component, 0) + _Cijkl[0](2, 2, 1, 1) * (ov(0) * _shape[1](0, 1) + ov(1) * _shape[1](1, 1)) * _dgrad[1](component, 1));

    // fill in the row corresponding to the scalar variable
    kne(0, 0) += (dEidUi * nv_i - dEjdUj * nv_j) * bond_status_ij; // node i
    kne(0, 1) += (dEjdUj * nv_j - dEidUi * nv_i) * bond_status_ij; // node j
  }
}

void
GPSPDOffDiag::computeDispPartialOffDiagJacobianScalar(unsigned int component, unsigned int jvar)
{
  if (jvar == _strain_zz_var)
  {
    DenseMatrix<Number> & ken = _assembly.jacobianBlock(_var.number(), jvar);
    DenseMatrix<Number> & kne = _assembly.jacobianBlock(jvar, _var.number());
    MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);

    const NumericVector<Number> & nsys_sln = *_nsys.currentSolution();
    dof_id_type node_i = _current_elem->get_node(0)->id();
    Node * nd_i = _mesh.nodePtr(node_i);
    dof_id_type node_j = _current_elem->get_node(1)->id();
    Node * nd_j = _mesh.nodePtr(node_j);

    RealGradient ori_ij(3);
    for (unsigned int i = 0; i < _pddim; ++i)
      ori_ij(i) = _pdmesh.coord(node_j)(i) + nsys_sln(nd_j->dof_number(_nsys.number(), _disp_var[i]->number(), 0)) - _pdmesh.coord(node_i)(i) - nsys_sln(nd_i->dof_number(_nsys.number(), _disp_var[i]->number(), 0));
    ori_ij /= ori_ij.size();

    dof_id_type bs_dof = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
    Number bond_status_ij = _aux_sln(bs_dof);

    // fill in the column corresponding to the scalar variable from bond ij
    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < jv.order(); _j++)
        ken(_i, _j) += (_i == _j ? - 1 : 1) * ori_ij(component) * _bond_dfdE_ij[0] * bond_status_ij;

    // off-diagonal jacobian entries on the row
    // nodal area for node i and j
    double nv_i = _pdmesh.volume(node_i);
    double nv_j = _pdmesh.volume(node_j);

    // horizon radius of node i and j
    double h_i = _pdmesh.horizon(node_i);
    double h_j = _pdmesh.horizon(node_j);

    RealGradient ov = *_current_elem->get_node(1) - *_current_elem->get_node(0);
    double dEidUi = - nv_j * h_i / ov.size() * (_Cijkl[0](2, 2, 0, 0) * (ov(0) * _shape[0](0, 0) + ov(1) * _shape[0](1, 0)) * _dgrad[0](component, 0) + _Cijkl[0](2, 2, 1, 1) * (ov(0) * _shape[0](0, 1) + ov(1) * _shape[0](1, 1)) * _dgrad[0](component, 1));
    double dEjdUj = nv_i * h_j / ov.size() * (_Cijkl[0](2, 2, 0, 0) * (ov(0) * _shape[1](0, 0) + ov(1) * _shape[1](1, 0)) * _dgrad[1](component, 0) + _Cijkl[0](2, 2, 1, 1) * (ov(0) * _shape[1](0, 1) + ov(1) * _shape[1](1, 1)) * _dgrad[1](component, 1));

    // fill in the row corresponding to the scalar variable
    kne(0, 0) += (dEidUi * nv_i - dEjdUj * nv_j) * bond_status_ij; // node i
    kne(0, 1) += (dEjdUj * nv_j - dEidUi * nv_i) * bond_status_ij; // node j
  }
}

void
GPSPDOffDiag::computeTempOffDiagJacobianScalar(unsigned int jvar)
{
  if (jvar == _strain_zz_var)
  {
    DenseMatrix<Number> & kne = _assembly.jacobianBlock(jvar, _var.number());

    // off-diagonal jacobian entries on the column
    dof_id_type node_i = _current_elem->get_node(0)->id();
    dof_id_type node_j = _current_elem->get_node(1)->id();

    // nodal area for node i and j
    double nv_i = _pdmesh.volume(node_i);
    double nv_j = _pdmesh.volume(node_j);

    // number of neighbors for node i and j
    unsigned int nn_i = _pdmesh.neighbors(node_i).size();
    unsigned int nn_j = _pdmesh.neighbors(node_j).size();

    /*
      one-way coupling between the strain_zz and temperature.
      fill in the row corresponding to the scalar variable strain_zz
    */
    kne(0, 0) += - _alpha[0] * (_Cijkl[0](2, 2, 0, 0) + _Cijkl[0](2, 2, 1, 1) + _Cijkl[0](2, 2, 2, 2)) * nv_i / nn_i; // node i
    kne(0, 1) += - _alpha[0] * (_Cijkl[0](2, 2, 0, 0) + _Cijkl[0](2, 2, 1, 1) + _Cijkl[0](2, 2, 2, 2)) * nv_j / nn_j; // node j
  }
}
