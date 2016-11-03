/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergenceSPD.h"
#include "Assembly.h"
#include "NonlinearSystem.h"

template<>
InputParameters validParams<StressDivergenceSPD>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("SPD Stress divergence kernel");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the variable this kernel acts on (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements", "Variables for the displacements");
  params.addCoupledVar("temp", "Variable for the temperature");
  params.addRequiredParam<NonlinearVariableName>("bond_status", "Auxiliary variable for the bond failure status");
  params.addParam<std::string>("full_jacobian", "whether to use full SPD jacobian or not");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

StressDivergenceSPD::StressDivergenceSPD(const InputParameters & parameters) :
  Kernel(parameters),
  _bond_force_ij(getMaterialProperty<Real>("bond_force_ij")),
  _bond_force_i_j(getMaterialProperty<Real>("bond_force_i_j")),
  _bond_dfdU_ij(getMaterialProperty<Real>("bond_dfdU_ij")),
  _bond_dfdU_i_j(getMaterialProperty<Real>("bond_dfdU_i_j")),
  _bond_dfdT_ij(getMaterialProperty<Real>("bond_dfdT_ij")),
  _bond_dfdT_i_j(getMaterialProperty<Real>("bond_dfdT_i_j")),
  _aux(_fe_problem.getAuxiliarySystem()),
  _aux_sln(*_aux.currentSolution()),
  _nsys(_fe_problem.getNonlinearSystem()),
  _component(getParam<unsigned int>("component")),
  _ndisp(coupledComponents("displacements")),
  _temp_coupled(isCoupled("temp")),
  _temp_var(_temp_coupled ? coupled("temp") : 0),
  _bond_status_var(&_fe_problem.getVariable(_tid, getParam<NonlinearVariableName>("bond_status"))),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh)),
  _pddim(_pdmesh.dim()),
  _orientation(NULL)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var.push_back(coupled("displacements", i));
}

void
StressDivergenceSPD::initialSetup()
{
  _orientation = &_assembly.getFE(FEType(), 1)->get_dxyzdxi();
}

void
StressDivergenceSPD::computeResidual()
{
  // LOCAL residual contribution
  // bond status for current element ij
  dof_id_type bs_dof_ij = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  Number bond_status_ij = _aux_sln(bs_dof_ij);

  // calculation of residual contribution to node_i and node_j from current element
  RealGradient ori_ij = (*_orientation)[0];
  ori_ij /= ori_ij.size();

  _local_re.resize(2);
  _local_re(0) = - _bond_force_ij[0] * ori_ij(_component) * bond_status_ij;
  _local_re(1) = - _local_re(0);

  std::vector<dof_id_type> dof_ij(2);
  dof_id_type node_i = _current_elem->get_node(0)->id();
  Node * nd_i = _mesh.nodePtr(node_i);
  dof_id_type node_j = _current_elem->get_node(1)->id();
  Node * nd_j = _mesh.nodePtr(node_j);
  dof_ij[0] = nd_i->dof_number(_sys.number(), _var.number(), 0);
  dof_ij[1] = nd_j->dof_number(_sys.number(), _var.number(), 0);

  // cache the residual contribution to node_i and node_j using their global dof indices
  _assembly.cacheResidualNodes(_local_re, dof_ij);

  // save in the displacement residuals
  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _save_in.size(); i++)
      _save_in[i]->sys().solution().add_vector(_local_re, dof_ij);
  }

  // NONLOCAL residual contribution
  const NumericVector<Number> & nsys_sln = *_nsys.currentSolution();
  std::vector<dof_id_type> dof(2), neighbors_i, bonds_i, neighbors_j, bonds_j;

  // calculation of residual contribution to node_i's neighbors
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
      // should use MooseVariable->number() to get the value of j in order to access the nodal solution
      ori_ik(j) = _pdmesh.coord(neighbors_i[k])(j) + nsys_sln(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_i)(j) - nsys_sln(nd_i->dof_number(_nsys.number(), j, 0));
    }
    origin_length_ik = std::sqrt(origin_length_ik);
    ori_ik /= ori_ik.size();

    // bond status for bond ik
    Elem * elem_k = _mesh.elemPtr(bonds_i[k]);
    dof_id_type bs_dof_ik = elem_k->dof_number(_aux.number(), _bond_status_var->number(), 0);
    Number bond_status_ik = _aux_sln(bs_dof_ik);

    _local_re(0) = - _bond_force_i_j[0] * vol_k / origin_length_ik * ori_ik(_component) * bond_status_ik * bond_status_ij;
    _local_re(1) = - _local_re(0);

    // cache the residual contribution to node_i and its neighbor k using their global dof indices
    _assembly.cacheResidualNodes(_local_re, dof);

    // save in the displacement residuals
    if (_has_save_in)
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (unsigned int i = 0; i < _save_in.size(); i++)
        _save_in[i]->sys().solution().add_vector(_local_re, dof);
    }
  }

  // calculation of residual contribution to node_j's neighbors
  dof[1] = dof_ij[1];
  neighbors_j = _pdmesh.neighbors(node_j);
  bonds_j = _pdmesh.bonds(node_j);
  for (unsigned int k = 0; k < neighbors_j.size(); ++k)
  {
    Node * nd_k = _mesh.nodePtr(neighbors_j[k]);
    dof[0] = nd_k->dof_number(_sys.number(), _var.number(), 0);
    double vol_k = _pdmesh.volume(neighbors_j[k]);
    // obtain bond jk's origin length and current orientation
    double origin_length_jk = 0;
    RealGradient ori_jk(3);
    for (unsigned int j = 0; j < _pddim; ++j)
    {
      origin_length_jk += std::pow(_pdmesh.coord(neighbors_j[k])(j) - _pdmesh.coord(node_j)(j), 2);
      // should use MooseVariable->number() to get the value of j in order to access the nodal solution
      ori_jk(j) = _pdmesh.coord(neighbors_j[k])(j) + nsys_sln(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_j)(j) - nsys_sln(nd_j->dof_number(_nsys.number(), j, 0));
    }
    origin_length_jk = std::sqrt(origin_length_jk);
    ori_jk /= ori_jk.size();

    // bond status for bond jk
    Elem * elem_k = _mesh.elemPtr(bonds_j[k]);
    dof_id_type bs_dof_jk = elem_k->dof_number(_aux.number(), _bond_status_var->number(), 0);
    Number bond_status_jk = _aux_sln(bs_dof_jk);

    _local_re(0) = _bond_force_i_j[1] * vol_k / origin_length_jk * ori_jk(_component) * bond_status_jk * bond_status_ij;
    _local_re(1) = - _local_re(0);

    // cache the residual contribution to node_j and its neighbor k using their global dof indices
    _assembly.cacheResidualNodes(_local_re, dof);

    // save in the displacement residuals
    if (_has_save_in)
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (unsigned int i = 0; i < _save_in.size(); i++)
        _save_in[i]->sys().solution().add_vector(_local_re, dof);
    }
  }
}

void
StressDivergenceSPD::computeJacobian()
{
  if (isParamValid("full_jacobian"))
    computeFullJacobian();
  else
    computePartialJacobian();
}

void
StressDivergenceSPD::computeOffDiagJacobian(unsigned int jvar)
{
  if (isParamValid("full_jacobian"))
    computeFullOffDiagJacobian(jvar);
  else
    computePartialOffDiagJacobian(jvar);
}

void
StressDivergenceSPD::computePartialJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  RealGradient ori = (*_orientation)[0]; // ori gives the direction from j to i, from i to j is needed
  double current_length_ij = 2.0 * ori.size(); // ori.size() only gives half current length
  ori /= ori.size();

  // bond status for current element ij
  dof_id_type bs_dof_ij = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  Number bond_status_ij = _aux_sln(bs_dof_ij);

  double stiff_elem = ori(_component) * ori(_component) * _bond_dfdU_ij[0] + _bond_force_ij[0] * (1.0 - ori(_component) * ori(_component)) / current_length_ij;

  for (unsigned int i = 0; i < _test.size(); ++i)
    for (unsigned int j = 0; j < _phi.size(); ++j)
       _local_ke(i, j) += (i == j ? 1 : -1) * stiff_elem * bond_status_ij;

  ke += _local_ke;

  if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Real> diag(rows);
    for (unsigned int i = 0; i < rows; ++i)
      diag(i) = _local_ke(i,i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
      _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
  }
}

void
StressDivergenceSPD::computePartialOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _var.number())
    computePartialJacobian();
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
      coupled_component = 3;
      active = true;
    }

    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);
    if (active)
    {
      RealGradient ori = (*_orientation)[0]; // ori gives the direction from j to i, from i to j is needed
      double current_length_ij = 2.0 * ori.size(); // ori.size() only gives half current length
      ori /= ori.size();

      // bond status for current element ij
      dof_id_type bs_dof_ij = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
      Number bond_status_ij = _aux_sln(bs_dof_ij);

      double off_stiff_elem;
      if (coupled_component == 3)
      {
        off_stiff_elem = ori(_component) * _bond_dfdT_ij[0];
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            ke(i, j) += (i == 1 ? 1 : -1) * off_stiff_elem * bond_status_ij;
      }
      else
      {
        off_stiff_elem = ori(_component) * ori(coupled_component) * _bond_dfdU_ij[0] - _bond_force_ij[0] * ori(_component) * ori(coupled_component) / current_length_ij;
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            ke(i, j) += (i == j ? 1 : -1) * off_stiff_elem * bond_status_ij;
      }
    }
  }
}

void
StressDivergenceSPD::computeFullJacobian()
{
  // LOCAL jacobian contribution
  _local_ke.resize(2, 2);
  _local_ke.zero();

  RealGradient ori_ij(3); // orientation of bond ij
  ori_ij = (*_orientation)[0];
  Real current_length_ij = 2.0 * ori_ij.size(); // current length for bond ij
  ori_ij /= ori_ij.size();

  // bond status for current element ij
  dof_id_type bs_dof_ij = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
  Number bond_status_ij = _aux_sln(bs_dof_ij);

  double stiff_elem = ori_ij(_component) * ori_ij(_component) * _bond_dfdU_ij[0] + _bond_force_ij[0] * (1.0 - ori_ij(_component) * ori_ij(_component)) / current_length_ij;

  for (unsigned int i = 0; i < _test.size(); ++i)
    for (unsigned int j = 0; j < _phi.size(); ++j)
       _local_ke(i, j) += (i == j ? 1 : -1) * stiff_elem * bond_status_ij;

  std::vector<dof_id_type> dof_ij(2);
  dof_id_type node_i = _current_elem->get_node(0)->id();
  Node * nd_i = _mesh.nodePtr(node_i);
  dof_id_type node_j = _current_elem->get_node(1)->id();
  Node * nd_j = _mesh.nodePtr(node_j);
  dof_ij[0] = nd_i->dof_number(_sys.number(), _var.number(), 0);
  dof_ij[1] = nd_j->dof_number(_sys.number(), _var.number(), 0);

  _assembly.cacheJacobianBlock(_local_ke, dof_ij, dof_ij, _var.scalingFactor());

  if (_has_diag_save_in)
  {
    unsigned int rows = _test.size();
    DenseVector<Real> diag(rows);
    for (unsigned int i = 0; i < rows; ++i)
      diag(i) = _local_ke(i,i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
    _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
  }

  // NONLOCAL jacobian contribution
  const NumericVector<Number> & nsys_sln = *_nsys.currentSolution();
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
    double origin_length_ik = 0, current_length_ik = 0;
    RealGradient ori_ik(3);
    for (unsigned int j = 0; j < _pddim; ++j)
    {
      origin_length_ik += std::pow(_pdmesh.coord(node_i)(j) - _pdmesh.coord(neighbors_i[k])(j), 2);
      // should use MooseVariable->number() to get the value of j in order to access the nodal solution
      ori_ik(j) = _pdmesh.coord(neighbors_i[k])(j) + nsys_sln(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_i)(j) - nsys_sln(nd_i->dof_number(_nsys.number(), j, 0));
    }
    origin_length_ik = std::sqrt(origin_length_ik);
    current_length_ik = ori_ik.size();
    ori_ik /= current_length_ik;

    // bond status for bond ik
    Elem * elem_k = _mesh.elemPtr(bonds_i[k]);
    dof_id_type bs_dof_ik = elem_k->dof_number(_aux.number(), _bond_status_var->number(), 0);
    Number bond_status_ik = _aux_sln(bs_dof_ik);

    double diag1 = ori_ik(_component) * ori_ij(_component) * _bond_dfdU_i_j[0];
    double diag2 =  _bond_force_i_j[0] * (1.0 - ori_ik(_component) * ori_ik(_component)) / current_length_ik;

    _local_ke.zero();
    for (unsigned int i = 0; i < _test.size(); ++i)
      for (unsigned int j = 0; j < _phi.size(); ++j)
        _local_ke(i, j) += (i == j ? 1 : -1) * diag1 / origin_length_ik * vol_k * bond_status_ik * bond_status_ij;

    _assembly.cacheJacobianBlock(_local_ke, dof, dof_ij, _var.scalingFactor());

    if (_has_diag_save_in)
    {
      unsigned int rows = _test.size();
      DenseVector<Real> diag(rows);
      for (unsigned int i = 0; i < rows; ++i)
        diag(i) = _local_ke(i,i);

      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
        _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
    }

    _local_ke.zero();
    for (unsigned int i = 0; i < _test.size(); ++i)
      for (unsigned int j = 0; j < _phi.size(); ++j)
        _local_ke(i, j) += (i == j ? 1 : -1) * diag2 / origin_length_ik * vol_k * bond_status_ik * bond_status_ij;

    _assembly.cacheJacobianBlock(_local_ke, dof, dof, _var.scalingFactor());

    if (_has_diag_save_in)
    {
      unsigned int rows = _test.size();
      DenseVector<Real> diag(rows);
      for (unsigned int i = 0; i < rows; ++i)
        diag(i) = _local_ke(i,i);

      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
        _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
    }
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
    double origin_length_jk = 0, current_length_jk = 0;
    RealGradient ori_jk(3);
    for (unsigned int j = 0; j < _pddim; ++j)
    {
      origin_length_jk += std::pow(_pdmesh.coord(node_j)(j) - _pdmesh.coord(neighbors_j[k])(j), 2);
      // should use MooseVariable->number() to get the value of j in order to access the nodal solution
      ori_jk(j) = _pdmesh.coord(neighbors_j[k])(j) + nsys_sln(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_j)(j) - nsys_sln(nd_j->dof_number(_nsys.number(), j, 0));
    }
    origin_length_jk = std::sqrt(origin_length_jk);
    current_length_jk = ori_jk.size();
    ori_jk /= current_length_jk;

    // bond status for bond jk
    Elem * elem_k = _mesh.elemPtr(bonds_j[k]);
    dof_id_type bs_dof_jk = elem_k->dof_number(_aux.number(), _bond_status_var->number(), 0);
    Number bond_status_jk = _aux_sln(bs_dof_jk);

    double diag1 = - ori_jk(_component) * ori_ij(_component) * _bond_dfdU_i_j[1];
    double diag2 = _bond_force_i_j[1] * (1.0 - ori_jk(_component) * ori_jk(_component)) / current_length_jk;

    _local_ke.zero();
    for (unsigned int i = 0; i < _test.size(); ++i)
      for (unsigned int j = 0; j < _phi.size(); ++j)
        _local_ke(i, j) += (i == j ? 1 : -1) * diag1 / origin_length_jk * vol_k * bond_status_jk * bond_status_ij;

    _assembly.cacheJacobianBlock(_local_ke, dof, dof_ij, _var.scalingFactor());

    if (_has_diag_save_in)
    {
      unsigned int rows = _test.size();
      DenseVector<Real> diag(rows);
      for (unsigned int i = 0; i < rows; ++i)
        diag(i) = _local_ke(i,i);

      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
        _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
    }

    _local_ke.zero();
    for (unsigned int i = 0; i < _test.size(); ++i)
      for (unsigned int j = 0; j < _phi.size(); ++j)
        _local_ke(i, j) += (i == j ? 1 : -1) * diag2 / origin_length_jk * vol_k * bond_status_jk * bond_status_ij;

    _assembly.cacheJacobianBlock(_local_ke, dof, dof, _var.scalingFactor());

    if (_has_diag_save_in)
    {
      unsigned int rows = _test.size();
      DenseVector<Real> diag(rows);
      for (unsigned int i = 0; i < rows; ++i)
        diag(i) = _local_ke(i,i);

      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
        _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
    }
  }
}

void
StressDivergenceSPD::computeFullOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _var.number())
    computeFullJacobian();
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
      coupled_component = 3;
      active = true;
    }

    if (active)
    {
      RealGradient ori_ij;
      ori_ij = (*_orientation)[0];
      Real current_length_ij = 2.0 * ori_ij.size();
      ori_ij /= ori_ij.size();

      // bond status for current element ij
      dof_id_type bs_dof_ij = _current_elem->dof_number(_aux.number(), _bond_status_var->number(), 0);
      Number bond_status_ij = _aux_sln(bs_dof_ij);

      //the effect of truss ori change has been accounted for
      _local_ke.resize(2, 2);
      _local_ke.zero();
      double off_stiff_elem;
      if(coupled_component == 3)
      {
        off_stiff_elem = ori_ij(_component) * _bond_dfdT_ij[0];
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            _local_ke(i, j) += (i == 1 ? 1 : -1) * off_stiff_elem * bond_status_ij;
      }
      else
      {
        off_stiff_elem = ori_ij(_component) * ori_ij(coupled_component) * _bond_dfdU_ij[0] - _bond_force_ij[0] * ori_ij(_component) * ori_ij(coupled_component) / current_length_ij;
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            _local_ke(i, j) += (i == j ? 1 : -1) * off_stiff_elem * bond_status_ij;
      }

      std::vector<dof_id_type> dof_ij(2), dof_ij_jvar(2);
      dof_id_type node_i = _current_elem->get_node(0)->id();
      Node * nd_i = _mesh.nodePtr(node_i);
      dof_id_type node_j = _current_elem->get_node(1)->id();
      Node * nd_j = _mesh.nodePtr(node_j);
      dof_ij[0] = nd_i->dof_number(_sys.number(), _var.number(), 0);
      dof_ij[1] = nd_j->dof_number(_sys.number(), _var.number(), 0);
      dof_ij_jvar[0] = nd_i->dof_number(_sys.number(), jvar, 0);
      dof_ij_jvar[1] = nd_j->dof_number(_sys.number(), jvar, 0);

      _assembly.cacheJacobianBlock(_local_ke, dof_ij, dof_ij_jvar, _var.scalingFactor());

      // NONLOCAL jacobian contribution
      const NumericVector<Number> & nsys_sln = *_nsys.currentSolution();
      std::vector<dof_id_type> dof(2), dof_jvar(2), neighbors_i, bonds_i, neighbors_j, bonds_j;

      // calculation of jacobian contribution to node_i's neighbors
      dof[0] = dof_ij[0];
      dof_jvar[0] = dof_ij_jvar[0];
      neighbors_i = _pdmesh.neighbors(node_i);
      bonds_i = _pdmesh.bonds(node_i);
      for (unsigned int k = 0; k < neighbors_i.size(); ++k)
      {
        Node * nd_k = _mesh.nodePtr(neighbors_i[k]);
        dof[1] = nd_k->dof_number(_sys.number(), _var.number(), 0);
        dof_jvar[1] = nd_k->dof_number(_sys.number(), jvar, 0);
        double vol_k = _pdmesh.volume(neighbors_i[k]);

        // obtain bond ik's origin length and current orientation
        double origin_length_ik = 0, current_length_ik = 0;
        RealGradient ori_ik(3);
        for (unsigned int j = 0; j < _pddim; ++j)
        {
          origin_length_ik += std::pow(_pdmesh.coord(node_i)(j) - _pdmesh.coord(neighbors_i[k])(j), 2);
          // should use MooseVariable->number() to get the value of j in order to access the nodal solution
          ori_ik(j) = _pdmesh.coord(neighbors_i[k])(j) + nsys_sln(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_i)(j) - nsys_sln(nd_i->dof_number(_nsys.number(), j, 0));
        }
        origin_length_ik = std::sqrt(origin_length_ik);
        current_length_ik = ori_ik.size();
        ori_ik /= current_length_ik;

        // bond status for bond ik
        Elem * elem_k = _mesh.elemPtr(bonds_i[k]);
        dof_id_type bs_dof_ik = elem_k->dof_number(_aux.number(), _bond_status_var->number(), 0);
        Number bond_status_ik = _aux_sln(bs_dof_ik);

        _local_ke.zero();
        double off_diag1, off_diag2;
        if (coupled_component == 3)
        {
          off_diag1 = ori_ik(_component) * _bond_dfdT_i_j[0];
          off_diag2 = 0;
          _local_ke(0, 0) += - off_diag1 / origin_length_ik * vol_k * bond_status_ik * bond_status_ij;
          _local_ke(1, 0) += off_diag1 / origin_length_ik * vol_k * bond_status_ik * bond_status_ij;
        }
        else
        {
          off_diag1 = ori_ik(_component) * ori_ij(coupled_component) * _bond_dfdU_i_j[0];
          off_diag2 = - _bond_force_i_j[0] * ori_ik(_component) * ori_ik(coupled_component) / current_length_ik;
          for (unsigned int i = 0; i < _test.size(); ++i)
            for (unsigned int j = 0; j < _phi.size(); ++j)
              _local_ke(i, j) += (i == j ? 1 : -1) * off_diag1 / origin_length_ik * vol_k * bond_status_ik * bond_status_ij;
        }

        _assembly.cacheJacobianBlock(_local_ke, dof, dof_ij_jvar, _var.scalingFactor());

        _local_ke.zero();
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            _local_ke(i, j) += (i == j ? 1 : -1) * off_diag2 / origin_length_ik * vol_k * bond_status_ik * bond_status_ij;

        _assembly.cacheJacobianBlock(_local_ke, dof, dof_jvar, _var.scalingFactor());
      }

      // calculation of jacobian contribution to node_j's neighbors
      dof[1] = dof_ij[1];
      dof_jvar[1] = dof_ij_jvar[1];
      neighbors_j = _pdmesh.neighbors(node_j);
      bonds_j = _pdmesh.bonds(node_j);
      for (unsigned int k = 0; k < neighbors_j.size(); ++k)
      {
        Node * nd_k = _mesh.nodePtr(neighbors_j[k]);
        dof[0] = nd_k->dof_number(_sys.number(), _var.number(), 0);
        dof_jvar[0] = nd_k->dof_number(_sys.number(), jvar, 0);
        double vol_k = _pdmesh.volume(neighbors_j[k]);

        // obtain bond ik's origin length and current orientation
        double origin_length_jk = 0, current_length_jk = 0;
        RealGradient ori_jk(3);
        for (unsigned int j = 0; j < _pddim; ++j)
        {
          origin_length_jk += std::pow(_pdmesh.coord(node_j)(j) - _pdmesh.coord(neighbors_j[k])(j), 2);
          // should use MooseVariable->number() to get the value of j in order to access the nodal solution
          ori_jk(j) = _pdmesh.coord(neighbors_j[k])(j) + nsys_sln(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_j)(j) - nsys_sln(nd_j->dof_number(_nsys.number(), j, 0));
        }
        origin_length_jk = std::sqrt(origin_length_jk);
        current_length_jk = ori_jk.size();
        ori_jk /= current_length_jk;

        // bond status for bond jk
        Elem * elem_k = _mesh.elemPtr(bonds_j[k]);
        dof_id_type bs_dof_jk = elem_k->dof_number(_aux.number(), _bond_status_var->number(), 0);
        Number bond_status_jk = _aux_sln(bs_dof_jk);

        _local_ke.zero();
        double off_diag1, off_diag2;
        if (coupled_component == 3)
        {
          off_diag1 = ori_jk(_component) * _bond_dfdT_i_j[1];
          off_diag2 = 0;
          _local_ke(0, 1) += off_diag1 / origin_length_jk * vol_k * bond_status_jk * bond_status_ij;
          _local_ke(1, 1) += - off_diag1 / origin_length_jk * vol_k * bond_status_jk * bond_status_ij;
        }
        else
        {
          off_diag1 = - ori_jk(_component) * ori_ij(coupled_component) * _bond_dfdU_i_j[1];
          off_diag2 = - _bond_force_i_j[1] * ori_jk(_component) * ori_jk(coupled_component) / current_length_jk;
          for (unsigned int i = 0; i < _test.size(); ++i)
            for (unsigned int j = 0; j < _phi.size(); ++j)
              _local_ke(i, j) += (i == j ? 1 : -1) * off_diag1 / origin_length_jk * vol_k * bond_status_jk * bond_status_ij;
        }

        _assembly.cacheJacobianBlock(_local_ke, dof, dof_ij_jvar, _var.scalingFactor());

        _local_ke.zero();
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            _local_ke(i, j) += (i == j ? 1 : -1) * off_diag2 / origin_length_jk * vol_k * bond_status_jk * bond_status_ij;

        _assembly.cacheJacobianBlock(_local_ke, dof, dof_jvar, _var.scalingFactor());
      }
    }
  }
}
