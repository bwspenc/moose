/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergenceSPD.h"
#include "MooseMesh.h"
#include "Material.h"
#include "Assembly.h"
#include "NonlinearSystem.h"

template<>
InputParameters validParams<StressDivergenceSPD>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("SPD Stress divergence kernel");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the variable this kernel acts on (0 for x, 1 for y, 2 for z)");
  params.addParam<std::string>("full_jacobian", "whether to use full SPD jacobian or not");
  params.addRequiredCoupledVar("displacements", "The displacement variables");
  params.addCoupledVar("temp", "The temperature variable");
  params.addCoupledVar("strain_zz", "The strain_zz variable");
  params.addCoupledVar("bond_status", "Auxiliary variable for failure status of each bond");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

StressDivergenceSPD::StressDivergenceSPD(const InputParameters & parameters) :
  Kernel(parameters),
  _bond_force_ij(getMaterialProperty<Real>("bond_force_ij")),
  _bond_force_i_j(getMaterialProperty<Real>("bond_force_i_j")),
  _bond_dfdU_ij(getMaterialProperty<Real>("bond_dfdU_ij")),
  _bond_dfdU_i_j(getMaterialProperty<Real>("bond_dfdU_i_j")),
  _bond_dfdE_ij(getMaterialProperty<Real>("bond_dfdE_ij")),
  _bond_dfdE_i_j(getMaterialProperty<Real>("bond_dfdE_i_j")),
  _bond_dfdT_ij(getMaterialProperty<Real>("bond_dfdT_ij")),
  _bond_dfdT_i_j(getMaterialProperty<Real>("bond_dfdT_i_j")),
  _aux(_fe_problem.getAuxiliarySystem()),
  _nsys(_fe_problem.getNonlinearSystem()),
  _component(getParam<unsigned int>("component")),
  _ndisp(coupledComponents("displacements")),
  _temp_coupled(isCoupled("temp")),
  _temp_var(_temp_coupled ? coupled("temp") : 0),
  _strain_zz_coupled(isCoupled("strain_zz")),
  _strain_zz_var(_strain_zz_coupled ? coupled("strain_zz") : 0),
  _bond_status_var(getVar("bond_status", 0)),
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
  // calculation of residual contribution to node_i and node_j from current element
  RealGradient ori_ij = (*_orientation)[0];
  ori_ij /= ori_ij.size();

  _local_re.resize(2);
  _local_re(0) = - _bond_force_ij[0] * ori_ij(_component);
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
  const NumericVector<Number> & sol = *_nsys.currentSolution();
  std::vector<dof_id_type> dof(2), neighbors;

  // calculation of residual contribution to node_i's neighbors
  dof[0] = dof_ij[0];
  neighbors = _pdmesh.neighbors(node_i);
  for (unsigned int k = 0; k < neighbors.size(); ++k)
  {
    Node * nd_k = _mesh.nodePtr(neighbors[k]);
    dof[1] = nd_k->dof_number(_sys.number(), _var.number(), 0);
    double vol_k = _pdmesh.volume(neighbors[k]);
    // obtain bond ik's origin length and current orientation
    double origin_length = 0;
    RealGradient ori_ik(3);
    for (unsigned int j = 0; j < _pddim; ++j)
    {
      origin_length += std::pow(_pdmesh.coord(node_i)(j) - _pdmesh.coord(neighbors[k])(j), 2);
      ori_ik(j) = _pdmesh.coord(neighbors[k])(j) + sol(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_i)(j) - sol(nd_i->dof_number(_nsys.number(), j, 0));
    }
    origin_length = std::sqrt(origin_length);
    ori_ik /= ori_ik.size();
    _local_re(0) = - _bond_force_i_j[0] * vol_k / origin_length * ori_ik(_component);
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

  // calculation of residual contribution to node_j's neighbor
  dof[1] = dof_ij[1];
  neighbors = _pdmesh.neighbors(node_j);
  for (unsigned int k = 0; k < neighbors.size(); ++k)
  {
    Node * nd_k = _mesh.nodePtr(neighbors[k]);
    dof[0] = nd_k->dof_number(_sys.number(), _var.number(), 0);
    double vol_k = _pdmesh.volume(neighbors[k]);
    // obtain bond jk's origin length and current orientation
    double origin_length = 0;
    RealGradient ori_jk(3);
    for (unsigned int j = 0; j < _pddim; ++j)
    {
      origin_length += std::pow(_pdmesh.coord(neighbors[k])(j) - _pdmesh.coord(node_j)(j), 2);
      ori_jk(j) = _pdmesh.coord(neighbors[k])(j) + sol(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_j)(j) - sol(nd_j->dof_number(_nsys.number(), j, 0));
    }
    origin_length = std::sqrt(origin_length);
    ori_jk /= ori_jk.size();
    _local_re(0) = _bond_force_i_j[1] * vol_k / origin_length * ori_jk(_component);
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
  double current_length = 2.0 * ori.size(); // ori.size() only gives half current length
  ori /= ori.size();

  double stiff_elem = ori(_component) * ori(_component) * _bond_dfdU_ij[0] + _bond_force_ij[0] * (1.0 - ori(_component) * ori(_component)) / current_length;

  for (unsigned int i = 0; i < _test.size(); ++i)
    for (unsigned int j = 0; j < _phi.size(); ++j)
       _local_ke(i, j) += (i == j ? 1 : -1) * stiff_elem;

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

    if (_strain_zz_coupled && jvar == _strain_zz_var)
    {
      coupled_component = 4;
      active = true;
    }

    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);
    if (active)
    {
      RealGradient ori = (*_orientation)[0]; // ori gives the direction from j to i, from i to j is needed
      double current_length = 2.0 * ori.size(); // ori.size() only gives half current length
      ori /= ori.size();

      double off_stiff_elem;
      if (coupled_component == 3)
      {
        off_stiff_elem = ori(_component) * _bond_dfdT_ij[0];
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            ke(i, j) += (i == 1 ? 1 : -1) * off_stiff_elem;
      }
      else if (coupled_component == 4)
      {
        off_stiff_elem = ori(_component) * _bond_dfdE_ij[0];
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            ke(i, j) += (i == 1 ? 1 : -1) * off_stiff_elem;
      }
      else
      {
        off_stiff_elem = ori(_component) * ori(coupled_component) * _bond_dfdU_ij[0] - _bond_force_ij[0] * ori(_component) * ori(coupled_component) / current_length;
      for (unsigned int i = 0; i < _test.size(); ++i)
        for (unsigned int j = 0; j < _phi.size(); ++j)
          ke(i, j) += (i == j ? 1 : -1) * off_stiff_elem;
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

  double stiff_elem = ori_ij(_component) * ori_ij(_component) * _bond_dfdU_ij[0] + _bond_force_ij[0] * (1.0 - ori_ij(_component) * ori_ij(_component)) / current_length_ij;

  for (unsigned int i = 0; i < _test.size(); ++i)
    for (unsigned int j = 0; j < _phi.size(); ++j)
       _local_ke(i, j) += (i == j ? 1 : -1) * stiff_elem;

  std::vector<dof_id_type> dof_ij(2);
  dof_id_type node_i = _current_elem->get_node(0)->id();
  Node * nd_i = _mesh.nodePtr(node_i);
  dof_id_type node_j = _current_elem->get_node(1)->id();
  Node * nd_j = _mesh.nodePtr(node_j);
  dof_ij[0] = nd_i->dof_number(_sys.number(), _var.number(), 0);
  dof_ij[1] = nd_j->dof_number(_sys.number(), _var.number(), 0);

  _assembly.cacheJacobianBlock(_local_ke, dof_ij, dof_ij, 1.0);

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
  const NumericVector<Number> & sol = *_nsys.currentSolution();
  std::vector<dof_id_type> dof(2), neighbors;

  // calculation of jacobian contribution to node_i's neighbors
  dof[0] = dof_ij[0];
  neighbors = _pdmesh.neighbors(node_i);
  for (unsigned int k = 0; k < neighbors.size(); ++k)
  {
    Node * nd_k = _mesh.nodePtr(neighbors[k]);
    dof[1] = nd_k->dof_number(_sys.number(), _var.number(), 0);
    double vol_k = _pdmesh.volume(neighbors[k]);
    // obtain bond ik's origin length and current orientation
    double origin_length = 0, current_length = 0;
    RealGradient ori_ik(3);
    for (unsigned int j = 0; j < _pddim; ++j)
    {
      origin_length += std::pow(_pdmesh.coord(node_i)(j) - _pdmesh.coord(neighbors[k])(j), 2);
      ori_ik(j) = _pdmesh.coord(neighbors[k])(j) + sol(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_i)(j) - sol(nd_i->dof_number(_nsys.number(), j, 0));
    }
    origin_length = std::sqrt(origin_length);
    current_length = ori_ik.size();
    ori_ik /= current_length;
    double stiff_elem1 = ori_ik(_component) * ori_ij(_component) * _bond_dfdU_i_j[0];
    double stiff_elem2 =  _bond_force_i_j[0] * (1.0 - ori_ik(_component) * ori_ik(_component)) / current_length;

    _local_ke.zero();
    for (unsigned int i = 0; i < _test.size(); ++i)
      for (unsigned int j = 0; j < _phi.size(); ++j)
        _local_ke(i, j) += (i == j ? 1 : -1) * stiff_elem1 / origin_length * vol_k;

    _assembly.cacheJacobianBlock(_local_ke, dof, dof_ij, 1.0);

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
        _local_ke(i, j) += (i == j ? 1 : -1) * stiff_elem2 / origin_length * vol_k;

    _assembly.cacheJacobianBlock(_local_ke, dof, dof, 1.0);

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
  neighbors = _pdmesh.neighbors(node_j);
  for (unsigned int k = 0; k < neighbors.size(); ++k)
  {
    Node * nd_k = _mesh.nodePtr(neighbors[k]);
    dof[0] = nd_k->dof_number(_sys.number(), _var.number(), 0);
    double vol_k = _pdmesh.volume(neighbors[k]);
    // obtain bond ik's origin length and current orientation
    double origin_length = 0, current_length = 0;
    RealGradient ori_jk(3);
    for (unsigned int j = 0; j < _pddim; ++j)
    {
      origin_length += std::pow(_pdmesh.coord(node_j)(j) - _pdmesh.coord(neighbors[k])(j), 2);
      ori_jk(j) = _pdmesh.coord(neighbors[k])(j) + sol(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_j)(j) - sol(nd_j->dof_number(_nsys.number(), j, 0));
    }
    origin_length = std::sqrt(origin_length);
    current_length = ori_jk.size();
    ori_jk /= current_length;
    double stiff_elem1 = - ori_jk(_component) * ori_ij(_component) * _bond_dfdU_i_j[1];
    double stiff_elem2 = _bond_force_i_j[1] * (1.0 - ori_jk(_component) * ori_jk(_component)) / current_length;

    _local_ke.zero();
    for (unsigned int i = 0; i < _test.size(); ++i)
      for (unsigned int j = 0; j < _phi.size(); ++j)
        _local_ke(i, j) += (i == j ? 1 : -1) * stiff_elem1 / origin_length * vol_k;

    _assembly.cacheJacobianBlock(_local_ke, dof, dof_ij, 1.0);

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
        _local_ke(i, j) += (i == j ? 1 : -1) * stiff_elem2 / origin_length * vol_k;

    _assembly.cacheJacobianBlock(_local_ke, dof, dof, 1.0);

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

    if (_strain_zz_coupled && jvar == _strain_zz_var)
    {
      coupled_component = 4;
      active = true;
    }

    if (active)
    {
      RealGradient ori_ij;
      ori_ij = (*_orientation)[0];
      Real current_length_ij = 2.0 * ori_ij.size();
      ori_ij /= ori_ij.size();
      //the effect of truss ori change has been accounted for
      _local_ke.resize(2, 2);
      _local_ke.zero();
      double off_stiff_elem;
      if(coupled_component == 3)
      {
        off_stiff_elem = ori_ij(_component) * _bond_dfdT_ij[0];
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            _local_ke(i, j) += (i == 1 ? 1 : -1) * off_stiff_elem;
      }
      else if(coupled_component == 4)
      {
        off_stiff_elem = ori_ij(_component) * _bond_dfdE_ij[0];
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            _local_ke(i, j) += (i == 1 ? 1 : -1) * off_stiff_elem;
      }
      else
      {
        off_stiff_elem = ori_ij(_component) * ori_ij(coupled_component) * _bond_dfdU_ij[0] - _bond_force_ij[0] * ori_ij(_component) * ori_ij(coupled_component) / current_length_ij;
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            _local_ke(i, j) += (i == j ? 1 : -1) * off_stiff_elem;
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

      _assembly.cacheJacobianBlock(_local_ke, dof_ij, dof_ij_jvar, 1.0);

      // NONLOCAL jacobian contribution
      const NumericVector<Number> & sol = *_nsys.currentSolution();
      std::vector<dof_id_type> dof(2), dof_jvar(2), neighbors;
      // calculation of jacobian contribution to node_i's neighbors
      dof[0] = dof_ij[0];
      dof_jvar[0] = dof_ij_jvar[0];
      neighbors = _pdmesh.neighbors(node_i);
      for (unsigned int k = 0; k < neighbors.size(); ++k)
      {
        Node * nd_k = _mesh.nodePtr(neighbors[k]);
        dof[1] = nd_k->dof_number(_sys.number(), _var.number(), 0);
        dof_jvar[1] = nd_k->dof_number(_sys.number(), jvar, 0);
        double vol_k = _pdmesh.volume(neighbors[k]);
        // obtain bond ik's origin length and current orientation
        double origin_length = 0, current_length = 0;
        RealGradient ori_ik(3);
        for (unsigned int j = 0; j < _pddim; ++j)
        {
          origin_length += std::pow(_pdmesh.coord(node_i)(j) - _pdmesh.coord(neighbors[k])(j), 2);
          ori_ik(j) = _pdmesh.coord(neighbors[k])(j) + sol(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_i)(j) - sol(nd_i->dof_number(_nsys.number(), j, 0));
        }
        origin_length = std::sqrt(origin_length);
        current_length = ori_ik.size();
        ori_ik /= current_length;
        _local_ke.zero();
        double off_stiff_elem1, off_stiff_elem2;
        if (coupled_component == 3)
        {
          off_stiff_elem1 = ori_ik(_component) * _bond_dfdT_i_j[0];
          off_stiff_elem2 = 0;
          _local_ke(0, 0) += - off_stiff_elem1 / origin_length * vol_k;
          _local_ke(1, 0) += off_stiff_elem1 / origin_length * vol_k;
        }
        else if (coupled_component == 4)
        {
          off_stiff_elem1 = ori_ik(_component) * _bond_dfdE_i_j[0];
          off_stiff_elem2 = 0;
          _local_ke(0, 0) += - off_stiff_elem1 / origin_length * vol_k;
          _local_ke(1, 0) += off_stiff_elem1 / origin_length * vol_k;
        }
        else
        {
          off_stiff_elem1 = ori_ik(_component) * ori_ij(coupled_component) * _bond_dfdU_i_j[0];
          off_stiff_elem2 = - _bond_force_i_j[0] * ori_ik(_component) * ori_ik(coupled_component) / current_length;
          for (unsigned int i = 0; i < _test.size(); ++i)
            for (unsigned int j = 0; j < _phi.size(); ++j)
              _local_ke(i, j) += (i == j ? 1 : -1) * off_stiff_elem1 / origin_length * vol_k;
        }

        _assembly.cacheJacobianBlock(_local_ke, dof, dof_ij_jvar, 1.0);

        _local_ke.zero();
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            _local_ke(i, j) += (i == j ? 1 : -1) * off_stiff_elem2 / origin_length * vol_k;

        _assembly.cacheJacobianBlock(_local_ke, dof, dof_jvar, 1.0);
      }

      // calculation of jacobian contribution to node_j's neighbors
      dof[1] = dof_ij[1];
      dof_jvar[1] = dof_ij_jvar[1];
      neighbors = _pdmesh.neighbors(node_j);
      for (unsigned int k = 0; k < neighbors.size(); ++k)
      {
        Node * nd_k = _mesh.nodePtr(neighbors[k]);
        dof[0] = nd_k->dof_number(_sys.number(), _var.number(), 0);
        dof_jvar[0] = nd_k->dof_number(_sys.number(), jvar, 0);
        double vol_k = _pdmesh.volume(neighbors[k]);
        // obtain bond ik's origin length and current orientation
        double origin_length = 0, current_length = 0;
        RealGradient ori_jk(3);
        for (unsigned int j = 0; j < _pddim; ++j)
        {
          origin_length += std::pow(_pdmesh.coord(node_j)(j) - _pdmesh.coord(neighbors[k])(j), 2);
          ori_jk(j) = _pdmesh.coord(neighbors[k])(j) + sol(nd_k->dof_number(_nsys.number(), j, 0)) - _pdmesh.coord(node_j)(j) - sol(nd_j->dof_number(_nsys.number(), j, 0));
        }
        origin_length = std::sqrt(origin_length);
        current_length = ori_jk.size();
        ori_jk /= current_length;
        _local_ke.zero();
        double off_stiff_elem1, off_stiff_elem2;
        if (coupled_component == 3)
        {
          off_stiff_elem1 = ori_jk(_component) * _bond_dfdT_i_j[1];
          off_stiff_elem2 = 0;
          _local_ke(0, 1) += off_stiff_elem1 / origin_length * vol_k;
          _local_ke(1, 1) += - off_stiff_elem1 / origin_length * vol_k;
        }
        else if (coupled_component == 4)
        {
          off_stiff_elem1 = ori_jk(_component) * _bond_dfdE_i_j[1];
          off_stiff_elem2 = 0;
          _local_ke(0, 1) += off_stiff_elem1 / origin_length * vol_k;
          _local_ke(1, 1) += - off_stiff_elem1 / origin_length * vol_k;
        }
        else
        {
          off_stiff_elem1 = - ori_jk(_component) * ori_ij(coupled_component) * _bond_dfdU_i_j[1];
          off_stiff_elem2 = - _bond_force_i_j[1] * ori_jk(_component) * ori_jk(coupled_component) / current_length;
          for (unsigned int i = 0; i < _test.size(); ++i)
            for (unsigned int j = 0; j < _phi.size(); ++j)
              _local_ke(i, j) += (i == j ? 1 : -1) * off_stiff_elem1 / origin_length * vol_k;
        }

        _assembly.cacheJacobianBlock(_local_ke, dof, dof_ij_jvar, 1.0);

        _local_ke.zero();
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            _local_ke(i, j) += (i == j ? 1 : -1) * off_stiff_elem2 / origin_length * vol_k;

        _assembly.cacheJacobianBlock(_local_ke, dof, dof_jvar, 1.0);
      }
    }
  }
}
