CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes      !   num_elem      
num_el_blk        num_node_sets         num_side_sets         num_el_in_blk1        num_nod_per_el1       num_side_ss1   
   num_side_ss2      num_side_ss3      num_side_ss4   
   num_nod_ns1       num_nod_ns2       num_nod_ns3       num_nod_ns4       num_nod_ns5       num_nod_ns6       num_nod_ns7       num_nod_var       num_info  �         api_version       @�
=   version       @�
=   floating_point_word_size            	file_size               int64_status             title         bending_jacobian_out.e     maximum_name_length                 #   
time_whole                            ��   	eb_status                             	�   eb_prop1               name      ID              	�   	ns_status         	                    	�   ns_prop1      	         name      ID              	�   	ss_status         
                    	�   ss_prop1      
         name      ID              
   coordx                           
   coordy                               eb_names                       $      (   ns_names      	                 �      L   ss_names      
                 �      4   
coor_names                         D      �   node_num_map                    �      �   connect1                  	elem_type         QUAD4        @      �   elem_num_map                    P      �   elem_ss1                    (         side_ss1                    (      8   elem_ss2                          `   side_ss2                          h   elem_ss3                          p   side_ss3                          x   elem_ss4                    (      �   side_ss4                    (      �   node_ns1                    ,      �   node_ns2                          �   node_ns3                             node_ns4                    ,         node_ns5                          8   node_ns6                          <   node_ns7                          H   vals_nod_var1                               ��   vals_nod_var2                               ��   name_nod_var                       D      L   info_records                      |      �                                      e      g      f                                 ?�      ?�              @       @       @      @      @      @      @      @      @      @      @      @      @       @       @"      @"      @$      @$      ?�              @       @      @      @      @      @      @       @"      @$                      ?�      ?�              ?�              ?�              ?�              ?�              ?�              ?�              ?�              ?�              ?�      @       @       @       @       @       @       @       @       @       @       @                                           bottom                           left                             101                              top                              103                              right                            102                               bottom                           left                             right                            top                                                                                                                             	   
                                                                      !                                          	   
      	         
                                                                                                      
         
                                                                   !                               	   
                                       	                     
                                             
                                                                                    	                                                             !         !         disp_x                           disp_y                             ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                   -i                                                                               bending_jacobian.i                                                                                                                                                ### Version Info ###                                                             Framework Information:                                                           MOOSE version:           git commit 5048397 on 2016-06-21                        PETSc Version:           3.5.3                                                   Current Time:            Tue Jun 21 17:41:56 2016                                Executable Timestamp:    Tue Jun 21 17:41:35 2016                                                                                                                                                                                                  ### Input File ###                                                                                                                                                []                                                                                 initial_from_file_timestep     = LATEST                                          initial_from_file_var          = INVALID                                         block                          = INVALID                                         coord_type                     = XYZ                                             fe_cache                       = 0                                               kernel_coverage_check          = 1                                               material_coverage_check        = 1                                               name                           = 'MOOSE Problem'                                 restart_file_base              = INVALID                                         rz_coord_axis                  = Y                                               type                           = FEProblem                                       use_legacy_uo_aux_computation  = INVALID                                         use_legacy_uo_initialization   = INVALID                                         element_order                  = AUTO                                            order                          = AUTO                                            side_order                     = AUTO                                            control_tags                   = INVALID                                         dimNearNullSpace               = 0                                               dimNullSpace                   = 0                                               enable                         = 1                                               error_on_jacobian_nonzero_reallocation = 0                                       force_restart                  = 0                                               petsc_inames                   =                                                 petsc_options                  = INVALID                                         petsc_values                   =                                                 solve                          = 1                                               use_nonlinear                  = 1                                             []                                                                                                                                                                [BCs]                                                                                                                                                               [./fix_corner_x]                                                                   boundary                     = 101                                               control_tags                 = INVALID                                           enable                       = 1                                                 implicit                     = 1                                                 type                         = PresetBC                                          use_displaced_mesh           = 0                                                 variable                     = disp_x                                            diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 0                                               [../]                                                                                                                                                             [./fix_corner_y]                                                                   boundary                     = 101                                               control_tags                 = INVALID                                           enable                       = 1                                                 implicit                     = 1                                                 type                         = PresetBC                                          use_displaced_mesh           = 0                                                 variable                     = disp_y                                            diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 0                                               [../]                                                                                                                                                             [./fix_y]                                                                          boundary                     = 102                                               control_tags                 = INVALID                                           enable                       = 1                                                 implicit                     = 1                                                 type                         = PresetBC                                          use_displaced_mesh           = 0                                                 variable                     = disp_y                                            diag_save_in                 = INVALID                                           save_in                      = INVALID                                           seed                         = 0                                                 value                        = 0                                               [../]                                                                                                                                                             [./move_y]                                                                         boundary                     = 103                                               control_tags                 = INVALID                                           enable                       = 1                                                 implicit                     = 1                                                 type                         = FunctionPresetBC                                  use_displaced_mesh           = 0                                                 variable                     = disp_y                                            diag_save_in                 = INVALID                                           function                     = -t                                                save_in                      = INVALID                                           seed                         = 0                                               [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      type                           = Transient                                       abort_on_solve_fail            = 0                                               compute_initial_residual_before_preset_bcs = 0                                   control_tags                   =                                                 dt                             = 0.1                                             dtmax                          = 1e+30                                           dtmin                          = 0.1                                             enable                         = 1                                               end_time                       = 1e+30                                           l_abs_step_tol                 = -1                                              l_max_its                      = 50                                              l_tol                          = 0.0001                                          line_search                    = default                                         max_xfem_update                = 4294967295                                      n_startup_steps                = 0                                               nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-50                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 10                                              nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-10                                           no_fe_reinit                   = 0                                               num_steps                      = 2                                               petsc_options                  = INVALID                                         petsc_options_iname            = _PC_TYPE                                        petsc_options_value            = lu                                              picard_abs_tol                 = 1e-50                                           picard_max_its                 = 1                                               picard_rel_tol                 = 1e-08                                           reset_dt                       = 0                                               restart_file_base              =                                                 scheme                         = INVALID                                         solve_type                     = NEWTON                                          splitting                      = INVALID                                         ss_check_tol                   = 1e-08                                           ss_tmin                        = 0                                               start_time                     = 0                                               time_period_ends               = INVALID                                         time_period_starts             = INVALID                                         time_periods                   = INVALID                                         timestep_tolerance             = 2e-14                                           trans_ss_check                 = 0                                               use_multiapp_dt                = 0                                               verbose                        = 0                                             []                                                                                                                                                                [Executioner]                                                                      _fe_problem                    = 0x7f968301ae00                                []                                                                                                                                                                [Kernels]                                                                                                                                                           [./TensorMechanics]                                                                base_name                    = INVALID                                           block                        = INVALID                                           diag_save_in                 = INVALID                                           displacements                = 'disp_x disp_y'                                   save_in                      = INVALID                                           temp                         = INVALID                                           use_displaced_mesh           = 1                                                 use_finite_deform_jacobian   = 1                                               [../]                                                                          []                                                                                                                                                                [Materials]                                                                                                                                                         [./elasticity_tensor]                                                              type                         = ComputeElasticityTensor                           C_ijkl                       = '168400 121400 121400 168400 121400 168400...  75400 75400 75400'                                                                  base_name                    = INVALID                                           block                        = INVALID                                           boundary                     = INVALID                                           compute                      = 1                                                 control_tags                 = Materials                                         elasticity_tensor_prefactor  = INVALID                                           enable                       = 1                                                 euler_angle_1                = 0                                                 euler_angle_2                = 0                                                 euler_angle_3                = 0                                                 fill_method                  = symmetric9                                        implicit                     = 1                                                 output_properties            = INVALID                                           outputs                      = none                                              seed                         = 0                                                 use_displaced_mesh           = 0                                               [../]                                                                                                                                                             [./strain]                                                                         type                         = ComputeFiniteStrain                               base_name                    = INVALID                                           block                        = INVALID                                           boundary                     = INVALID                                           compute                      = 1                                                 control_tags                 = Materials                                         displacements                = 'disp_x disp_y'                                   enable                       = 1                                                 implicit                     = 1                                                 output_properties            = INVALID                                           outputs                      = none                                              seed                         = 0                                                 temperature                  = 273                                               temperature_ref              = 273                                               thermal_expansion_coeff      = 0                                                 use_displaced_mesh           = 0                                                 volumetric_locking_correction = 0                                              [../]                                                                                                                                                             [./stress]                                                                         type                         = ComputeFiniteStrainElasticStress                  base_name                    = INVALID                                           block                        = INVALID                                           boundary                     = INVALID                                           compute                      = 1                                                 control_tags                 = Materials                                         enable                       = 1                                                 implicit                     = 1                                                 initial_stress               = INVALID                                           output_properties            = INVALID                                           outputs                      = none                                              seed                         = 0                                                 store_stress_old             = 0                                                 use_displaced_mesh           = 0                                               [../]                                                                          []                                                                                                                                                                [Mesh]                                                                             displacements                  = 'disp_x disp_y'                                 block_id                       = INVALID                                         block_name                     = INVALID                                         boundary_id                    = INVALID                                         boundary_name                  = INVALID                                         construct_side_list_from_node_list = 0                                           ghosted_boundaries             = INVALID                                         ghosted_boundaries_inflation   = INVALID                                         patch_size                     = 40                                              second_order                   = 0                                               skip_partitioning              = 0                                               type                           = GeneratedMesh                                   uniform_refine                 = 0                                               bias_x                         = 1                                               bias_y                         = 1                                               bias_z                         = 1                                               centroid_partitioner_direction = INVALID                                         construct_node_list_from_side_list = 1                                           control_tags                   =                                                 dim                            = 2                                               distribution                   = DEFAULT                                         elem_type                      = QUAD4                                           enable                         = 1                                               nemesis                        = 0                                               nx                             = 10                                              ny                             = 2                                               nz                             = 1                                               parallel_type                  = DEFAULT                                         partitioner                    = default                                         patch_update_strategy          = never                                           xmax                           = 10                                              xmin                           = 0                                               ymax                           = 2                                               ymin                           = 0                                               zmax                           = 1                                               zmin                           = 0                                             []                                                                                                                                                                [Mesh]                                                                           []                                                                                                                                                                [MeshModifiers]                                                                                                                                                     [./corner]                                                                         type                         = AddExtraNodeset                                   _mesh                        = 0x7f9683017e00                                    control_tags                 = MeshModifiers                                     coord                        = '0 0'                                             depends_on                   = INVALID                                           enable                       = 1                                                 force_prepare                = 0                                                 new_boundary                 = 101                                               nodes                        = INVALID                                           tolerance                    = 1e-06                                           [../]                                                                                                                                                             [./mid]                                                                            type                         = AddExtraNodeset                                   _mesh                        = 0x7f9683017e00                                    control_tags                 = MeshModifiers                                     coord                        = '5 2'                                             depends_on                   = INVALID                                           enable                       = 1                                                 force_prepare                = 0                                                 new_boundary                 = 103                                               nodes                        = INVALID                                           tolerance                    = 1e-06                                           [../]                                                                                                                                                             [./side]                                                                           type                         = AddExtraNodeset                                   _mesh                        = 0x7f9683017e00                                    control_tags                 = MeshModifiers                                     coord                        = '10 0'                                            depends_on                   = INVALID                                           enable                       = 1                                                 force_prepare                = 0                                                 new_boundary                 = 102                                               nodes                        = INVALID                                           tolerance                    = 1e-06                                           [../]                                                                          []                                                                                                                                                                [Outputs]                                                                          append_date                    = 0                                               append_date_format             = INVALID                                         checkpoint                     = 0                                               color                          = 1                                               console                        = 1                                               controls                       = 0                                               csv                            = 0                                               dofmap                         = 0                                               execute_on                     = 'INITIAL TIMESTEP_END'                          exodus                         = 1                                               file_base                      = INVALID                                         gmv                            = 0                                               gnuplot                        = 0                                               hide                           = INVALID                                         interval                       = 1                                               nemesis                        = 0                                               output_if_base_contains        = INVALID                                         print_linear_residuals         = 1                                               print_mesh_changed_info        = 0                                               print_perf_log                 = 0                                               show                           = INVALID                                         solution_history               = 0                                               sync_times                     =                                                 tecplot                        = 0                                               vtk                            = 0                                               xda                            = 0                                               xdr                            = 0                                             []                                                                                                                                                                [Preconditioning]                                                                                                                                                   [./smp]                                                                            line_search                  = default                                           petsc_options                = INVALID                                           petsc_options_iname          = INVALID                                           petsc_options_value          = INVALID                                           solve_type                   = INVALID                                           type                         = SMP                                               control_tags                 = Preconditioning                                   coupled_groups               = INVALID                                           enable                       = 1                                                 full                         = 1                                                 off_diag_column              = INVALID                                           off_diag_row                 = INVALID                                           pc_side                      = right                                           [../]                                                                          []                                                                                                                                                                [Variables]                                                                                                                                                         [./disp_x]                                                                         block                        = INVALID                                           eigen                        = 0                                                 family                       = LAGRANGE                                          initial_condition            = INVALID                                           order                        = FIRST                                             outputs                      = INVALID                                           scaling                      = 1                                                 initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                         [../]                                                                                                                                                             [./disp_y]                                                                         block                        = INVALID                                           eigen                        = 0                                                 family                       = LAGRANGE                                          initial_condition            = INVALID                                           order                        = FIRST                                             outputs                      = INVALID                                           scaling                      = 1                                                 initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                         [../]                                                                          []                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ?�������        ?c3���)�?�����?��� ��?w@LQI��?�l� �?�٧"�]�?�����?��]7x�k?��!�H�L?��,���?��,��r?�+��π�?����M�?�Rsc�_�?���#]�?� Ӣb��?���j�S�?�՞���?�!�V/)?��,��T?�[���!�?�8(�e�?��ǈ��?�~�ԳZ?���
��?��lS�e?��,��o?����?����Q�{?|V��� `?m�9Y�?d1ZH���        ��s��ygJ��b�Bd<��o8��χ���� ��ѿ��۪8��ݍ�sk�����tt^q��qC�-�(��c�A��忸��[ډ���\�8�¿�qC�-�"��c�A��߿�ݍ�sk�����tt^k��� ��Ŀ��۪"��s��yg'��b�Bd<�        �o8��χm��H��҉��t+���PT��'&:j⫿���t�s�����^ؿ������������^տ���t�o��'&:j⠿�H��҉��t+���P<?ə�����        ?o|^*ڈ�?�p�l� �?��݁��?�;1t��?��*c�u[?��N�kN?�hPk�v?���D?�?���qL?�?��y�W�Z?��y�W�[?��L���?�3c2�?���>���?�܋TC�9?�p��C�?�X�@��V?���� �
?�~7�Q�?��y�W�V?�L#��?��	>n�?��Z�L|�?�D슔P?��(��J?��tQ�w?��y�W�L?��ܼk�?�C�I?���.i(?p͜���?S����        ��eh����ü7��4�,��5������B���ĺ�������K�Ċ+�4����m.���b�����Ȫ�R�u���m���G��m.�	��b���������K�Ċ+�4��������:���ĺ����eh����ü7        ��4�,��-����򧿇'w�?d��``/��� �
ǐt�Ǩ���x�ə������Ǩ���x�� �
ǐn��``/�����򟿇'w�?8