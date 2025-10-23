import sympy as sp
import numpy as np
import scipy as sc
import sympy.physics.mechanics as me
from scipy.spatial.transform import Rotation as spat
import pickle

def objective_activation_diff(num_nodes, interval_value):
    node_val = sp.Matrix(sp.symbols('act0:' + str(num_nodes)))
    obj_list = []

    for i in range(num_nodes-1):
        act_diff = node_val[i+1] - node_val[i]
        obj_list.append(interval_value * (act_diff**2))

    obj = sum(obj_list)
    obj_jac = sp.Matrix([obj]).jacobian(node_val)
    obj_np = sp.lambdify([node_val],obj)
    obj_jac_np = sp.lambdify([node_val],obj_jac)

    return obj_np,obj_jac_np


def max_GH_stab(num_reactions,tilt_y,tilt_z,w_rx):
    # reactions = sp.symbols('r1:' + str(num_reactions+1))
    # reactions_mag = sp.sqrt(sum([reactions[i]**2 for i in range(num_reactions)]))
    # 
    reactions = sp.symbols('rx ry rz')
    reactions_rotated = R_y(tilt_y*sp.pi/180).T * R_z(tilt_z*sp.pi/180).T * sp.Matrix([*reactions,1])

    obj = (w_rx * reactions_rotated[0]**2 + (1-w_rx) * (reactions_rotated[1]**2 + reactions_rotated[2]**2)) # + 0.02 * sp.sqrt(reactions[0]**2 + reactions[1]**2 + reactions[2]**2)
    # y_comp = sp.atan2(-reactions_rotated[1],-reactions_rotated[0])
    # z_comp = sp.atan2(-reactions_rotated[2],-reactions_rotated[0])
    # obj = 0.01 * (y_comp**2 + z_comp**2)
    # print(obj)

    obj_jac = sp.Matrix([obj]).jacobian(reactions)[:]
    # print(obj_jac)
    obj_np = sp.lambdify([reactions],obj)
    obj_jac_np = sp.lambdify([reactions],obj_jac)
    return obj_np, obj_jac_np

def optimize_parameters(num_params,mus_group,optimize,w_fmax,w_lceopt):

    ipar = 0
    mus_dict = {}
    full_obj = 0

    # do fmax first
    if optimize == 'all':
        for igroup in range(len(mus_group)):
            group_obj = 0
            group_pars = []
            for imus in range(len(mus_group[igroup])):
                current_par = sp.Symbol('p' + str(ipar+1))
                group_pars.append(current_par)
                ipar += 1
                mus_dict.update({'fmax_' + mus_group[igroup][imus]: current_par})
                group_obj += (group_pars[imus])**2
                if imus > 0:
                    group_obj += (group_pars[imus] - group_pars[imus-1])**2

            full_obj += group_obj

    elif optimize == 'one':
        fmax_scaler = sp.Symbol('fmax_scaler')
        mus_dict.update({'fmax_scaler': fmax_scaler})
        full_obj += w_fmax * (fmax_scaler)**2

    elif optimize == 'group':
        for igroup in range(len(mus_group)):
            current_group = mus_group[igroup][0][0:-1]
            fmax_scaler = sp.Symbol('fmax_scaler_' + current_group)
            mus_dict.update({'fmax_scaler_' + current_group: fmax_scaler})
            full_obj += w_fmax * ((1-fmax_scaler)**2)

        # then do lceopt
        for igroup in range(len(mus_group)):
            current_group = mus_group[igroup][0][0:-1]
            lceopt_scaler = sp.Symbol('lceopt_scaler_' + current_group)
            mus_dict.update({'lceopt_scaler_' + current_group: lceopt_scaler})
            full_obj += w_lceopt * (1-lceopt_scaler)**2

    # for igroup in range(len(mus_group)):
    #     group_obj = 0
    #     group_pars = []
    #     for imus in range(len(mus_group[igroup])):
    #         current_par = sp.Symbol('p' + str(ipar+1))
    #         group_pars.append(current_par)
    #         ipar += 1
    #         mus_dict.update({'lceopt_' + mus_group[igroup][imus]: current_par})
    #         group_obj += w_lceopt * (1 - group_pars[imus])**2
    #         if imus > 0:
    #             group_obj += (group_pars[imus] - group_pars[imus-1])**2

    #     full_obj += group_obj


    
    


    myKeys = list(mus_dict.keys())
    myKeys.sort()
    sorted = {i: mus_dict[i] for i in myKeys}
    sorted_params = [*sorted.values()]
    print(myKeys)
    print(sorted_params)

    # params = sp.symbols('p1:' + str(num_params+1))
    # obj = sum([((params[i] - init_values[i]))**2 for i in range(num_params)])
    print(full_obj)
    obj_jac = sp.Matrix([full_obj]).jacobian(sorted_params)
    obj_np = sp.lambdify(sorted_params,full_obj)
    obj_jac_np = sp.lambdify(sorted_params,obj_jac)

    return obj_np, obj_jac_np

        
    # params = sp.symbols('p1:' + str(num_params+1))
    # obj = sum([((params[i] - init_values[i]))**2 for i in range(num_params)])
    # print(obj)
    # obj_jac = sp.Matrix([obj]).jacobian(params)
    # obj_np = sp.lambdify(params,obj)
    # obj_jac_np = sp.lambdify(params,obj_jac)

    # return obj_np, obj_jac_np


def min_activation(acts,interval_value):

    a_str = sorted([str(acts[x]).replace('(t)','') for x in range(len(acts))])
    muscles_emg_optim = [['act_47','act_48'],['act_44','act_45','act_46'],['act_38','act_39','act_40'],['act_56','act_57'],['act_1','act_2'],['act_6','act_7','act_8'],['act_12'],['act_26','act_27','act_28']]
    a_emg = sp.symbols(f'a_emg1:{len(muscles_emg_optim)+1}')
    weight_index = sp.Symbol('weight_index')
    w_min_squared = sp.Symbol('w_min_squared')
    w_min_emg = sp.Symbol('w_min_emg')
    obj_emg = 0
    a = sp.symbols(a_str)
    for j,imus in enumerate(muscles_emg_optim):
        mus_indeces = [a_str.index(mus) for mus in imus]
        # obj_emg += 10 * interval_value * ((a[mus_indeces[0]] + a[mus_indeces[1]])/2 - a_emg)**2
        # obj_emg += 10 * interval_value * ((sum([a[mus_indeces[i]] for i in range(len(mus_indeces))]))/len(mus_indeces) - a_emg)**2
        mus_current = 0
        for i in range(len(mus_indeces)):
            mus_current += a[mus_indeces[i]]
        # mus_current += [a[mus_indeces[i]] for i in range(len(mus_indeces))]
        obj_emg += (mus_current/len(mus_indeces)-a_emg[j])**2 * interval_value
    
    obj_a = interval_value*sum(sp.Matrix(a).applyfunc(lambda x:x**2))
    
    obj = sp.Matrix([obj_a * w_min_squared + weight_index * obj_emg * w_min_emg])
    # print(obj)
    obj_jac = (obj).jacobian(a)[:]
    obj_np = sp.lambdify((a,a_emg, weight_index, w_min_squared, w_min_emg),obj)
    obj_jac_np = sp.lambdify((a,a_emg, weight_index, w_min_squared, w_min_emg),obj_jac)

    return obj_np, obj_jac_np

# def activation_constraint(acts, interval_value, num_nodes, EMG, EMG_indexes):
#     muscles_emg_optim = [[47,48],[44,45,46],[38,39,40],[56,57],[1,2],[6,7,8],[12],[26,27,28]]
#     constraints = []
#     for j,imus in enumerate(muscles_emg_optim):
#         for inode in range(num_nodes):
#             constraint = 0
#             for i in range(len(imus)):
#                 constraint += acts[imus][i].func(interval_value*inode)
#             constraints.append(constraint-EMG[inode,])

def custom_objective_quat(num_coords,interval_value, clav_pos,include_quat_norm = False):
    x = sp.Matrix(sp.symbols(f'x1:{num_coords+1}'))
    x_traj = sp.Matrix(sp.symbols(f'x_traj1:{num_coords+1}'))

    SCrot = sp.Matrix(x[0:4])
    scapula_thorax = sp.Matrix([mulQuat_sp(x[0:4],x[4:8])])
    humerus_thorax = sp.Matrix([mulQuat_sp(scapula_thorax,x[8:12])])
    AC_pos = Qrm(x[0:4]) * position([clav_pos,0,0])
    # GH_rot = sp.Matrix(x[8:12])
    elrot = sp.Matrix(x[12:13])

    # obj_SCrot = sp.Matrix([interval_value*sum((SCrot-sp.Matrix(x_traj[0:4])).applyfunc(lambda x: x**2))])
    obj_SC_wo_x = sp.Matrix([interval_value*sum((sp.Matrix(AC_pos[0:3])-sp.Matrix(x_traj[0:3])).applyfunc(lambda x: x**2))])
    # obj_scapula_thorax = sp.Matrix([interval_value*sum((sp.Matrix(scapula_thorax[0:4])-sp.Matrix(x_traj[4:8])).applyfunc(lambda x: x**2))])
    obj_scapula_thorax = sp.Matrix([interval_value*sum(((mulQuat_sp(scapula_thorax,invQuat_sp(x_traj[4:8])))).applyfunc(lambda x: x**2)[1:4])])
    # obj_GH_rot = sp.Matrix([interval_value*sum((GH_rot-sp.Matrix(x_traj[8:12])).applyfunc(lambda x: x**2))])
    # obj_humerus_thorax = sp.Matrix([interval_value*sum((sp.Matrix(humerus_thorax[0:4])-sp.Matrix(x_traj[8:12])).applyfunc(lambda x: x**2))])
    obj_humerus_thorax = sp.Matrix([interval_value*sum(((mulQuat_sp(humerus_thorax,invQuat_sp(x_traj[8:12])))).applyfunc(lambda x: x**2)[1:4])])


    obj_elrot = sp.Matrix([interval_value*sum((elrot-sp.Matrix(x_traj[12:13])).applyfunc(lambda x: x**2))])
    # obj = obj_SC_wo_x + obj_scapula_thorax + obj_humerus_thorax + obj_elrot
    obj = obj_humerus_thorax + obj_elrot
    if include_quat_norm:
        SC_norm = sp.Matrix([x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2-1])**2
        AC_norm = sp.Matrix([x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2-1])**2
        GH_norm = sp.Matrix([x[8]**2 + x[9]**2 + x[10]**2 + x[11]**2-1])**2
        obj += (SC_norm + AC_norm + GH_norm)*interval_value

    # obj = obj_SCrot + obj_scapula_thorax + obj_GH_rot + obj_elrot
    
    obj_jac = (obj).jacobian(x)[:]
    obj_np = sp.lambdify((x,x_traj),obj)
    obj_jac_np = sp.lambdify((x,x_traj),obj_jac)

    return obj_np,obj_jac_np

def custom_objective_eul(num_coords,interval_value,GH_seq = 'YZY'):
    x = sp.symbols(f'x1:{num_coords+1}')
    x_traj = sp.symbols(f'x_traj1:{num_coords+1}')

    SCrot = sp.Matrix(x[0:2])
    scapula_thorax_RM = R_y(x[0]) * R_z(x[1]) * R_x(x[2]) * R_y(x[3]) * R_z(x[4]) * R_x(x[5])
    scapula_thorax_z = sp.asin(scapula_thorax_RM[1,0])
    scapula_thorax_x = sp.atan2(-scapula_thorax_RM[1,2],scapula_thorax_RM[1,1])
    scapula_thorax_y = sp.atan2(-scapula_thorax_RM[2,0],scapula_thorax_RM[0,0])
    scapula_thorax = sp.Matrix([scapula_thorax_x,scapula_thorax_y,scapula_thorax_z])

    if GH_seq == 'YZX':
        humerus_thorax_RM = scapula_thorax_RM * R_y(x[6]) * R_z(x[7]) * R_x(x[8])
        humerus_thorax_z = sp.asin(humerus_thorax_RM[1,0])
        humerus_thorax_x = sp.atan2(-humerus_thorax_RM[1,2],humerus_thorax_RM[1,1])
        humerus_thorax_y = sp.atan2(-humerus_thorax_RM[2,0],humerus_thorax_RM[0,0])
        humerus_thorax = sp.Matrix([humerus_thorax_y,humerus_thorax_z,humerus_thorax_x])
    elif GH_seq == 'YZY':
        humerus_thorax_RM = scapula_thorax_RM * R_y(x[6]) * R_z(x[7]) * R_y(x[8])
        humerus_thorax_z = sp.acos(humerus_thorax_RM[1,1])
        humerus_thorax_yy = sp.atan2(humerus_thorax_RM[1,2],humerus_thorax_RM[1,0])
        humerus_thorax_y = sp.atan2(humerus_thorax_RM[2,1],-humerus_thorax_RM[0,1])
        humerus_thorax = sp.Matrix([humerus_thorax_y,humerus_thorax_z,humerus_thorax_yy])
    
    elrot = sp.Matrix(x[9:10])

    obj_SCrot = 1 * sp.Matrix([interval_value*sum((SCrot-sp.Matrix(x_traj[0:2])).applyfunc(lambda x: x**2))])
    obj_scapula_thorax = sp.Matrix([interval_value*sum((scapula_thorax-sp.Matrix(x_traj[3:6])).applyfunc(lambda x: x**2))])
    # obj_GH_rot = sp.Matrix([interval_value*sum((GH_rot-sp.Matrix(x_traj[6:9])).applyfunc(lambda x: x**2))])
    obj_humerus_thorax = sp.Matrix([interval_value*sum((humerus_thorax-sp.Matrix(x_traj[6:9])).applyfunc(lambda x: x**2))])
    obj_elrot = sp.Matrix([interval_value*sum((elrot-sp.Matrix(x_traj[9:10])).applyfunc(lambda x: x**2))])
    # obj = obj_SCrot + obj_scapula_thorax + obj_GH_rot + obj_elrot

    obj = obj_SCrot + obj_scapula_thorax + obj_humerus_thorax + obj_elrot
    obj_jac = (obj).jacobian(x)[:]
    obj_np = sp.lambdify((x,x_traj),obj)
    obj_jac_np = sp.lambdify((x,x_traj),obj_jac)

    return obj_np,obj_jac_np

def polynomials_euler(model_struct,q,u,derive, motion_folder = None, gen_matlab_functions = None):

    q_thorax1 = me.dynamicsymbols('q_thorax1')
    q_thorax2 = me.dynamicsymbols('q_thorax2')
    q_thorax3 = me.dynamicsymbols('q_thorax3')
    q_radius = (me.dynamicsymbols('q_radius'))

    q_new = q.copy()
    q_new.append(q_radius)
    q_new.insert(0,q_thorax3)
    q_new.insert(0,q_thorax2)
    q_new.insert(0,q_thorax1)

    qpol = q_new
    nmus = len(model_struct['model']['muscles'][0,0][0])
    
    actSym = sp.symbols('actSym_1:' + str(nmus + 1))
    
    if derive == 'symbolic':
        fmax = sp.symbols('fmax_1:' + str(nmus + 1))
        lceopt = sp.symbols('lceopt_1:' + str(nmus + 1))
        lslack = sp.symbols('lslack_1:' + str(nmus + 1))
        
        muscle_constants = {'fmax': fmax,
                           'lceopt': lceopt,
                           'lslack': lslack}
    elif derive == 'numeric':
        fmax = []
        lceopt = []
        lslack = []
        vmax = []
        muscle_constants = {}
        
        for i in range(nmus):
            # fmax.append(model_params_struct['params'][initCond_name][0,0]['fmax'][0,0][i,0].item())
            # lceopt.append(model_params_struct['params'][initCond_name][0,0]['lceopt'][0,0][i,0].item())
            # lslack.append(model_params_struct['params'][initCond_name][0,0]['lslack'][0,0][i,0].item())
            muscle = model_struct['model']['muscles'].item()[0,i]
            fmax.append(muscle['fmax'][0,0].item())
            lceopt.append(muscle['lceopt'][0,0].item())
            lslack.append(muscle['lslack'][0,0].item())
            vmax.append(muscle['vmax'][0,0].item())

        # print(lceopt)
            
    mus_lengths_full = sp.zeros(nmus,1)
    mus_forces = sp.zeros(nmus,1)
    jacobian_full = sp.zeros(11,nmus)
    
    for imus in range(nmus):
        muscle = model_struct['model']['muscles'].item()[0,imus]
        isWrapped = muscle['isWrapped'].item()
        
        if isWrapped == 0:
            origin = muscle['origin_frame'].item()
            insertion = muscle['insertion_frame'].item()
            O_pos = muscle['origin_position'].item()[0]
            I_pos = muscle['insertion_position'].item()[0]
            # print('analytic')
            L = analytic_length_eul(origin, insertion, O_pos, I_pos, qpol[3:], model_struct)
            jacobian_full[:,imus] = -sp.Matrix([L]).jacobian(qpol[3:]).T
            vce = -jacobian_full[:-1,imus].T * sp.Matrix(u)
            
        
        elif isWrapped == 1:
            
            name = muscle['name'].item()
            npolterms = muscle['Euler'][0,0]['lparam_count'].item()
            polcoeff_np = muscle['Euler'][0,0]['lcoefs'].item()
            polcoeff = sp.Matrix(polcoeff_np)
            expon_np = muscle['Euler'][0,0]['lparams'].item()
            expon = sp.Matrix(expon_np)
            musdof = muscle['dof_indeces'].item()
            nmusdof = muscle['dof_count'].item()
            L = 0

            for i in range(npolterms.item()):
             # Add this term's contribution to the muscle length
                term = polcoeff[i]
                for j in range(nmusdof.item()):
                    mdof = musdof[j]
                    for k in range(expon_np[i,j].item()):
                        term = term * qpol[int(mdof-1)]

                L = L + term

            jacobian_full[:,imus] = -sp.Matrix([L]).jacobian(qpol[3:]).T
            vce = -jacobian_full[:-1,imus].T * sp.Matrix(u)

        mus_lengths_full[imus] = L
        mus_forces[imus] = muscle_force(actSym[imus],L,vce,fmax[imus],lceopt[imus],lslack[imus],vmax[imus])

    # jacobian_full = -mus_lengths_full.jacobian(qpol[3:]).T
    FQ = jacobian_full * mus_forces
    TE = sp.Matrix(FQ[:-1])
    
    TE = TE.subs(q_new[13],120*np.pi/180)
    jacobian = jacobian_full.subs(q_new[13],120*np.pi/180)
    mus_lengths = mus_lengths_full.subs(q_new[13],120*np.pi/180)
    mus_forces = mus_forces.subs(q_new[13],120*np.pi/180)
    
    symbols_list = TE.free_symbols
    t = sp.Symbol('t')
    actSym_list = [i for i in symbols_list if str(i).startswith('actSym')]
    
    act = []
    for i in actSym_list:
        act.append(me.dynamicsymbols(str(i).replace('Sym','')))
    act_subs = dict(zip(actSym_list,act))
    TE_act_subbed = me.msubs(TE, act_subs)
    
    # conoid_lopt = model_params_struct['params']['model'].item()['conoid_length'].item()[0][0]
    # conoid_lopt = model_struct['model']['conoid_length'].item()[0][0]
    conoid_lopt = 0.018
    # conoid_k = model_params_struct['params']['model'].item()['conoid_stiffness'].item()[0][0]
    # conoid_k = model_struct['model']['conoid_stiffness'].item()[0][0]
    conoid_k = 1e4
    # conoid_eps = model_params_struct['params']['model'].item()['conoid_eps'].item()[0][0]
    # conoid_eps = model_struct['model']['conoid_eps'].item()[0][0]
    conoid_eps = 1e-3
    # conoid_origin = model_params_struct['params']['model'].item()['conoid_origin'].item()[0]
    conoid_origin = model_struct['model']['conoid_origin'].item()[0]
    print(conoid_origin)
    # conoid_insertion = model_params_struct['params']['model'].item()['conoid_insertion'].item()[0]
    conoid_insertion = model_struct['model']['conoid_insertion'].item()[0]
    print(conoid_insertion)
    conoid_length = analytic_length_eul('clavicle_r','scapula_r', conoid_origin, conoid_insertion, qpol[3:], model_struct)
    F_conoid = conoid_force(conoid_length, conoid_lopt, conoid_k, conoid_eps)
    jac_conoid = -sp.Matrix([conoid_length]).jacobian(qpol[3:]).T
    TE_conoid = F_conoid * jac_conoid

    # jnt_spring = joint_spring(q,model_params_struct,initCond_name)
    
    if gen_matlab_functions == 1:
        qsubs = sp.symbols('qsubs1:12')
        q_subs_dict = dict(zip(qpol[3:],qsubs))
        segments = []
        TE_conoid_subbed = me.msubs(TE_conoid,q_subs_dict)
        TE_subbed = me.msubs(TE,q_subs_dict)
        jacobian_subbed = me.msubs(jacobian_full,q_subs_dict)
        mus_lengths_subbed = me.msubs(mus_lengths_full,q_subs_dict)
        mus_forces_subbed = me.msubs(mus_forces,q_subs_dict)
        # jnt_spring_subbed = me.msubs(jnt_spring,q_subs_dict)
        
        # MatlabFunction(function = TE_conoid_subbed[:-1],
        #                fun_name = 'TE_conoid_eul', assignto = 'TE_conoid',
        #                coordinates = qsubs,
        #                speeds = [],
        #                inputs = [],
        #                body_constants = {},
        #                segments = [],
        #                other_constants = {},
        #                muscle_constants = {},
        #                parameters = [],
        #                folder = 'euler')
        # MatlabFunction(function = TE_subbed,
        #                fun_name = 'TE_eul', assignto = 'TE',
        #                coordinates = qsubs,
        #                speeds = [],
        #                inputs = actSym,
        #                body_constants = {},
        #                segments = [],
        #                other_constants={},
        #                muscle_constants = muscle_constants,
        #                parameters = [],
        #                folder = 'euler')
        # MatlabFunction(function = mus_forces_subbed,
        #                fun_name = 'mus_forces_eul', assignto = 'mus_forces',
        #                coordinates = qsubs,
        #                speeds = [],
        #                inputs = actSym,
        #                body_constants = {},
        #                segments = [],
        #                other_constants={},
        #                muscle_constants = muscle_constants,
        #                parameters = [],
        #                folder = 'euler')
        MatlabFunction(function = jacobian_subbed,
                       fun_name = 'jacobiannoi_eul',assignto = 'jacobian',
                       coordinates = qsubs,
                       speeds = [],
                       inputs = [],
                       body_constants = {},
                       segments = [],
                       other_constants={},
                       muscle_constants = {},
                       parameters = [],
                       folder = '../Motions/' + motion_folder + '/Poly_functions')
        MatlabFunction(function = mus_lengths_subbed,
                       fun_name = 'mus_lengths_eul',assignto = 'mus_lengths',
                       coordinates = qsubs,
                       speeds = [],
                       inputs = [],
                       body_constants = {},
                       segments = [],
                       other_constants={},
                       muscle_constants = {},
                       parameters = [],
                       folder = '../Motions/' + motion_folder + '/Poly_functions')
        # MatlabFunction(function = jnt_spring_subbed,
        #                fun_name = 'jnt_spring_eul',assignto = 'jnt_spring_eul',
        #                coordinates = qsubs,
        #                speeds = [],
        #                inputs = [],
        #                body_constants = {},
        #                segments = [],
        #                other_constants={},
        #                muscle_constants = {},
        #                parameters = [],
        #                folder = 'euler')
    
    return TE_act_subbed, act, TE_conoid[:-1]


def polynomials_quat(model_struct,q,u,derive, optim_params = None, gen_matlab_functions = None):
    
    q_thorax0 = me.dynamicsymbols('q_thorax0')
    q_thorax1 = me.dynamicsymbols('q_thorax1')
    q_thorax2 = me.dynamicsymbols('q_thorax2')
    q_thorax3 = me.dynamicsymbols('q_thorax3')
    q_radius = (me.dynamicsymbols('q_radius'))

    q_new = q.copy()
    q_new.append(q_radius)
    q_new.insert(0,q_thorax3)
    q_new.insert(0,q_thorax2)
    q_new.insert(0,q_thorax1)
    q_new.insert(0,q_thorax0)

    qpol = q_new[1:4]+q_new[5:8]+q_new[9:12]+q_new[13:16]+q_new[16:18] #exclude q0 from quaternions
    nmus = len(model_struct['model']['muscles'][0,0][0])

    dquatdSC = dquatdt(q_new[4:8],u[0:3])
    dquatdAC = dquatdt(q_new[8:12],u[3:6])
    dquatdGH = dquatdt(q_new[12:16],u[6:9])
    dqdtEL = u[9]
    
    actSym = sp.symbols('actSym_1:' + str(nmus + 1))
    
    if derive == 'symbolic':
        fmax = sp.symbols('fmax_1:' + str(nmus + 1))
        lceopt = sp.symbols('lceopt_1:' + str(nmus + 1))
        lslack = sp.symbols('lslack_1:' + str(nmus + 1))
        muscle_constants = {'fmax': fmax,
                           'lceopt': lceopt,
                           'lslack': lslack}
    elif derive == 'numeric':
        fmax = []
        fmax_init = {}
        fmax_range = {}
        lceopt = []
        lceopt_init = {}
        lceopt_range = {}
        lslack = []
        vmax = []
        # optimized_muscles = ['deltscap','deltclav','trapscap','trapclav','serr','rhomb','lev','pec']
        optimized_muscles = ['deltscap','deltclav','trapscap','trapclav','serr','infra','rhomb','lev','pect','tric','termin','subscap','supra']
        mus_groups = [[] for _ in range(len(optimized_muscles))]
        muscle_constants = {}
        fmax_scaler_one = sp.Symbol('fmax_scaler')
        fmax_scaler_group = [sp.Symbol('fmax_scaler_' + i) for i in optimized_muscles]
        lceopt_scaler_group = [sp.Symbol('lceopt_scaler_' + i) for i in optimized_muscles]
        
        for i in range(nmus):
            # fmax.append(model_params_struct['params'][initCond_name][0,0]['fmax'][0,0][i,0].item())
            # lceopt.append(model_params_struct['params'][initCond_name][0,0]['lceopt'][0,0][i,0].item())
            # lslack.append(model_params_struct['params'][initCond_name][0,0]['lslack'][0,0][i,0].item())
            muscle = model_struct['model']['muscles'].item()[0,i]
            
            
            lslack.append(muscle['lslack'][0,0].item())
            vmax.append(muscle['vmax'][0,0].item())
            mus_name = muscle['name'][0,0].item()
            # lceopt.append(muscle['lceopt'][0,0].item())
            
            

            

            
            if optim_params == 'all' or optim_params == 'one' or optim_params == 'group':
                if 'deltscap' in mus_name or 'deltclav' in mus_name or 'trapscap' in mus_name or 'trapclav' in mus_name or 'serr' in mus_name or 'infra' in mus_name or 'rhomb' in mus_name or 'lev' in mus_name or 'pect' in mus_name or 'tric' in mus_name:
                    mus_index = next((i for i, item in enumerate(optimized_muscles) if mus_name.startswith(item)), None)
                    mus_groups[mus_index].append(mus_name)

                    fmax_current = muscle['fmax'][0,0].item()
                    if optim_params == 'all':
                        fmax_scaler = sp.Symbol('fmax_scaler_' + mus_name)
                        fmax_init.update({('fmax_scaler_' + mus_name): 1})
                        fmax.append(fmax_scaler*fmax_current)
                        fmax_range.update({fmax_scaler: (0.3, 2)})
                    elif optim_params == 'one':
                        fmax.append(fmax_scaler_one*fmax_current)
                        fmax_init.update({('fmax_scaler_one'): 1})
                        fmax_range.update({fmax_scaler_one: (0.5, 2)})
                    elif optim_params == 'group':
                        fmax_scaler = fmax_scaler_group[mus_index]
                        fmax.append(fmax_scaler*fmax_current)
                        fmax_init.update({'fmax_scaler_' + optimized_muscles[mus_index]: 1})
                        fmax_range.update({fmax_scaler: (0.3, 1.4)})

                        lceopt_scaler = lceopt_scaler_group[mus_index]
                        lceopt.append(lceopt_scaler*muscle['lceopt'][0,0].item())
                        lceopt_init.update({'lceopt_scaler_' + optimized_muscles[mus_index]: 1})
                        lceopt_range.update({lceopt_scaler: (0.6, 1.4)})

                    # lceopt_current = muscle['lceopt'][0,0].item()
                    # lceopt_init.update({('lceopt_scaler_' + mus_name): 1})
                    # lceopt_symbol = sp.Symbol('lceopt_scaler_' + mus_name)
                    # lceopt.append(lceopt_symbol*lceopt_current)
                    # lceopt_range.update({lceopt_symbol: (0.2, 2)})

                else:
                    fmax.append(muscle['fmax'][0,0].item())
                    lceopt.append(muscle['lceopt'][0,0].item())

            else:
                try:
                    try:
                        mus_index = next((i for i, item in enumerate(optimized_muscles) if mus_name.startswith(item)), None)
                        mus_groups[mus_index].append(mus_name)
                        fmax.append(optim_params['optimized_params']['fmax_scaler_' + optimized_muscles[mus_index]].item()[0][0] * muscle['fmax'][0,0].item())
                        lceopt.append(optim_params['optimized_params']['lceopt_scaler_' + optimized_muscles[mus_index]].item()[0][0] * muscle['lceopt'][0,0].item())
                        # print(mus_name)
                        # fmax.append(optim_params['optimized_params']['fmax_scaler_one'].item()[0][0] * muscle['fmax'][0,0].item())
                    except:
                        fmax.append(optim_params['optimized_params']['fmax_scaler_' + optimized_muscles[mus_index]].item()[0][0] * muscle['fmax'][0,0].item())
                        lceopt.append(optim_params['optimized_params']['lceopt_scaler_' + optimized_muscles[mus_index]].item()[0][0] * muscle['lceopt'][0,0].item())
                    # print('mus')
                    # except:
                    #     fmax.append(optim_params['optimized_params']['fmax_scaler_' + mus_name].item()[0][0] * muscle['fmax'][0,0].item())
                    # lceopt.append(optim_params['optimized_params']['lceopt_scaler_' + mus_name].item()[0][0] * muscle['lceopt'][0,0].item())
                except:
                    fmax.append(muscle['fmax'][0,0].item()*0.7)
                    lceopt.append(muscle['lceopt'][0,0].item())
                


            # print(mus_name)
            # print(fmax[-1])
                    
    print(lceopt)
    print(fmax)
    mus_lengths_full = sp.zeros(nmus,1)
    mus_forces = sp.zeros(nmus,1)
    JacInSpat_full = sp.zeros(11,nmus)
    xforce = 0
    yforce = 0
    zforce = 0
    
    for imus in range(nmus):
    # for imus in range(44,45):
        muscle = model_struct['model']['muscles'].item()[0,imus]
        isWrapped = muscle['isWrapped'].item()

        if isWrapped == 0:
            # if not wrapped, calculate length and moment arms analytically
            origin = muscle['origin_frame'].item()
            insertion = muscle['insertion_frame'].item()
            O_pos = muscle['origin_position'].item()[0]
            I_pos = muscle['insertion_position'].item()[0]
            L = analytic_length_quat(origin, insertion, O_pos, I_pos, q_new[4:], model_struct)
            Jac = -sp.Matrix([L]).jacobian(q_new).T
            vce = -sp.Matrix(Jac[4:-1]).T * sp.Matrix([dquatdSC,dquatdAC,dquatdGH,dqdtEL])
            # vce = 0
                        
            TEsc = 1/2 * G(q_new[4:8])*sp.Matrix(Jac[4:8])
            TEac = 1/2 * G(q_new[8:12])*sp.Matrix(Jac[8:12])
            TEgh = 1/2 * G(q_new[12:16])*sp.Matrix(Jac[12:16])
            TEel = 1/2 * sp.Matrix(Jac[16:18])
            
            iJacInSpat = sp.Matrix(TEsc).col_join(TEac).col_join(TEgh).col_join(TEel)
            JacInSpat_full[:,imus] = iJacInSpat

            # print('analytic')
            
        elif isWrapped == 1:
            
            name = muscle['name'].item()
            # print(name)
            npolterms = muscle['Quaternion'][0,0]['lparam_count'].item()
            polcoeff_np = muscle['Quaternion'][0,0]['lcoefs'].item()
            polcoeff = sp.Matrix(polcoeff_np)
            expon_np = muscle['Quaternion'][0,0]['lparams'].item()
            expon = sp.Matrix(expon_np)
            musdof = muscle['dof_indeces'].item()
            nmusdof = muscle['dof_count'].item()
            L = 0
            

            for i in range(npolterms.item()):
             # Add this term's contribution to the muscle length
                term = polcoeff[i]
                for j in range(nmusdof.item()):
                    mdof = musdof[j];
                    for k in range(expon_np[i,j].item()):
                        term = term * qpol[int(mdof-1)];

                L = L + term;
            jac = -sp.Matrix([L]).jacobian(qpol).T
            vce = -sp.Matrix(jac[3:-1]).T * sp.Matrix([sp.Matrix(dquatdSC[1:4]),sp.Matrix(dquatdAC[1:4]),sp.Matrix(dquatdGH[1:4]),dqdtEL])
            # vce = 0
            TEsc = invJtrans(q_new[4:8])*sp.Matrix(jac[3:6])
            TEac = invJtrans(q_new[8:12])*sp.Matrix(jac[6:9])
            TEgh = invJtrans(q_new[12:16])*sp.Matrix(jac[9:12])
            TEel = sp.Matrix(jac[12:14])
            iJacInSpat = sp.Matrix(TEsc).col_join(TEac).col_join(TEgh).col_join(TEel)
            JacInSpat_full[:,imus] = iJacInSpat
            # print('poly')
            
        mus_lengths_full[imus] = L     
        mus_forces[imus] = muscle_force(actSym[imus],L,vce,fmax[imus],lceopt[imus],lslack[imus],vmax[imus])

        if muscle['crossesGH'] == 1:
            xparam_count = muscle['xparam_count'].item()
            xparams = muscle['xparams'].item()
            xcoefs = muscle['xcoefs'].item()
            yparam_count = muscle['yparam_count'].item()
            yparams = muscle['yparams'].item()
            ycoefs = muscle['ycoefs'].item()
            zparam_count = muscle['zparam_count'].item()
            zparams = muscle['zparams'].item()
            zcoefs = muscle['zcoefs'].item()
            xvec = 0
            yvec = 0
            zvec = 0
            

            for i in range(xparam_count.item()):
             # Add this term's contribution to the muscle length
                term = xcoefs[i]
                for j in range(nmusdof.item()):
                    mdof = musdof[j];
                    for k in range(xparams[i,j].item()):
                        term = term * qpol[int(mdof-1)];

                xvec = xvec + term
    
            for i in range(yparam_count.item()):
             # Add this term's contribution to the muscle length
                term = ycoefs[i]
                for j in range(nmusdof.item()):
                    mdof = musdof[j];
                    for k in range(yparams[i,j].item()):
                        term = term * qpol[int(mdof-1)];

                yvec = yvec + term
    
            for i in range(zparam_count.item()):
             # Add this term's contribution to the muscle length
                term = zcoefs[i]
                for j in range(nmusdof.item()):
                    mdof = musdof[j];
                    for k in range(zparams[i,j].item()):
                        term = term * qpol[int(mdof-1)];

                zvec = zvec + term;
            
            # print(xvec,yvec,zvec)
            # vecs = sp.Matrix([xvec,yvec,zvec])
            # vecs_np = sp.lambdify(q_new,vecs)
            # print(vecs_np(*[0,0,0,0,0.919602359906538,-0.167184664249180,-0.352129042135548,-0.0488459351319762,0.850386214641430,0.251785446200559,0.456800642893191,0.0691415047693658,0.959918850208316,0.0301071839771582,-0.254775233279737,0.112866908325537,120*np.pi/180,0]))

            xforce += xvec * mus_forces[imus]
            yforce += yvec * mus_forces[imus]
            zforce += zvec * mus_forces[imus]
            

    FQ = JacInSpat_full * mus_forces
    
    TE = sp.Matrix(FQ[:-1])
    TE = TE.subs(q_new[17],120*np.pi/180)
    JacInSpat = JacInSpat_full.subs(q_new[17],120*np.pi/180)
    mus_lengths = mus_lengths_full.subs(q_new[17],120*np.pi/180)
    mus_forces = mus_forces.subs(q_new[17],120*np.pi/180)
    GH_mus_forces = sp.Matrix([xforce, yforce, zforce])
    GH_mus_forces = GH_mus_forces.subs(q_new[17],120*np.pi/180)
    
    
    symbols_list = TE.free_symbols
    t = sp.Symbol('t')
    actSym_list = [i for i in symbols_list if str(i).startswith('actSym')]
    
    act = []
    for i in actSym_list:
        act.append(me.dynamicsymbols(str(i).replace('Sym','')))
    act_subs = dict(zip(actSym_list,act))
    TE_subs = me.msubs(TE, act_subs)
    GH_mus_forces_subbed = me.msubs(GH_mus_forces, act_subs)

    # conoid_lopt = model_params_struct['params']['model'].item()['conoid_length'].item()[0][0]
    # conoid_lopt = model_struct['model']['conoid_length'].item()[0][0]
    # conoid_k = model_params_struct['params']['model'].item()['conoid_stiffness'].item()[0][0]
    conoid_lopt = 0.018
    # conoid_k = model_params_struct['params']['model'].item()['conoid_stiffness'].item()[0][0]
    # conoid_k = model_struct['model']['conoid_stiffness'].item()[0][0]
    conoid_k = 5e3
    # conoid_eps = model_params_struct['params']['model'].item()['conoid_eps'].item()[0][0]
    # conoid_eps = model_struct['model']['conoid_eps'].item()[0][0]
    conoid_eps = 1e-3
    # conoid_insertion = model_params_struct['params']['model'].item()['conoid_insertion'].item()[0]
    conoid_origin = model_struct['model']['conoid_origin'].item()[0]
    print(conoid_origin)
    conoid_insertion = model_struct['model']['conoid_insertion'].item()[0]
    conoid_length = analytic_length_quat('clavicle_r','scapula_r', conoid_origin, conoid_insertion, q_new[4:], model_struct)
    F_conoid = conoid_force(conoid_length, conoid_lopt, conoid_k, conoid_eps)
    Jac_conoid = -sp.Matrix([conoid_length]).jacobian(q_new).T                
    TEsc_conoid = 1/2 * G(q_new[4:8])*sp.Matrix(Jac_conoid[4:8])
    TEac_conoid = 1/2 * G(q_new[8:12])*sp.Matrix(Jac_conoid[8:12])
    TEgh_conoid = 1/2 * G(q_new[12:16])*sp.Matrix(Jac_conoid[12:16])
    TEel_conoid = 1/2 * sp.Matrix(Jac_conoid[16:18])
    JacInSpat_conoid = sp.Matrix(TEsc_conoid).col_join(TEac_conoid).col_join(TEgh_conoid).col_join(TEel_conoid)
    TE_conoid = F_conoid * JacInSpat_conoid
    
    if gen_matlab_functions == 1:
        qsubs = sp.symbols('qsubs0:13')
        q_subs_dict = dict(zip(q_new[4:-1],qsubs))
        segments = []
        # TE_conoid_subbed = me.msubs(TE_conoid,q_subs_dict)
        TE_subbed = me.msubs(TE,q_subs_dict)
        JacInSpat_subbed = me.msubs(JacInSpat_full,q_subs_dict)
        mus_lengths_subbed = me.msubs(mus_lengths_full,q_subs_dict)
        mus_forces_subbed = me.msubs(mus_forces,q_subs_dict)
        TE_conoid_subbed = me.msubs(TE_conoid,q_subs_dict)
        
        MatlabFunction(function = TE_conoid_subbed[:-1],
                       fun_name = 'TE_conoid_quat', assignto = 'TE_conoid',
                       coordinates = qsubs,
                       speeds = [],
                       inputs = [],
                       body_constants = {},
                       segments = [],
                       other_constants = {},
                       muscle_constants = {},
                       parameters = [],
                       folder = '../equations_of_motion/quaternion')
        # MatlabFunction(function = TE_subbed,
        #                fun_name = 'TE_quat', assignto = 'TE',
        #                coordinates = qsubs,
        #                speeds = [],
        #                inputs = actSym,
        #                body_constants = {},
        #                segments = [],
        #                other_constants={},
        #                muscle_constants = muscle_constants,
        #                parameters = [],
        #                folder = folder)
        # MatlabFunction(function = mus_forces_subbed,
        #                fun_name = 'mus_forces_quat', assignto = 'mus_forces',
        #                coordinates = qsubs,
        #                speeds = [],
        #                inputs = actSym,
        #                body_constants = {},
        #                segments = [],
        #                other_constants={},
        #                muscle_constants = muscle_constants,
        #                parameters = [])
        # MatlabFunction(function = JacInSpat_subbed,
        #                fun_name = 'JacInSpatnoi_quat',assignto = 'JacInSpat',
        #                coordinates = qsubs,
        #                speeds = [],
        #                inputs = [],
        #                body_constants = {},
        #                segments = [],
        #                other_constants={},
        #                muscle_constants = {},
        #                parameters = [],
        #                folder = '../Motions/' + motion_folder + '/Poly_functions')
        # MatlabFunction(function = mus_lengths_subbed,
        #                fun_name = 'mus_lengths_quat',assignto = 'mus_lengths',
        #                coordinates = qsubs,
        #                speeds = [],
        #                inputs = [],
        #                body_constants = {},
        #                segments = [],
        #                other_constants={},
        #                muscle_constants = {},
        #                parameters = [],
        #                folder = '../Motions/' + motion_folder + '/Poly_functions')
    
    return TE_subs, act, TE_conoid[:-1], fmax_init, fmax_range, lceopt_init, lceopt_range, mus_groups, GH_mus_forces_subbed

    
def create_parameters_dict(model_params_struct, initCond_name):
    data_struct = model_params_struct['params']['model'][0,0]
    data_names = data_struct.dtype.names
    # data = data_struct[data_names[1]]
    element_names = []
    element_values = []
    symbolic_list = []
    # element = data_struct['offset_scapula']
    
    for name in data_names:
        element = data_struct[name].item()
        element_length = element.shape[1]
        
        if element_length == 1:
            element_names.append(name)
            element_values.append(element.item())
            
        else:
            for k in range(element_length):
                element_names.append((name+'_'+str(k+1)))
                element_values.append(element[0,k])
                
    muscles_params_struct = model_params_struct['params'][initCond_name][0,0]
    muscles_params_names = ['fmax','lceopt','lslack']
    
    muscle_names = []
    muscle_values = []
    missing = (20,31,32,33)
    
    for name in muscles_params_names:
        muscles_param = muscles_params_struct[name].item()
        muscles_params_len = muscles_param.shape[0]
        
        for k in range(muscles_params_len):
            if k == 20:
                pass
            elif k == 31:
                pass
            elif k == 32:
                pass
            elif k == 33:
                pass
            else:
                element_names.append((name+'_'+str(k+1)))
                element_values.append(muscles_param[k,0])
                
    symbolic_list = sp.symbols(element_names)
    dict_vals = dict(zip(symbolic_list, element_values))
        
    return dict_vals, symbolic_list, element_values
    
def save_eoms_explicit(MM,FO,q,w,name):
    # Sym2Dym,Dym2Sym = create_sym_dym_dict(q,w)
    qs = sp.symbols('q1:14')
    us = sp.symbols('u1:11')
    
    states = [qs,us]
    subs_q = {q[i]: qs[i] for i in range(len(q))}
    subs_u = {w[i]: us[i] for i in range(len(w))}
    
    MM_sym = me.msubs(MM,subs_q,subs_u)
    FO_sym = me.msubs(FO,subs_q,subs_u)
        
    with open(name,'wb') as f:
        pickle.dump([MM_sym,FO_sym,qs,us],f)
        
def save_eoms_implicit(eoms_implicit,q,w,name):
    subs_q, subs_u, subs_dq, subs_du, _, _, _, _ = create_sym_dym_dict(q,w)
    eoms_subs = me.msubs(eoms_implicit,subs_q,subs_u,subs_dq,subs_du)
    
    with open(name,'wb') as f:
        pickle.dump([eoms_subs],f)
        
def load_eoms_implicit(q,w,name):
    with open(name,'rb') as f:
        eoms_impl = pickle.load(f)
        
    eoms_impl_sym = sp.Matrix(eoms_impl)
        
    _, _, _, _, subs_qback, subs_uback, subs_dqback, subs_duback = create_sym_dym_dict(q,w)
    eoms_back = me.msubs(eoms_impl_sym,subs_qback,subs_uback,subs_dqback,subs_duback)
        
    return eoms_impl_sym
    
        
def load_eoms(name):

    with open(name,'rb') as f:
        MM,FO,q,w = pickle.load(f)
        
    return MM,FO,list(q),list(w)
    
    
def create_sym_dym_dict(q,w):
    qs = sp.symbols('q0:13')
    us = sp.symbols('u1:11')
    dqs = sp.symbols('dq0:13')
    dus = sp.symbols('du1:11')
    
    subs_q = {q[i]: qs[i] for i in range(len(q))}
    subs_u = {w[i]: us[i] for i in range(len(w))}
    subs_dq = {q[i].diff(): dqs[i] for i in range(len(q))}
    subs_du = {w[i].diff(): dus[i] for i in range(len(w))}
    
    subs_qback = swtich_dict(subs_q)
    subs_uback = swtich_dict(subs_u)
    subs_dqback = swtich_dict(subs_dq)
    subs_duback = swtich_dict(subs_du)
    
    return subs_q, subs_u, subs_dq, subs_du, subs_qback, subs_uback, subs_dqback, subs_duback
    
    for i,seg in enumerate(segment):
        if joints[i] == 'quat':
            for j in ('0','1','2','3'):
                qsym.append(sp.symbols('q' + j+ '_' + seg))
        elif joints[i] == 'rotaxis':
            qsym.append(sp.symbols('q_' + seg))
        else:
            pass
        
        if joints[i] == 'quat':
            for j in ('1','2','3'):
                wsym.append(sp.symbols('w' + j + '_' + seg))
        elif joints[i] == 'rotaxis':
            wsym.append(sp.symbols('w_'+ seg))
        else:
            pass
        
    symbolic_list = qsym+wsym
    dynamic_list = qdyn+wdyn
    sym_dym_dict = dict(zip(symbolic_list,dynamic_list))
    dym_sym_dict = {y: x for x, y in sym_dym_dict.items()}
        
        
    return sym_dym_dict, dym_sym_dict
    
    
#############################################
                #QUAT#
#############################################
    
    
def create_eoms_quat(OS_struct, derive = 'symbolic',gen_matlab_functions = None):
    
    # symbols
    t = sp.symbols('t')

    states = ['q','w','u0']
    segment = ['clavicula','scapula','humerus','ulna','radius','hand']
    segment_index = [3,6,9,10,11,12]
    joints = ['quat','quat','quat','rotaxis','weld','weld']
    joint_index = [4,7,10,11,12,13]
    inertia = []
    mass = []
    com = []
    offset = []
    rot_offset = []
    q = []
    w = []
    u0 = []
    frame = []
    point_offset = []
    masscenter = []
    inertia_elem = []

    if derive == 'symbolic':
        g,c = sp.symbols('g,c')  # 
    elif derive == 'numeric':
        g = 9.81
        c = 0.5


    for i,seg in enumerate(segment):
        for idat,j in enumerate(('1','2','3')):
            if derive == 'symbolic':
                inertia.append(sp.symbols('I_' + seg + '_' + j))
                com.append(sp.symbols('com_' + seg + '_' + j))
                offset.append(sp.symbols('offset_' + seg + '_' + j))
            elif derive == 'numeric':
                com.append(OS_struct['model']['segments'][0,0][0,segment_index[i]]['mass_center'][0,0][0,idat].item())
                inertia_diag = np.diag(OS_struct['model']['segments'][0,0][0,segment_index[i]]['inertia'][0,0])
                inertia.append(inertia_diag[idat].item())
                try:
                    offset.append(OS_struct['model']['joints'][0,0][0,joint_index[i]]['location'][0,0][0,idat].item())
                except:
                    offset.append(0)

        if joints[i] == 'quat':
            for j in ('1','2','3'):
                w.append(me.dynamicsymbols('w' + j + '_' + seg))
            for j in ('0','1','2','3'):
                q.append(me.dynamicsymbols('q' + j+ '_' + seg))
            u0.append(me.dynamicsymbols('u0'+ '_' + seg))
        elif joints[i] == 'rotaxis':
            w.append(me.dynamicsymbols('w_'+ seg))
            q.append(me.dynamicsymbols('q_' + seg))
        else:
            pass

        if derive == 'symbolic':
            mass.append(sp.symbols('mass_'+seg))
        elif derive == 'numeric':
            mass.append(OS_struct['model']['segments'][0,0][0,segment_index[i]]['mass'][0,0].item())
        print(mass)
        frame.append(me.ReferenceFrame('frame_' + str(seg)))
        point_offset.append(me.Point('point_offset_' + str(seg)))
        masscenter.append(me.Point('masscenter_' + str(seg)))

    # inertial frame and point
    frame_ground = me.ReferenceFrame('frame_ground')
    point_ground = me.Point('point_ground')
    point_ground.set_vel(frame_ground,0)

    offset_ground = me.Point('offset_ground')
    if derive == 'symbolic':
        offset_thorax = sp.symbols('offset_thorax_1:4')
    elif derive == 'numeric':
        offset_thorax = OS_struct['model']['joints'][0,0][0,1]['location'][0,0][0]

    offset_ground.set_pos(point_ground, offset_thorax[0]*frame_ground.x 
                          + offset_thorax[1]*frame_ground.y + offset_thorax[2]*frame_ground.z)
    offset_ground.set_vel(frame_ground,0)

    #rotate first body
    frame[0].orient(frame_ground, 'Quaternion', q[0:4])

    for i in range(1,3):
        frame[i].orient(frame[i-1], 'Quaternion', q[0+i*4:4+i*4])
    kinematical =[]
    for i in range(0,3):
        kinematical.append(q[0+i*4].diff(t) - 0.5 * (-w[0+i*3]*q[1+i*4] - w[1+i*3]*q[2+i*4] - w[2+i*3]*q[3+i*4]))
        kinematical.append(q[1+i*4].diff(t) - 0.5 * (w[0+i*3]*q[0+i*4] + w[2+i*3]*q[2+i*4] - w[1+i*3]*q[3+i*4]))
        kinematical.append(q[2+i*4].diff(t) - 0.5 * (w[1+i*3]*q[0+i*4] - w[2+i*3]*q[1+i*4] + w[0+i*3]*q[3+i*4]))
        kinematical.append(q[3+i*4].diff(t) - 0.5 * (w[2+i*3]*q[0+i*4] + w[1+i*3]*q[1+i*4] - w[0+i*3]*q[2+i*4]))

    frame[0].set_ang_vel(frame_ground,w[0]*frame[0].x + w[1]*frame[0].y + w[2]*frame[0].z)

    for i in range(1,3):
        frame[i].set_ang_vel(frame[i-1],w[0+i*3]*frame[i].x + w[1+i*3]*frame[i].y + w[2+i*3]*frame[i].z)

    # set masscenter of first body
    masscenter[0].set_pos(offset_ground,com[0]*frame[0].x + com[1]*frame[0].y + com[2]*frame[0].z)
    masscenter[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set offset of first joint in first body
    point_offset[0].set_pos(offset_ground,offset[0]*frame[0].x + offset[1]*frame[0].y + offset[2]*frame[0].z)
    point_offset[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set gravity force and damping of first body
    FG = [(masscenter[0], -mass[0] * g * frame_ground.y)]
    DAMP = [(frame[0], -c*(w[0]*frame[0].x+w[1]*frame[0].y+w[2]*frame[0].z))]
    # iterate over segments 2:end (first body is already done)
    for i in range(1,3):

        # set masscenter points 
        masscenter[i].set_pos(point_offset[i-1],com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
        masscenter[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set gravity force in masscenter (-y direction in frame_ground)
        FG.append((masscenter[i], -mass[i] * g * frame_ground.y))

        # set offsent points (where the next joint is)
        point_offset[i].set_pos(point_offset[i-1],offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
        point_offset[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set damping in joints (c * angular_velocity)
        damping = -c*(w[0+i*3]*frame[i].x+w[1+i*3]*frame[i].y+w[2+i*3]*frame[i].z)

        # apply damping in frame, opposite moment is applied in previous frame (action and reaction)
        DAMP.append((frame[i], damping))
        DAMP.append((frame[i-1], -damping))
        # set kinematic differential equations (q_dot = f(q,u))

    # symbols for ulna
    ulna_rot_frame = me.ReferenceFrame('ulna_rot_frame')

    if derive == 'symbolic':
        offset_humerus_rot = sp.symbols('offset_humerus_rot_1:4')
        EL_rot_axis = sp.symbols('EL_rot_axis_1:4')
        PSY_rot_axis = sp.symbols('PSY_rot_axis_1:4')
    if derive == 'numeric':
        offset_humerus_rot = OS_struct['model']['joints'][0,0][0,10]['orientation'][0,0][0]
        EL_rot_axis = OS_struct['model']['joints'][0,0][0,10]['r1_axis'][0,0][0]
        PSY_rot_axis = OS_struct['model']['joints'][0,0][0,11]['r1_axis'][0,0][0]


    # offset frame rotated in humerus frame (frame[2])
    ulna_rot_frame.orient_axis(frame[2],frame[2].z,offset_humerus_rot[2])

    # ulna and elbow joint
    frame[3].orient_axis(ulna_rot_frame,ulna_rot_frame.x*EL_rot_axis[0]
                         +ulna_rot_frame.y*EL_rot_axis[1]+ulna_rot_frame.z*EL_rot_axis[2],
                         q[12])

    frame[3].set_ang_vel(ulna_rot_frame,
                          w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2]))

    masscenter[3].set_pos(point_offset[2],com[0+3*3]*frame[3].x + com[1+3*3]*frame[3].y + com[2+3*3]*frame[3].z)
    masscenter[3].v2pt_theory(point_offset[2],frame_ground,frame[3])
    FG.append((masscenter[3], -mass[3] * g * frame_ground.y))

    DAMP.append(((frame[3]),-c*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                            +ulna_rot_frame.z*EL_rot_axis[2])))
    DAMP.append((ulna_rot_frame,c*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2])))
    kinematical.append(w[9]-q[12].diff(t))

    point_offset[3].set_pos(point_offset[2],offset[0+3*3]*frame[3].x + offset[1+3*3]*frame[3].y + offset[2+3*3]*frame[3].z)
    point_offset[3].v2pt_theory(point_offset[2],frame_ground,frame[3]);

    # radius and PSY joint (weld joint for now)
    frame[4].orient_axis(frame[3],frame[3].z,0)
    frame[4].set_ang_vel(frame[3],0)
    masscenter[4].set_pos(point_offset[3],com[0+4*3]*frame[4].x + com[1+4*3]*frame[4].y + com[2+4*3]*frame[4].z)
    masscenter[4].v2pt_theory(point_offset[3],frame_ground,frame[4])
    FG.append((masscenter[4], -mass[4] * g * frame_ground.y))
    frame[4].ang_vel_in(frame[3])

    point_offset[4].set_pos(point_offset[3],offset[0+4*3]*frame[3].x + offset[1+4*3]*frame[3].y + offset[2+4*3]*frame[3].z)
    point_offset[4].v2pt_theory(point_offset[3],frame_ground,frame[3])

    # hand
    frame[5].orient_axis(frame[4],frame[4].z,0)
    frame[5].set_ang_vel(frame[4],0)
    masscenter[5].set_pos(point_offset[4],com[0+5*3]*frame[5].x + com[1+5*3]*frame[5].y + com[2+5*3]*frame[5].z)
    masscenter[5].v2pt_theory(point_offset[4],frame_ground,frame[5])
    FG.append((masscenter[5], -mass[5] * g * frame_ground.y))

    BODY = []

    for i in range(len(segment)):
        # set inertias of each body and create RigidBodies
        I = me.inertia(frame[i], inertia[0+i*3], inertia[1+i*3], inertia[2+i*3])
        BODY.append(me.RigidBody('body' + str(i), masscenter[i], frame[i], mass[i], (I, masscenter[i])))

    # Contact between scapula and thorax
    if derive == 'symbolic':
        contTS = sp.symbols('contTS_1:4')
        contAI = sp.symbols('contAI_1:4')
        elips_trans = sp.symbols('elips_trans_1:4')
        elips_dim = sp.symbols('elips_dim_1:4')
        k_contact_in, eps_in = sp.symbols('k_contact_in eps_in')
        k_contact_out, eps_out = sp.symbols('k_contact_out eps_out')
        # first_elips_scale = sp.Symbol('first_elips_scale')
        second_elips_scale = sp.Symbol('second_elips_scale')
    elif derive == 'numeric':
        # k_contact_in = data_struct['k_contact_in'][0,0].item()
        k_contact_out = 0
        # eps_in = data_struct['eps_in'][0,0].item()
        
        k_contact_in = OS_struct['model']['scap_thorax_k'][0,0].item()
        # eps_in = OS_struct['model']['scap_thorax_eps'][0,0].item()
        eps_in = 0.01
        eps_out = 1e-3
        # first_elips_scale = model_params_struct['params'][initCond_name][0,0]['first_elips_scale'][0,0].item()
        # first_elips_scale = sp.symbols('first_elips_scale1:4')
        # elips_dim = OS_struct['model']['thorax_dim_optimized'].item()[0]
        elips_trans = OS_struct['model']['thoracic_wall'][0,0].item()[0][0]
        elips_dim = OS_struct['model']['thoracic_wall'][0,0].item()[1][0]
        contTS = OS_struct['model']['TScontact'][0,0].item()[0][0]
        contAI = OS_struct['model']['AIcontact'][0,0].item()[0][0]
        first_elips_scale = [1.0,1.0,1.0]
        # elips_trans = sp.symbols('elips_trans1:4')
        second_elips_scale = [1.0,1.0,1.0]

    # contact points 
    contact_point1 = me.Point('CP1')
    contact_point1.set_pos(point_offset[0],contTS[0]*frame[1].x+contTS[1]*frame[1].y  +contTS[2]*frame[1].z)
    # contact_point1.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point1.v2pt_theory(point_offset[0],frame_ground,frame[1])

    contact_point2 = me.Point('CP2')
    contact_point2.set_pos(point_offset[0],contAI[0]*frame[1].x+contAI[1]*frame[1].y  +contAI[2]*frame[1].z)
    # contact_point2.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point2.v2pt_theory(point_offset[0],frame_ground,frame[1])

    ## contact forces

    # Distances between contact points and thorax frame
    x_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.x)
    y_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.y)
    z_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.z)
    x_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.x)
    y_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.y)
    z_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.z)

    # Contact forces
    elips_dim_scaled = []
    elips_dim_scaled.append(elips_dim[0]*first_elips_scale[0])
    elips_dim_scaled.append(elips_dim[1]*first_elips_scale[1])
    elips_dim_scaled.append(elips_dim[2]*first_elips_scale[2])
    # print(elips_dim_scaled)
    f1_in = ((x_pos1-elips_trans[0])/(elips_dim_scaled[0]))**2+((y_pos1-elips_trans[1])/(elips_dim_scaled[1]))**2+((z_pos1-elips_trans[2])/(elips_dim_scaled[2]))**2-1
    # f1_out = ((x_pos1-elips_trans[0])/(second_elips_scale*elips_dim[0]))**2+((y_pos1-elips_trans[1])/(second_elips_scale*elips_dim[1]))**2+((z_pos1-elips_trans[2])/(second_elips_scale*elips_dim[2]))**2-1
    F1_in = 1/2*(f1_in-sp.sqrt(f1_in**2+eps_in**2))
    # F1_out = 1/2*(f1_out+sp.sqrt(f1_out**2+eps_out**2))
    Fx1 = -(k_contact_in*F1_in)*(x_pos1-elips_trans[0])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[0]**2)
    Fy1 = -(k_contact_in*F1_in)*(y_pos1-elips_trans[1])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[1]**2)
    Fz1 = -(k_contact_in*F1_in)*(z_pos1-elips_trans[2])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[2]**2)

    f2_in = ((x_pos2-elips_trans[0])/(elips_dim_scaled[0]))**2+((y_pos2-elips_trans[1])/(elips_dim_scaled[1]))**2+((z_pos2-elips_trans[2])/(elips_dim_scaled[2]))**2-1
    # f2_out = ((x_pos2-elips_trans[0])/(second_elips_scale*elips_dim[0]))**2+((y_pos2-elips_trans[1])/(second_elips_scale*elips_dim[1]))**2+((z_pos2-elips_trans[2])/(second_elips_scale*elips_dim[2]))**2-1
    F2_in = 1/2*(f2_in-sp.sqrt(f2_in**2+eps_in**2))
    # F2_out = 1/2*(f2_out+sp.sqrt(f2_out**2+eps_out**2))
    Fx2 = -(k_contact_in*F2_in)*(x_pos2-elips_trans[0])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[0]**2)
    Fy2 = -(k_contact_in*F2_in)*(y_pos2-elips_trans[1])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[1]**2)
    Fz2 = -(k_contact_in*F2_in)*(z_pos2-elips_trans[2])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[2]**2)

    # applying contact forces to contact points in thorax frame
    cont_force1 = [(contact_point1,frame_ground.x*Fx1+frame_ground.y*Fy1+frame_ground.z*Fz1)]
    cont_force2 = [(contact_point2,frame_ground.x*Fx2+frame_ground.y*Fy2+frame_ground.z*Fz2)]
    CONT = cont_force1+cont_force2

    # print(kinematical)
    KM = me.KanesMethod(frame_ground, q_ind=q, u_ind=w, kd_eqs = kinematical)

    (fr, frstar) = KM.kanes_equations(BODY, (FG+DAMP+CONT))
    MM = KM.mass_matrix_full
    FO = KM.forcing_full
    xdot = (KM.q.col_join(KM.u)).diff()
    
    return MM,FO,q,w,fr,frstar,kinematical,xdot,first_elips_scale,elips_trans



#############################################
                #EUL#
#############################################

def create_eoms_eul(OS_struct, derive = 'symbolic',gen_matlab_functions = None,GH_seq = 'YZY'):
    # symbols
    t = sp.symbols('t')
    g,c = sp.symbols('g,c')  # 
    states = ['q','u']
    segment = ['clavicula','scapula','humerus','ulna','radius','hand']
    segment_index = [3,6,9,10,11,12]
    joints = ['YZX','YZX',GH_seq,'rotaxis','weld','weld']
    joint_index = [4,7,10,11,12,13]
    inertia = []
    mass = []
    com = []
    offset = []
    rot_offset = []
    q = []
    u = []
    rot = [] # for angular velocities
    frame = []
    point_offset = []
    masscenter = []
    inertia_elem = []
        
    if derive == 'symbolic':
        g,c = sp.symbols('g,c')  # 
    elif derive == 'numeric':
        g = 9.81
        c = 0.5

    for i,seg in enumerate(segment):
        for idat,j in enumerate(('1','2','3')):
            if derive == 'symbolic':
                inertia.append(sp.symbols('I_' + seg + '_' + j))
                com.append(sp.symbols('com_' + seg + '_' + j))
                offset.append(sp.symbols('offset_' + seg + '_' + j))
            elif derive == 'numeric':
                # inertia.append(data_struct['I_' + seg][0,0][0,idat])
                # com.append(data_struct['com_' + seg][0,0][0,idat])
                # offset.append(data_struct['offset_' + seg][0,0][0,idat])
                com.append(OS_struct['model']['segments'][0,0][0,segment_index[i]]['mass_center'][0,0][0,idat].item())
                inertia_diag = np.diag(OS_struct['model']['segments'][0,0][0,segment_index[i]]['inertia'][0,0])
                inertia.append(inertia_diag[idat].item())

                try:
                    offset.append(OS_struct['model']['joints'][0,0][0,joint_index[i]]['location'][0,0][0,idat].item())
                except:
                    offset.append(0)

        if joints[i] == 'YZX' or joints[i] == 'YZY':
            for j in ('1','2','3'):
                q.append(me.dynamicsymbols('q' + j + '_' + seg))
                u.append(me.dynamicsymbols('u' + j + '_' + seg))
        elif joints[i] == 'rotaxis':
            q.append(me.dynamicsymbols('q_' + seg))
            u.append(me.dynamicsymbols('u_' + seg))
        else:
            pass

        if derive == 'symbolic':
            mass.append(sp.symbols('mass_'+seg))
        elif derive == 'numeric':
            # mass.append(data_struct['mass_' + seg][0,0].item())
            mass.append(OS_struct['model']['segments'][0,0][0,segment_index[i]]['mass'][0,0].item()/2)


        frame.append(me.ReferenceFrame('frame_' + str(seg)))
        point_offset.append(me.Point('point_offset_' + str(seg)))
        masscenter.append(me.Point('masscenter_' + str(seg)))

    # inertial frame and point
    frame_ground = me.ReferenceFrame('frame_ground')
    point_ground = me.Point('point_ground')
    point_ground.set_vel(frame_ground,0)

    offset_ground = me.Point('offset_ground')
    if derive == 'symbolic':
        offset_thorax = sp.symbols('offset_thorax_1:4')
    elif derive == 'numeric':
        offset_thorax = OS_struct['model']['joints'][0,0][0,1]['location'][0,0][0]

    offset_ground.set_pos(point_ground, offset_thorax[0]*frame_ground.x 
                          + offset_thorax[1]*frame_ground.y + offset_thorax[2]*frame_ground.z)
    offset_ground.set_vel(frame_ground,0)

    #rotate first body
    frame[0].orient_body_fixed(frame_ground, q[0:3], rotation_order = joints[0])
    rot.append(frame[0].ang_vel_in(frame_ground))

    # set masscenter of first body
    masscenter[0].set_pos(offset_ground,com[0]*frame[0].x + com[1]*frame[0].y + com[2]*frame[0].z)
    masscenter[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set offset of first joint in first body
    point_offset[0].set_pos(offset_ground,offset[0]*frame[0].x + offset[1]*frame[0].y + offset[2]*frame[0].z)
    point_offset[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set gravity force and damping of first body
    FG = [(masscenter[0], -mass[0] * g * frame_ground.y)]
    DAMP = [(frame[0], -c*(rot[0].dot(frame[0].x)*frame[0].x+rot[0].dot(frame[0].y)*frame[0].y+rot[0].dot(frame[0].z)*frame[0].z))]
    kindeq = []

    # iterate over segments 2:end (first body is already done)
    for i in range(1,3):

        # orient frame w.r.t. child's frame
        frame[i].orient_body_fixed(frame[i-1], q[0+i*3:3+i*3], rotation_order = joints[i])
        rot.append(frame[i].ang_vel_in(frame[i-1]))

        # set masscenter points 
        masscenter[i].set_pos(point_offset[i-1],com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
        masscenter[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set gravity force in masscenter (-y direction in frame_ground)
        FG.append((masscenter[i], -mass[i] * g * frame_ground.y))

        # set offsent points (where the next joint is)
        point_offset[i].set_pos(point_offset[i-1],offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
        point_offset[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set damping in joints (c * angular_velocity)
        damping = -c*(rot[i].dot(frame[i].x)*frame[i].x+rot[i].dot(frame[i].y)*frame[i].y+rot[i].dot(frame[i].z)*frame[i].z)

        # apply damping in frame, opposite moment is applied in previous frame (action and reaction)
        DAMP.append((frame[i], damping))
        DAMP.append((frame[i-1], -damping))

    for i in range(9):
        kindeq.append(q[i].diff()-u[i])

    # symbols (or values) for ulna and radius
    ulna_rot_frame = me.ReferenceFrame('ulna_rot_frame')
    
    if derive == 'symbolic':
        offset_humerus_rot = sp.symbols('offset_humerus_rot_1:4')
        EL_rot_axis = sp.symbols('EL_rot_axis_1:4')
        PSY_rot_axis = sp.symbols('PSY_rot_axis_1:4')
    if derive == 'numeric':
        offset_humerus_rot = OS_struct['model']['joints'][0,0][0,10]['orientation'][0,0][0]
        EL_rot_axis = OS_struct['model']['joints'][0,0][0,10]['r1_axis'][0,0][0]
        PSY_rot_axis = OS_struct['model']['joints'][0,0][0,11]['r1_axis'][0,0][0]

        # for i in range(3):
        #     offset_humerus_rot.append(data_struct['offset_humerus_rot'][0,0][0,i].item())
        #     EL_rot_axis.append(data_struct['EL_rot_axis'][0,0][0,i].item())
        #     PSY_rot_axis.append(data_struct['PSY_rot_axis'][0,0][0,i].item())

    # offset frame rotated in humerus frame (frame[2])
    ulna_rot_frame.orient_axis(frame[2],frame[2].z,offset_humerus_rot[2])

    # ulna and elbow joint
    frame[3].orient_axis(ulna_rot_frame,ulna_rot_frame.x*EL_rot_axis[0]
                         +ulna_rot_frame.y*EL_rot_axis[1]+ulna_rot_frame.z*EL_rot_axis[2],
                         q[9])

    frame[3].set_ang_vel(ulna_rot_frame,
                          u[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2]))

    masscenter[3].set_pos(point_offset[2],com[0+3*3]*frame[3].x + com[1+3*3]*frame[3].y + com[2+3*3]*frame[3].z)
    masscenter[3].v2pt_theory(point_offset[2],frame_ground,frame[3])
    FG.append((masscenter[3], -mass[3] * g * frame_ground.y))

    DAMP.append(((frame[3]),-c*u[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                            +ulna_rot_frame.z*EL_rot_axis[2])))
    DAMP.append((ulna_rot_frame,c*u[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2])))
    kindeq.append(u[9]-q[9].diff(t))

    point_offset[3].set_pos(point_offset[2],offset[0+3*3]*frame[3].x + offset[1+3*3]*frame[3].y + offset[2+3*3]*frame[3].z)
    point_offset[3].v2pt_theory(point_offset[2],frame_ground,frame[3]);

    # radius and PSY joint (weld joint for now)
    frame[4].orient_axis(frame[3],frame[3].x*PSY_rot_axis[0]
                         +frame[3].y*PSY_rot_axis[1]+frame[3].z*PSY_rot_axis[2],
                         0) #q[10]
    masscenter[4].set_pos(point_offset[3],com[0+4*3]*frame[4].x + com[1+4*3]*frame[4].y + com[2+4*3]*frame[4].z)
    masscenter[4].v2pt_theory(point_offset[3],frame_ground,frame[4])
    FG.append((masscenter[4], -mass[4] * g * frame_ground.y))
    # DAMP.append(((frame[4]),-c*u[10]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
    #                         +ulna_rot_frame.z*EL_rot_axis[2])))
    # DAMP.append((ulna_rot_frame,c*u[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
    #                             +ulna_rot_frame.z*EL_rot_axis[2])))
    # kindeq.append(u[10]-q[10].diff(t))

    point_offset[4].set_pos(point_offset[3],offset[0+4*3]*frame[4].x + offset[1+4*3]*frame[4].y + offset[2+4*3]*frame[4].z)
    point_offset[4].v2pt_theory(point_offset[3],frame_ground,frame[4])

    # hand
    frame[5].orient_axis(frame[4],frame[4].z,0)
    frame[5].set_ang_vel(frame[4],0)
    masscenter[5].set_pos(point_offset[4],com[0+5*3]*frame[5].x + com[1+5*3]*frame[5].y + com[2+5*3]*frame[5].z)
    masscenter[5].v2pt_theory(point_offset[4],frame_ground,frame[5])
    FG.append((masscenter[5], -mass[5] * g * frame_ground.y))

    BODY = []

    for i in range(len(segment)):
        # set inertias of each body and create RigidBodies
        I = me.inertia(frame[i], inertia[0+i*3], inertia[1+i*3], inertia[2+i*3])
        BODY.append(me.RigidBody('body' + str(i), masscenter[i], frame[i], mass[i], (I, masscenter[i])))

    # Contact between scapula and thorax
    if derive == 'symbolic':
        contTS = sp.symbols('contTS_1:4')
        contAI = sp.symbols('contAI_1:4')
        elips_trans = sp.symbols('elips_trans_1:4')
        elips_dim = sp.symbols('elips_dim_1:4')
        k_contact_in, eps_in = sp.symbols('k_contact_in eps_in')
        k_contact_out, eps_out = sp.symbols('k_contact_out eps_out')
        first_elips_scale = sp.symbols('first_elips_scale')
        second_elips_scale = sp.Symbol('second_elips_scale')
    elif derive == 'numeric':
        # k_contact_in = data_struct['k_contact_in'][0,0].item()
        k_contact_out = 0
        # eps_in = data_struct['eps_in'][0,0].item()
        
        k_contact_in = OS_struct['model']['scap_thorax_k'][0,0].item()
        # eps_in = OS_struct['model']['scap_thorax_eps'][0,0].item()
        eps_in = 0.01
        eps_out = 1e-3
        # first_elips_scale = model_params_struct['params'][initCond_name][0,0]['first_elips_scale'][0,0].item()
        # first_elips_scale = sp.symbols('first_elips_scale1:4')
        # elips_dim = OS_struct['model']['thorax_dim_optimized'].item()[0]
        elips_trans = OS_struct['model']['thoracic_wall'][0,0].item()[0][0]
        elips_dim = OS_struct['model']['thoracic_wall'][0,0].item()[1][0]
        # elips_trans = OS_struct['model']['thorax_center'].item()[0]
        contTS = OS_struct['model']['TScontact'][0,0].item()[0][0]
        contAI = OS_struct['model']['AIcontact'][0,0].item()[0][0]
        first_elips_scale = [1.0,1.0,1.0]
        # elips_trans = sp.symbols('elips_trans1:4')
        second_elips_scale = [1.0,1.0,1.0]
        print('dimensions:',elips_dim)
        print('trans',elips_trans)

    # contact points 
    contact_point1 = me.Point('CP1')
    contact_point1.set_pos(point_offset[0],contTS[0]*frame[1].x+contTS[1]*frame[1].y  +contTS[2]*frame[1].z)
    # contact_point1.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point1.v2pt_theory(point_offset[0],frame_ground,frame[1])

    contact_point2 = me.Point('CP2')
    contact_point2.set_pos(point_offset[0],contAI[0]*frame[1].x+contAI[1]*frame[1].y  +contAI[2]*frame[1].z)
    # contact_point2.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point2.v2pt_theory(point_offset[0],frame_ground,frame[1])

    # Distances between contact points and thorax frame
    x_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.x)
    y_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.y)
    z_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.z)
    x_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.x)
    y_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.y)
    z_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.z)

    # Contact forces
    elips_dim_scaled = []
    elips_dim_scaled.append(elips_dim[0]*first_elips_scale[0])
    elips_dim_scaled.append(elips_dim[1]*first_elips_scale[1])
    elips_dim_scaled.append(elips_dim[2]*first_elips_scale[2])
    # print(elips_dim_scaled)
    f1_in = ((x_pos1-elips_trans[0])/(elips_dim_scaled[0]))**2+((y_pos1-elips_trans[1])/(elips_dim_scaled[1]))**2+((z_pos1-elips_trans[2])/(elips_dim_scaled[2]))**2-1
    # f1_out = ((x_pos1-elips_trans[0])/(second_elips_scale*elips_dim[0]))**2+((y_pos1-elips_trans[1])/(second_elips_scale*elips_dim[1]))**2+((z_pos1-elips_trans[2])/(second_elips_scale*elips_dim[2]))**2-1
    F1_in = 1/2*(f1_in-sp.sqrt(f1_in**2+eps_in**2))
    # F1_out = 1/2*(f1_out+sp.sqrt(f1_out**2+eps_out**2))
    Fx1 = -(k_contact_in*F1_in)*(x_pos1-elips_trans[0])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[0]**2)
    Fy1 = -(k_contact_in*F1_in)*(y_pos1-elips_trans[1])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[1]**2)
    Fz1 = -(k_contact_in*F1_in)*(z_pos1-elips_trans[2])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[2]**2)

    f2_in = ((x_pos2-elips_trans[0])/(elips_dim_scaled[0]))**2+((y_pos2-elips_trans[1])/(elips_dim_scaled[1]))**2+((z_pos2-elips_trans[2])/(elips_dim_scaled[2]))**2-1
    # f2_out = ((x_pos2-elips_trans[0])/(second_elips_scale*elips_dim[0]))**2+((y_pos2-elips_trans[1])/(second_elips_scale*elips_dim[1]))**2+((z_pos2-elips_trans[2])/(second_elips_scale*elips_dim[2]))**2-1
    F2_in = 1/2*(f2_in-sp.sqrt(f2_in**2+eps_in**2))
    # F2_out = 1/2*(f2_out+sp.sqrt(f2_out**2+eps_out**2))
    Fx2 = -(k_contact_in*F2_in)*(x_pos2-elips_trans[0])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[0]**2)
    Fy2 = -(k_contact_in*F2_in)*(y_pos2-elips_trans[1])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[1]**2)
    Fz2 = -(k_contact_in*F2_in)*(z_pos2-elips_trans[2])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[2]**2)

    # applying contact forces to contact points in thorax frame
    cont_force1 = [(contact_point1,frame_ground.x*Fx1+frame_ground.y*Fy1+frame_ground.z*Fz1)]
    cont_force2 = [(contact_point2,frame_ground.x*Fx2+frame_ground.y*Fy2+frame_ground.z*Fz2)]
    CONT = cont_force1+cont_force2
    
    KM = me.KanesMethod(frame_ground, q_ind=q, u_ind=u, kd_eqs=kindeq)
    (fr, frstar) = KM.kanes_equations(BODY, (FG+DAMP+CONT))
    MM = KM.mass_matrix_full
    FO = KM.forcing_full
    xdot = (KM.q.col_join(KM.u)).diff()
    print('equations created')
    
    # if gen_matlab_functions == 1:
    
    #     body_constants = {'I_': inertia,'mass_':mass,'com_':com,'offset_':offset,'c': c,'g': g}
    #     other_constants = {'offset_humerus_rot':list(offset_humerus_rot),'EL_rot_axis': list(EL_rot_axis),'PSY_rot_axis': list(PSY_rot_axis),
    #                         'k_contact_in': k_contact_in,'eps_in': eps_in,'contTS': list(contTS),
    #                         'contAI': list(contAI), 'elips_trans':list(elips_trans), 'elips_dim':list(elips_dim),
    #                         'k_contact_out': k_contact_out,'eps_out': eps_out,'second_elips_scale':second_elips_scale, 'offset_thorax': list(offset_thorax)}
    #     usubs = sp.symbols('u1:11')
    #     qsubs = sp.symbols('q1:11')

    # # sympy dynamicsymbols has to be substituted with symbols (so it can be printed in octave_code)

    #     subs_q = {q[i]: qsubs[i] for i in range(len(q))}
    #     subs_u = {u[i]: usubs[i] for i in range(len(u))}
    #     mm = me.msubs(KM.mass_matrix_full,subs_q,subs_u)
    #     fo = me.msubs(KM.forcing_full,subs_q,subs_u)

    #     MatlabFunction(function = mm,
    #                    fun_name = 'mm_eul',assignto = 'mm',
    #                    coordinates = qsubs,
    #                    speeds = usubs,
    #                    inputs = [],
    #                    body_constants = body_constants,
    #                    segments = segment,
    #                    other_constants=other_constants,
    #                    muscle_constants = {},
    #                    parameters = [first_elips_scale],
    #                    folder = 'euler')
    #     MatlabFunction(function = fo,
    #                    fun_name = 'fo_eul',assignto = 'fo',
    #                    coordinates = qsubs,
    #                    speeds = usubs,
    #                    inputs = [],
    #                    body_constants = body_constants,
    #                    segments = segment,
    #                    other_constants=other_constants,
    #                    muscle_constants = {},
    #                    parameters = [first_elips_scale],
    #                    folder = 'euler')

    #     print('matlab functions generated')

    return MM,FO,q,u,fr,frstar,kindeq,xdot,first_elips_scale,elips_trans
    # return q

def muscle_force(act, lmt, vce, fmax, lceopt, lslack, vmax):
    lm = lmt - lslack
    
    f_gauss = 0.25
    kpe = 5
    epsm0 = 0.6
    fpe = (sp.exp(kpe*(lm / lceopt - 1)/epsm0)-1)/(sp.exp(kpe)-1)
    flce = (sp.exp(-(lm / lceopt - 1)**2 / f_gauss))

    d1 = -0.318
    d2 = -8.149
    d3 = -0.374
    d4 = 0.886
    vmax_norm = vmax * lceopt
    vnorm = vce[0,0]/vmax_norm

    fvce = d1 * sp.log(d2 * vnorm + d3 + sp.sqrt((d2 * vnorm + d3)**2 + 1)) + d4

    force = (flce * fvce * act +  fpe) * fmax
    # force = (flce * act +  fpe) * fmax
    
    return force

def act_dynamics(a,u,t_act = 0.05,t_deact = 0.06):

    da = (u/t_act + (1-u)/t_deact)*(u-a)

    return da

def conoid_force(lmt, lopt, k, eps):
    d = lmt-lopt
    F = k / 2 * (d + sp.sqrt(d**2 + eps))
    
    return F

def dquatdt(quat,w):
    # print(sp.shape(sp.Matrix([w])))
    res = 1/2 * G(quat).T * sp.Matrix([w]).T
    return res

def invQuat_sp(q):
    return sp.Matrix([q[0],-q[1],-q[2],-q[3]])

def invJtrans(quat):
    q1 = quat[0];
    q2 = quat[1];
    q3 = quat[2];
    q4 = quat[3];
    res = sp.Matrix([[ q1/2,  q4/2, -q3/2],
                     [-q4/2,  q1/2,  q2/2],
                     [ q3/2, -q2/2,  q1/2]])

    return res

def G(quat):
    Q0 = quat[0];
    Q1 = quat[1];
    Q2 = quat[2];
    Q3 = quat[3];
    res = np.array([[-Q1, Q0, Q3, -Q2],
                    [-Q2,-Q3, Q0, Q1],
                    [-Q3, Q2, -Q1, Q0]])
    return res

def switch_dict(my_dict):
    my_new_dict = dict(zip(my_dict.values(), my_dict.keys()))
    return my_new_dict

def T_trans(vec):
    trans_y = sp.Matrix([[1,0,0,vec[0]],
                         [0,1,0,vec[1]],
                         [0,0,1,vec[2]],
                         [0,0,0,1]])
    return trans_y

def R_x(phix):
    rot_phix = sp.Matrix([[1,0          ,0           ,0],
                         [0,sp.cos(phix),-sp.sin(phix),0],
                         [0,sp.sin(phix), sp.cos(phix),0],
                         [0,0          ,0           ,1]])
    return rot_phix

def R_y(phiy):
    rot_phiy = sp.Matrix([[sp.cos(phiy) ,0,sp.sin(phiy),0],
                         [0           ,1,0          ,0],
                         [-sp.sin(phiy),0,sp.cos(phiy),0],
                         [0           ,0,0           ,1]])
    return rot_phiy

def R_z(phiz):
    rot_phiz = sp.Matrix([[sp.cos(phiz),-sp.sin(phiz),0,0],
                        [sp.sin(phiz), sp.cos(phiz),0,0],
                        [0           ,0            ,1,0],
                        [0           ,0            ,0,1]])
    return rot_phiz

def position(vec):
    r = sp.Matrix([[vec[0]],[vec[1]],[vec[2]],[1]])
    return r

def YZX_seq(phi_vec):
    
    res = R_y(phi_vec[0]) * R_z(phi_vec[1]) * R_x(phi_vec[2])
    
    return res

def analytic_length_eul(origin, insertion, O_pos, I_pos, q, model):
    jnts = model['model']['joints'].item()
    offset_thorax = jnts[0,1]['location'].item()[0]
    offset_clavicle = jnts[0,4]['location'].item()[0]

    if origin == 'thorax' and insertion == 'clavicle_r':
        O = position(O_pos);
        I = T_trans(offset_thorax) * YZX_seq(q[0:3]) * position(I_pos);
        
    elif origin == 'thorax' and insertion == 'scapula_r':
        O = position(O_pos);
        RW_C = T_trans(offset_thorax) * YZX_seq(q[0:3]);
        TC_S = T_trans(offset_clavicle);
        RC_S = YZX_seq(q[3:6]);
        I = RW_C * TC_S * RC_S * position(I_pos);
        
    elif origin == 'clavicle_r' and insertion == 'scapula_r':
        O = position(O_pos);
        TC_S = T_trans(offset_clavicle);
        RC_S = YZX_seq(q[3:6]);
        I = TC_S * RC_S * position(I_pos);

    muscle_length = sp.sqrt((O[0] - I[0])**2 + (O[1] - I[1])**2 + (O[2] - I[2])**2);
    
    return muscle_length

def analytic_length_quat(origin, insertion, O_pos, I_pos, q, model):
    jnts = model['model']['joints'].item()
    offset_thorax = jnts[0,1]['location'].item()[0]
    offset_clavicle = jnts[0,4]['location'].item()[0]

    if origin == 'thorax' and insertion == 'clavicle_r':
        O = position(O_pos);
        I = T_trans(offset_thorax) * Qrm(q[0:4]) * position(I_pos);
        
    elif origin == 'thorax' and insertion == 'scapula_r':
        O = position(O_pos);
        RW_C = T_trans(offset_thorax) * Qrm(q[0:4]);
        TC_S = T_trans(offset_clavicle);
        RC_S = Qrm(q[4:8]);
        I = RW_C * TC_S * RC_S * position(I_pos);
        
    elif origin == 'clavicle_r' and insertion == 'scapula_r':
        O = position(O_pos);
        TC_S = T_trans(offset_clavicle);
        RC_S = Qrm(q[4:8]);
        I = TC_S * RC_S * position(I_pos);

    muscle_length = sp.sqrt((O[0] - I[0])**2 + (O[1] - I[1])**2 + (O[2] - I[2])**2)
    
    return muscle_length

def Qrm(q):
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]
    res =  sp.Matrix([[1-2*(y**2+z**2), 2*(x*y-z*w), 2*(x*z+y*w),0],
                      [2*(x*y+z*w), 1-2*(x**2+z**2), 2*(y*z-x*w),0],
                      [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x**2+y**2),0],
                      [0           ,0          ,0             ,1]])

    return res

def joint_spring_euler(q,model_params_struct,initCond_name):
    q0_joint = model_params_struct['params'][initCond_name][0,0]['initCondEul'].item()
    # joint_stiffnes = model_params_struct['params']['model'].item()['joint_stiffnes'][0,0].item()
    joint_stiffness = sp.symbols('joint_stiffness')
    # q0_joint = sp.symbols('q0joint1:7')
    SCx_spring = joint_stiffness*(-q[0]+q0_joint[0].item())
    SCy_spring = joint_stiffness*(-q[1]+q0_joint[1].item())
    SCz_spring = joint_stiffness*(-q[2]+q0_joint[2].item())
    ACx_spring = joint_stiffness*(-q[3]+q0_joint[3].item())
    ACy_spring = joint_stiffness*(-q[4]+q0_joint[4].item())
    ACz_spring = joint_stiffness*(-q[5]+q0_joint[5].item())

    res = sp.Matrix([SCx_spring-ACx_spring,SCy_spring-ACy_spring,SCz_spring-ACz_spring,ACx_spring,ACy_spring,ACz_spring,0,0,0,0])
    # res = sp.Matrix([SCx_spring,SCy_spring,SCz_spring,0,0,0,0,0,0,0])

    return res

def joint_spring_quat(q,Qeq):
    # Qeq = model_params_struct['params'][initCond_name][0,0]['initCondQuat'].item()
    # Qeq = sp.Matrix([1,0,0,0])
    # joint_stiffness = sp.symbols('joint_stiffness')
    joint_stiffness = 2.0
    Qdif = mulQuat_sp(Qinv_sp(Qeq),q)

    angle = 2*sp.atan2(sp.sqrt(Qdif[1]**2 + Qdif[2]**2 + Qdif[3]**2),Qdif[0])
    scale = 1/sp.sqrt(1-Qdif[0]**2 + 1e-3)
    axis = scale * sp.Matrix([Qdif[1],Qdif[2],Qdif[3]])
    res = -axis * (sp.exp(angle) - 1) * joint_stiffness

    return res




def Qinv_sp(q):
    res = sp.Matrix([q[0],-q[1],-q[2],-q[3]])

    return res

def MatlabFunction(function,fun_name,assignto,coordinates,speeds,inputs,body_constants,segments,other_constants,muscle_constants,parameters,folder):
    list_of_variables = [speeds,inputs, body_constants, muscle_constants, parameters]
    list_of_variables_names = [',u',',inputs',',model',',fmax, lceopt, lslack',',parameters']
    
    text_file = open(f"{folder}/{fun_name}.m","w")
    with open(f"{folder}/{fun_name}.m","w") as text_file:
        # function .. = .. ()
        header = f"function {assignto} = {fun_name}(t,q" #,u,act,model,opt_var)
        
        for i, current_list in enumerate(list_of_variables):
            if current_list:
                header += list_of_variables_names[i]
        header += ")"
        
        print(header,file=text_file)
        
        # acces to  constants in model struct

        for var_name, var_tuple in body_constants.items():
            n = 0
            k = 0
            if not type(var_tuple) == list:
                print(f'{var_name} = model.{var_name};',file=text_file)
            
            elif len(var_tuple) == len(segments): 
                for i, elem in enumerate(var_tuple):
                    str_res = f'{var_name}{segments[i]} = model.{var_name}{segments[i]};'
                    print(str_res,file=text_file)
            else:
                for i, elem in enumerate(var_tuple):
                    str_res = f'{elem} = model.{var_name}{segments[n]}(:,{k+1});'
                    print(str_res,file=text_file)
                    k += 1
                    if (i + 1) % 3 == 0:
                        k = 0
                        n += 1
                    else:
                        pass
                    
        for var_name, var_list in muscle_constants.items():
            for i, elem in enumerate(var_list):
                str_res = f'{elem} = {var_name}({i+1});'
                print(str_res,file=text_file)
        
        for var_name, var_tuple in other_constants.items():
            if not type(var_tuple) == list:
                print(f'{var_name} = model.{var_name};',file=text_file)
            else:
                for i, elem in enumerate(var_tuple):
                    str_res = f'{elem} = model.{var_name}(:,{i+1});'
                    print(str_res,file=text_file)
        
        for i, coord in enumerate(coordinates):
            print(f"{str(coord)} = q({i+1},:);", file=text_file)
        for i, speed in enumerate(speeds):
            print(f"{str(speed)} = u({i+1},:);", file=text_file)
                    
        for i, act in enumerate(inputs):
            print(f"{str(act)} = inputs({i+1},:);", file=text_file)
            
        for i, param in enumerate(parameters):
            print(f"{param} = parameters({i+1});", file=text_file)

                
        sub_exprs, simplified_rhs = sp.cse(function)
        for var, expr in sub_exprs:
            print('%s = %s;' % (sp.octave_code(var),sp.octave_code(expr)),file=text_file)
        print('%s' % sp.octave_code(sp.Matrix([simplified_rhs]), assign_to = assignto),file=text_file)
        
        
        
        
        
###################################################################################################
                #QUAT_U0_STATE#
###################################################################################################


def create_eoms_u0state(OS_struct, derive = 'symbolic',gen_matlab_functions = None):
    # symbols
    t = sp.symbols('t')

    states = ['q','w','u0']
    segment = ['clavicula','scapula','humerus','ulna','radius','hand']
    segment_index = [3,6,9,10,11,12]
    joints = ['quat','quat','quat','rotaxis','weld','weld']
    joint_index = [4,7,10,11,12,13]
    inertia = []
    mass = []
    com = []
    offset = []
    rot_offset = []
    q = []
    w = []
    u0 = []
    frame = []
    point_offset = []
    masscenter = []
    inertia_elem = []

    if derive == 'symbolic':
        g,c = sp.symbols('g,c')  # 
    elif derive == 'numeric':
        g = 9.81
        c = 0.5


    for i,seg in enumerate(segment):
        for idat,j in enumerate(('1','2','3')):
            if derive == 'symbolic':
                inertia.append(sp.symbols('I_' + seg + '_' + j))
                com.append(sp.symbols('com_' + seg + '_' + j))
                offset.append(sp.symbols('offset_' + seg + '_' + j))
            elif derive == 'numeric':
                com.append(OS_struct['model']['segments'][0,0][0,segment_index[i]]['mass_center'][0,0][0,idat].item())
                inertia_diag = np.diag(OS_struct['model']['segments'][0,0][0,segment_index[i]]['inertia'][0,0])
                inertia.append(inertia_diag[idat].item())
                try:
                    offset.append(OS_struct['model']['joints'][0,0][0,joint_index[i]]['location'][0,0][0,idat].item())
                except:
                    offset.append(0)

        if joints[i] == 'quat':
            for j in ('1','2','3'):
                w.append(me.dynamicsymbols('w' + j + '_' + seg))
            for j in ('0','1','2','3'):
                q.append(me.dynamicsymbols('q' + j+ '_' + seg))
            u0.append(me.dynamicsymbols('u0'+ '_' + seg))
        elif joints[i] == 'rotaxis':
            w.append(me.dynamicsymbols('w_'+ seg))
            q.append(me.dynamicsymbols('q_' + seg))
        else:
            pass

        if derive == 'symbolic':
            mass.append(sp.symbols('mass_'+seg))
        elif derive == 'numeric':
            mass.append(OS_struct['model']['segments'][0,0][0,segment_index[i]]['mass'][0,0].item())

        frame.append(me.ReferenceFrame('frame_' + str(seg)))
        point_offset.append(me.Point('point_offset_' + str(seg)))
        masscenter.append(me.Point('masscenter_' + str(seg)))

    # inertial frame and point
    frame_ground = me.ReferenceFrame('frame_ground')
    point_ground = me.Point('point_ground')
    point_ground.set_vel(frame_ground,0)

    offset_ground = me.Point('offset_ground')
    if derive == 'symbolic':
        offset_thorax = sp.symbols('offset_thorax_1:4')
    elif derive == 'numeric':
        offset_thorax = OS_struct['model']['joints'][0,0][0,1]['location'][0,0][0]

    offset_ground.set_pos(point_ground, offset_thorax[0]*frame_ground.x 
                          + offset_thorax[1]*frame_ground.y + offset_thorax[2]*frame_ground.z)
    offset_ground.set_vel(frame_ground,0)

    #rotate first body
    frame[0].orient(frame_ground, 'Quaternion', q[0:4])
    N_w_A = (frame[0].ang_vel_in(frame_ground))
    kinematical = sp.Matrix([
        u0[0] - q[0].diff(t),
        w[0] - N_w_A.dot(frame[0].x),
        w[1] - N_w_A.dot(frame[0].y),
        w[2] - N_w_A.dot(frame[0].z),
    ])

    for i in range(1,3):
        frame[i].orient(frame[i-1], 'Quaternion', q[0+i*4:4+i*4])
        N_w_A = (frame[i].ang_vel_in(frame[i-1]))
        kinematical = kinematical.col_join(sp.Matrix([
                            u0[i] - q[i*4].diff(t),
                             w[i*3] - N_w_A.dot(frame[i].x),
                             w[i*3+1] - N_w_A.dot(frame[i].y),
                             w[i*3+2] - N_w_A.dot(frame[i].z),
                                ]))

    frame[0].set_ang_vel(frame_ground,w[0]*frame[0].x + w[1]*frame[0].y + w[2]*frame[0].z)

    for i in range(1,3):
        frame[i].set_ang_vel(frame[i-1],w[0+i*3]*frame[i].x + w[1+i*3]*frame[i].y + w[2+i*3]*frame[i].z)

    # set masscenter of first body
    masscenter[0].set_pos(offset_ground,com[0]*frame[0].x + com[1]*frame[0].y + com[2]*frame[0].z)
    masscenter[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set offset of first joint in first body
    point_offset[0].set_pos(offset_ground,offset[0]*frame[0].x + offset[1]*frame[0].y + offset[2]*frame[0].z)
    point_offset[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set gravity force and damping of first body
    FG = [(masscenter[0], -mass[0] * g * frame_ground.y)]
    DAMP = [(frame[0], -c*(w[0]*frame[0].x+w[1]*frame[0].y+w[2]*frame[0].z))]
    # iterate over segments 2:end (first body is already done)
    for i in range(1,3):

        # set masscenter points 
        masscenter[i].set_pos(point_offset[i-1],com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
        masscenter[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set gravity force in masscenter (-y direction in frame_ground)
        FG.append((masscenter[i], -mass[i] * g * frame_ground.y))

        # set offsent points (where the next joint is)
        point_offset[i].set_pos(point_offset[i-1],offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
        point_offset[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set damping in joints (c * angular_velocity)
        damping = -c*(w[0+i*3]*frame[i].x+w[1+i*3]*frame[i].y+w[2+i*3]*frame[i].z)

        # apply damping in frame, opposite moment is applied in previous frame (action and reaction)
        DAMP.append((frame[i], damping))
        DAMP.append((frame[i-1], -damping))
        # set kinematic differential equations (q_dot = f(q,u))

    # symbols for ulna
    ulna_rot_frame = me.ReferenceFrame('ulna_rot_frame')

    if derive == 'symbolic':
        offset_humerus_rot = sp.symbols('offset_humerus_rot_1:4')
        EL_rot_axis = sp.symbols('EL_rot_axis_1:4')
        PSY_rot_axis = sp.symbols('PSY_rot_axis_1:4')
    if derive == 'numeric':
        offset_humerus_rot = OS_struct['model']['joints'][0,0][0,10]['orientation'][0,0][0]
        EL_rot_axis = OS_struct['model']['joints'][0,0][0,10]['r1_axis'][0,0][0]
        PSY_rot_axis = OS_struct['model']['joints'][0,0][0,11]['r1_axis'][0,0][0]


    # offset frame rotated in humerus frame (frame[2])
    ulna_rot_frame.orient_axis(frame[2],frame[2].z,offset_humerus_rot[2])

    # ulna and elbow joint
    frame[3].orient_axis(ulna_rot_frame,ulna_rot_frame.x*EL_rot_axis[0]
                         +ulna_rot_frame.y*EL_rot_axis[1]+ulna_rot_frame.z*EL_rot_axis[2],
                         q[12])

    frame[3].set_ang_vel(ulna_rot_frame,
                          w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2]))

    masscenter[3].set_pos(point_offset[2],com[0+3*3]*frame[3].x + com[1+3*3]*frame[3].y + com[2+3*3]*frame[3].z)
    masscenter[3].v2pt_theory(point_offset[2],frame_ground,frame[3])
    FG.append((masscenter[3], -mass[3] * g * frame_ground.y))

    DAMP.append(((frame[3]),-c*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                            +ulna_rot_frame.z*EL_rot_axis[2])))
    DAMP.append((ulna_rot_frame,c*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2])))
    kinematical = kinematical.col_join(sp.Matrix([w[9]-q[12].diff(t)]))

    point_offset[3].set_pos(point_offset[2],offset[0+3*3]*frame[3].x + offset[1+3*3]*frame[3].y + offset[2+3*3]*frame[3].z)
    point_offset[3].v2pt_theory(point_offset[2],frame_ground,frame[3]);

    # radius and PSY joint (weld joint for now)
    frame[4].orient_axis(frame[3],frame[3].z,0)
    frame[4].set_ang_vel(frame[3],0)
    masscenter[4].set_pos(point_offset[3],com[0+4*3]*frame[4].x + com[1+4*3]*frame[4].y + com[2+4*3]*frame[4].z)
    masscenter[4].v2pt_theory(point_offset[3],frame_ground,frame[4])
    FG.append((masscenter[4], -mass[4] * g * frame_ground.y))
    frame[4].ang_vel_in(frame[3])

    point_offset[4].set_pos(point_offset[3],offset[0+4*3]*frame[3].x + offset[1+4*3]*frame[3].y + offset[2+4*3]*frame[3].z)
    point_offset[4].v2pt_theory(point_offset[3],frame_ground,frame[3])

    # hand
    frame[5].orient_axis(frame[4],frame[4].z,0)
    frame[5].set_ang_vel(frame[4],0)
    masscenter[5].set_pos(point_offset[4],com[0+5*3]*frame[5].x + com[1+5*3]*frame[5].y + com[2+5*3]*frame[5].z)
    masscenter[5].v2pt_theory(point_offset[4],frame_ground,frame[5])
    FG.append((masscenter[5], -mass[5] * g * frame_ground.y))

    BODY = []

    for i in range(len(segment)):
        # set inertias of each body and create RigidBodies
        I = me.inertia(frame[i], inertia[0+i*3], inertia[1+i*3], inertia[2+i*3])
        BODY.append(me.RigidBody('body' + str(i), masscenter[i], frame[i], mass[i], (I, masscenter[i])))

    # Contact between scapula and thorax
    if derive == 'symbolic':
        contTS = sp.symbols('contTS_1:4')
        contAI = sp.symbols('contAI_1:4')
        elips_trans = sp.symbols('elips_trans_1:4')
        elips_dim = sp.symbols('elips_dim_1:4')
        k_contact_in, eps_in = sp.symbols('k_contact_in eps_in')
        k_contact_out, eps_out = sp.symbols('k_contact_out eps_out')
        # first_elips_scale = sp.Symbol('first_elips_scale')
        second_elips_scale = sp.Symbol('second_elips_scale')
    elif derive == 'numeric':
        # k_contact_in = data_struct['k_contact_in'][0,0].item()
        k_contact_out = 0
        # eps_in = data_struct['eps_in'][0,0].item()
        
        k_contact_in = OS_struct['model']['scap_thorax_k'][0,0].item()
        # eps_in = OS_struct['model']['scap_thorax_eps'][0,0].item()
        eps_in = 0.01
        eps_out = 1e-3
        # first_elips_scale = model_params_struct['params'][initCond_name][0,0]['first_elips_scale'][0,0].item()
        # first_elips_scale = sp.symbols('first_elips_scale1:4')
        # elips_dim = OS_struct['model']['thorax_dim_optimized'].item()[0]
        elips_trans = OS_struct['model']['thoracic_wall'][0,0].item()[0][0]
        elips_dim = OS_struct['model']['thoracic_wall'][0,0].item()[1][0]
        contTS = OS_struct['model']['TScontact'][0,0].item()[0][0]
        contAI = OS_struct['model']['AIcontact'][0,0].item()[0][0]
        first_elips_scale = [1.0,1.0,1.0]
        # elips_trans = sp.symbols('elips_trans1:4')
        second_elips_scale = [1.0,1.0,1.0]

    # contact points 
    contact_point1 = me.Point('CP1')
    contact_point1.set_pos(point_offset[0],contTS[0]*frame[1].x+contTS[1]*frame[1].y  +contTS[2]*frame[1].z)
    # contact_point1.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point1.v2pt_theory(point_offset[0],frame_ground,frame[1])

    contact_point2 = me.Point('CP2')
    contact_point2.set_pos(point_offset[0],contAI[0]*frame[1].x+contAI[1]*frame[1].y  +contAI[2]*frame[1].z)
    # contact_point2.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point2.v2pt_theory(point_offset[0],frame_ground,frame[1])

    ## contact forces

    # Distances between contact points and thorax frame
    x_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.x)
    y_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.y)
    z_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.z)
    x_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.x)
    y_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.y)
    z_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.z)

    # Contact forces
    elips_dim_scaled = []
    elips_dim_scaled.append(elips_dim[0]*first_elips_scale[0])
    elips_dim_scaled.append(elips_dim[1]*first_elips_scale[1])
    elips_dim_scaled.append(elips_dim[2]*first_elips_scale[2])
    # print(elips_dim_scaled)
    f1_in = ((x_pos1-elips_trans[0])/(elips_dim_scaled[0]))**2+((y_pos1-elips_trans[1])/(elips_dim_scaled[1]))**2+((z_pos1-elips_trans[2])/(elips_dim_scaled[2]))**2-1
    # f1_out = ((x_pos1-elips_trans[0])/(second_elips_scale*elips_dim[0]))**2+((y_pos1-elips_trans[1])/(second_elips_scale*elips_dim[1]))**2+((z_pos1-elips_trans[2])/(second_elips_scale*elips_dim[2]))**2-1
    F1_in = 1/2*(f1_in-sp.sqrt(f1_in**2+eps_in**2))
    # F1_out = 1/2*(f1_out+sp.sqrt(f1_out**2+eps_out**2))
    Fx1 = -(k_contact_in*F1_in)*(x_pos1-elips_trans[0])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[0]**2)
    Fy1 = -(k_contact_in*F1_in)*(y_pos1-elips_trans[1])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[1]**2)
    Fz1 = -(k_contact_in*F1_in)*(z_pos1-elips_trans[2])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[2]**2)

    f2_in = ((x_pos2-elips_trans[0])/(elips_dim_scaled[0]))**2+((y_pos2-elips_trans[1])/(elips_dim_scaled[1]))**2+((z_pos2-elips_trans[2])/(elips_dim_scaled[2]))**2-1
    # f2_out = ((x_pos2-elips_trans[0])/(second_elips_scale*elips_dim[0]))**2+((y_pos2-elips_trans[1])/(second_elips_scale*elips_dim[1]))**2+((z_pos2-elips_trans[2])/(second_elips_scale*elips_dim[2]))**2-1
    F2_in = 1/2*(f2_in-sp.sqrt(f2_in**2+eps_in**2))
    # F2_out = 1/2*(f2_out+sp.sqrt(f2_out**2+eps_out**2))
    Fx2 = -(k_contact_in*F2_in)*(x_pos2-elips_trans[0])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[0]**2)
    Fy2 = -(k_contact_in*F2_in)*(y_pos2-elips_trans[1])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[1]**2)
    Fz2 = -(k_contact_in*F2_in)*(z_pos2-elips_trans[2])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[2]**2)

    # applying contact forces to contact points in thorax frame
    cont_force1 = [(contact_point1,frame_ground.x*Fx1+frame_ground.y*Fy1+frame_ground.z*Fz1)]
    cont_force2 = [(contact_point2,frame_ground.x*Fx2+frame_ground.y*Fy2+frame_ground.z*Fz2)]
    CONT = cont_force1+cont_force2


    q_dep = [q[0],q[4],q[8]]
    q_ind = q[1:4]+q[5:8]+q[9:]
    holonomic = sp.Matrix([[q[0]**2+q[1]**2+q[2]**2+q[3]**2-1],
                             [q[4]**2+q[5]**2+q[6]**2+q[7]**2-1],
                             [q[8]**2+q[9]**2+q[10]**2+q[11]**2-1]])
    KM = me.KanesMethod(frame_ground, q_ind=q_ind, u_ind=w, kd_eqs=kinematical,
                        q_dependent = q_dep, u_dependent = u0, 
                        configuration_constraints=holonomic,
                        velocity_constraints=holonomic.diff(t))

    (fr, frstar) = KM.kanes_equations(BODY, (FG+DAMP+CONT))
    MM = KM.mass_matrix_full
    FO = KM.forcing_full
    xdot = (KM.q.col_join(KM.u)).diff()
    
    return MM,FO,q,w,u0,fr,frstar,kinematical,xdot,holonomic,first_elips_scale,elips_trans


def mulQuat_sp(qa,qb):
    res = sp.Matrix([[qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3]],
                     [qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2]],
                     [qa[0]*qb[2] - qa[1]*qb[3] + qa[2]*qb[0] + qa[3]*qb[1]],
                     [qa[0]*qb[3] + qa[1]*qb[2] - qa[2]*qb[1] + qa[3]*qb[0]]])
    return res

def create_eoms_eul_torqueDriven(model_params_struct, derive = 'symbolic',gen_matlab_functions = None):
    # symbols
    t = sp.symbols('t')
    g,c = sp.symbols('g,c')  # 
    states = ['q','u']
    segment = ['clavicula','scapula','humerus','ulna','radius','hand']
    joints = ['YZX','YZX','YZY','rotaxis','weld','weld']
    inertia = []
    mass = []
    com = []
    offset = []
    rot_offset = []
    q = []
    u = []
    rot = [] # for angular velocities
    frame = []
    point_offset = []
    masscenter = []
    inertia_elem = []
    Tau = me.dynamicsymbols('Tau1:11')
        
    if derive == 'symbolic':
        g,c = sp.symbols('g,c')  # 
    elif derive == 'numeric':
        data_struct = model_params_struct['params']['model'][0,0]
        g = data_struct['g'][0,0].item()
        c = data_struct['c'][0,0].item()

    for i,seg in enumerate(segment):
        for idat,j in enumerate(('1','2','3')):
            if derive == 'symbolic':
                inertia.append(sp.symbols('I_' + seg + '_' + j))
                com.append(sp.symbols('com_' + seg + '_' + j))
                offset.append(sp.symbols('offset_' + seg + '_' + j))
            elif derive == 'numeric':
                inertia.append(data_struct['I_' + seg][0,0][0,idat])
                com.append(data_struct['com_' + seg][0,0][0,idat])
                offset.append(data_struct['offset_' + seg][0,0][0,idat])

        if joints[i] == 'YZX' or joints[i] == 'YZY':
            for j in ('1','2','3'):
                q.append(me.dynamicsymbols('q' + j + '_' + seg))
                u.append(me.dynamicsymbols('u' + j + '_' + seg))
        elif joints[i] == 'rotaxis':
            q.append(me.dynamicsymbols('q_' + seg))
            u.append(me.dynamicsymbols('u_' + seg))
        else:
            pass

        if derive == 'symbolic':
            mass.append(sp.symbols('mass_'+seg))
        elif derive == 'numeric':
            mass.append(data_struct['mass_' + seg][0,0].item())

        frame.append(me.ReferenceFrame('frame_' + str(seg)))
        point_offset.append(me.Point('point_offset_' + str(seg)))
        masscenter.append(me.Point('masscenter_' + str(seg)))


    # inertial frame and point
    frame_ground = me.ReferenceFrame('frame_ground')
    point_ground = me.Point('point_ground')
    point_ground.set_vel(frame_ground,0)

    offset_ground = me.Point('offset_ground')
    if derive == 'symbolic':
        offset_thorax = sp.symbols('offset_thorax_1:4')
    elif derive == 'numeric':
        offset_thorax = []
        for i in range(3):
            offset_thorax.append(data_struct['offset_thorax'][0,0][0,i].item())

    offset_ground.set_pos(point_ground, offset_thorax[0]*frame_ground.x 
                          + offset_thorax[1]*frame_ground.y + offset_thorax[2]*frame_ground.z)
    offset_ground.set_vel(frame_ground,0)

    #rotate first body
    frame[0].orient_body_fixed(frame_ground, q[0:3], rotation_order = joints[0])
    rot.append(frame[0].ang_vel_in(frame_ground))

    # set masscenter of first body
    masscenter[0].set_pos(offset_ground,com[0]*frame[0].x + com[1]*frame[0].y + com[2]*frame[0].z)
    masscenter[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set offset of first joint in first body
    point_offset[0].set_pos(offset_ground,offset[0]*frame[0].x + offset[1]*frame[0].y + offset[2]*frame[0].z)
    point_offset[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set gravity force and damping of first body
    FG = [(masscenter[0], -mass[0] * g * frame_ground.y)]
    DAMP = [(frame[0], -c*(rot[0].dot(frame[0].x)*frame[0].x+rot[0].dot(frame[0].y)*frame[0].y+rot[0].dot(frame[0].z)*frame[0].z))]
    Torque = [(frame[0],10*Tau[0]*frame[0].x+10*Tau[1]*frame[0].y+10*Tau[2]*frame[0].z)]
    kindeq = []

    # iterate over segments 2:end (first body is already done)
    for i in range(1,3):

        # orient frame w.r.t. child's frame
        frame[i].orient_body_fixed(frame[i-1], q[0+i*3:3+i*3], rotation_order = joints[i])
        rot.append(frame[i].ang_vel_in(frame[i-1]))

        # set masscenter points 
        masscenter[i].set_pos(point_offset[i-1],com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
        masscenter[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set gravity force in masscenter (-y direction in frame_ground)
        FG.append((masscenter[i], -mass[i] * g * frame_ground.y))

        # set offsent points (where the next joint is)
        point_offset[i].set_pos(point_offset[i-1],offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
        point_offset[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set damping in joints (c * angular_velocity)
        damping = -c*(rot[i].dot(frame[i].x)*frame[i].x+rot[i].dot(frame[i].y)*frame[i].y+rot[i].dot(frame[i].z)*frame[i].z)

        # apply damping in frame, opposite moment is applied in previous frame (action and reaction)
        DAMP.append((frame[i], damping))
        DAMP.append((frame[i-1], -damping))
        Torque.append((frame[i],10*Tau[i*3]*frame[i].x+10*Tau[i*3+1]*frame[i].y+10*Tau[i*3+2]*frame[i].z))

    for i in range(9):
        kindeq.append(q[i].diff()-u[i])

    # symbols (or values) for ulna and radius
    ulna_rot_frame = me.ReferenceFrame('ulna_rot_frame')
    
    if derive == 'symbolic':
        offset_humerus_rot = sp.symbols('offset_humerus_rot_1:4')
        EL_rot_axis = sp.symbols('EL_rot_axis_1:4')
        PSY_rot_axis = sp.symbols('PSY_rot_axis_1:4')
    if derive == 'numeric':
        offset_humerus_rot = []
        EL_rot_axis = []
        PSY_rot_axis = []

        for i in range(3):
            offset_humerus_rot.append(data_struct['offset_humerus_rot'][0,0][0,i].item())
            EL_rot_axis.append(data_struct['EL_rot_axis'][0,0][0,i].item())
            PSY_rot_axis.append(data_struct['PSY_rot_axis'][0,0][0,i].item())

    # offset frame rotated in humerus frame (frame[2])
    ulna_rot_frame.orient_axis(frame[2],frame[2].z,offset_humerus_rot[2])

    # ulna and elbow joint
    frame[3].orient_axis(ulna_rot_frame,ulna_rot_frame.x*EL_rot_axis[0]
                         +ulna_rot_frame.y*EL_rot_axis[1]+ulna_rot_frame.z*EL_rot_axis[2],
                         q[9])

    frame[3].set_ang_vel(ulna_rot_frame,
                          u[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2]))

    masscenter[3].set_pos(point_offset[2],com[0+3*3]*frame[3].x + com[1+3*3]*frame[3].y + com[2+3*3]*frame[3].z)
    masscenter[3].v2pt_theory(point_offset[2],frame_ground,frame[3])
    FG.append((masscenter[3], -mass[3] * g * frame_ground.y))

    DAMP.append(((frame[3]),-c*u[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                            +ulna_rot_frame.z*EL_rot_axis[2])))
    DAMP.append((ulna_rot_frame,c*u[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2])))
    Torque.append((frame[3],10*Tau[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                            +ulna_rot_frame.z*EL_rot_axis[2])))
    kindeq.append(u[9]-q[9].diff(t))

    point_offset[3].set_pos(point_offset[2],offset[0+3*3]*frame[3].x + offset[1+3*3]*frame[3].y + offset[2+3*3]*frame[3].z)
    point_offset[3].v2pt_theory(point_offset[2],frame_ground,frame[3]);

    # radius and PSY joint (weld joint for now)
    frame[4].orient_axis(frame[3],frame[3].x*PSY_rot_axis[0]
                         +frame[3].y*PSY_rot_axis[1]+frame[3].z*PSY_rot_axis[2],
                         0) #q[10]
    masscenter[4].set_pos(point_offset[3],com[0+4*3]*frame[4].x + com[1+4*3]*frame[4].y + com[2+4*3]*frame[4].z)
    masscenter[4].v2pt_theory(point_offset[3],frame_ground,frame[4])
    FG.append((masscenter[4], -mass[4] * g * frame_ground.y))
    # DAMP.append(((frame[4]),-c*u[10]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
    #                         +ulna_rot_frame.z*EL_rot_axis[2])))
    # DAMP.append((ulna_rot_frame,c*u[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
    #                             +ulna_rot_frame.z*EL_rot_axis[2])))
    # kindeq.append(u[10]-q[10].diff(t))

    point_offset[4].set_pos(point_offset[3],offset[0+4*3]*frame[4].x + offset[1+4*3]*frame[4].y + offset[2+4*3]*frame[4].z)
    point_offset[4].v2pt_theory(point_offset[3],frame_ground,frame[4])

    # hand
    frame[5].orient_axis(frame[4],frame[4].z,0)
    frame[5].set_ang_vel(frame[4],0)
    masscenter[5].set_pos(point_offset[4],com[0+5*3]*frame[5].x + com[1+5*3]*frame[5].y + com[2+5*3]*frame[5].z)
    masscenter[5].v2pt_theory(point_offset[4],frame_ground,frame[5])
    FG.append((masscenter[5], -mass[5] * g * frame_ground.y))

    BODY = []

    for i in range(len(segment)):
        # set inertias of each body and create RigidBodies
        I = me.inertia(frame[i], inertia[0+i*3], inertia[1+i*3], inertia[2+i*3])
        BODY.append(me.RigidBody('body' + str(i), masscenter[i], frame[i], mass[i], (I, masscenter[i])))
    
    KM = me.KanesMethod(frame_ground, q_ind=q, u_ind=u, kd_eqs=kindeq)
    (fr, frstar) = KM.kanes_equations(BODY, (FG+DAMP+Torque))
    MM = KM.mass_matrix_full
    FO = KM.forcing_full
    xdot = (KM.q.col_join(KM.u)).diff()
    print('equations created')

    return MM,FO,q,u,fr,frstar,kindeq,xdot,Tau, Torque
    # return Torque

def create_eoms_quat_w_RF(OS_struct,weight = 0, derive = 'symbolic',gen_matlab_functions = None):
    
    # symbols
    t = sp.symbols('t')

    states = ['q','w','u0']
    segment = ['clavicula','scapula','humerus','ulna','radius','hand']
    segment_index = [3,6,9,10,11,12]
    joints = ['quat','quat','quat','rotaxis','weld','weld']
    joint_index = [4,7,10,11,12,13]
    inertia = []
    mass = []
    com = []
    offset = []
    rot_offset = []
    q = []
    w = []
    u0 = []
    frame = []
    point_offset = []
    masscenter = []
    inertia_elem = []
    

    if derive == 'symbolic':
        g,c = sp.symbols('g,c')  # 
    elif derive == 'numeric':
        g = 9.81
        c = 1.0


    for i,seg in enumerate(segment):
        for idat,j in enumerate(('1','2','3')):
            if derive == 'symbolic':
                inertia.append(sp.symbols('I_' + seg + '_' + j))
                com.append(sp.symbols('com_' + seg + '_' + j))
                offset.append(sp.symbols('offset_' + seg + '_' + j))
            elif derive == 'numeric':
                com.append(OS_struct['model']['segments'][0,0][0,segment_index[i]]['mass_center'][0,0][0,idat].item())
                inertia_diag = np.diag(OS_struct['model']['segments'][0,0][0,segment_index[i]]['inertia'][0,0])
                inertia.append(inertia_diag[idat].item())
                try:
                    offset.append(OS_struct['model']['joints'][0,0][0,joint_index[i]]['location'][0,0][0,idat].item())
                except:
                    offset.append(0)

        if joints[i] == 'quat':
            for j in ('1','2','3'):
                w.append(me.dynamicsymbols('w' + j + '_' + seg))
            for j in ('0','1','2','3'):
                q.append(me.dynamicsymbols('q' + j+ '_' + seg))
            u0.append(me.dynamicsymbols('u0'+ '_' + seg))
        elif joints[i] == 'rotaxis':
            w.append(me.dynamicsymbols('w_'+ seg))
            q.append(me.dynamicsymbols('q_' + seg))
        else:
            pass

        if derive == 'symbolic':
            mass.append(sp.symbols('mass_'+seg))
        elif derive == 'numeric':
            mass.append(OS_struct['model']['segments'][0,0][0,segment_index[i]]['mass'][0,0].item()/1.4) # why 1.5? just because dude
        
        frame.append(me.ReferenceFrame('frame_' + str(seg)))
        point_offset.append(me.Point('point_offset_' + str(seg)))
        masscenter.append(me.Point('masscenter_' + str(seg)))
    mass[-1] += weight
    print(mass)
    # inertial frame and point
    frame_ground = me.ReferenceFrame('frame_ground')
    point_ground = me.Point('point_ground')
    point_ground.set_vel(frame_ground,0)

    offset_ground = me.Point('offset_ground')
    if derive == 'symbolic':
        offset_thorax = sp.symbols('offset_thorax_1:4')
    elif derive == 'numeric':
        offset_thorax = OS_struct['model']['joints'][0,0][0,1]['location'][0,0][0]

    offset_ground.set_pos(point_ground, offset_thorax[0]*frame_ground.x 
                          + offset_thorax[1]*frame_ground.y + offset_thorax[2]*frame_ground.z)
    offset_ground.set_vel(frame_ground,0)

    #rotate first body
    frame[0].orient(frame_ground, 'Quaternion', q[0:4])

    for i in range(1,3):
        frame[i].orient(frame[i-1], 'Quaternion', q[0+i*4:4+i*4])
    kinematical =[]
    for i in range(0,3):
        kinematical.append(q[0+i*4].diff(t) - 0.5 * (-w[0+i*3]*q[1+i*4] - w[1+i*3]*q[2+i*4] - w[2+i*3]*q[3+i*4]))
        kinematical.append(q[1+i*4].diff(t) - 0.5 * (w[0+i*3]*q[0+i*4] + w[2+i*3]*q[2+i*4] - w[1+i*3]*q[3+i*4]))
        kinematical.append(q[2+i*4].diff(t) - 0.5 * (w[1+i*3]*q[0+i*4] - w[2+i*3]*q[1+i*4] + w[0+i*3]*q[3+i*4]))
        kinematical.append(q[3+i*4].diff(t) - 0.5 * (w[2+i*3]*q[0+i*4] + w[1+i*3]*q[1+i*4] - w[0+i*3]*q[2+i*4]))

    frame[0].set_ang_vel(frame_ground,w[0]*frame[0].x + w[1]*frame[0].y + w[2]*frame[0].z)

    for i in range(1,3):
        frame[i].set_ang_vel(frame[i-1],w[0+i*3]*frame[i].x + w[1+i*3]*frame[i].y + w[2+i*3]*frame[i].z)

    # set masscenter of first body
    masscenter[0].set_pos(offset_ground,com[0]*frame[0].x + com[1]*frame[0].y + com[2]*frame[0].z)
    masscenter[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set offset of first joint in first body
    point_offset[0].set_pos(offset_ground,offset[0]*frame[0].x + offset[1]*frame[0].y + offset[2]*frame[0].z)
    point_offset[0].v2pt_theory(offset_ground,frame_ground,frame[0])

    # set gravity force and damping of first body
    FG = [(masscenter[0], -mass[0] * g * frame_ground.y)]
    DAMP = [(frame[0], -c*(w[0]*frame[0].x+w[1]*frame[0].y+w[2]*frame[0].z))]
    # iterate over segments 2:end (first body is already done)
    for i in range(1,3):

        # set masscenter points 
        if i == 1:
            # calculate reaction forces in GH joint
            uaux = me.dynamicsymbols('uaux1 uaux2 uaux3')
            faux = me.dynamicsymbols('faux1 faux2 faux3')
            
            N_GH = point_offset[i-1].locatenew('N_GH', 0)
            N_GH.set_vel(frame_ground, N_GH.vel(frame_ground) + uaux[0]*frame[i].x + uaux[1]*frame[i].y + uaux[2]*frame[i].z)

            masscenter[i].set_pos(N_GH,com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
            masscenter[i].v2pt_theory(N_GH,frame_ground,frame[i])

            point_offset[i].set_pos(N_GH,offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
            point_offset[i].v2pt_theory(N_GH,frame_ground,frame[i])

            f_aux = [(point_offset[i-1], 800*(faux[0]*frame[i].x + faux[1]*frame[i].y + faux[2]*frame[i].z)),(N_GH, -800*(faux[0]*frame[i].x + faux[1]*frame[i].y + faux[2]*frame[i].z))]

            GH_spring_torque = joint_spring_quat(q[8:12],sp.Matrix([1,0,0,0]))
            spring = [(frame[i], GH_spring_torque[0]*frame[i].x + GH_spring_torque[1]*frame[i].y + GH_spring_torque[2]*frame[i].z),
                      (frame[i-1], -GH_spring_torque[0]*frame[i].x - GH_spring_torque[1]*frame[i].y - GH_spring_torque[2]*frame[i].z)]
        else:
            masscenter[i].set_pos(point_offset[i-1],com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
            masscenter[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

            # set offsent points (where the next joint is)
            point_offset[i].set_pos(point_offset[i-1],offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
            point_offset[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set gravity force in masscenter (-y direction in frame_ground)
        FG.append((masscenter[i], -mass[i] * g * frame_ground.y))
        # set damping in joints (c * angular_velocity)
        damping = -c*(w[0+i*3]*frame[i].x+w[1+i*3]*frame[i].y+w[2+i*3]*frame[i].z)

        # apply damping in frame, opposite moment is applied in previous frame (action and reaction)
        DAMP.append((frame[i], damping))
        DAMP.append((frame[i-1], -damping))
        # set kinematic differential equations (q_dot = f(q,u))

    # symbols for ulna
    ulna_rot_frame = me.ReferenceFrame('ulna_rot_frame')

    if derive == 'symbolic':
        offset_humerus_rot = sp.symbols('offset_humerus_rot_1:4')
        EL_rot_axis = sp.symbols('EL_rot_axis_1:4')
        PSY_rot_axis = sp.symbols('PSY_rot_axis_1:4')
    if derive == 'numeric':
        offset_humerus_rot = OS_struct['model']['joints'][0,0][0,10]['orientation'][0,0][0]
        EL_rot_axis = OS_struct['model']['joints'][0,0][0,10]['r1_axis'][0,0][0]
        PSY_rot_axis = OS_struct['model']['joints'][0,0][0,11]['r1_axis'][0,0][0]


    # offset frame rotated in humerus frame (frame[2])
    ulna_rot_frame.orient_axis(frame[2],frame[2].z,offset_humerus_rot[2])

    # ulna and elbow joint
    frame[3].orient_axis(ulna_rot_frame,ulna_rot_frame.x*EL_rot_axis[0]
                         +ulna_rot_frame.y*EL_rot_axis[1]+ulna_rot_frame.z*EL_rot_axis[2],
                         q[12])

    frame[3].set_ang_vel(ulna_rot_frame,
                          w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2]))

    masscenter[3].set_pos(point_offset[2],com[0+3*3]*frame[3].x + com[1+3*3]*frame[3].y + com[2+3*3]*frame[3].z)
    masscenter[3].v2pt_theory(point_offset[2],frame_ground,frame[3])
    FG.append((masscenter[3], -mass[3] * g * frame_ground.y))

    DAMP.append(((frame[3]),-c*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                            +ulna_rot_frame.z*EL_rot_axis[2])))
    DAMP.append((ulna_rot_frame,c*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                                +ulna_rot_frame.z*EL_rot_axis[2])))
    kinematical.append(w[9]-q[12].diff(t))

    point_offset[3].set_pos(point_offset[2],offset[0+3*3]*frame[3].x + offset[1+3*3]*frame[3].y + offset[2+3*3]*frame[3].z)
    point_offset[3].v2pt_theory(point_offset[2],frame_ground,frame[3]);

    # radius and PSY joint (weld joint for now)
    frame[4].orient_axis(frame[3],frame[3].z,0)
    frame[4].set_ang_vel(frame[3],0)
    masscenter[4].set_pos(point_offset[3],com[0+4*3]*frame[4].x + com[1+4*3]*frame[4].y + com[2+4*3]*frame[4].z)
    masscenter[4].v2pt_theory(point_offset[3],frame_ground,frame[4])
    FG.append((masscenter[4], -mass[4] * g * frame_ground.y))
    frame[4].ang_vel_in(frame[3])

    point_offset[4].set_pos(point_offset[3],offset[0+4*3]*frame[3].x + offset[1+4*3]*frame[3].y + offset[2+4*3]*frame[3].z)
    point_offset[4].v2pt_theory(point_offset[3],frame_ground,frame[3])

    # hand
    frame[5].orient_axis(frame[4],frame[4].z,0)
    frame[5].set_ang_vel(frame[4],0)
    masscenter[5].set_pos(point_offset[4],com[0+5*3]*frame[5].x + com[1+5*3]*frame[5].y + com[2+5*3]*frame[5].z)
    masscenter[5].v2pt_theory(point_offset[4],frame_ground,frame[5])
    FG.append((masscenter[5], -mass[5] * g * frame_ground.y))

    BODY = []

    for i in range(len(segment)):
        # set inertias of each body and create RigidBodies
        I = me.inertia(frame[i], inertia[0+i*3], inertia[1+i*3], inertia[2+i*3])
        BODY.append(me.RigidBody('body' + str(i), masscenter[i], frame[i], mass[i], (I, masscenter[i])))

    # Contact between scapula and thorax
    if derive == 'symbolic':
        contTS = sp.symbols('contTS_1:4')
        contAI = sp.symbols('contAI_1:4')
        elips_trans = sp.symbols('elips_trans_1:4')
        elips_dim = sp.symbols('elips_dim_1:4')
        k_contact_in, eps_in = sp.symbols('k_contact_in eps_in')
        k_contact_out, eps_out = sp.symbols('k_contact_out eps_out')
        # first_elips_scale = sp.Symbol('first_elips_scale')
        second_elips_scale = sp.Symbol('second_elips_scale')
    elif derive == 'numeric':
        # k_contact_in = data_struct['k_contact_in'][0,0].item()
        k_contact_out = 0
        # eps_in = data_struct['eps_in'][0,0].item()
        
        k_contact_in = OS_struct['model']['scap_thorax_k'][0,0].item()
        # eps_in = OS_struct['model']['scap_thorax_eps'][0,0].item()
        eps_in = 0.01
        eps_out = 1e-3
        # first_elips_scale = model_params_struct['params'][initCond_name][0,0]['first_elips_scale'][0,0].item()
        # first_elips_scale = sp.symbols('first_elips_scale1:4')
        # elips_dim = OS_struct['model']['thorax_dim_optimized'].item()[0]
        elips_trans = OS_struct['model']['thoracic_wall'][0,0].item()[0][0]
        elips_dim = OS_struct['model']['thoracic_wall'][0,0].item()[1][0]
        contTS = OS_struct['model']['TScontact'][0,0].item()[0][0]
        contAI = OS_struct['model']['AIcontact'][0,0].item()[0][0]
        first_elips_scale = [1.0,1.0,1.0]
        # elips_trans = sp.symbols('elips_trans1:4')
        second_elips_scale = [1.0,1.0,1.0]

    # contact points 
    contact_point1 = me.Point('CP1')
    contact_point1.set_pos(point_offset[0],contTS[0]*frame[1].x+contTS[1]*frame[1].y  +contTS[2]*frame[1].z)
    # contact_point1.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point1.v2pt_theory(point_offset[0],frame_ground,frame[1])

    contact_point2 = me.Point('CP2')
    contact_point2.set_pos(point_offset[0],contAI[0]*frame[1].x+contAI[1]*frame[1].y  +contAI[2]*frame[1].z)
    # contact_point2.set_vel(scapula.frame,0) # point is fixed in scapula
    contact_point2.v2pt_theory(point_offset[0],frame_ground,frame[1])

    ## contact forces

    # Distances between contact points and thorax frame
    x_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.x)
    y_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.y)
    z_pos1 = contact_point1.pos_from(point_ground).dot(frame_ground.z)
    x_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.x)
    y_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.y)
    z_pos2 = contact_point2.pos_from(point_ground).dot(frame_ground.z)

    # Contact forces
    elips_dim_scaled = []
    elips_dim_scaled.append(elips_dim[0]*first_elips_scale[0])
    elips_dim_scaled.append(elips_dim[1]*first_elips_scale[1])
    elips_dim_scaled.append(elips_dim[2]*first_elips_scale[2])
    # print(elips_dim_scaled)
    f1_in = ((x_pos1-elips_trans[0])/(elips_dim_scaled[0]))**2+((y_pos1-elips_trans[1])/(elips_dim_scaled[1]))**2+((z_pos1-elips_trans[2])/(elips_dim_scaled[2]))**2-1
    # f1_out = ((x_pos1-elips_trans[0])/(second_elips_scale*elips_dim[0]))**2+((y_pos1-elips_trans[1])/(second_elips_scale*elips_dim[1]))**2+((z_pos1-elips_trans[2])/(second_elips_scale*elips_dim[2]))**2-1
    F1_in = 1/2*(f1_in-sp.sqrt(f1_in**2+eps_in**2))
    # F1_out = 1/2*(f1_out+sp.sqrt(f1_out**2+eps_out**2))
    Fx1 = -(k_contact_in*F1_in)*(x_pos1-elips_trans[0])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[0]**2)
    Fy1 = -(k_contact_in*F1_in)*(y_pos1-elips_trans[1])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[1]**2)
    Fz1 = -(k_contact_in*F1_in)*(z_pos1-elips_trans[2])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[2]**2)

    f2_in = ((x_pos2-elips_trans[0])/(elips_dim_scaled[0]))**2+((y_pos2-elips_trans[1])/(elips_dim_scaled[1]))**2+((z_pos2-elips_trans[2])/(elips_dim_scaled[2]))**2-1
    # f2_out = ((x_pos2-elips_trans[0])/(second_elips_scale*elips_dim[0]))**2+((y_pos2-elips_trans[1])/(second_elips_scale*elips_dim[1]))**2+((z_pos2-elips_trans[2])/(second_elips_scale*elips_dim[2]))**2-1
    F2_in = 1/2*(f2_in-sp.sqrt(f2_in**2+eps_in**2))
    # F2_out = 1/2*(f2_out+sp.sqrt(f2_out**2+eps_out**2))
    Fx2 = -(k_contact_in*F2_in)*(x_pos2-elips_trans[0])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[0]**2)
    Fy2 = -(k_contact_in*F2_in)*(y_pos2-elips_trans[1])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[1]**2)
    Fz2 = -(k_contact_in*F2_in)*(z_pos2-elips_trans[2])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[2]**2)

    # applying contact forces to contact points in thorax frame
    cont_force1 = [(contact_point1,frame_ground.x*Fx1+frame_ground.y*Fy1+frame_ground.z*Fz1)]
    cont_force2 = [(contact_point2,frame_ground.x*Fx2+frame_ground.y*Fy2+frame_ground.z*Fz2)]
    CONT = cont_force1+cont_force2

    # print(kinematical)
    KM = me.KanesMethod(frame_ground, q_ind=q, u_ind=w, kd_eqs = kinematical, u_auxiliary=uaux)

    (fr, frstar) = KM.kanes_equations(BODY, (FG+DAMP+CONT+f_aux+spring))
    MM = KM.mass_matrix_full
    FO = KM.forcing_full
    xdot = (KM.q.col_join(KM.u)).diff()
    
    return MM,FO,q,w,faux,fr,frstar,kinematical,xdot,first_elips_scale,elips_trans