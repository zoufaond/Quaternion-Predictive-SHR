import sympy as sp
import numpy as np
import scipy as sc
import sympy.physics.mechanics as me
from scipy.spatial.transform import Rotation as spat
import pickle

def objective_state_diff(num_nodes,interval_value):
    
    """Build objective and Jacobian terms that penalize state rate of change squared."""

    node_val = sp.Matrix(sp.symbols('state_0:' + str(num_nodes)))
    obj_list = []

    for i in range(num_nodes-1):
        state_diff = node_val[i+1] - node_val[i]
        obj_list.append((state_diff**2) / interval_value)

    obj = sum(obj_list)
    obj_jac = sp.Matrix([obj]).jacobian(node_val)

    # Lambda functions for numerical evaluation
    obj_np = sp.lambdify([node_val],obj)
    obj_jac_np = sp.lambdify([node_val],obj_jac)
    
    return obj_np,obj_jac_np



def objective_max_GH_stab(tilt_y,tilt_z,interval_value):

    """Create objective and Jacobian functions for GH joint stability.

    The reaction force vector is rotated into the GH frame using the provided
    version and inclination angles (degrees). The objective penalizes shear-force components
    relative to the compressive component.
    """

    reactions = sp.symbols('rx ry rz')
    # Rotate reaction forces into GH coordinates using y- and z-tilt angles.
    reactions_rotated = R_y(tilt_y*sp.pi/180).T * R_z(tilt_z*sp.pi/180).T * sp.Matrix([*reactions,1])

    # Penalize normalized non-compressive GH reaction components.
    # 0.31 scales the y-shear relative to z-shear based on GH fossa geometry (Lippitt et al. 1993)
    obj = interval_value*(0.31 * (reactions_rotated[1]/reactions_rotated[0])**2 + (reactions_rotated[2]/reactions_rotated[0])**2)
    # Symbolic gradient of the objective with respect to [rx, ry, rz].
    obj_jac = sp.Matrix([obj]).jacobian(reactions)[:]

    # Convert symbolic expressions to callable numerical functions.
    obj_np = sp.lambdify([reactions],obj)
    obj_jac_np = sp.lambdify([reactions],obj_jac)
    return obj_np, obj_jac_np



def objective_scalers_restriction(mus_group,w_fmax,w_lceopt):

    """Create objective and Jacobian functions for muscle scalers restriction.
    
    Extracts muscle group names from mus_group and creates symbolic variables for
    fmax and lceopt scalers. The objective penalizes deviations of each scaler from 1
    using weighted squared-error terms (1-scaler)^2.
    
    Parameters:
    -----------
        mus_group: List of muscle groups; each group name is extracted from index [0][0:-1].
        w_fmax: Weight for fmax scaler objective term.
        w_lceopt: Weight for lceopt scaler objective term.
    
    Returns:
        obj_np: Lambdified objective function callable with sorted scaler parameters.
        obj_jac_np: Lambdified Jacobian function callable with sorted scaler parameters.
    """
    
    mus_dict = {}
    obj = 0

    # Create fmax_scaler symbolic variables for each muscle group and add penalty term.
    for igroup in range(len(mus_group)):
        current_group = mus_group[igroup][0][0:-1]
        fmax_scaler = sp.Symbol('fmax_scaler_' + current_group)
        mus_dict.update({'fmax_scaler_' + current_group: fmax_scaler})
        obj += w_fmax * ((1-fmax_scaler)**2)

    # Create lceopt_scaler symbolic variables for each muscle group and add penalty term.
    for igroup in range(len(mus_group)):
        current_group = mus_group[igroup][0][0:-1]
        lceopt_scaler = sp.Symbol('lceopt_scaler_' + current_group)
        mus_dict.update({'lceopt_scaler_' + current_group: lceopt_scaler})
        obj += w_lceopt * (1-lceopt_scaler)**2

    # Sort parameters alphabetically and extract as list for consistent ordering.
    myKeys = list(mus_dict.keys())
    myKeys.sort()
    sorted = {i: mus_dict[i] for i in myKeys}
    sorted_params = [*sorted.values()]

    obj_jac = sp.Matrix([obj]).jacobian(sorted_params)
    obj_np = sp.lambdify(sorted_params,obj)
    obj_jac_np = sp.lambdify(sorted_params,obj_jac)

    return obj_np, obj_jac_np


def objective_min_activation(acts,interval_value,include_EMG = False, muscles_emg_tracked = []):
    """
    Objective function for activation minimization with optional EMG-based constraints.
    
    Parameters:
    -----------
        acts: List of activation symbols.
        interval_value: Interval value.
        include_EMG: Boolean flag to enable EMG-based minimization (True) or activation-only 
                     minimization (False). When True, muscle group activations are minimized 
                     against corresponding EMG signals. When False, only activation minimization 
                     is performed.
        muscles_emg_tracked: List of muscle elements matching EMG signal. Each group is a list of corresponding activation symbols.
    
    Returns:
        obj_np: Lambdified objective function.
        obj_jac_np: Lambdified Jacobian function.
    """
    
    a_str = sorted([str(acts[x]).replace('(t)','') for x in range(len(acts))])
    a_emg = sp.symbols(f'a_emg1:{len(muscles_emg_tracked)+1}')
    weight_index = sp.Symbol('weight_index')
    w_min_squared = sp.Symbol('w_min_squared')
    w_min_emg = sp.Symbol('w_min_emg')
    obj_emg = 0
    a = sp.symbols(a_str)
    for j,imus in enumerate(muscles_emg_tracked):
        mus_indeces = [a_str.index(mus) for mus in imus]
        mus_current = 0
        for i in range(len(mus_indeces)):
            mus_current += a[mus_indeces[i]]
        obj_emg += (mus_current/len(mus_indeces)-a_emg[j])**2
    
    if include_EMG:
        obj_a = sum(sp.Matrix(a).applyfunc(lambda x:x**2))
        
        obj = sp.Matrix([interval_value * (obj_a * w_min_squared + weight_index * obj_emg * w_min_emg)])
        obj_jac = (obj).jacobian(a)[:]
        obj_np = sp.lambdify((a,a_emg, weight_index, w_min_squared, w_min_emg),obj)
        obj_jac_np = sp.lambdify((a,a_emg, weight_index, w_min_squared, w_min_emg),obj_jac)
    else:
        obj_a = sum(sp.Matrix(a).applyfunc(lambda x:x**2))
        obj = sp.Matrix([interval_value * obj_a])
        obj_jac = (obj).jacobian(a)[:]
        obj_np = sp.lambdify((a,),obj)
        obj_jac_np = sp.lambdify((a,),obj_jac)

    return obj_np, obj_jac_np



def objective_traj_quat(num_coords,interval_value, clav_pos,include_quat_norm = False):
    """
    Constructs a quaternion-based objective function for trajectory tracking and its Jacobian.
    
    Parameters:
    -----------
        num_coords: Number of coordinate variables (13 for full shoulder model)
        interval_value: Scaling factor for objective function terms
        clav_pos: Position along clavicle x-axis for AC joint position calculation
        include_quat_norm: If True, adds quaternion norm constraints to keep quaternions valid
    
    Returns:
        obj_np: Lambdified objective function
        obj_jac_np: Lambdified Jacobian of objective with respect to state variables
        obj_SC_np: Lambdified SC joint rotation objective function
        obj_jac_SC_np: Lambdified Jacobian of SC objective with respect to SC rotation
    """
    
    # Define symbolic variables for the optimization problem (state variables)
    # x[0:4] - SC joint rotation (quaternion), x[4:8] - AC joint rotation (quaternion)
    # x[8:12] - GH joint rotation (quaternion), x[12] - elbow flexion/extension
    x = sp.Matrix(sp.symbols(f'x1:{num_coords+1}'))
    
    # Define symbolic variables for the trajectory reference points
    x_traj = sp.Matrix(sp.symbols(f'x_traj1:{num_coords+1}'))

    # Index parameters to selectively enable/disable objective terms during optimization
    index_clav_scap = sp.Symbol('index_clav_scap')  # Enables SC and scapula-thorax objectives
    index_hum = sp.Symbol('index_hum')              # Enables humerus-thorax objective

    # Extract and compute rotations for each joint segment
    SCrot = sp.Matrix(x[0:4])                                              # SC joint rotation (quaternion)
    scapula_thorax = sp.Matrix([mulQuat_sp(x[0:4],x[4:8])])              # Scapula orientation in thorax frame
    humerus_thorax = sp.Matrix([mulQuat_sp(scapula_thorax,x[8:12])])     # Humerus orientation in thorax frame
    elrot = sp.Matrix(x[12:13])                                            # Elbow rotation (flexion/extension angle)

    # Calculate AC joint position in thorax frame (point on clavicle x-axis, axial rotation free)
    AC_pos = Qrm(x[0:4]) * position([clav_pos,0,0])
    
    # Objective function components - minimize squared differences between current and reference trajectory
    # AC position objective: tracks AC joint position (translation component)
    obj_SC_wo_x = sp.Matrix([interval_value*sum((sp.Matrix(AC_pos[0:3])-sp.Matrix(x_traj[0:3])).applyfunc(lambda x: x**2))])
    
    # Scapula-thorax orientation objective: tracks scapula orientation using vector part of quaternion error
    obj_scapula_thorax = sp.Matrix([interval_value*sum(((mulQuat_sp(scapula_thorax,invQuat_sp(x_traj[4:8])))).applyfunc(lambda x: x**2)[1:4])])
    
    # Humerus-thorax orientation objective: tracks humerus orientation using vector part of quaternion error
    obj_humerus_thorax = sp.Matrix([interval_value*sum(((mulQuat_sp(humerus_thorax,invQuat_sp(x_traj[8:12])))).applyfunc(lambda x: x**2)[1:4])])
    
    # Elbow rotation objective: tracks elbow flexion/extension angle
    obj_elrot = sp.Matrix([interval_value*sum((elrot-sp.Matrix(x_traj[12:13])).applyfunc(lambda x: x**2))])
    
    # SC joint full rotation objective: used at initial node to minimize conoid ligament lengthening
    # Minimizes the vector part of quaternion error (rotation angle without axial component)
    obj_SC = sp.Matrix([interval_value*sum(((mulQuat_sp(SCrot,invQuat_sp(x_traj[0:4])))).applyfunc(lambda x: x**2)[1:4])])

    # Combine objective terms with selective weighting indices
    obj = index_clav_scap*(obj_SC_wo_x + obj_scapula_thorax) + index_hum*(obj_humerus_thorax) + obj_elrot

    # Optional constraint: minimize quaternion norm deviation from 1 to maintain quaternion validity
    # This term penalizes deviation from unit norm for all three rotation quaternions
    if include_quat_norm:
        SC_norm = sp.Matrix([x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2-1])**2
        AC_norm = sp.Matrix([x[4]**2 + x[5]**2 + x[6]**2 + x[7]**2-1])**2
        GH_norm = sp.Matrix([x[8]**2 + x[9]**2 + x[10]**2 + x[11]**2-1])**2
        obj += (SC_norm + AC_norm + GH_norm)*interval_value
    
    # Lambda functions for numerical evaluation of the objective and its Jacobian with respect to the optimization variables.
    obj_jac = (obj).jacobian(x)[:]
    obj_np = sp.lambdify((x,x_traj,index_clav_scap,index_hum),obj)
    obj_jac_np = sp.lambdify((x,x_traj,index_clav_scap,index_hum),obj_jac)

    # Lambda functions for numerical evaluation of the full SC rotation objective and its Jacobian with respect to the SC rotation variables.
    obj_SC_jac = (obj_SC).jacobian(x[0:4])[:]
    obj_SC_np = sp.lambdify((x[0:4],x_traj[0:4]),obj_SC)
    obj_jac_SC_np = sp.lambdify((x[0:4],x_traj[0:4]),obj_SC_jac)

    return obj_np,obj_jac_np, obj_SC_np, obj_jac_SC_np



def objective_traj_eul(num_coords,interval_value,GH_seq = 'YZY'):

    """"Building the objective for Euler-angle-based model follows the same logic as for quaternion-based model"""

    x = sp.symbols(f'x1:{num_coords+1}')
    x_traj = sp.symbols(f'x_traj1:{num_coords+1}')

    SCrot = sp.Matrix(x[0:2])
    SCrot_t0 = sp.Matrix(x[0:3])
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

    obj_SCrot = sp.Matrix([interval_value*sum((SCrot-sp.Matrix(x_traj[0:2])).applyfunc(lambda x: x**2))]) # We exclude axial rotation of clavicle from the objective function
    obj_SC_t0 = sp.Matrix([interval_value*sum((SCrot_t0-sp.Matrix(x_traj[0:3])).applyfunc(lambda x: x**2))]) # This is used at the inital node to prevent excessive conoid ligmament lengthening.

    obj_scapula_thorax = sp.Matrix([interval_value*sum((scapula_thorax-sp.Matrix(x_traj[3:6])).applyfunc(lambda x: x**2))])
    obj_humerus_thorax = sp.Matrix([interval_value*sum((humerus_thorax-sp.Matrix(x_traj[6:9])).applyfunc(lambda x: x**2))])
    obj_elrot = sp.Matrix([interval_value*sum((elrot-sp.Matrix(x_traj[9:10])).applyfunc(lambda x: x**2))])

    obj = obj_SCrot + obj_scapula_thorax + obj_humerus_thorax + obj_elrot
    obj_jac = (obj).jacobian(x)[:]
    obj_np = sp.lambdify((x,x_traj),obj)
    obj_jac_np = sp.lambdify((x,x_traj),obj_jac)
    obj_SC_t0_jac = (obj_SC_t0).jacobian(x[0:3])[:]
    obj_SC_np = sp.lambdify((x[0:3],x_traj[0:3]),obj_SC_t0)
    obj_jac_SC_np = sp.lambdify((x[0:3],x_traj[0:3]),obj_SC_t0_jac)

    return obj_np,obj_jac_np, obj_SC_np, obj_jac_SC_np


def polynomials_quat(model_struct,q,u,derive,RC_lim, calibrated_params = None):
    """
    Calculates muscle-tendon forces, moment arms, and torques for a shoulder musculoskeletal model using quaternion-based kinematics and polynomial approximations for wrapped muscles.

    Parameters
    ----------
    model_struct : Matlab Structure containing shoulder model data from readosim.m function.
    q : List of generalized coordinates representing joint positions (excluding thorax quaternion and radius).
    u : List of generalized velocities.
    derive : Derivation method, either 'symbolic' or 'numeric'. 'Numeric' uses the corresponding values instead of symbols and is quicker.
    RC_lim : Rotator cuff muscle force limiter, scaling factor applied to rotator cuff muscles (supra, infra). 1.0 means no scaling, values <1 reduce force, values >1 increase force.
    calibrated_params : Parameter calibration mode. If None, uses default muscle parameters. If 'calibrate_params', enables parameter optimization with scalers. If dict, applies pre-calibrated parameters from previous optimization. Default is None.
    
    Returns
    -------
    TE_subs : Generalized torques/forces (10x1) excluding pronation/supination.
    act :  List of muscle activation symbols as dynamic symbols.
    TE_conoid : Torque contributions from conoid ligament (10x1) (to be consistent with EoMs shape, but it affects only the AC joint rotation).
    fmax_init :  Initial values for maximum muscle force scalers.
    fmax_range : Parameter ranges for maximum muscle force scalers.
    lceopt_init : Initial values for optimal muscle length scalers.
    lceopt_range : Parameter ranges for optimal muscle length scalers.
    mus_groups : Grouped muscle elements that are scaled by the same scaler.
    GH_mus_forces_subbed : Glenohumeral joint reaction forces (3x1: x, y, z components).
    """

    # polynomials were build with coordinates for thorax tilting and for pronation and supination, so we need to artifically add these coordinates
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

    qpol = q_new[1:4]+q_new[5:8]+q_new[9:12]+q_new[13:16]+q_new[16:18] #exclude q0 from quaternions, polynomials are build from vector part only (Zoufaly et al., 2025)
    nmus = len(model_struct['model']['muscles'][0,0][0])

    # calculate derivative of quaternions using angular velocities and quaternions.
    dquatdSC = dquatdt(q_new[4:8],u[0:3])
    dquatdAC = dquatdt(q_new[8:12],u[3:6])
    dquatdGH = dquatdt(q_new[12:16],u[6:9])
    dqdtEL = u[9]
    
    # symbolically define activations, these will be later substituted by dynamicsymbols for optimization.
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
        # these muscle groups will be calibrated, deltscap group will be calibrated by fmax_scaler and lceopt_scaler etc.
        calibrated_muscles = ['deltscap','deltclav','trapscap','trapclav','serr','infra']
        mus_groups = [[] for _ in range(len(calibrated_muscles))]
        muscle_constants = {}
        fmax_scaler_one = sp.Symbol('fmax_scaler')
        fmax_scaler_group = [sp.Symbol('fmax_scaler_' + i) for i in calibrated_muscles]
        lceopt_scaler_group = [sp.Symbol('lceopt_scaler_' + i) for i in calibrated_muscles]
        
        for i in range(nmus):
            muscle = model_struct['model']['muscles'].item()[0,i]
            
            lslack.append(muscle['lslack'][0,0].item())
            vmax.append(muscle['vmax'][0,0].item())
            mus_name = muscle['name'][0,0].item()

            if calibrated_params == None:
                fmax.append(muscle['fmax'][0,0].item())
                lceopt.append(muscle['lceopt'][0,0].item())

            elif calibrated_params == 'calibrate_params':
                if 'deltscap' in mus_name or 'deltclav' in mus_name or 'trapscap' in mus_name or 'trapclav' in mus_name or 'serr' in mus_name or 'infra' in mus_name:
                    mus_index = next((i for i, item in enumerate(calibrated_muscles) if mus_name.startswith(item)), None)
                    mus_groups[mus_index].append(mus_name)
                    fmax_current = muscle['fmax'][0,0].item()
                    fmax_scaler = fmax_scaler_group[mus_index]
                    fmax.append(fmax_scaler*fmax_current)
                    fmax_init.update({'fmax_scaler_' + calibrated_muscles[mus_index]: 1})
                    fmax_range.update({fmax_scaler: (0.5, 1.5)})

                    lceopt_scaler = lceopt_scaler_group[mus_index]
                    lceopt.append(lceopt_scaler*muscle['lceopt'][0,0].item())
                    lceopt_init.update({'lceopt_scaler_' + calibrated_muscles[mus_index]: 1})
                    lceopt_range.update({lceopt_scaler: (0.7, 1.3)})
                else:
                    fmax.append(muscle['fmax'][0,0].item())
                    lceopt.append(muscle['lceopt'][0,0].item())

            else:
                try:
                    mus_index = next((i for i, item in enumerate(calibrated_muscles) if mus_name.startswith(item)), None)
                    mus_groups[mus_index].append(mus_name)
                    fmax.append(calibrated_params['calibrated_params']['fmax_scaler_' + calibrated_muscles[mus_index]].item()[0][0] * muscle['fmax'][0,0].item())
                    lceopt.append(calibrated_params['calibrated_params']['lceopt_scaler_' + calibrated_muscles[mus_index]].item()[0][0] * muscle['lceopt'][0,0].item())
                except:
                    fmax.append(muscle['fmax'][0,0].item())
                    lceopt.append(muscle['lceopt'][0,0].item())

            if 'supra' in mus_name or 'infra' in mus_name:
                fmax[i] = fmax[i] * RC_lim

    # preallocate muscle length, force, and Jacobian arrays
    mus_lengths_full = sp.zeros(nmus,1)
    mus_forces = sp.zeros(nmus,1)
    JacInSpat_full = sp.zeros(11,nmus)
    # preallocate GH force vector components
    xforce = 0
    yforce = 0
    zforce = 0
    
    for imus in range(nmus):
        muscle = model_struct['model']['muscles'].item()[0,imus]
        isWrapped = muscle['isWrapped'].item()
        mus_name = muscle['name'][0,0].item()

        if isWrapped == 0:
            # if not wrapped, calculate length and moment arms analytically
            # This is based on section 3.2 in Zoufaly et al., 2025
            origin = muscle['origin_frame'].item()
            insertion = muscle['insertion_frame'].item()
            O_pos = muscle['origin_position'].item()[0]
            I_pos = muscle['insertion_position'].item()[0]
            # Calculate muscle length using analytical function for quaternion representation
            L = analytic_length_quat(origin, insertion, O_pos, I_pos, q_new[4:], model_struct)
            # Calculate Jacobian of muscle length and velocity of muscle length)
            Jac = -sp.Matrix([L]).jacobian(q_new).T
            vce = -sp.Matrix(Jac[4:-1]).T * sp.Matrix([dquatdSC,dquatdAC,dquatdGH,dqdtEL])

            # Map the jacobian from the quaternion representation to spatial coordinates using the G matrix for each joint segment
            # Elbow is a revolute joint, so no mapping is needed         
            TEsc = 1/2 * G(q_new[4:8])*sp.Matrix(Jac[4:8])
            TEac = 1/2 * G(q_new[8:12])*sp.Matrix(Jac[8:12])
            TEgh = 1/2 * G(q_new[12:16])*sp.Matrix(Jac[12:16])
            TEel =  sp.Matrix(Jac[16:18])
            
            iJacInSpat = sp.Matrix(TEsc).col_join(TEac).col_join(TEgh).col_join(TEel)
            JacInSpat_full[:,imus] = iJacInSpat
            
        elif isWrapped == 1:
            
            # If muscle is wrapped, calculate the polynomial approximaiton of the muscle-tendon path
            # The mapping used here is based on section 3.3 in Zoufaly et al., 2025.
            name = muscle['name'].item()
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
                    mdof = musdof[j]
                    for k in range(expon_np[i,j].item()):
                        term = term * qpol[int(mdof-1)]

                L = L + term
            jac = -sp.Matrix([L]).jacobian(qpol).T
            vce = -sp.Matrix(jac[3:-1]).T * sp.Matrix([sp.Matrix(dquatdSC[1:4]),sp.Matrix(dquatdAC[1:4]),sp.Matrix(dquatdGH[1:4]),dqdtEL])

            # inEtrans introduced in eq. 20 in Zoufaly et al., 2025
            TEsc = invEtrans(q_new[4:8])*sp.Matrix(jac[3:6])
            TEac = invEtrans(q_new[8:12])*sp.Matrix(jac[6:9])
            TEgh = invEtrans(q_new[12:16])*sp.Matrix(jac[9:12])
            TEel = sp.Matrix(jac[12:14])
            iJacInSpat = sp.Matrix(TEsc).col_join(TEac).col_join(TEgh).col_join(TEel)
            JacInSpat_full[:,imus] = iJacInSpat
            # print('poly')
            
        mus_lengths_full[imus] = L
        mus_forces[imus] = muscle_force(actSym[imus],L,vce,fmax[imus],lceopt[imus],lslack[imus],vmax[imus])

        # Not all OS_model.m files include polynomial coefficients for GH force vectors (only OS_model_prediction for par2 does).
        try:
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
                        term = term * qpol[int(mdof-1)]

                zvec = zvec + term
            
            xforce += xvec * mus_forces[imus]
            yforce += yvec * mus_forces[imus]
            zforce += zvec * mus_forces[imus]

        except:
            pass
            
    # Calculate full torque from muscle forces
    FQ = JacInSpat_full * mus_forces
    # Exclude the pronation and supination row and set pronation to 120 degrees.
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

    # Conoid ligament moment arms are calculated the same way as non-wrapped muscles
    conoid_lopt = model_struct['model']['conoid_length'].item()[0][0]
    conoid_k = model_struct['model']['conoid_stiffness'].item()[0][0]
    conoid_eps = model_struct['model']['conoid_eps'].item()[0][0]
    conoid_origin = model_struct['model']['conoid_origin'].item()[0]
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
    
    return TE_subs, act, TE_conoid[:-1], fmax_init, fmax_range, lceopt_init, lceopt_range, mus_groups, GH_mus_forces_subbed


def polynomials_euler(model_struct,q,u,derive, motion_folder = None, gen_matlab_functions = None):

    """ Polynomials for Euler angle version of the model is applied with the same logic as the quaternion version. """

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
    
    conoid_lopt = model_struct['model']['conoid_length'].item()[0][0]
    conoid_k = model_struct['model']['conoid_stiffness'].item()[0][0]
    conoid_eps = model_struct['model']['conoid_eps'].item()[0][0]
    conoid_origin = model_struct['model']['conoid_origin'].item()[0]
    conoid_insertion = model_struct['model']['conoid_insertion'].item()[0]
    conoid_length = analytic_length_eul('clavicle_r','scapula_r', conoid_origin, conoid_insertion, qpol[3:], model_struct)
    F_conoid = conoid_force(conoid_length, conoid_lopt, conoid_k, conoid_eps)
    jac_conoid = -sp.Matrix([conoid_length]).jacobian(qpol[3:]).T
    TE_conoid = F_conoid * jac_conoid
    
    return TE_act_subbed, act, TE_conoid[:-1]

def muscle_force(act, lmt, vce, fmax, lceopt, lslack, vmax):
    """Compute muscle force using a Hill-type muscle model.

        force = (flce * act * fvce + fpe) * fmax

    Parameters
    ----------
    act : Activation symbol.
    lmt : Muscle-tendon length - symbolic expression from polynomial approximation.
    vce : Contractile element velocity.
    fmax : Maximum isometric muscle force.
    lceopt : Optimal contractile element (muscle fiber) length.
    lslack : Tendon slack length.
    vmax : Maximum normalized shortening velocity coefficient.

    Returns
    -------
    sympy expression
        Muscle force produced by the model.
    """

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

    force = (flce * act * fvce +  fpe) * fmax
    
    return force

def act_dynamics(a,u,t_act = 0.015,t_deact = 0.05):
    """First-order muscle activation dynamics (Thelen 2003).

    Parameters
    ----------
    a : Current muscle activation state.
    u : Neural excitation input.
    t_act : Activation time constant in seconds (default 0.015).
    t_deact : Deactivation time constant in seconds (default 0.05).

    Returns
    -------
    Time derivative of activation da/dt as a symbolic expression.
    """

    da = (u/t_act + (1-u)/t_deact)*(u-a)

    return da

def conoid_force(lmt, lopt, k, eps):
    """Conoid ligament force (Chadwick, 2014) with smooth tension-only behavior.

    Parameters
    ----------
    lmt : Symbolic ligament length.
    lopt : Slack/reference ligament length.
    k : Ligament stiffness.
    eps : Smoothing constant at the zero-force transition.

    Returns
    -------
        Ligament tensile force expression.
    """

    d = lmt-lopt
    F = k / 2 * (d + sp.sqrt(d**2 + eps))
    
    return F

def dquatdt(quat,w):
    """Compute quaternion derivative from quaternion and angular velocity.

    Parameters
    ----------
        quat : Quaternion components as [q0, q1, q2, q3].
        w : Angular velocity components as [wx, wy, wz].

    Returns:
        Quaternion time derivative dq/dt.
    """
    res = 1/2 * G(quat).T * sp.Matrix([w]).T
    return res

def invQuat_sp(q):
    """Return the quaternion inverse (conjugate for unit quaternions)."""

    return sp.Matrix([q[0],-q[1],-q[2],-q[3]])


def invEtrans(quat):
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

def analytic_length_quat(origin, insertion, O_pos, I_pos, q, model):
    
    """Symbolic length of a non-wrapped muscle using quaternion rotations.

    Parameters
    ----------
    origin : Origin body name in the OpenSim model.
    insertion : Insertion body name in the OpenSim model.
    O_pos : Origin value point position.
    I_pos : Insertion value point position.
    q : Quaternion coordinates used for the body rotations.
    model : OpenSim model structure.

    Returns
    -------
    Muscle-tendon length.
    """

    jnts = model['model']['joints'].item()
    offset_thorax = jnts[0,1]['location'].item()[0]
    offset_clavicle = jnts[0,4]['location'].item()[0]

    if origin == 'thorax' and insertion == 'clavicle_r':
        O = position(O_pos)
        I = T_trans(offset_thorax) * Qrm(q[0:4]) * position(I_pos)
        
    elif origin == 'thorax' and insertion == 'scapula_r':
        O = position(O_pos)
        RW_C = T_trans(offset_thorax) * Qrm(q[0:4])
        TC_S = T_trans(offset_clavicle)
        RC_S = Qrm(q[4:8])
        I = RW_C * TC_S * RC_S * position(I_pos)
        
    elif origin == 'clavicle_r' and insertion == 'scapula_r':
        O = position(O_pos)
        TC_S = T_trans(offset_clavicle)
        RC_S = Qrm(q[4:8])
        I = TC_S * RC_S * position(I_pos)

    muscle_length = sp.sqrt((O[0] - I[0])**2 + (O[1] - I[1])**2 + (O[2] - I[2])**2)
    
    return muscle_length

def analytic_length_eul(origin, insertion, O_pos, I_pos, q, model):
    """This follows the same logic as analytic_length_quat but uses Euler angles."""

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

def Qrm(q):
    """Quaternion-based rotation matrix.

    Parameters
    ----------
    q : list or sympy.Matrix
        Quaternion components [q0, q1, q2, q3].

    Returns
    -------
    sympy.Matrix
        Rotation matrix.
    """
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]
    res =  sp.Matrix([[1-2*(y**2+z**2), 2*(x*y-z*w), 2*(x*z+y*w),0],
                      [2*(x*y+z*w), 1-2*(x**2+z**2), 2*(y*z-x*w),0],
                      [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x**2+y**2),0],
                      [0           ,0          ,0               ,1]])

    return res

def joint_spring_quat(q,Qeq,kpe = 2.0, epsm0 = 0.6):
    joint_stiffness = 1.0
    Qdif = mulQuat_sp(invQuat_sp(Qeq),q)

    angle = 2*sp.atan2(sp.sqrt(Qdif[1]**2 + Qdif[2]**2 + Qdif[3]**2),Qdif[0])
    scale = 1/sp.sqrt(1-Qdif[0]**2 + 1e-2)
    axis = scale * sp.Matrix([Qdif[1],Qdif[2],Qdif[3]])
    nonlin_spring = (sp.exp(kpe*angle/epsm0)-1)/(sp.exp(kpe)-1)
    res = -axis * nonlin_spring * joint_stiffness * 1

    return res





def mulQuat_sp(qa,qb):
    """Quaternion multiplication.

    Parameters
    ----------
    qa : First quaternion [q0, q1, q2, q3].
    qb : Second quaternion [q0, q1, q2, q3].

    Returns
    -------
    Quaternion product qa * qb.
    """

    res = sp.Matrix([[qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3]],
                     [qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2]],
                     [qa[0]*qb[2] - qa[1]*qb[3] + qa[2]*qb[0] + qa[3]*qb[1]],
                     [qa[0]*qb[3] + qa[1]*qb[2] - qa[2]*qb[1] + qa[3]*qb[0]]])
    return res



def create_eoms_quat_w_RF(OS_struct,hand_weight = 0, derive = 'symbolic',gen_matlab_functions = None):
    
    """Create equations of motion with reaction forces in GH joint.

    Parameters
    ----------
    OS_struct : OpenSim model structure.
    hand_weight : Weight in hand segment to include in the model (default 0).
    derive : 'symbolic' or 'numeric' to specify whether to use symbolic parameters or numeric values from the OpenSim model.

    Returns
    -------
    q : list of sympy dynamic symbols for generalized coordinates.
    w : list of sympy dynamic symbols for angular velocities.
    faux : list of sympy dynamic symbols for GH reaction forces (in scapular frame).
    fr,frstar : fr + frstar = 0 are implicitly defined Equations of motion. The last 3 elements correspond to implcitly defined reaction forces in the GH joint.
    kinematical : list of kinematic equations relating quaternion derivatives to angular velocities.
    """

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
        c = [0.5,0.5,0.5]


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
    mass[-1] += hand_weight

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
    DAMP = [(frame[0], -c[0]*(w[0]*frame[0].x+w[1]*frame[0].y+w[2]*frame[0].z))]

    # axis-angle springs have the equilibrium position defined in OpenSim, recalculated to quaternions
    SC_spring_torque = joint_spring_quat(q[0:4],sp.Matrix([0.981546017184242,-0.005342836650565,-0.186060194671936,0.043823215364912]))
    spring = [(frame[0], (SC_spring_torque[0]*frame[0].x + SC_spring_torque[1]*frame[0].y + SC_spring_torque[2]*frame[0].z))]
    # spring = []
             
    # iterate over segments 2:end (first body is already done)
    for i in range(1,3):

        # set masscenter points 
        if i == 2:
            # calculate reaction forces in GH joint in scapular frame
            uaux = me.dynamicsymbols('uaux1 uaux2 uaux3')
            faux = me.dynamicsymbols('faux1 faux2 faux3')
            
            N_GH = point_offset[i-1].locatenew('N_GH', 0)
            N_GH.set_vel(frame_ground, N_GH.vel(frame_ground) + uaux[0]*frame[i-1].x + uaux[1]*frame[i-1].y + uaux[2]*frame[i-1].z)

            masscenter[i].set_pos(N_GH,com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
            masscenter[i].v2pt_theory(N_GH,frame_ground,frame[i])

            point_offset[i].set_pos(N_GH,offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
            point_offset[i].v2pt_theory(N_GH,frame_ground,frame[i])

            f_aux = [(point_offset[i-1], 800*(faux[0]*frame[i-1].x + faux[1]*frame[i-1].y + faux[2]*frame[i-1].z)),(N_GH, -800*(faux[0]*frame[i-1].x + faux[1]*frame[i-1].y + faux[2]*frame[i-1].z))]

            GH_spring_torque = joint_spring_quat(q[8:12],sp.Matrix([1,0,0,0])) * 0.25
            spring.append((frame[i], GH_spring_torque[0]*frame[i].x + GH_spring_torque[1]*frame[i].y + GH_spring_torque[2]*frame[i].z))
            spring.append((frame[i-1], -GH_spring_torque[0]*frame[i].x - GH_spring_torque[1]*frame[i].y - GH_spring_torque[2]*frame[i].z))

            
        else:
            masscenter[i].set_pos(point_offset[i-1],com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
            masscenter[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

            # set offsent points (where the next joint is)
            point_offset[i].set_pos(point_offset[i-1],offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
            point_offset[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

            # 
            AC_spring_torque = joint_spring_quat(q[4:8],sp.Matrix([0.894427818265390,0.034323247153231,0.443000552929629,0.050708014375610]))
            spring.append((frame[i], AC_spring_torque[0]*frame[i].x + AC_spring_torque[1]*frame[i].y + AC_spring_torque[2]*frame[i].z))
            spring.append((frame[i-1], -AC_spring_torque[0]*frame[i].x - AC_spring_torque[1]*frame[i].y - AC_spring_torque[2]*frame[i].z))

            # GH_spring_torque = joint_spring_quat(q[8:12],sp.Matrix([0.934216430050072,-0.004829341953100,0.356581663437207,-0.008115206784339]))
            

        # set gravity force in masscenter (-y direction in frame_ground)
        FG.append((masscenter[i], -mass[i] * g * frame_ground.y))
        # set damping in joints (c * angular_velocity)
        damping = -c[i]*(w[0+i*3]*frame[i].x+w[1+i*3]*frame[i].y+w[2+i*3]*frame[i].z)

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

    DAMP.append(((frame[3]),-c[2]*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                            +ulna_rot_frame.z*EL_rot_axis[2])))
    DAMP.append((ulna_rot_frame,c[2]*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
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
    contact_point1.v2pt_theory(point_offset[0],frame_ground,frame[1])

    contact_point2 = me.Point('CP2')
    contact_point2.set_pos(point_offset[0],contAI[0]*frame[1].x+contAI[1]*frame[1].y  +contAI[2]*frame[1].z)
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

    f1_in = ((x_pos1-elips_trans[0])/(elips_dim_scaled[0]))**2+((y_pos1-elips_trans[1])/(elips_dim_scaled[1]))**2+((z_pos1-elips_trans[2])/(elips_dim_scaled[2]))**2-1
    F1_in = 1/2*(f1_in-sp.sqrt(f1_in**2+eps_in**2))
    Fx1 = -(k_contact_in*F1_in)*(x_pos1-elips_trans[0])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[0]**2)
    Fy1 = -(k_contact_in*F1_in)*(y_pos1-elips_trans[1])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[1]**2)
    Fz1 = -(k_contact_in*F1_in)*(z_pos1-elips_trans[2])*(elips_dim_scaled[0]**2+elips_dim_scaled[1]**2+elips_dim_scaled[2]**2)/(elips_dim_scaled[2]**2)

    f2_in = ((x_pos2-elips_trans[0])/(elips_dim_scaled[0]))**2+((y_pos2-elips_trans[1])/(elips_dim_scaled[1]))**2+((z_pos2-elips_trans[2])/(elips_dim_scaled[2]))**2-1
    F2_in = 1/2*(f2_in-sp.sqrt(f2_in**2+eps_in**2))
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
    
    return q,w,faux,fr,frstar,kinematical


def create_eoms_quat_no_RF(OS_struct,hand_weight = 0, derive = 'symbolic',gen_matlab_functions = None):
    """Create equations of motion without reaction forces in GH joint. It follows the same logic as create_eoms_quat_w_RF but without the reaction forces in the GH joint."""
    
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
        c = [0.5,0.5,0.5]


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
    mass[-1] += hand_weight

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
    DAMP = [(frame[0], -c[0]*(w[0]*frame[0].x+w[1]*frame[0].y+w[2]*frame[0].z))]
             
    # iterate over segments 2:end (first body is already done)
    for i in range(1,3):

        # set masscenter points 

        masscenter[i].set_pos(point_offset[i-1],com[0+i*3]*frame[i].x + com[1+i*3]*frame[i].y + com[2+i*3]*frame[i].z)
        masscenter[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

        # set offsent points (where the next joint is)
        point_offset[i].set_pos(point_offset[i-1],offset[0+i*3]*frame[i].x + offset[1+i*3]*frame[i].y + offset[2+i*3]*frame[i].z)
        point_offset[i].v2pt_theory(point_offset[i-1],frame_ground,frame[i])

            # GH_spring_torque = joint_spring_quat(q[8:12],sp.Matrix([0.934216430050072,-0.004829341953100,0.356581663437207,-0.008115206784339]))
            

        # set gravity force in masscenter (-y direction in frame_ground)
        FG.append((masscenter[i], -mass[i] * g * frame_ground.y))
        # set damping in joints (c * angular_velocity)
        damping = -c[i]*(w[0+i*3]*frame[i].x+w[1+i*3]*frame[i].y+w[2+i*3]*frame[i].z)

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

    DAMP.append(((frame[3]),-c[2]*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
                            +ulna_rot_frame.z*EL_rot_axis[2])))
    DAMP.append((ulna_rot_frame,c[2]*w[9]*(ulna_rot_frame.x*EL_rot_axis[0]+ulna_rot_frame.y*EL_rot_axis[1]
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
    
    return q,w,fr,frstar,kinematical

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
            mass.append(OS_struct['model']['segments'][0,0][0,segment_index[i]]['mass'][0,0].item()/1.4)


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
    kinematical = []

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
        kinematical.append(q[i].diff()-u[i])

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
    kinematical.append(u[9]-q[9].diff(t))

    point_offset[3].set_pos(point_offset[2],offset[0+3*3]*frame[3].x + offset[1+3*3]*frame[3].y + offset[2+3*3]*frame[3].z)
    point_offset[3].v2pt_theory(point_offset[2],frame_ground,frame[3]);

    # radius and PSY joint (weld joint for now)
    frame[4].orient_axis(frame[3],frame[3].x*PSY_rot_axis[0]
                         +frame[3].y*PSY_rot_axis[1]+frame[3].z*PSY_rot_axis[2],
                         0) #q[10]
    masscenter[4].set_pos(point_offset[3],com[0+4*3]*frame[4].x + com[1+4*3]*frame[4].y + com[2+4*3]*frame[4].z)
    masscenter[4].v2pt_theory(point_offset[3],frame_ground,frame[4])
    FG.append((masscenter[4], -mass[4] * g * frame_ground.y))

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
        eps_in = 0.01
        eps_out = 1e-3
        elips_trans = OS_struct['model']['thoracic_wall'][0,0].item()[0][0]
        elips_dim = OS_struct['model']['thoracic_wall'][0,0].item()[1][0]
        contTS = OS_struct['model']['TScontact'][0,0].item()[0][0]
        contAI = OS_struct['model']['AIcontact'][0,0].item()[0][0]
        first_elips_scale = [1.0,1.0,1.0]
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
    
    KM = me.KanesMethod(frame_ground, q_ind=q, u_ind=u, kd_eqs=kinematical)
    (fr, frstar) = KM.kanes_equations(BODY, (FG+DAMP+CONT))
    MM = KM.mass_matrix_full
    FO = KM.forcing_full
    xdot = (KM.q.col_join(KM.u)).diff()

    return q,u,fr,frstar,kinematical
        