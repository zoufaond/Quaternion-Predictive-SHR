import sympy as sp
import numpy as np
import scipy as sc
import sympy.physics.mechanics as me
from scipy.spatial.transform import Rotation as spat
import matplotlib.pyplot as plt
import pickle

def plot_results(solution, vals , time, num_nodes):
    solution_mat = solution[:13*num_nodes].reshape(13,num_nodes)
    sol_glob = solution_mat.copy()
    sol_glob[0:4,:] = solution_mat[0:4,:]
    for i in range(num_nodes):
        sol_glob[4:8,i] = mulQuat_np(solution_mat[0:4,i],solution_mat[4:8,i])[:,0]
        # print(mulQuat_np(solution_mat[0:4,i],solution_mat[4:8,i])[:,0])

    vals_mat = vals.reshape(13,num_nodes)

    fig, axs = plt.subplots(13)
    for j in range(13):
        axs[j].plot(time,vals_mat[j,:],marker = 'o')
        axs[j].plot(time,sol_glob[j,:])
        fig.set_figheight(10)

    # print(sol_glob-vals_mat)

def initial_guess_from_solution(solution_file,num_free):

    solution = sc.io.loadmat(solution_file)['data'][0,0]
    initial_guess = solution['solution'][0]

    # initial_guess = np.zeros(num_free)
    # trajectories = solution['trajectories'][0].transpose().flatten()
    # initial_guess[0:len(trajectories)] = trajectories
    return initial_guess

def exp_emg(emg_struct_name,num_nodes,EMG_num):
    emg_struct = sc.io.loadmat(emg_struct_name)
    data = emg_struct['data']['num_'+str(EMG_num)]
    emg_names = ('AnteriorDelt','IntermediateDelt','PosteriorDelt','Infrasp','Suprasp','MiddleTrap','UpperTrap','Serrupper')
    # print(data[0,0][emg_names[0]])
    time = np.linspace(0,1,len(data[0,0][emg_names[0]].item()[0]))
    time_new = np.linspace(0,1,num_nodes)
    emg_exp = np.zeros([len(emg_names),num_nodes])

    for i in range(len(emg_names)):
        cs = sc.interpolate.CubicSpline(time,data[0,0][emg_names[i]].item()[0])
        emg_exp[i,:] = cs(time_new)

    # indexes = (emg_exp[0,:]>1e-4)*1
    indexes = np.ones(num_nodes,dtype=int)
    # indexes[30:40] = int(1)
    # indexes[140:150] = int(1)
    # indexes[240:250] = int(1)
    # indexes = np.z

    return emg_exp, indexes

def exp_trajectory_quat(mot_struct_name,num_nodes):
    mot_struct = sc.io.loadmat(mot_struct_name)
    time = mot_struct['mot_struct']['time'][0,0][:,0]
    duration = time[-1]
    time_new = np.linspace(0,duration,num_nodes)
    
    quat_coords = mot_struct['mot_struct']['quat'][0,0]
    num_coords = np.shape(quat_coords)[1]
    trajectory = np.zeros([num_coords,num_nodes])
    interval_value = duration/(num_nodes - 1)
    # x0 = mot_struct['mot_struct']['mot_quaternion_IC'][0,0][0]
    
    for i in range(num_coords):
        cs = sc.interpolate.CubicSpline(time,quat_coords[:,i])
        trajectory[i,:] = cs(time_new)
    
    return trajectory, interval_value, time_new #, x0

def exp_trajectory_quat_myobj(trajectory, clav_pos):
    new_traj = trajectory.copy()
    num_nodes = np.shape(trajectory)[1]

    for i in range(num_nodes):
        SC_wo_x = Qrm_np(trajectory[0:4,i]) @ position(clav_pos,0,0)
        new_traj[0:3,i] = SC_wo_x[0:3]
        scapula_thorax = mulQuat_np(trajectory[0:4,i],trajectory[4:8,i]).reshape(4,)
        new_traj[4:8,i] = scapula_thorax[0:4]
        humerus_thorax = mulQuat_np(scapula_thorax,trajectory[8:12,i]).reshape(4,)
        new_traj[8:12,i] = humerus_thorax[0:4]

    return new_traj

def exp_trajectory_eul(mot_struct_name,num_nodes):
    mot_struct = sc.io.loadmat(mot_struct_name)
    time = mot_struct['mot_struct']['time'][0,0][:,0]
    duration = time[-1]
    time_new = np.linspace(0,duration,num_nodes)
    
    eul_coords = mot_struct['mot_struct']['euler'][0,0] #mot_euler_mod
    num_coords = np.shape(eul_coords)[1]
    eul_new = np.zeros([num_coords,num_nodes])
    interval_value = duration/(num_nodes - 1)
    # x0 = mot_struct['mot_struct']['mot_euler_IC'][0,0][0]
    
    for i in range(num_coords):
        cs = sc.interpolate.CubicSpline(time,eul_coords[:,i])
        eul_new[i,:] = cs(time_new)
    trajectory = eul_new
    
    return trajectory, interval_value, time_new #, x0

def exp_trajectory_eul_myobj(trajectory, GH_seq = 'YZY'):
    new_traj = trajectory.copy()
    num_nodes = np.shape(trajectory)[1]

    for i in range(num_nodes):
        scapula_thorax_RM = R_y_np(trajectory[0,i]) @ R_z_np(trajectory[1,i]) @ R_x_np(trajectory[2,i]) @ R_y_np(trajectory[3,i]) @ R_z_np(trajectory[4,i]) @ R_x_np(trajectory[5,i])
        scapula_thorax_z = np.arcsin(scapula_thorax_RM[1,0])
        scapula_thorax_x = np.arctan2(-scapula_thorax_RM[1,2],scapula_thorax_RM[1,1])
        scapula_thorax_y = np.arctan2(-scapula_thorax_RM[2,0],scapula_thorax_RM[0,0])
        new_traj[3,i] = scapula_thorax_x
        new_traj[4,i] = scapula_thorax_y
        new_traj[5,i] = scapula_thorax_z

        if GH_seq == 'YZX':
            humerus_thorax_RM = scapula_thorax_RM @ R_y_np(trajectory[6,i]) @ R_z_np(trajectory[7,i]) @ R_x_np(trajectory[8,i])
            humerus_thorax_z = np.arcsin(humerus_thorax_RM[1,0])
            humerus_thorax_x = np.arctan2(-humerus_thorax_RM[1,2],humerus_thorax_RM[1,1])
            humerus_thorax_y = np.arctan2(-humerus_thorax_RM[2,0],humerus_thorax_RM[0,0])
            new_traj[6,i] = humerus_thorax_y
            new_traj[7,i] = humerus_thorax_z
            new_traj[8,i] = humerus_thorax_x
        elif GH_seq == 'YZY':
            humerus_thorax_RM = scapula_thorax_RM @ R_y_np(trajectory[6,i]) @ R_z_np(trajectory[7,i]) @ R_y_np(trajectory[8,i])
            humerus_thorax_z = np.arccos(humerus_thorax_RM[1,1])
            humerus_thorax_yy = np.arctan2(humerus_thorax_RM[1,2],humerus_thorax_RM[1,0])
            humerus_thorax_y = np.arctan2(humerus_thorax_RM[2,1],-humerus_thorax_RM[0,1])
            new_traj[6,i] = humerus_thorax_y
            new_traj[7,i] = humerus_thorax_z
            new_traj[8,i] = humerus_thorax_yy
            
    return new_traj

def das_trajectory(data_struct,num_nodes,duration,weight, coords):
    time = np.linspace(0.0, duration, num=num_nodes)
    interval_value = duration/(num_nodes - 1)
    x0 = data_struct['params']['InitPosOptQuat'][0,0]['initCondQuat'].item()
    x0eul = data_struct['params']['InitPosOptQuat'][0,0]['initCondEul'].item()
    GH_motion_Eul = np.array([np.ones(num_nodes)*x0eul[6],
                          (-np.cos(time*np.pi)+1)*weight+x0eul[7],
                          np.ones(num_nodes)*x0eul[8]]).T
    GH_motion_R = spat.from_euler('YZY',GH_motion_Eul)
    GH_motion_Q = GH_motion_R.as_quat(scalar_first=True).T
    
    
    
    if coords == 'quaternion':
        traj = np.zeros(13*num_nodes)
        for i in range(8):
            traj[i*num_nodes:(i+1)*num_nodes] = x0[i]

            traj[8*num_nodes:(8+4)*num_nodes] = GH_motion_Q.flatten()
            traj[(8+4)*num_nodes:(8+5)*num_nodes] = x0[12]

            traj_split = np.vstack(np.split(traj,13))
            d_traj = np.concatenate((np.zeros((13,1)),np.diff(traj_split)),axis=1)/interval_value
            init_guess = np.concatenate((traj,d_traj.flatten()))
        
    elif coords == 'euler':
        traj = np.zeros(10*num_nodes)
        for i in range(6):
            traj[i*num_nodes:(i+1)*num_nodes] = x0eul[i]

            traj[6*num_nodes:(6+3)*num_nodes] = GH_motion_Eul.T.flatten()
            traj[(6+3)*num_nodes:(6+4)*num_nodes] = x0eul[9]

            traj_split = np.vstack(np.split(traj,10))
            d_traj = np.concatenate((np.zeros((10,1)),np.diff(traj_split)),axis=1)/interval_value
            init_guess = np.concatenate((traj,d_traj.flatten()))
    
    return traj, init_guess

def sol2mot_quat(solution, num_nodes, num_q, time, file_name = 'traj_opt.mot', GH_seq = 'YZY'):
    traj_quat = solution[:(num_nodes*num_q)]
    traj_splitted = np.vstack(np.split(traj_quat,num_q)).T
    joints = ('YZX','YZX',GH_seq,'rev')
    traj_eul = np.zeros([num_nodes,10])

    # print(np.shape(traj_splitted))
    
    for i,jnt in enumerate(joints):
        if jnt != 'rev':
            for j in range(num_nodes):
                traj_eul[j,i*3:(i+1)*3] = quat2eul(traj_splitted[j,i*4:(i+1)*4],jnt)
        else:
            traj_eul[:,i*3] = traj_splitted[:,i*4]
            
    traj_eul = traj_eul*180/np.pi
    
    text_file = open(f"{file_name}","w")
    with open(f'{file_name}',"w") as text_file:
        print(f'Simulation',file=text_file)
        print(f'nRows={num_nodes}',file=text_file)
        print(f'nColumns={10+5}',file = text_file)
        print(f'endheader',file = text_file)
        if GH_seq == 'YZY':
            print(f'time  TH_x TH_y TH_z SC_y  SC_z  SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_yy  EL_x  PS_y',file = text_file)
        elif GH_seq == 'YZX':
            print(f'time  TH_x TH_y TH_z SC_y  SC_z  SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_x  EL_x  PS_y',file = text_file)
        
        for i in range(num_nodes):
            for j in range(10):
                if j == 0 :
                    print(f'{time[i]}  0.000000  0.000000  0.000000  {traj_eul[i,j]}', end = '  ',file = text_file)
                elif j == 9:
                    print(f'{traj_eul[i,j]}  120.000000', file = text_file)
                else:
                    print(f'{traj_eul[i,j]}', end = '  ',file = text_file)
                        
        
    text_file.close()

    print('Saved to .mot file')

def sol2mot_eul(solution, num_nodes, num_q, time, file_name = 'traj_opt.mot',GH_seq = 'YZY'):
    traj_eul = np.vstack(np.split(solution[:(num_nodes*num_q)],num_q)).T
    traj_eul = traj_eul*180/np.pi
    
    text_file = open(f"{file_name}","w")
    with open(f'{file_name}',"w") as text_file:
        print(f'Simulation',file=text_file)
        print(f'nRows={num_nodes}',file=text_file)
        print(f'nColumns={10+5}',file = text_file)
        print(f'endheader',file = text_file)
        if GH_seq == 'YZY':
            print(f'time  TH_x TH_y TH_z SC_y  SC_z  SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_yy  EL_x  PS_y',file = text_file)
        elif GH_seq == 'YZX':
            print(f'time  TH_x TH_y TH_z SC_y  SC_z  SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_x  EL_x  PS_y',file = text_file)
        
        for i in range(num_nodes):
            for j in range(10):
                if j == 0:
                    print(f'{time[i]}  0.000000  0.000000  0.000000  {traj_eul[i,j]}', end = '  ',file = text_file)
                elif j == 9:
                    print(f'{traj_eul[i,j]}  120.000000', file = text_file)
                else:
                    print(f'{traj_eul[i,j]}', end = '  ',file = text_file)
                        
        
    text_file.close()

    print('Saved to .mot file')

def sol2struct(solution,activations_list,num_q,num_u,num_inputs,num_nodes,time,objective_value,time2sol,file_name,act_dyn = False,torqueDriven = False):
    
    trajectories = np.zeros([num_nodes, num_q])
    speeds = np.zeros([num_nodes, num_u])
    reactions = np.zeros([num_nodes, 3])
    for i in range(num_q):
        trajectories[:,i] = solution[(i)*num_nodes:(i+1)*num_nodes]

    for i in range(num_u):
        speeds[:,i] = solution[(num_q + i)*num_nodes:(num_q + i +1)*num_nodes]

    for i in range(3):
        reactions[:,i] = solution[(num_q + num_u + i)*num_nodes:(num_q + num_u + i +1)*num_nodes]*800

    if torqueDriven is not True:
        activations, excitations = input2mat(solution, num_nodes,num_q,num_u, num_inputs, activations_list, act_dyn)
    else:
        activations = np.zeros([num_nodes, num_inputs])
        excitations = np.zeros([num_nodes, num_inputs])
        for i in range(num_inputs):
            activations[:,i] = solution[(num_q + num_u + i)*num_nodes:(num_q + num_u + i +1)*num_nodes]

    data = {
                'tout': time,
                'trajectories': trajectories,
                'activations': activations,
                'excitations': excitations,
                'speeds': speeds,
                'time2sol': time2sol,
                'objective_value': objective_value,
                'solution': solution,
                'reactions': reactions,
                    }
    sc.io.savemat(f'{file_name}', {'data': data})
    print('Saved to .mat file')

def input2mat(solution, num_nodes, num_q, num_u, num_inputs, activations, act_dyn):
    activations_str = [str(activations[x]).replace('(t)','') for x in range(len(activations))]
    act_index = np.linspace(0,137,138,dtype=int)

    if act_dyn:
        dict_act_index = dict(zip((activations_str), act_index))
    else:
        dict_act_index = dict(zip(sorted(activations_str), act_index))

    dict_exc_index = dict(zip(sorted(activations_str), act_index))

    data_act = np.zeros([num_nodes, 137])
    data_exc = np.zeros([num_nodes, 137])
    for i in range(137):
        try:
            ind = dict_act_index.get(f'act_{i+1}')
            data_act[:,i] = solution[(num_q + num_u + 3 +ind)*num_nodes:(num_q + num_u + 3 + ind +1)*num_nodes]
        except:
            continue

    for i in range(137):
        try:
            ind = dict_exc_index.get(f'act_{i+1}')
            data_exc[:,i] = solution[(num_q + num_u + num_inputs + 3 + ind)*num_nodes:(num_q + num_u + num_inputs + 3 + ind + 1)*num_nodes]
        except:
            continue

    return data_act, data_exc

def quat2eul(quat,seq):
    rotm = Qrm_np(quat)

    if seq == 'YZY':
        z = np.arccos(rotm[1,1])
        yy = np.arctan2(rotm[1,2],rotm[1,0])
        y = np.arctan2(rotm[2,1],-rotm[0,1])
        eul = np.array([y,z,yy])
    elif seq == 'YZX':
        z = np.arcsin(rotm[1,0])
        x = np.arctan2(-rotm[1,2],rotm[1,1])
        y = np.arctan2(-rotm[2,0],rotm[0,0])
        eul = np.array([y,z,x])

    return eul

def eul2quat(eul,seq):
    rot = spat.from_euler(seq,eul)
    quat = rot.as_quat(scalar_first=True)
    
    return quat

def Qrm_np(q):
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]
    res =  np.array([[1-2*(y**2+z**2), 2*(x*y-z*w), 2*(x*z+y*w),0],
                      [2*(x*y+z*w), 1-2*(x**2+z**2), 2*(y*z-x*w),0],
                      [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x**2+y**2),0],
                      [0           ,0          ,0             ,1]])

    return res

def mulQuat_sp(qa,qb):
    res = sp.Matrix([[qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3]],
                     [qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2]],
                     [qa[0]*qb[2] - qa[1]*qb[3] + qa[2]*qb[0] + qa[3]*qb[1]],
                     [qa[0]*qb[3] + qa[1]*qb[2] - qa[2]*qb[1] + qa[3]*qb[0]]])
    return res

def mulQuat_np(qa,qb):
    res = np.array([[qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3]],
                     [qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2]],
                     [qa[0]*qb[2] - qa[1]*qb[3] + qa[2]*qb[0] + qa[3]*qb[1]],
                     [qa[0]*qb[3] + qa[1]*qb[2] - qa[2]*qb[1] + qa[3]*qb[0]]])
    return res


def Qrm_sp(q):
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]
    res =  sp.Matrix([[1-2*(y**2+z**2), 2*(x*y-z*w), 2*(x*z+y*w),0],
                      [2*(x*y+z*w), 1-2*(x**2+z**2), 2*(y*z-x*w),0],
                      [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x**2+y**2),0],
                      [0           ,0          ,0             ,1]])

    return res

def YZY_seq(phi1,phi2,phi3):
    res = R_y(phi1) * R_z(phi2) * R_y(phi3)

    return res

def T_x(x):
    trans_x = sp.Matrix([[1,0,0,x],
                         [0,1,0,0],
                         [0,0,1,0],
                         [0,0,0,1]])
    return trans_x

def T_y(y):
    trans_y = sp.Matrix([[1,0,0,0],
                         [0,1,0,y],
                         [0,0,1,0],
                         [0,0,0,1]])
    return trans_y

def T_y(z):
    trans_y = sp.Matrix([[1,0,0,0],
                         [0,1,0,0],
                         [0,0,1,z],
                         [0,0,0,1]])
    return trans_y

def R_x_sp(phix):
    rot_phix = sp.Matrix([[1,0          ,0           ,0],
                         [0,sp.cos(phix),-sp.sin(phix),0],
                         [0,sp.sin(phix), sp.cos(phix),0],
                         [0,0          ,0           ,1]])
    return rot_phix

def R_y_sp(phiy):
    rot_phiy = sp.Matrix([[sp.cos(phiy) ,0,sp.sin(phiy),0],
                         [0           ,1,0          ,0],
                         [-sp.sin(phiy),0,sp.cos(phiy),0],
                         [0           ,0,0           ,1]])
    return rot_phiy

def R_z_sp(phiz):
    rot_phiz = sp.Matrix([[sp.cos(phiz),-sp.sin(phiz),0,0],
                        [sp.sin(phiz), sp.cos(phiz),0,0],
                        [0           ,0            ,1,0],
                        [0           ,0            ,0,1]])
    return rot_phiz

def R_x_np(phix):
    rot_phix = np.array([[1,0          ,0           ,0],
                         [0,np.cos(phix),-np.sin(phix),0],
                         [0,np.sin(phix), np.cos(phix),0],
                         [0,0          ,0           ,1]])
    return rot_phix

def R_y_np(phiy):
    rot_phiy = np.array([[np.cos(phiy) ,0,np.sin(phiy),0],
                         [0           ,1,0          ,0],
                         [-np.sin(phiy),0,np.cos(phiy),0],
                         [0           ,0,0           ,1]])
    return rot_phiy

def R_z_np(phiz):
    rot_phiz = np.array([[np.cos(phiz),-np.sin(phiz),0,0],
                        [np.sin(phiz), np.cos(phiz),0,0],
                        [0           ,0            ,1,0],
                        [0           ,0            ,0,1]])
    return rot_phiz

def position(x,y,z):
    r = np.array([x,y,z,1])
    return r

def G(quat):
    Q0 = quat[0];
    Q1 = quat[1];
    Q2 = quat[2];
    Q3 = quat[3];
    res = np.array([[-Q1, Q0, Q3, -Q2],
                    [-Q2,-Q3, Q0, Q1],
                    [-Q3, Q2, -Q1, Q0]])
    return res