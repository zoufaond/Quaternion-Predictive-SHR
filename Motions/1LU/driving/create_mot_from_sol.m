addpath ..\..\..\Matlab_functions\
load('res_quat_driving_init.mat')
traj_quat = data.trajectories;
time = data.tout;

mot_struct.quat = traj_quat;
mot_struct.euler = quat2eul_motion(traj_quat,'YZY');
mot_struct.time = time';
save('driving.mat','mot_struct');