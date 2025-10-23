addpath ..\..\..\Matlab_functions\
load('Elevation_prediction.mat')
traj_quat = data.trajectories;
time = data.tout;

mot_struct.quat = traj_quat;
mot_struct.euler = quat2eul_motion(traj_quat,'YZY');
mot_struct.time = time';
save('Elevation_prediction.mat','mot_struct');