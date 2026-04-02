
addpath ..\..\..\Matlab_functions\
% elevation = load('res_quat_Elevation.mat');
% scabduction = load('res_quat_Scabduction.mat');
% flexion = load('res_quat_Flexion.mat');


names = ["res_quat_Elevation.mat","res_quat_Scabduction","res_quat_Flexion.mat"];


traj_all_eul = [];
scalings = [0.3,0.3,0.3];
timestep = 0.04;

for i=1:3
    mot_struct = load(names(i));
    traj_quat = mot_struct.data.trajectories;
    eul_traj = quat2eul_motion(traj_quat,'YZY');
    eul_traj(:,8) = scale_dof(eul_traj(:,8),scalings(i));
    if i < 3
        traj_all_eul = [traj_all_eul; eul_traj; zeros(40,10)];
    else
        traj_all_eul = [traj_all_eul; eul_traj];
    end
end

idx = traj_all_eul(:,1)~=0;
x = linspace(0,100,length(traj_all_eul(:,1)));
x1 = linspace(0,100,301);

traj_all_interp = [];
for i=1:10
    yn = interp1(x(idx),traj_all_eul(idx,i),x);
    xinter = x(1:8:length(x));
    yinter = yn(1:8:length(yn));
    yn1 = interp1(xinter,yinter,x1,'spline');
    traj_all_interp = [traj_all_interp,yn1'];
    % figure
    % plot(x,yn,xinter,yinter,x1,yn1)
    % legend('orig','inter','new')
end

time = linspace(0,0.028*length(x),length(yn1));
kinematics = eul2quat_motion(traj_all_interp,'YZY');
% 
% % 
mot_struct.euler = traj_all_interp;
mot_struct.quat = kinematics;
mot_struct.time = time';
% % 
save('All_motions.mat','mot_struct');