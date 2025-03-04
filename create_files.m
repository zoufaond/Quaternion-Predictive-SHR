clearvars
clear all
addpath Matlab_functions\
motion_name = 'Scabduction_yzx';
IK_file = 'IK.mot';
GH_seq = 'YZX';
osim_file = ['das3_',GH_seq,'.osim'];
OS_model = das3_readosim(osim_file);
polyfiles_dir = 'Polyfiles';
musclepoly_file = 'musclepoly';

%% Load and create training/simulation data from inverse kinematics
IK_data = loadmot(['Motions\',motion_name,'\',IK_file]);
IK_data(:,15) = scale_dof(IK_data(:,15),0.17);

%%
[time_training, trajectory_training] = motion_polyfit(IK_data,40,40,0);
[time_simulation, trajectory_simulation] = motion_polyfit(IK_data,30,100,0);
thorax_dim_optimized = optimize_thorax_dim(trajectory_simulation,OS_model);

%% Recalculate to YZX
if strcmp(GH_seq,'YZX')
    trajectory_training = GH_yzy2yzx(trajectory_training);
    trajectory_simulation = GH_yzy2yzx(trajectory_simulation);
end

%% Add noise to training data
clavicle_xrot_vals = deg2rad(linspace(-10,40,20));
num_noised = 200;
noise_upper = deg2rad( [30,20,0,20,20,20,20,15,20,10,10]);
noise_lower = deg2rad(-[10,20,0,20,20,20,20,15,20,10,10]);
noise_nearGH = deg2rad([70,-70]);
[time_training_noised,data_training_noised] = noise2data(trajectory_training, clavicle_xrot_vals,noise_upper,noise_lower,noise_nearGH,num_noised, OS_model,GH_seq, thorax_dim_optimized);
data2mot(data_training_noised,time_training_noised,['Motions\',motion_name,'\',motion_name,'_training.mot'], 'euler', 'struct', GH_seq);

%% Create .mot files of simulation
data2mot(trajectory_simulation,time_simulation,['Motions\',motion_name,'\',motion_name,'_simulation.mot'], 'euler', 'struct', GH_seq);
% create struct with simulation
data2struct(trajectory_simulation,time_simulation,['Motions\',motion_name,'\',motion_name],GH_seq)

%%
training_data = ['Motions\',motion_name,'\',motion_name,'_training.mot'];
das3_polynomials(osim_file,polyfiles_dir,training_data,GH_seq,0,musclepoly_file); %'gen_polyvalues'
model = das3_readosim(osim_file,['Motions\',motion_name,'\',polyfiles_dir,'\',musclepoly_file]);
model.thorax_dim_optimized = thorax_dim_optimized;
save(['Motions\',motion_name,'\','OS_model'],'model')