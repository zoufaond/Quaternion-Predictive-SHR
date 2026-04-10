% This is an example how to create the input files used in optimal control
% The example folder is in Motions/par2/Polynomials example,
% Scabduction_IK.mot motion from inverse kinematics was used as an example.
%

addpath Matlab_functions\
participant = 'par2';
motion_folder = 'Polynomials_example';
IK_file = 'Scabduction_IK.mot';
motion_name = 'Scabduction'; % with this name the simulation file used in optimal control will be created
GH_seq = 'YZY';
osim_file = ['Motions\',participant,'\scaled_par2.osim']; % path to scaled opensim model
OS_model = das3_readosim(osim_file); % read OS model
polyfiles_dir = 'Polyfiles'; % Folder that will contain generated muscle-tendon paths data. Here it is located in Motions/par2/Polynomials_example/Polyfiles
musclepoly_file = 'musclepoly';
%%
% % Load and create training/simulation data from inverse kinematics
% IK file contains rotation of thorax, simulation data do not.
IK_data = loadmot(['Motions\',participant,'\',motion_folder,'\',IK_file]);
[time_training, trajectory_training] = motion_polyfit(IK_data,40,40,0);
[time_simulation, trajectory_simulation] = motion_polyfit(IK_data,20,100,0);

%% We calculate axial rotation of clavicle based on minimal conoid length.
for i=1:size(trajectory_simulation,1)
    trajectory_min_con(i,:) = min_conoid_length(trajectory_simulation(i,:),OS_model);
end

% % Create .mot files of simulation
data2mot(trajectory_min_con,time_simulation,['Motions\',participant,'\',motion_folder,'\',motion_name,'_simulation.mot'], 'euler', 'struct', GH_seq);
% create .mat struct with simulation to be tracked
data2struct(trajectory_min_con,time_simulation,['Motions\',participant,'\',motion_folder,'\',motion_name],GH_seq)
%% Create sample data for polynomials by adding noise to map wide range of configurations 
clavicle_xrot_vals = deg2rad(linspace(-10,50,8));
num_noised = 3;
noise_upper = deg2rad( [35,35,0,35,35,35,25,25,25,115,10])/1;
noise_lower = deg2rad(-[35,35,0,35,35,35,25,25,25,10,10])/1;
noise_nearGH = deg2rad([25,-25]);
[time_training_noised,data_training_noised] = noise2data(trajectory_training, clavicle_xrot_vals,noise_upper,noise_lower,noise_nearGH,num_noised, OS_model,osim_file,GH_seq);

% Save training data
data2mot(data_training_noised,time_training_noised,['Motions\',participant,'\',motion_folder,'\',motion_name,'_training.mot'], 'euler', 'struct', GH_seq);


%%
% Load training data
training_data = ['Motions\',participant,'\',motion_folder,'\',motion_name,'_training.mot'];

% Generate lengths, moment arms and caluclate polynomial approximation.
das3_polynomials(osim_file,polyfiles_dir,training_data,'YZY','ngen_polyvalues',musclepoly_file);

% Read OS model with polynomials and save it to OS_model.mat
model = das3_readosim(osim_file,['Motions\',participant,'\',motion_folder,'\',polyfiles_dir,'\',musclepoly_file],['Motions\',participant,'\',motion_folder,'\',polyfiles_dir,'\',musclepoly_file]);
save(['Motions\',participant,'\',motion_folder,'\','OS_model'],'model')