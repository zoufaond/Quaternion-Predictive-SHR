clear all
addpath Matlab_functions\
participant = '3NA';
motion_name = 'All_motions';
IK_file = 'All_motions.mot';
GH_seq = 'YZY';
osim_file = ['Motions\',participant,'\scaled_3NA_GH.osim'];
OS_model = das3_readosim(osim_file);
polyfiles_dir = 'Polyfiles';
musclepoly_file = 'musclepoly';
%%
% % Load and create training/simulation data from inverse kinematics
IK_data = loadmot(['Motions\',participant,'\',motion_name,'\',IK_file]);
IK_data(:,15) = scale_dof(IK_data(:,15),0.0);
[time_training, trajectory_training] = motion_polyfit(IK_data,40,40,0);
[time_simulation, trajectory_simulation] = motion_polyfit(IK_data,20,100,0);
% %%
% IK_struct1 = load(['Motions\11NI\Elevation_yzy\res_euler_Elevation_yzy_80.mat']);
% IK_struct2 = load(['Motions\11NI\Scabduction_yzy\res_euler_Scabduction_yzy_80.mat']);
% IK_struct3 = load(['Motions\11NI\Flexion_yzy\res_euler_Flexion_yzy_80.mat']);
% trajectory_training = [IK_struct1.data.trajectories,ones(101,1)*120*pi/180;
%                        IK_struct2.data.trajectories,ones(101,1)*120*pi/180;
%                        IK_struct3.data.trajectories,ones(101,1)*120*pi/180];
% IK_struct1 = load(['Motions\3NA\Elevation_yzy\Elevation_yzy.mat']);
% IK_struct2 = load(['Motions\3NA\Scabduction_yzy\Scabduction_yzy.mat']);
% IK_struct3 = load(['Motions\3NA\Flexion_yzy\Flexion_yzy.mat']);
% trajectory_training = [IK_struct1.mot_struct.euler,ones(100,1)*120*pi/180;
%                        IK_struct2.mot_struct.euler,ones(100,1)*120*pi/180;
%                        IK_struct3.mot_struct.euler,ones(100,1)*120*pi/180];

%%
% % trajectory_min_con = zeros(size(trajectory_simulation));
for i=1:size(trajectory_simulation,1)
    trajectory_min_con(i,:) = min_conoid_length(trajectory_simulation(i,:),OS_model);
    % trajectory_min_con(i,:) = change_clavx(trajectory_simulation(i,:),0);
    % conoid_length(i) = get_conoid_length(trajectory_min_con(i,4:6),OS_model);
    % conoid_force(i) = get_conoid_force(conoid_length(i),1);
end
% figure
% plot(conoid_length)
% figure
% plot(conoid_force)
% figure
% plot(trajectory_min_con(:,3)*180/pi)
% % Create .mot files of simulation
data2mot(trajectory_min_con,time_simulation,['Motions\',participant,'\',motion_name,'\',motion_name,'_simulation.mot'], 'euler', 'struct', GH_seq);
% % % % % create struct with simulation
data2struct(trajectory_min_con,time_simulation,['Motions\',participant,'\',motion_name,'\',motion_name],GH_seq)
%%
% correct = 0;
% while correct == 0
%     try
clavicle_xrot_vals = deg2rad(linspace(-10,50,8));
num_noised = 3;
noise_upper = deg2rad( [35,35,0,35,35,35,25,25,25,115,10])/1;
noise_lower = deg2rad(-[35,35,0,35,35,35,25,25,25,10,10])/1;
noise_nearGH = deg2rad([25,-25]);
[time_training_noised,data_training_noised] = noise2data(trajectory_training, clavicle_xrot_vals,noise_upper,noise_lower,noise_nearGH,num_noised, OS_model,osim_file,GH_seq);

if strcmp(GH_seq,'YZY')
    data_training_noised = GH_yzx2yzy(data_training_noised);
end

data2mot(data_training_noised,time_training_noised,['Motions\',participant,'\',motion_name,'\',motion_name,'_training.mot'], 'euler', 'struct', GH_seq);


%%
training_data = ['Motions\',participant,'\',motion_name,'\',motion_name,'_training.mot'];

das3_polynomials(osim_file,polyfiles_dir,training_data,'YZY','gen_polyvalues',musclepoly_file); %'gen_polyvalues'
model = das3_readosim(osim_file,['Motions\',participant,'\',motion_name,'\',polyfiles_dir,'\',musclepoly_file],['Motions\',participant,'\',motion_name,'\',polyfiles_dir,'\',musclepoly_file]);
save(['Motions\',participant,'\',motion_name,'\','OS_model'],'model')
disp("jsme tu")
%         correct = 1;
%     catch
%         correct = 0;
%         continue
%     end
% end
%%
% training_data = ['Motions\',participant,'\',motion_name,'\',motion_name,'_training.mot'];
% 
% das3_polynomials(osim_file,polyfiles_dir,training_data,'YZY','gen_polyvalues',musclepoly_file); %'gen_polyvalues'
% model = das3_readosim(osim_file,['Motions\',participant,'\',motion_name,'\',polyfiles_dir,'\',musclepoly_file],['Motions\',participant,'\',motion_name,'\',polyfiles_dir,'\',musclepoly_file]);
% save(['Motions\',participant,'\',motion_name,'\','OS_model'],'model')

% %% Recalculate to YZX 
% 
% % if strcmp(GH_seq,'YZY')
% %     trajectory_training = GH_yzy2yzx(trajectory_training);
% % end
% 
% % Add noise to training data
% clavicle_xrot_vals = deg2rad(linspace(-20,50,3));
% num_noised = 4;
% noise_upper = deg2rad( [35,35,0,35,35,30,20,55,30,110,10]);
% noise_lower = deg2rad(-[35,35,0,35,35,30,20,15,20,10,10]);
% noise_nearGH = deg2rad([50,-50]);
% [time_training_noised,data_training_noised] = noise2data(trajectory_training, clavicle_xrot_vals,noise_upper,noise_lower,noise_nearGH,num_noised, OS_model,osim_file,GH_seq);
% 
% % if strcmp(GH_seq,'YZY')
% %     data_training_noised = GH_yzx2yzy(data_training_noised);
% % end
% %%
% data2mot(data_training_noised,time_training_noised,['Motions\',participant,'\',motion_name,'\',motion_name,'_training.mot'], 'euler', 'struct', GH_seq);
% 
% %%
training_data = ['Motions\',participant,'\',motion_name,'\',motion_name,'_training.mot'];
das3_polynomials(osim_file,polyfiles_dir,training_data,'YZY','nogen_polyvalues',musclepoly_file); %'gen_polyvalues'
model = das3_readosim(osim_file,['Motions\',participant,'\',motion_name,'\',polyfiles_dir,'\',musclepoly_file],['Motions\',participant,'\',motion_name,'\',polyfiles_dir,'\',musclepoly_file]);
save(['Motions\',participant,'\',motion_name,'\','OS_model'],'model')