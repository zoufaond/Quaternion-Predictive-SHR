clear all; clc;
names = ["EMG_1LU_Elevation","EMG_1LU_Scabduction","EMG_1LU_Flexion.mat"];
muscles = ["Pect","AnteriorDelt","Biceps","IntermediateDelt","Infrasp","Suprasp","MiddleTrap","UpperTrap","Triceps","PosteriorDelt","Serrupper","Serrlower"];
trials = {'3','3','2'};

for imus=1:length(muscles)
    % current_mus = EMG_current.([muscles(imus)]);
    EMG_all.num_1.(muscles(imus)) = [];
end

time2add = 40*0.03;

for i=1:3
    EMG_struct = load(names(i));
    EMG_current = EMG_struct.data.(['num_',trials{i}]);
    time = linspace(0,3,length(EMG_current.Pect));
    for imus=1:length(muscles)
        current_mus = EMG_current.([muscles(imus)]);
        if i < 3
            EMG_all.num_1.(muscles(imus)) = [EMG_all.num_1.(muscles(imus)), current_mus, zeros(1,ceil(time2add/time(2)))];
        else
            EMG_all.num_1.(muscles(imus)) = [EMG_all.num_1.(muscles(imus)), current_mus];
        end
    end
    % eul_traj = quat2eul_motion(traj_quat,'YZY');
    % eul_traj(:,8) = scale_dof(eul_traj(:,8),scalings(i));
    % if i < 3
    %     traj_all_eul = [traj_all_eul; eul_traj; zeros(40,10)];
    % else
    %     traj_all_eul = [traj_all_eul; eul_traj];
    % end
end

indexes = EMG_all.num_1.Pect~=0;

data = EMG_all;
data.indexes = indexes;
save('EMG_1LU_all.mat','data');