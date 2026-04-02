clear all; clc;
names = ["EMG_3NA_Elevation","EMG_3NA_Scabduction","EMG_3NA_Flexion.mat"];
muscles = ["Pect","AnteriorDelt","Biceps","IntermediateDelt","Infrasp","Suprasp","MiddleTrap","UpperTrap","Triceps","PosteriorDelt","Serrupper","Serrlower"];
trials = {'4','3','1'};

for imus=1:length(muscles)
    % current_mus = EMG_current.([muscles(imus)]);
    EMG_all.num_1.(muscles(imus)) = [];
end

time2add = 100*0.028;

for i=1:3
    EMG_struct = load(names(i));
    EMG_current = EMG_struct.data.(['num_',trials{i}]);
    time = linspace(0,2.8,length(EMG_current.Pect));
    
    for imus=1:length(muscles)
        current_mus = EMG_current.([muscles(imus)]);
        time_new = linspace(0,2.8,100);
        EMG_new = interp1(time,current_mus,time_new);
        if i < 3
            EMG_all.num_1.(muscles(imus)) = [EMG_all.num_1.(muscles(imus)), EMG_new, nan(1,40)];
        else
            EMG_all.num_1.(muscles(imus)) = [EMG_all.num_1.(muscles(imus)), EMG_new];
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
save('EMG_3NA_All_motions_fig.mat','data');