addpath Matlab_functions/

%choose two results struct to compare
motion_name = 'Flexion';
GH_seq = 'YZY';
weighteul = '102';
weightquat = '1002';
muscle_group = {'delt'};
% folders = {[folder_path,'results_quat_QuatInit_100.mat']};


folder_path = ['Motions/',motion_name,'/'];
OS_model = [folder_path,'OS_model.mat'];
OS_struct = load([folder_path,motion_name,'.mat']);
folders = {[folder_path,'res_euler_',motion_name,'_',weighteul,'.mat'],[folder_path,'res_quat_',motion_name,'_',weightquat,'.mat']};
% folders = {[folder_path,'results_euler_QuatInit_',weighteul,'.mat'],[folder_path,'results_quat_QuatInit_',weightquat,'.mat']};
% plot_kinematics(folders,OS_struct,GH_seq);
% % plot_conoid_length(folders);
% plot_SCx(folders);
% plot_ESB(folders, OS_model, OS_struct)
% plot_elipsoid_eq([folder_path,'res_euler_',motion_name,'_',weighteul,'.mat'],OS_model)

% plot_activations(folders,muscle_group,OS_model,OS_struct);
% plot_muscles(folders,muscle_group,OS_model, GH_seq);
% 
% 
% 
% EMG_struct = 'Experimental_data/EMG_struct.mat';
% plot_activations_EMG(folders, EMG_struct, OS_model);
% plot_activations(folders,muscle_group,OS_model)
% mus_group = 'ter';
% mus_group = ["delt_scap10","delt_scap11","delt_clav_1","delt_clav_2","delt_clav_3","serr_ant_1","serr_ant_2","serr_ant_3","trap_clav_1","trap_clav_2","trap_scap10","trap_scap11"];
% % 
res_eul_yzx = {'Scabduction_yzx','euler','102','YZX'};
res_quat_yzx = {'Scabduction_yzx','quat','252','YZX'};
res_eul_yzy = {'Scabduction_yzy', 'euler','103','YZY'};
plot_GH_seq = 'YZX';
motion_name = 'Scabduction_yzx';
OS_struct = load(['Motions/',motion_name,'/',motion_name,'.mat']);
% plot_paper_kinematics(OS_struct,plot_GH_seq,res_eul_yzx,res_quat_yzx,res_eul_yzy)
plot_paper_activations(OS_struct,OS_model,mus_group,res_eul_yzx,res_quat_yzx,res_eul_yzy)
% plot_objective_values(res_quat_yzx,res_eul_yzy)

% res_eul_yzy = {'Flexion','euler','102','YZY'};
% res_quat_yzy = {'Flexion','quat','602','YZY'};
% plot_GH_seq = 'YZY';
% motion_name = 'Flexion';
% OS_struct = load(['Motions/',motion_name,'/',motion_name,'.mat']);
% plot_paper_kinematics(OS_struct,plot_GH_seq,res_quat_yzy)
% plot_paper_activations(OS_struct,OS_model,mus_group,res_quat_yzy)

function plot_objective_values(varargin)
    num_res = length(varargin);
    figure
    for ires = 1:num_res
        iresult = varargin{ires};
        motion_name = iresult{1};
        rot_type = iresult{2};
        weight = iresult{3};
        GH_seq = iresult{4};
        ires = load(['Motions\',motion_name,'\res_',rot_type,'_',motion_name,'_',weight,'.mat']);
        obj_values = ires.data.objective_value;
        yscale log
        plot(obj_values); hold on
    end
    hold off

    
end

function plot_paper_activations(OS_struct,OS_model,mus_group,varargin)
    model = load(OS_model);
    muscles = model.model.muscles;
    for i = 1:length(muscles)
        muscle_names{i} = muscles{i}.osim_name;
    end

    colors_act = [0, 0, 1, 1; 1, 0, 0, 1; 0.4660, 0.6740, 0.1880, 1];
    colors_exc = [0, 0, 1, 0.5; 1, 0, 0, 0.5;0.4660, 0.6740, 0.1880, 0.5];


    mask = startsWith(muscle_names,mus_group);
    num_in_group = nnz(mask);
    plot_rows = ceil(num_in_group/3);
    if rem(num_in_group,3) == 0
        plot_rows = plot_rows + 1;
    end
    current_names = muscle_names(mask);
    legend_names = {};
    num_res = length(varargin);
    figure
    for ires = 1:num_res
        iresult = varargin{ires};
        motion_name = iresult{1};
        rot_type = iresult{2};
        weight = iresult{3};
        GH_seq = iresult{4};
        result = load(['Motions\',motion_name,'\res_',rot_type,'_',motion_name,'_',weight,'.mat']);
    
        activations = result.data.activations(:,mask);
        excitations = result.data.excitations(:,mask);
        time = result.data.tout;
        

        for imus = 1:num_in_group
            subplot(plot_rows,3,imus)
            plot(time, activations(:,imus),'Color',colors_act(ires,:),'LineWidth',1); hold on
            plot(time, excitations(:,imus),'Color',colors_exc(ires,:),'LineWidth',1); hold on
            axis([-inf inf -inf inf])  
            title(current_names{imus},'Interpreter','none');
            xlabel('Time [s]')
            ylabel('Activation [-]')
        end
        
    if strcmp(rot_type,'euler')
        legend_names{end+1} = ['Euler ',GH_seq];
        legend_names{end+1} = [''];
    elseif strcmp(rot_type,'quat')
        legend_names{end+1} = ['Quaternion'];
        legend_names{end+1} = [''];
    end

    fig = gcf;
    
    % fig.Position(3:4) = [1200,900];
    fig.Position(3) = fig.Position(3) + 100;
    fig.Position(4) = fig.Position(4) + 100;
    Lgnd = legend(legend_names,'FontSize',15);
    Lgnd.Position(1) = 0.7;
    Lgnd.Position(2) = 0.1;
    Lgnd.Direction = "reverse";
    
    % sgtitle('Activations', 'Interpreter', 'none')
    end
    hold off
    movegui(fig,'center');
end

function plot_paper_kinematics(OS_struct,plot_GH_seq,varargin)
    num_res = length(varargin);
    legend_names = {};
    line_styles = {'--','-.','-'};
    line_colors = [0, 0, 1, 1; 1, 0, 0, 1; 0.4660, 0.6740, 0.1880, 1];
    figure
    for ires = 1:num_res
        iresult = varargin{ires};
        motion_name = iresult{1};
        rot_type = iresult{2};
        weight = iresult{3};
        GH_seq = iresult{4};
        result = load(['Motions\',motion_name,'\res_',rot_type,'_',motion_name,'_',weight,'.mat']);
        time = result.data.tout;
        traj = result.data.trajectories;
        dofs_names = {'SCy','SCz','SCx (not tracked)','Yscap in global','Zscap in global','Xscap in global','Yhum in global','Zhum in global',[GH_seq(end),'hum in global'],'ELx'};
    
        if strcmp(rot_type,'euler')
            traj_obj = create_objective_traj_eul(traj,GH_seq);
        elseif strcmp(rot_type,'quat')
            traj_in_eul = quat2eul_motion(traj,GH_seq);
            traj_obj = create_objective_traj_eul(traj_in_eul,GH_seq);
        end
        traj_obj = rad2deg(traj_obj);

        for i = 1:length(dofs_names)
        subplot(4,3,i)
        plot(time,traj_obj(:,i),'Color',line_colors(ires,:),'Linestyle',line_styles{ires},'LineWidth',1.5)
        hold on
        title(dofs_names{i})
        


        end
     if strcmp(rot_type,'euler')
        legend_names{end+1} = ['Euler ',GH_seq];
    elseif strcmp(rot_type,'quat')
        legend_names{end+1} = 'Quaternion';
     end
    
    end

    OS_interp = interp1(OS_struct.mot_struct.time,OS_struct.mot_struct.euler,time',"spline");
    OS_interp_obj = rad2deg(create_objective_traj_eul(OS_interp,plot_GH_seq));
    legend_names{end+1} = 'Inverse Kinematics';

    for i = 1:length(dofs_names)
        subplot(4,3,i)
        if i == 3
            continue
        else
            plot(time, OS_interp_obj(:,i),'g','LineWidth',1.5)
            % axis([-inf inf -inf inf])  
            
            xlabel('Time [s]')
            ylabel('Angle [deg]')

        end
        % OS_interp_quat = eul2quat_motion(OS_interp,GH_seq);
        % OS_interp_quat_obj = create_objective_traj_quat(OS_interp_quat);
    end

    hold off
    fig = gcf;
    Lgnd = legend(legend_names,'FontSize',15);
    Lgnd.Position(1) = 0.45;
    Lgnd.Position(2) = 0.08;
    Lgnd.Direction = "reverse";
    fig.Position(3) = fig.Position(3) + 250;

    % exportgraphics(fig,'kinematics_flexion.png','Resolution',600);

end

function plot_elipsoid_eq(OS_struct,OS_model)
    struct = load(OS_struct);
    traj_eul = struct.data.trajectories;
    model = load(OS_model);
    thorax_dim = model.model.thorax_dim_optimized;
    thorax_pos = model.model.thorax_center;
    numdata = size(traj_eul,1);
    for i = 1:numdata
        [iTS,iAI] = contact_points_position(traj_eul(i,:),model.model);
        eq_TS(i) = elips_eq(iTS,thorax_pos,thorax_dim);
        eq_AI(i) = elips_eq(iAI,thorax_pos,thorax_dim);
    end

    figure
    plot(eq_TS); hold on
    plot(eq_AI)

    figure
    plot(sqrt(eq_TS.^2+eq_AI.^2))

end

function plot_SCx(folders)
    GH_names = {'GHy','GHz','GHyy'};
    num_coords = 10;
    euler_struct = load(folders{1});
    time = euler_struct.data.tout;
    traj_eul = euler_struct.data.trajectories;
    quat_struct = load(folders{2});
    traj_quat_orig = quat_struct.data.trajectories;
    traj_quat_eul = quat2eul_motion(traj_quat_orig);

    for i = 1:length(time)
        newSCx_eul(i,:) = min_conoid_length(traj_eul(i,:));
        length_eul(i,:) = conoid_length(traj_eul(i,4:6));
        newSCx_quat(i,:) = min_conoid_length(traj_quat_eul(i,:));
    end

    figure
    subplot(2,1,1)
    plot(time,rad2deg(newSCx_eul(:,3)),time,rad2deg(traj_eul(:,3)),'*')
    legend('Min length','simulation')
    subplot(2,1,2)
    plot(time,length_eul)
    title('eul')

    figure
    plot(time,rad2deg(newSCx_quat(:,3)),time,rad2deg(traj_quat_eul(:,3)),'*')
    legend('Min length','simulation')
    title('quat')
    

end

function plot_ESB(folders, OS_model, OS_struct)

euler_struct = load(folders{1});
quat_struct = load(folders{2});
time = euler_struct.data.tout;
time_prcnt = time/time(end)*100;

muscle = 42;
euler_act_all = euler_struct.data.inputs;
euler_act_cur = euler_act_all(:,muscle);
quat_act_all = quat_struct.data.inputs;
quat_act_cur = quat_act_all(:,muscle);
model = load(OS_model);
OS_interp = interp1(OS_struct.mot_struct.time,OS_struct.mot_struct.euler,time',"spline");
gimbal_locks = find_gimbal_lock(time,OS_interp(:,8));

figure

plot(time_prcnt,euler_act_cur,time_prcnt,quat_act_cur,'LineWidth',1)
xline(gimbal_locks([2,3])/time(end)*100,'k-.','LineWidth',1)
legend('Euler model','Quaternion model','Gimbal lock in GH')
ylabel('Activation [-]')
xlabel('Humeral elevation cycle [%]')
axis([-inf inf 0 0.25])
fig = gcf;
fig.Position(3:4)=[400,200];
title('Lateral deltoid muscle element')
exportgraphics(fig,'Delt_scap_6.png','Resolution',600)


end

function plot_activations_EMG(folders, EMG_struct, OS_model)

EMG_muscles = {'AnteriorDelt','IntermediateDelt','PosteriorDelt','Infrasp','MiddleTrap','UpperTrap','Serrupper', 'Serrlower'};
model_names = {'delt_scap11','delt_scap10','delt_scap_3','infra_1','trap_scap_7','trap_scap10','serr_ant_5','serr_ant_2'};

simulation = load(folders{2});
ending_val = 0;
time_simulation = simulation.data.tout(1:end-ending_val);

emg_data = load(EMG_struct);
model = load(OS_model);
muscles = model.model.muscles;
num_muscles = length(muscles);

for i = 1:num_muscles
    muscle_names{i} = muscles{i}.osim_name;
end


figure
tiledlayout(4,2);

current_emg_rsmpld = zeros(length(time_simulation),size(emg_data.data.(EMG_muscles{1}),2));
for i = 1:length(EMG_muscles)
    nexttile
    current_emg_full = emg_data.data.(EMG_muscles{i});
    current_emg = current_emg_full(1:end/2,:);
    time = linspace(0,1.5,size(current_emg,1));
    for ipar = 1:size(current_emg,2)
        current_emg_rsmpld(:,ipar) = spline(time,current_emg(:,ipar),time_simulation);
    end
    [S,M] = std(current_emg_rsmpld,0,2);
    upper_bound = M+S;
    lower_bound = M-S;
    plot(time_simulation,M,'k')
    hold on
    plot(time_simulation,upper_bound,'k',time_simulation,lower_bound,'k')
    hold on
    mus_index = find(strcmp(muscle_names,model_names{i}));
    plot(time_simulation,simulation.data.inputs(1:end-ending_val,mus_index),'r','LineWidth',1.5)
    hold on
    % patch([time' fliplr(time')], [lower_bound fliplr(upper_bound)], 'g')
    fill([time_simulation'; flip(time_simulation')],[lower_bound; flip(upper_bound)], 'b', 'edgecolor', 'none', 'facealpha', 0.1);
    title(EMG_muscles{i})
    axis([-inf inf 0 inf])
    xlabel('Time[s]')
    ylabel('Activation[-]')
    legend('Exp','','','Sim');
end

% fig = gcf;
% fig.Position(3) = fig.Position(3) + 250;
% Lgnd = legend('Exp','','','Sim');
% Lgnd.Position(1) = 0.9;
% Lgnd.Position(2) = 0.9;
% sgtitle('Activations', 'Interpreter', 'none')


end

function plot_activations(folders, muscle_group, OS_model, OS_struct)

model = load(OS_model);
muscles = model.model.muscles;
num_muscles = length(muscles);
euler_struct = load(folders{1});
time = euler_struct.data.tout;
OS_interp = interp1(OS_struct.mot_struct.time,OS_struct.mot_struct.euler,time',"spline");
gimbal_locks = find_gimbal_lock(time,OS_interp(:,8));

for i = 1:num_muscles
    muscle_names{i} = muscles{i}.osim_name;
end

colors_act = [0, 0, 1, 1; 1, 0, 0, 1];
colors_exc = [0, 0, 1, 0.5; 1, 0, 0, 0.5];

for imus_group = 1:length(muscle_group)
    current_mus_group = muscle_group{imus_group};
    mask = startsWith(muscle_names,current_mus_group);
    num_in_group = nnz(mask);
    plot_rows = ceil(num_in_group/3);
    current_names = muscle_names(mask);


    figure
    tiledlayout(plot_rows,3);
    for j = 1:num_in_group
        nexttile
        for i = 1:length(folders)
        current_struct = load(folders{i});
        activations = current_struct.data.activations(:,mask);
        excitations = current_struct.data.excitations(:,mask);
        if i == 1
            activations_eul_obj = sum(activations.^2,'all');
        end
        if i == 2
            activations_quat_obj = sum(activations.^2,'all');
        end
        time = current_struct.data.tout;
        plot(time,activations(:,j),'Color',colors_act(i,:),'LineWidth',1);hold on
        plot(time,excitations(:,j),'Color',colors_exc(i,:),'LineWidth',1)
        xlabel('time')
        ylabel('a[-]')
        
        end
        title(current_names(j), 'Interpreter', 'none')
        axis([-inf,inf,0,inf])
        for GL = 1:length(gimbal_locks)
            xline(gimbal_locks(GL))
        end
        

    end
    
    fig = gcf;
    fig.Position(3) = fig.Position(3) + 250;
    Lgnd = legend('eul_{act}','eul_{exc}','quat_{act}','q_{exc}');
    Lgnd.Position(1) = 0.01;
    Lgnd.Position(2) = 0.5;
    text1 = annotation('textbox',[0.7,0.12,0.1,0.1],'string',{['Act Eul = ', num2str(activations_eul_obj),newline,'Act Quat = ', num2str(activations_quat_obj)]});
    sgtitle('Activations', 'Interpreter', 'none')
end

end

function plot_kinematics(folders, OS_struct, GH_seq)
dofs_names = {'SCy','SCz','SCx','Yscap in global','Zscap in global','Xscap in global','Yhum in global','Zhum in global','YYhum in global','ELx','PSy'};
GH_names = {'GHy','GHz','GHyy'};
num_coords = 10;
euler_struct = load(folders{1});
time = euler_struct.data.tout;
traj_eul = euler_struct.data.trajectories;
quat_struct = load(folders{2});
traj_quat_orig = quat_struct.data.trajectories;
traj_quat_eul = quat2eul_motion(traj_quat_orig,GH_seq);
traj_eul_obj = create_objective_traj_eul(traj_eul,GH_seq);
traj_quat_eul_obj = create_objective_traj_eul(traj_quat_eul,GH_seq);
OS_interp = interp1(OS_struct.mot_struct.time,OS_struct.mot_struct.euler,time',"spline");
OS_interp_obj = create_objective_traj_eul(OS_interp,GH_seq);
OS_interp_quat = eul2quat_motion(OS_interp,GH_seq);
OS_interp_quat_obj = create_objective_traj_quat(OS_interp_quat);
traj_quat_obj = create_objective_traj_quat(traj_quat_orig);
gimbal_locks = find_gimbal_lock(time,OS_interp(:,8));

euler_act_obj = sum(euler_struct.data.activations.^2,'all');
quat_act_obj = sum(quat_struct.data.activations.^2,'all');

figure
idof = 1;
iplot = 1;
for i = 1:10
    subplot(5,3,iplot)
    plot(time,rad2deg(traj_eul_obj(:,idof)),'--','LineWidth',1.5)
    hold on
    plot(time,rad2deg(traj_quat_eul_obj(:,idof)),'-.','LineWidth',1.5)
    hold on
    for GL = 1:length(gimbal_locks)
        xline(gimbal_locks(GL))
    end
    title(dofs_names{idof})
    if i ~= 3
        plot(time,rad2deg(OS_interp_obj(:,idof)),'g','LineWidth',1.5)
    end
    axis([-inf inf -inf inf])  

    xlabel('Time [s]')
    ylabel('Angle [deg]')

    if idof == 9
        for GHdof = 1:3
        GHidof = GHdof+6;
        GHiplot = GHdof+9;
        subplot(5,3,GHiplot)
        plot(time,rad2deg(traj_eul(:,GHidof)),'--','LineWidth',1.5)
        hold on
        plot(time,rad2deg(traj_quat_eul(:,GHidof)),'-.','LineWidth',1.5)
        hold on
        title(GH_names{GHdof})
        plot(time,rad2deg(OS_interp(:,GHidof)),'g','LineWidth',1.5)
        for GL = 1:length(gimbal_locks)
            xline(gimbal_locks(GL))
        end

        if GHdof == 2
            yline(0,'r')
        end
        axis([-inf inf min(rad2deg(OS_interp(:,GHidof)))-10 max(rad2deg(OS_interp(:,GHidof)))+10])
    
        xlabel('Time [s]')
        ylabel('Angle [deg]')
        end
        iplot = iplot+3;
    end
    idof = idof+1;
    iplot = iplot+1;

end


fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('Euler','Quat','Experiment');
Lgnd.FontSize = 13;
Lgnd.Position(1) = 0.4;
Lgnd.Position(2) = 0.1;
text1 = annotation('textbox',[0.7,0.12,0.1,0.1],'string',{['Eul Act = ', num2str(euler_act_obj),newline,'Quat Act = ', num2str(quat_act_obj),]});

%%% QUAT KINEMATICS %%%

figure
idof = 1;
iplot = 1;
for i = 1:13
    subplot(5,4,i)
    plot(time,(traj_quat_obj(:,i)),'--','LineWidth',1.5)
    hold on
    plot(time,OS_interp_quat_obj(:,i))
    % title(dofs_names{idof})
    axis([-inf inf -inf inf])  

    xlabel('Time [s]')
    ylabel('Angle [deg]')

end

fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('Quat','Experiment');
Lgnd.FontSize = 13;
Lgnd.Position(1) = 0.4;
Lgnd.Position(2) = 0.1;

%%% END QUAT KINEMATICS %%%

end

function plot_muscles(folders,muscle_group, OS_model, GH_seq)

struct_euler = load(folders{1});
motion_euler = struct_euler.data.trajectories;
time = struct_euler.data.tout';
numdata = size(motion_euler,1);
struct_quat = load(folders{2});
motion_quat_orig = struct_quat.data.trajectories;
motion_quat_euler = quat2eul_motion(motion_quat_orig,GH_seq);



model = load(OS_model);
num_mus = length(model.model.muscles);
for i=1:num_mus
    all_muscle_names{i} = model.model.muscles{i}.osim_name;
end

mask = startsWith(all_muscle_names, muscle_group);
muscle_names = all_muscle_names(1,mask);
alljoints = {'YZX','YZX',GH_seq};
joints_names = {'SCy', 'SCz', 'SCx', 'ACy', 'ACz', 'ACx'};
indexes = find(mask);

for imus = 1:length(muscle_names)
    current_mus = model.model.muscles{indexes(imus)};
    current_name = current_mus.osim_name;
    dof_names = current_mus.dof_names;
    JQuatInJEul = zeros(numdata,11);
    dof_indeces = current_mus.dof_indeces-3;
    [lengths_eul,jacobian_eul] = momarms(current_mus.Euler, dof_indeces, motion_euler);
    motion_quat_WO_real = motion_quat_orig(:,[2:4,6:8,10:12,13]);
    [lengths_quat,jacobian_quat_current] = momarms(current_mus.Quaternion, dof_indeces, motion_quat_WO_real);
    jacobian_quat = zeros(numdata,11);
    jacobian_quat(:,dof_indeces) = jacobian_quat_current;

    for iframe=1:numdata
        for j = 1:3
            ind3 = ((j-1)*3+1:(j-1)*3+3);
            ind4 = ((j-1)*4+1:(j-1)*4+4);
            JQuatInSpat = invJtrans(motion_quat_orig(iframe,ind4)) * jacobian_quat(iframe,ind3)';
            JQuatInJEul(iframe,ind3) = GeomJ(motion_quat_euler(iframe,ind3),alljoints{j})*(JQuatInSpat);
        end
         JQuatInJEul(iframe,10) = jacobian_quat(10);
         JQuatInJEul(iframe,11) = jacobian_quat(11);
    end

    [lengths_OS_euler, dLdq_euler] = opensim_get_polyvalues(motion_euler, indexes(imus), current_mus.dof_indeces,GH_seq);
    [lengths_OS_quat, dLdq_quat] = opensim_get_polyvalues(motion_quat_euler, indexes(imus), current_mus.dof_indeces,GH_seq);
    max_momarm = max([max(abs(dLdq_euler),[],'all'),max(abs(dLdq_quat),[],'all'),max(abs(JQuatInJEul),[],'all'),max(abs(jacobian_eul),[],'all')]) + 0.01;
    figure
    
    for j=1:length(current_mus.dof_indeces)
        subplot(length(current_mus.dof_indeces)+1,2,j*2-1)
        plot(time,dLdq_euler(:,j),time,-jacobian_eul(:,j))
        if j==1
            title(['Eul, Number of params:', string(current_mus.Euler.lparam_count)])
        end
        % plot(time,dLdq_euler(:,j),time,imomarms(current_mus.dof_indeces(j)-3,:))
        axis([-inf inf -max_momarm max_momarm])
        subplot(length(current_mus.dof_indeces)+1,2,j*2)
        plot(time,dLdq_quat(:,j),time,-JQuatInJEul(:,dof_indeces(j)))
        if j==1
            title(['Quat, Number of params:', string(current_mus.Quaternion.lparam_count)'])
        end
        axis([-inf inf -max_momarm max_momarm])
    end

    subplot(length(current_mus.dof_indeces)+1,2,2*length(current_mus.dof_indeces)+1)
    plot(time, lengths_OS_euler, time, lengths_eul)
    subplot(length(current_mus.dof_indeces)+1,2,2*length(current_mus.dof_indeces)+2)
    plot(time, lengths_OS_quat, time, lengths_quat)
    % title(current_name)
    legend('osim','my')
    sgtitle(current_name)

    if strcmp(current_name,'serr_ant_s')

        figure
        tiledlayout(3,3);
        nexttile([1,3])
        plot(time,lengths_eul,'--','LineWidth',1.5);hold on
        plot(time,lengths_quat,'-.','LineWidth',1.5);hold on
        plot(time,lengths_OS_euler,'g','LineWidth',1.5);hold off
        title('Muscle length and moment arms approximation')
        ylabel('length [m]')
        xlabel('time [s]')
        axis([-inf inf 0 0.2])

        for j=1:length(current_mus.dof_indeces)
            nexttile
            plot(time,-jacobian_eul(:,j),'--','LineWidth',1.5);hold on
            plot(time,-JQuatInJEul(:,dof_indeces(j)),'-.','LineWidth',1.5);hold on
            plot(time,dLdq_euler(:,j),'g','LineWidth',1.5);
            axis([-inf inf -max_momarm max_momarm])
            ylabel('MA [m]')
            xlabel('time [s]')
            title(joints_names{j})
        end
        fig = gcf;
        fig.Position(3:4)=[700,350];
        Lgnd = legend('Euler approx','Quat approx','OpenSim');
        Lgnd.FontSize = 11;
        Lgnd.Position(1) = 0.7;
        Lgnd.Position(2) = 0.65;

        LRMS_eul = sqrt(sum((lengths_OS_euler-lengths_eul').^2,'all')/numdata)*1000
        LRMS_quat = sqrt(sum((lengths_OS_euler-lengths_quat').^2,'all')/numdata)*1000
        RRMS_eul = sqrt(sum((dLdq_euler+jacobian_eul(:,1:6)).^2,'all')/numdata)*1000
        RRMS_quat = sqrt(sum((dLdq_euler+JQuatInJEul(:,1:6)).^2,'all')/numdata)*1000

        exportgraphics(fig,'Serr_ant_12_approx.png','Resolution',600);
        
        

    end
end

% %%%%%%%%%%%%%%%%% plot with one trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for imus = 1:3 %length(muscle_names)
%     current_mus = model.model.muscles{indexes(imus)};
%     current_name = current_mus.osim_name;
%     dof_names = current_mus.dof_names;
%     JQuatInJEul = zeros(numdata,11);
%     dof_indeces = current_mus.dof_indeces-3;
%     jacobian_eul = momarms(current_mus.Euler, dof_indeces, motion_quat_euler);
%     motion_quat_WO_real = motion_quat_orig(:,[2:4,6:8,10:12,13]);
%     jacobian_quat_current = momarms(current_mus.Quaternion, dof_indeces, motion_quat_WO_real);
%     jacobian_quat = zeros(numdata,11);
%     jacobian_quat(:,dof_indeces) = jacobian_quat_current;
% 
%     for iframe=1:numdata
%         for j = 1:3
%             ind3 = ((j-1)*3+1:(j-1)*3+3);
%             ind4 = ((j-1)*4+1:(j-1)*4+4);
%             JQuatInSpat = invJtrans(motion_quat_orig(iframe,ind4)) * jacobian_quat(iframe,ind3)';
%             JQuatInJEul(iframe,ind3) = GeomJ(motion_quat_euler(iframe,ind3),alljoints{j})*(JQuatInSpat);
%         end
%          JQuatInJEul(iframe,10) = jacobian_quat(10);
%          JQuatInJEul(iframe,11) = jacobian_quat(11);
%     end
% 
%     [lengths_osim, dLdq_osim] = opensim_get_polyvalues(motion_quat_euler, indexes(imus), current_mus.dof_indeces);
%     max_momarm = max([max(abs(dLdq_euler),[],'all'),max(abs(dLdq_quat),[],'all'),max(abs(JQuatInJEul),[],'all'),max(abs(jacobian_eul),[],'all')]) + 0.01;
% 
%     figure
%     for j=1:length(current_mus.dof_indeces)
%         subplot(floor(length(current_mus.dof_indeces)/3),3,j)
%         plot(time,-jacobian_eul(:,j))
%         hold on
%         plot(time,-JQuatInJEul(:,dof_indeces(j)))
%         hold on
%         plot(time,dLdq_osim(:,j))
% 
%         axis([-inf inf -max_momarm max_momarm])
%     end
%     legend('Euler','Quaternion','Osim')
    % sgtitle([current_name,' Eul, Quat nterms:',string(current_mus.Euler.lparam_count),string(current_mus.Quaternion.lparam_count)]) 

% end

end

function time_positions = find_gimbal_lock(time,angles)
    numdata = length(angles);
    time_positions = [];
    for i = 1:numdata-1
        if angles(i) < 0 && angles(i+1) > 0 || angles(i) > 0 && angles(i+1) < 0 
            zero_crossing = -angles(i) * (time(i+1) - time(i))/(angles(i+1)-angles(i));
            time_positions = [time_positions time(i)+zero_crossing];
        end
    end
end


function index = find_mus_index(model,mus_names)
    muscles = model.model_full_eul.muscles;
    nmus = model.model_full_eul.nMus;

    for i = 1:length(mus_names)
        for j = 1:nmus
            if strcmp(mus_names(i),muscles{j}.osim_name)
                index(i) = j;
                break
            end
        end
    end

end

function res = invJtrans(quat)
    q1 = quat(1);
    q2 = quat(2);
    q3 = quat(3);
    q4 = quat(4);
    res = [ q1/2,  q4/2, -q3/2;
            -q4/2,  q1/2,  q2/2;
            q3/2, -q2/2,  q1/2];
end

function res = GeomJ(phi,seq)
    s2 = sin(phi(2));
    s3 = sin(phi(3));
    c2 = cos(phi(2));
    c3 = cos(phi(3));
    if seq == 'YZX'
        res = [s2, c2*c3, -s3*c2;0,s3,c3;1,0,0];

    elseif seq == 'YZY'
        res = [s2*c3, c2, s2*s3; -s3, 0 ,c3; 0, 1, 0];
    end
end

function res = mulQuat(qa,qb)
    res = [ qa(1)*qb(1) - qa(2)*qb(2) - qa(3)*qb(3) - qa(4)*qb(4);
            qa(1)*qb(2) + qa(2)*qb(1) + qa(3)*qb(4) - qa(4)*qb(3);
            qa(1)*qb(3) - qa(2)*qb(4) + qa(3)*qb(1) + qa(4)*qb(2);
            qa(1)*qb(4) + qa(2)*qb(3) - qa(3)*qb(2) + qa(4)*qb(1)];
end

function res = G(Q)
    Q0 = Q(1);
    Q1 = Q(2);
    Q2 = Q(3);
    Q3 = Q(4);
    res = [-Q1, Q0, Q3, -Q2;
            -Q2,-Q3, Q0, Q1;
            -Q3, Q2, -Q1, Q0];
end

function rot_phix = R_x(phix)
    rot_phix = [1,0        , 0        ,0;
                0,cos(phix),-sin(phix),0;
                0,sin(phix), cos(phix),0;
                0,0        , 0        ,1];
end

function rot_phiy = R_y(phiy)
    rot_phiy = [cos(phiy),0,sin(phiy),0;
                0        ,1,0        ,0;
               -sin(phiy),0,cos(phiy),0;
                0        ,0,0        ,1];
end

function rot_phiz = R_z(phiz)
    rot_phiz = [cos(phiz),-sin(phiz),0,0;
                sin(phiz), cos(phiz),0,0;
                0           ,0      ,1,0;
                0           ,0      ,0,1];
end

function res = YZX_seq(angles)
    res = R_y(angles(1)) * R_z(angles(2)) * R_x(angles(3));
end

function res = YZY_seq(angles)
    res = R_y(angles(1)) * R_z(angles(2)) * R_y(angles(3));
end

function r = position(vec)
    r = [vec(1);vec(2);vec(3);1];
end

function trans = T_trans(vec)
    trans = [1,0,0,vec(1);
               0,1,0,vec(2);
               0,0,1,vec(3);
               0,0,0,1];
end

function res = Qrm(q)
    % rotation matrix from quaternion
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);
    Rq =  [1-2*(y^2+z^2), 2*(x*y-z*w), 2*(x*z+y*w);
     2*(x*y+z*w), 1-2*(x^2+z^2), 2*(y*z-x*w);
     2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x^2+y^2)];
    res = [Rq,zeros(3,1);
            zeros(1,3),1];
end

function res = create_objective_traj_eul(trajectory,GH_seq)
    res = zeros(size(trajectory));
    for istep = 1:size(trajectory,1)
        scapula_thorax = YZX_seq(trajectory(istep,1:3)) * YZX_seq(trajectory(istep,4:6));
        humerus_thorax = scapula_thorax * YZY_seq (trajectory(istep,7:9));
        res(istep,4:6) = rotm2eul(scapula_thorax(1:3,1:3),'YZX');
        res(istep,7:9) = rotm2eul(humerus_thorax(1:3,1:3),GH_seq);
    end
    res(:,[1:3,10]) = trajectory(:,[1:3,10]);
end

function res = create_objective_traj_quat(trajectory)
    res = zeros(size(trajectory));
    for istep = 1:size(trajectory,1)
        AC_pos = Qrm(trajectory(istep,1:4)) * [1;0;0;1];
        scapula_thorax = mulQuat(trajectory(istep,1:4),trajectory(istep,5:8));
        humerus_thorax = mulQuat(scapula_thorax,trajectory(istep,(9:12)));
        res(istep,1:3) = AC_pos(1:3);
        res(istep,5:8) = scapula_thorax;
        res(istep,9:12) = humerus_thorax;
    end
    res(:,13) = trajectory(:,13);
end

% function res = hum_in_glob_quat(trajectory)
%     res = zeros(size(trajectory));
%     for istep = 1:size(trajectory,1)
%         scapula_thorax = mulQuat(trajectory(1:4),trajectory(5:8));
%         humerus_thorax = mulQuat(scapula_thorax,trajectory(9:12));
%         res(istep,4:6) = rotm2eul(scapula_thorax(1:3,1:3),"YZX");
%         res(istep,7:9) = rotm2eul(humerus_thorax(1:3,1:3),"YZY");
%     end
%     res(:,[1:3,10]) = trajectory(:,[1:3,10]);
% end


function [lengths, minusdLdq] = opensim_get_polyvalues(angles, iMus, Dofs, GH_seq)
% This function calculates the length and -dL/dq of muscle "Mus" about dof set 
% "Dofs" at a given angle matrix "angles" of opensim model "Mod"
%
% "Angles" can be a vector (one hand position) or a matrix (one hand
% position per row)
% "Dofs" is a vector of indeces
%
% Adapted from opensim_get_momarm.m 
% Dimitra Blana, March 2012
%
% 28/3/2012: Use setValue (with active constraints) instead of set state
% 1/10/2014: Simplified how the muscles and dofs are accessed
%
% 11/25/19: Update by Derek Wolf: iMus is used to access the muscle and a
% check is used to determine if the name (with no underscores) matches the
% name in Mus

import org.opensim.modeling.*
osimfile = ['das3_',GH_seq,'.osim'];
Mod = Model(osimfile);

% this is needed to get the GH lines of action in the scapular frame
SimEn = SimbodyEngine();
SimEn.connectSimbodyEngineToModel(Mod);
groundbody = Mod.getBodySet().get('thorax');
scapulabody = Mod.getBodySet().get('scapula_r');
% initialize the system to get the initial state
state = Mod.initSystem;

Mod.equilibrateMuscles(state);

% If we only want one moment arm, put it in a cell
CoordSet = Mod.getCoordinateSet();
num_request_dofs = length(Dofs);

% get the muscle
% currentMuscle = Mod.getMuscles().get(Mus);
currentMuscle = Mod.getMuscles().get(iMus-1);
% if ~strcmp(fixname(char(currentMuscle.getName())),Mus)
%     error('Current muscle name is incorrect')
% end

% angles matrix: one position per row
angles = [zeros(size(angles,1),3) angles zeros(size(angles,1),1)];
[nrows,ncols] = size(angles);
nDofs = CoordSet.getSize;

if ncols~=nDofs
    if nrows~=nDofs
        errordlg('Angle matrix not the right size','Input Error');
        minusdLdq = [];
        return;
    else
        angles = angles';
    end
end

% initialise matrices for output
lengths = zeros(size(angles,1),1);
minusdLdq = zeros(size(angles,1),num_request_dofs);
quat_J = zeros(size(angles,1),num_request_dofs);

input_vector = Dofs';
num_elements = numel(input_vector);
remainder = mod(num_elements, 3);
if remainder ~= 0
    num_to_add = 3 - remainder;
    input_vector = [input_vector, zeros(1, num_to_add)];
end

dofs_block = reshape(input_vector, 3, []).';

joints = [];
ndofsblock = size(dofs_block,1);
for i=1:ndofsblock
    currentdofs = dofs_block(i,:);
    if currentdofs == [4,5,6]
        joints = [joints 2];
    elseif currentdofs == [7,8,9]
        joints = [joints 3];
    elseif currentdofs == [10,11,12]
        joints = [joints 4];
    end
end

if nargout>3
    GHfvecs = zeros(size(angles,1),3); 
end

for istep = 1:size(angles,1)
    if ~mod(istep,50)
        disp(['Muscle ',char(currentMuscle.getName()), ' - step ',...
            num2str(istep),' of ', num2str(size(angles,1))]);
    end
    
    % set dof values for this step
    for idof = 1:nDofs
        currentDof = CoordSet.get(idof-1);    
        currentDof.setValue(state,angles(istep,idof),1);
    end

    % get GH force line of action (if the muscle crosses GH)
    if nargout>3
        muspath = currentMuscle.getGeometryPath();
        fdarray = ArrayPointForceDirection();
        fvec2 = Vec3();
        muspath.getPointForceDirections(state,fdarray);
        scap_pt_index = -1;
        % find the "effective" muscle attachment on the scapula
        for ipt=1:fdarray.getSize-1
            body1 = char(fdarray.get(ipt-1).frame); %4.0 uses frames not bodies
            body2 = char(fdarray.get(ipt).frame);
            if strcmp(body1,'scapula_r')&&strcmp(body2,'humerus_r')
                scap_pt_index=ipt;
                break;
            elseif strcmp(body2,'scapula_r')&&strcmp(body1,'humerus_r')
                scap_pt_index=ipt-1;
                break;
            end
        end
        
        if scap_pt_index==-1
            % find the "effective" muscle attachment on the clavicle
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'clavicle_r')&&strcmp(body2,'humerus_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'clavicle_r')&&strcmp(body1,'humerus_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
        end

         if scap_pt_index==-1
            % find the "effective" muscle attachment on the thorax
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'thorax')&&strcmp(body2,'humerus_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'thorax')&&strcmp(body1,'humerus_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
         end

         if scap_pt_index==-1
            % find the "effective" muscle attachment on the ulna
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'scapula_r')&&strcmp(body2,'ulna_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'scapula_r')&&strcmp(body1,'ulna_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
         end
         
         if scap_pt_index==-1
            % find the "effective" muscle attachment on the radius
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'scapula_r')&&strcmp(body2,'radius_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'scapula_r')&&strcmp(body1,'radius_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
        end

        % calculate muscle force direction at that point (in global frame)
        scap_pt = fdarray.get(scap_pt_index);
        fvec = scap_pt.direction();
        
        % transform to the scapular coordinate frame
        SimEn.transform(state,groundbody,fvec,scapulabody,fvec2);
        GHfvecs(istep,1)=fvec2.get(0);
        GHfvecs(istep,2)=fvec2.get(1);
        GHfvecs(istep,3)=fvec2.get(2);
    end    
    
    % get length
    lengths(istep) = currentMuscle.getLength(state);
    % currentMuscle.computeMomentArm(state,'SCx')

    % get moment arms
    for idof = 1:num_request_dofs
        dofindex  = Dofs(idof)-1;
        currentDof = CoordSet.get(Dofs(idof)-1);
        minusdLdq(istep,idof) = currentMuscle.computeMomentArm(state,currentDof);
    end
end

end

function [L,pmoment_arms] = momarms(musmodel, dof_indeces, angles)
% ...or use all angles
indeces = 1:size(angles,1);
sangles = angles(:,dof_indeces);

% calculate moment arms from polynomial
pmoment_arms = zeros(length(indeces),length(dof_indeces));
for iframe = 1:length(indeces)
    for i=1:musmodel.lparam_count

        % add this term's contribution to the muscle length 
        term = musmodel.lcoefs(i);

        for j=1:length(dof_indeces)
            for k=1:musmodel.lparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end

        % first derivatives of length with respect to all q's
        for  k=1:length(dof_indeces)
            % derivative with respect to q_k is zero unless exponent is 1 or higher and q is not zero
            if ((musmodel.lparams(i,k) > 0) && (sangles(iframe,k)))	
                dterm = musmodel.lparams(i,k)*term/sangles(iframe,k);
                pmoment_arms(iframe,k) = pmoment_arms(iframe,k) + dterm;
            end
        end
    end
end

L = zeros(1,length(indeces)); % Initialize the muscle length
for iframe = 1:length(indeces)
    Lterm = 0;
    for i=1:musmodel.lparam_count
        % Add this term's contribution to the muscle length
        term = musmodel.lcoefs(i);
        for j = 1:length(dof_indeces)
            for k = 1:musmodel.lparams(i, j)
                term = term * sangles(iframe,j);
            end
        end
        Lterm = Lterm + term;
    end

    L(iframe) = Lterm;
end

end

function res = min_conoid_length(traj)
    fun = @(x) conoid_length(x(2:4));
    % fun = @(x) sum((Rscap - R_scap_glob([traj(1),traj(2),x(1)],[x(2),x(3),x(4)])).^2,"all");
    nonlcon = @(x) mycon(x,traj);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = ones(1,4)*(-pi);
    ub = ones(1,4)*pi;
    x0 = traj(3:6);
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
    
    mot_eul_modified = traj;
    mot_eul_modified(3:6) = x;
    res = mot_eul_modified;
end

function res = conoid_length(AC)
    O = [0.1165,-0.0041,0.0143];
    I = T_trans([0.1575,0,0]) * YZX_seq(AC) * position([-0.0536, -0.0009, -0.0266]);
    res = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);
end

function [c,ceq] = mycon(x,traj)
    Rscap = R_scap_glob(traj(1:3),traj(4:6));
    ceq = Rscap - R_scap_glob([traj(1),traj(2),x(1)],[x(2),x(3),x(4)]);
    c = [];
end

function res = R_scap_glob(jnt1,jnt2)
    res = YZX_seq(jnt1) * YZX_seq(jnt2);
end