%% Run this first
addpath Matlab_functions/

%% Plot the computational performance and display the RMS errors
% computational_performance_eul_vs_quat
% computational_performance_RMS_angles
% %%  Plot the activations of two elements affected by gimbal lock. 
% % This generates Fig. 3.
% participant = 'par2';
% OS_model = ['Motions/',participant,'/OS_model.mat'];
% motion_name = 'Elevation';
% res_q1 = {['res_euler_',motion_name,'_50'],'euler',motion_name,'YZY',participant,'Euler angles'};
% res_q2 = {['res_quat_',motion_name,'_200'],'quat',motion_name,'YZY',participant,'Quaternions'};
% plot_IEEE_activations(OS_model,["pect_maj_c_2","delt_scap_6"],0,res_q1,res_q2)

%% This generates fig.5,6 and 7
%  - healthy SHR prediction (res_SHR_0) - this simulation is computed from
%  healthy_SHR_prediction_example.ipynb file
% - res_SHR_5 is taken from the sensitivity analysis study and correspond
% to the GH_weight of 10

% participant = 'par2';
% OS_model = ['Motions/',participant,'/OS_model.mat'];
% motion_name = 'All_elevations';
% OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
% res_q1 = {['res_SHR_Healthy_Elevation_Scabduction_wGH0'],'quat',motion_name,'YZY',participant,'Healthy Elev-Scab wGH=0'};
% res_q2 = {['res_SHR_Healthy_Elevation_Flexion_wGH0'],'quat',motion_name,'YZY',participant,'Healthy Elev-Flex wGH=0'};
% res_q3 = {['res_SHR_Healthy_Scabduction_Flexion_wGH0'],'quat',motion_name,'YZY',participant,'Healthy Scab-Flex wGH=0'};
% res_q4 = {['res_SHR_Healthy_drinking_shelf_reaching_wGH0'],'quat',motion_name,'YZY',participant,'Healthy ADLs wGH=0'};
% plot_IEEE_kinematics(OS_struct,res_q1,res_q2,res_q3,res_q4);

% participant = 'par2';
% OS_model = ['Motions/',participant,'/OS_model.mat'];
% motion_name = 'All_elevations';
% OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
% res_q1 = {['res_SHR_Healthy_Elevation_Scabduction_wGH2'],'quat',motion_name,'YZY',participant,'Healthy Elev-Scab wGH=2'};
% res_q2 = {['res_SHR_Healthy_Elevation_Flexion_wGH2'],'quat',motion_name,'YZY',participant,'Healthy Elev-Flex wGH=2'};
% res_q3 = {['res_SHR_Healthy_Scabduction_Flexion_wGH2'],'quat',motion_name,'YZY',participant,'Healthy Scab-Flex wGH=2'};
% res_q4 = {['res_SHR_Healthy_drinking_shelf_reaching_wGH2'],'quat',motion_name,'YZY',participant,'Healthy ADLs wGH=2'};
% plot_IEEE_kinematics(OS_struct,res_q1,res_q2,res_q3,res_q4);
clc
participant = 'par2';
OS_model = ['Motions/',participant,'/OS_model.mat'];
motion_name = 'All_elevations';
OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
res_q1 = {['res_SHR_5'],'quat',motion_name,'YZY',participant,'infra75'};
res_q2 = {['res_SHR_infra50_All_elevations_wGH2'],'quat',motion_name,'YZY',participant,'infra50'};
res_q3 = {['res_SHR_RClim0_Scabduction_Flexion_wGH10'],'quat',motion_name,'YZY',participant,'RClim0 Scab-Flex wGH=10'};
res_q4 = {['res_SHR_RClim0_drinking_shelf_reaching_wGH10'],'quat',motion_name,'YZY',participant,'RClim0 ADLs wGH=10'};
plot_IEEE_kinematics(OS_struct,res_q1);

% plotGHStabilityAnglesIEEE(res_q1)
% % res_q1 = {['res_SHR_0'],'quat',motion_name,'YZY',participant,'Healthy calib wGH=2'};
% % res_q2 = {['res_SHR_RClim050_calib_wgh2'],'quat',motion_name,'YZY',participant,'RClim calib wGH=2 RC 50%'};
% % plot_IEEE_kinematics(OS_struct,res_q1,res_q2)
% % plotGHStabilityAnglesIEEE(res_q1,res_q2)
% res_q1 = {['res_SHR_0'],'quat',motion_name,'YZY',participant,'Healthy calib wGH=2'};
% res_q2 = {['res_SHR_RClim050_calib_wgh4'],'quat',motion_name,'YZY',participant,'RClim calib wGH=4 RC 50%'};
% plot_IEEE_kinematics(OS_struct,res_q1,res_q2)
% plotGHStabilityAnglesIEEE(res_q1,res_q2)

% plot_EMG_healthy_RClim_IEEE(['Motions\',participant,'\',motion_name,'\EMG_',participant,'_',motion_name,'.mat'], OS_model,res_q1,res_q2)

%% This generates sensitivity analysis figures (Supplementary)
% plotGHStabilityAnglesIEEE_sanalysis()
% plot_IEEE_kinematics_sanalysis(OS_struct)
% %% This generates fig.4
%  % - res_"motion"_0 is the simulation with non-calibrated muscle parameters
%  % - res_"motion"_1 is the simulation with calibrated muscle parameters
% participant = 'par2';
% plot_GH_seq = 'YZY';
% motion_name = 'shelf_reaching';
% OS_model = ['Motions/',participant,'/OS_model.mat'];
% OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
% res_q1 = {['res_',motion_name,'_0'],'quat',motion_name,'YZY',participant,'Quat'};
% res_q2 = {['res_',motion_name,'_1'],'quat',motion_name,'YZY',participant,'Eul'};
% plot_EMG_optim_IEEE(['Motions\',participant,'\',motion_name,'\EMG_',participant,'_',motion_name,'.mat'], OS_model,res_q1,res_q2)
% 
% %% Print the RMSE between calibrated and generic model activation
% calibration_activation_RMSE



function plot_IEEE_activations(OS_model,mus_group,plot_excitation,varargin)
    alphabet = {'a','b'};
    model = load(OS_model);
    muscles = model.model.muscles;
    for i = 1:length(muscles)
        muscle_names{i} = muscles{i}.osim_name;
    end

    colors_act = [0, 0, 0, 1; 0.314, 0.784, 0.471, 1; 0.4660, 0.6740, 0.1880, 1];
    colors_exc = [0, 0, 0, 0.5; 0.314, 0.784, 0.471, 0.5;0.4660, 0.6740, 0.1880, 0.5];


    mask = startsWith(muscle_names,mus_group);
    num_in_group = nnz(mask);
    plot_rows = ceil(num_in_group/3);
    if rem(num_in_group,2) == 0
        plot_rows = plot_rows;
    end
    current_names = muscle_names(mask);
    legend_names = {};
    num_res = length(varargin);
    figure('Color','w',"Units","inches",'Position',[1 1 3.5 2])
    tiledlayout(1,2,"TileSpacing","compact","Padding","compact")

    for imus = 1:num_in_group
    nexttile; hold on; box on
        for ires = 1:num_res
        iresult = varargin{ires};
        file_name = iresult{1};
        rot_type = iresult{2};
        motion_name = iresult{3};
        GH_seq = iresult{4};
        participant = iresult{5};
        plot_name = iresult{6};
        % result = load(['Motions\',participant,'\',motion_name,'\res_',rot_type,'_',motion_name,'_',weight,'.mat']);
        result = load(['Motions\',participant,'\',motion_name,'\',file_name,'.mat']);
    
        activations = result.data.activations(:,mask);
        excitations = result.data.excitations(:,mask);
        time = result.data.tout;
        percent_of_motion = linspace(0,100,length(time));
        GL_pos_prcnt = 0;
    
        plot(percent_of_motion, activations(:,imus),'Color',colors_act(ires,:),'LineWidth',1.5); hold on
            if strcmp(rot_type,'euler')
                trajectory = result.data.trajectories;
                time = result.data.tout;
                GL_pos = find_gimbal_lock(time,trajectory(:,8));
                GL_pos_prcnt = GL_pos/time(end)*100;
            end
            if plot_excitation == 1
                plot(percent_of_motion, excitations(:,imus),'Color',colors_exc(ires,:),'LineWidth',1); hold on
            end
            axis([-inf inf -inf inf])  
            title(current_names{imus},'Interpreter','none','FontSize',10);
            % title('Lateral deltoid muscle element','FontSize',14)
            xlabel(['% of motion'])
            ylabel('Activation [-]')
            text(0.5, -0.25, ['(',alphabet{imus},')'], 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
            'FontName', 'Times New Roman', 'FontSize', 8);



            if strcmp(rot_type,'euler')
                legend_names{end+1} = plot_name;
                if plot_excitation==1
                    legend_names{end+1} = [''];
                end
            elseif strcmp(rot_type,'quat')
                legend_names{end+1} = plot_name;
                if plot_excitation==1
                    legend_names{end+1} = [''];
                end
            end

            if ~isempty(GL_pos_prcnt)
                % plot(plot_rows,2,imus)
                xline(GL_pos_prcnt,'--','LineWidth',1.0);
                for iGL = 1:length(GL_pos_prcnt)
                    legend_names{end+1} = [''];
                end
            end

        end
    end % end num_res
    legend_names{5} = ['Gimbal lock'];

    fig = gcf;
    lg = legend(legend_names,'FontSize',8); %
    lg.Box = 'off';
    lg.Layout.Tile = 'south';
    lg.Orientation = "horizontal";
    lg.ItemTokenSize = 20;

    % exportgraphics(gcf,'IEEE_GL_activation.png','Resolution',600);

end

function plot_IEEE_kinematics_sanalysis(kinematics)
weights = {'res_SHR_0','res_SHR_1','res_SHR_2','res_SHR_3','res_SHR_4','res_SHR_5'}; %,'200143'
wGHs = {'2','2','4','6','8','10'};
colors = {"blue","green","cyan","magenta","black","red"};
motion_name = 'All_motions';
legend_names = {};

figure('Color','w','Units','normalized','Position',[0.05 0.05 0.7 0.75]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

for i = 1:9
nexttile; hold on; box on;
for iweight = 1:length(weights)
RClim_struct   = load(['Motions\par2\All_elevations\',weights{iweight},'.mat']);

t            = RClim_struct.data.tout;

kin_exp      = kinematics.mot_struct.euler;
kin_exp = interp1(kinematics.mot_struct.time,kin_exp,t,"spline");
kin_exp = create_objective_traj_eul(kin_exp,'YZY',0);
% t = linspace(0,length(t),length(t));

kin_RC_loc     = quat2eul_motion(RClim_struct.data.trajectories,'YZY');
kin_RC  = create_objective_traj_eul(kin_RC_loc,'YZY',1);



% ---------- Labels ----------
labels = { ...
    'Clavicle protraction/retraction', 'Clavicle elevation','Clavicle axial rotation', ...
    'Scapula internal/external rotation','Scapula upward/downward rotation','Scapula anterior/posterior tilting', ...
    'Humerus plane of elevation','Humerus elevation','Humerus axial rotation'};

% ---------- Figure setup ----------

% 
% % ---------- Line styles ----------
lw = 1.6;
RC_style      = {'-','Color',colors{iweight},'LineWidth',lw};
RC_interp_style = {'--','Color',colors{iweight},'LineWidth',0.8};

    

    if i == 7 || i == 9

        kin_RC_interp = fillmissing(kin_RC,'linear');
        plot(t, rad2deg(kin_RC_interp(:,i)), RC_interp_style{:}); hold on
        plot(t, rad2deg(kin_RC(:,i)), RC_style{:}); hold on
    else
        plot(t, rad2deg(kin_RC(:,i)), RC_style{:});hold on
    end
    


    title(labels{i});
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    legend_names{end+1} = '';
    if iweight == 1
        legend_names{end+1} = ['$w_{GH}^{Healthy} = $',wGHs{iweight}];
    else
        legend_names{end+1} = ['$w_{GH}^{RC-lim} = $',wGHs{iweight}];
    end
    

    xlim([t(1) t(end)]);
    


end



end
set(gca,'FontSize',9,'LineWidth',0.8);
% "predicted" annotation
annotation('line', [0.96 0.96], [0.43 0.97], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);
annotation('line', [0.95 0.96], [0.97 0.97], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);
annotation('line', [0.95 0.96], [0.43 0.43], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);

annotation('textbox', [0.96 0.67 0.0 0.0], ...
    'String', 'Predicted', ...
    'EdgeColor', 'none', ...
    'Rotation', 90, ...
    'FontSize', 15, ...
    'FontAngle', 'italic', ...
    'HorizontalAlignment', 'center');


% "tracked" annotation
annotation('line', [0.96 0.96], [0.12 0.35], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);
annotation('line', [0.95 0.96], [0.12 0.12], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);
annotation('line', [0.95 0.96], [0.35 0.35], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);

annotation('textbox', [0.96 0.235 0.0 0.0], ...
    'String', 'Tracked', ...
    'EdgeColor', 'none', ...
    'Rotation', 90, ...
    'FontSize', 15, ...
    'FontAngle', 'italic', ...
    'HorizontalAlignment', 'center');


% ---------- Legend ----------
lg = legend(legend_names,'Interpreter','latex'); %, ...
             % 'Orientation','horizontal', ...
             % 'Location','southoutside');
lg.FontSize = 9;
lg.Box = 'off';
lg.Layout.Tile = 'south';
lg.Orientation = 'horizontal';
lg.FontSize = 9;
lg.Box = 'off';

% exportgraphics(gcf,'IEEE_kinematics_RClim_sanalyis.png','Resolution',600);
    
end

function plotGHStabilityAnglesIEEE_sanalysis()

weights = {'res_SHR_0','res_SHR_1','res_SHR_2','res_SHR_3','res_SHR_4','res_SHR_5'}; %,'200143'
wGHs = {'2','2','4','6','8','10'};
colors = {"blue","green","cyan","magenta","black","red"};

% --- Figure ---
figure('Color','w');
tiledlayout(2,length(weights),'TileSpacing','compact','Padding','compact');

for iweight = 1:length(weights)
res   = load(['Motions\par2\All_elevations\',weights{iweight},'.mat']);

Rx_r = res.data.reactions(:,1);
Ry_r = res.data.reactions(:,2);
Rz_r = res.data.reactions(:,3);
t            = res.data.tout;
t_norm = (t - t(1)) / (t(end) - t(1));
R_r_rot = zeros(size(res.data.reactions));
for i = 1:size(Rx_r,1)
    R_r_rot_i = R_y(13*pi/180)' * R_z(-6.5*pi/180)' * [Rx_r(i);Ry_r(i);Rz_r(i);1];
    R_r_rot(i,:) = R_r_rot_i(1:3)';
    GH_force_r(i) = norm(R_r_rot_i(1:3));
end

AP_lim = 14.84;   % anterior/posterior
SI_lim = 23.74;   % superior/inferior

theta_AP_r = -atan2d(R_r_rot(:,3), -R_r_rot(:,1));
theta_SI_r = atan2d(R_r_rot(:,2), -R_r_rot(:,1));

phi = linspace(0,2*pi,300);
ellipse_x = AP_lim * cos(phi);
ellipse_y = SI_lim * sin(phi);

cmap = parula(256);

% ========= Panel B: RC-limited =========
nexttile;
hold on; box off; axis equal;

plot(ellipse_x, ellipse_y, 'k','LineWidth',1.2);

scatter(theta_AP_r, theta_SI_r, 25, t, 'filled');

if iweight == 1
xlabel('Ant-Pos angle (deg)');
ylabel('Sup-Inf angle (deg)');
end

if iweight == 1
    title(['$w^{Healthy}_{GH} = $',wGHs{iweight}],'Interpreter','latex');
else
    title(['$w^{RC-lim}_{GH} = $',wGHs{iweight}],'Interpreter','latex');
end

xlim([-AP_lim-2 AP_lim+2]);
ylim([-SI_lim-2 SI_lim+2]);

set(gca,'FontSize',9,'LineWidth',0.8);
% colormap(cmap);
% cb = colorbar;

end
colormap(cmap);
cb = colorbar;

% ========= Panel C: Compression force =========
nexttile([1 length(weights)]);
hold on; box off;
legend_names = {};
for iweight = 1:length(weights)
res   = load(['Motions\par2\All_elevations\',weights{iweight},'.mat']);

Rx_r = res.data.reactions(:,1);
Ry_r = res.data.reactions(:,2);
Rz_r = res.data.reactions(:,3);
t            = res.data.tout;
t_norm = (t - t(1)) / (t(end) - t(1));
R_r_rot = zeros(size(res.data.reactions));
for i = 1:size(Rx_r,1)
    R_r_rot_i = R_y(13*pi/180)' * R_z(-6.5*pi/180)' * [Rx_r(i);Ry_r(i);Rz_r(i);1];
    R_r_rot(i,:) = R_r_rot_i(1:3)';
    GH_force_r(i) = norm(R_r_rot_i(1:3));
end
plot(t, GH_force_r/650*100,'Color',colors{iweight}, 'LineWidth',1.5);
xlim([-inf inf])

xlabel('Time (s)');
ylabel('GH force (%BW)');
if iweight == 1
    legend_names{end+1} = ['$w^{Healthy}_{GH} = $',wGHs{iweight}];
else
    legend_names{end+1} = ['$w^{RC-lim}_{GH} = $',wGHs{iweight}];
end

end
lg = legend(legend_names,'Location','northeast','Box','off','Interpreter','latex');
lg.Layout.Tile = 'south';
lg.Orientation = 'horizontal';
lg.FontSize = 9;
lg.Box = 'off';

hold off
set(gca,'FontSize',9,'LineWidth',0.8);

% ========= Colorbar =========
cb.Layout.Tile = 'east';
cb.Label.String = 'Time (s)';
cb.FontSize = 9;

% set(gcf, 'Units', 'inches', ' Position', [1 1 3.5 3])
% exportgraphics(gcf,'IEEE_GH_stability_sensitivity_analysis.png','Resolution',600);


end

function computational_performance_RMS_angles()
motions = {'Elevation', 'Scabduction', 'Flexion'};
participants = {'par1','par2','par3'};
rot_type = {'euler','quat'};
rot_type_weights = {'50','200'};
qSC_IK = [];
qSC_eul = [];
qSC_quat = [];
qAC_IK = [];
qAC_eul = [];
qAC_quat = [];
qGH_IK = [];
qGH_eul = [];
qGH_quat = [];
for ipar = 1:length(participants)
    for imot = 1:length(motions)
        for irot = 1:2
            struct_path = ['Motions\',participants{ipar},'\',motions{imot},'\res_',rot_type{irot},'_',motions{imot},'_',rot_type_weights{irot},'.mat'];
            load(struct_path);
            OS_struct = load(['Motions/',participants{ipar},'/',motions{imot},'/',motions{imot},'.mat']);
            IK = OS_struct.mot_struct.euler;
            IK_glob = create_objective_traj_eul(IK,'YZY',0);
            time = data.tout;
            time_new = linspace(0,time(end),100);
            if irot == 1
                qSC_IK = [qSC_IK; reshape(IK_glob(:,1:2),[200,1])];
                qAC_IK = [qAC_IK; reshape(IK_glob(:,4:6),[300,1])];
                qGH_IK = [qGH_IK; reshape(IK_glob(:,7:9),[300,1])];
            end
            if strcmp(rot_type{irot},'euler')
                trajectory = data.trajectories;
                trajectory = interp1(time,trajectory,time_new);
                trajectory_glob = create_objective_traj_eul(trajectory,'YZY',0);
                qSC_eul = [qSC_eul; reshape(trajectory_glob(:,1:2),[200,1])];
                qAC_eul = [qAC_eul; reshape(trajectory_glob(:,4:6),[300,1])];
                qGH_eul = [qGH_eul; reshape(trajectory_glob(:,7:9),[300,1])];

            elseif strcmp(rot_type{irot},'quat')
                trajectory = data.trajectories;
                trajectory = interp1(time,trajectory,time_new);
                trajectory = quat2eul_motion(trajectory,'YZY');
                
                trajectory_glob = create_objective_traj_eul(trajectory,'YZY',0);
                qSC_quat = [qSC_quat; reshape(trajectory_glob(:,1:2),[200,1])];
                qAC_quat = [qAC_quat; reshape(trajectory_glob(:,4:6),[300,1])];
                qGH_quat = [qGH_quat; reshape(trajectory_glob(:,7:9),[300,1])];
            end
        end
    end
end

fprintf('RMSE SC euler-angles: %.2f \n',rmse(rad2deg(qSC_eul),rad2deg(qSC_IK)))
fprintf('RMSE SC quaternions: %.2f \n',rmse(rad2deg(qSC_quat),rad2deg(qSC_IK)))

fprintf('RMSE AC euler-angles: %.2f \n',rmse(rad2deg(qAC_eul),rad2deg(qAC_IK)))
fprintf('RMSE AC quaternions: %.2f \n',rmse(rad2deg(qAC_quat),rad2deg(qAC_IK)))

fprintf('RMSE GH euler-angles: %.2f \n',rmse(rad2deg(qGH_eul),rad2deg(qGH_IK)))
fprintf('RMSE GH quaternions: %.2f \n',rmse(rad2deg(qGH_quat),rad2deg(qGH_IK)))
end

function calibration_activation_RMSE()
% motions = {'shelf_reaching','lifting_5kg','driving','drinking'};
motions = {'driving','lifting_5kg'};
calibrations = {'0','Elevation_Scabduction','Elevation_Flexion','Scabduction_Flexion','drinking_shelf_reaching'};
participant = 'par2';
OS_model = 'Motions/par2/OS_model_prediction.mat';
EMG_muscles = {'Infrasp','UpperTrap','Serrupper','IntermediateDelt','PosteriorDelt','AnteriorDelt','MiddleTrap'};
calib_file_names = {'infra','trapclav','serr','deltscap','deltscap','deltclav','trapscap'};
model_names = {["infra_3","infra_4","infra_5"],["trap_clav_1"],["serr_ant_2","serr_ant_3","serr_ant_4"],["delt_scap11","delt_scap10","delt_scap_9","delt_scap_8"],["delt_scap_3","delt_scap_4","delt_scap_5"],["delt_clav_1","delt_clav_2"],["trap_scap_6","trap_scap_7","trap_scap_8"]};


for icalib = 1:numel(calibrations)

for imus = 1:length(EMG_muscles)

model = load(OS_model);
muscles = model.model.muscles;
num_muscles = length(muscles);
results = {calibrations{icalib}};
activation_healthy = [];
activation_RClim = [];
activation_EMG = [];


for i = 1:num_muscles
    muscle_names{i} = muscles{i}.osim_name;
end
mus_index = {};
for igroup = 1:length(model_names)
    current_group = model_names{igroup};
    group_indeces = [];
        for imus_in_group = 1:length(current_group)
            group_indeces = [group_indeces,find(strcmp(muscle_names,current_group(imus_in_group)))];
        end
    mus_index{end+1} = int16(group_indeces);
end

for imot = 1:length(motions)

    
     for ires = 1:1
        emg_data = load(['Motions/par2/',motions{imot},'/EMG_',participant,'_',motions{imot},'.mat']);
        result = load(['Motions/par2/',motions{imot},'/res_',motions{imot},'_',results{ires},'.mat']);
        if icalib>1
            calibration_file = load(['Motions/par2/',calibrations{icalib},'/',calibrations{icalib},'_calibrated_params.mat']);
        end
        current_index = mus_index{imus};
        activations = zeros(size(result.data.activations(:,1)));
        excitations = zeros(size(result.data.excitations(:,1)));
        for ielement = 1:length(current_index)
            activations = activations+result.data.activations(:,current_index(ielement));
            excitations = excitations+result.data.excitations(:,current_index(ielement));
        end
        activations = excitations/length(current_index);
        excitations = excitations/length(current_index);
        time = linspace(0,100,length(result.data.tout));

        if ires == 1
            rsmpl_simulation = linspace(0,100,length(excitations));
            current_emg_rsmpld = zeros(length(excitations),6);
            for icase = 1:6
                try
                    current_emg = emg_data.data.(['num_',num2str(icase)]).(EMG_muscles{imus});
                end
                time_emg = linspace(0,100,length(current_emg));
                current_emg_rsmpld = spline(time_emg,current_emg,rsmpl_simulation);
            end

            activation_EMG = [activation_EMG;current_emg_rsmpld'];
            activation_healthy = [activation_healthy; excitations];
        
        else
            activation_RClim = [activation_RClim;excitations];
        end

        

     end
end

% if icalib == 1
%     fprintf('RMSE original, %s = %f6.2\n',EMG_muscles{imus},(rmse(activation_healthy,activation_EMG)))
% end

if icalib>1
    scaler_lceopt = calibration_file.calibrated_params.(['lceopt_scaler_',calib_file_names{imus}]);
    scaler_fmax = calibration_file.calibrated_params.(['fmax_scaler_',calib_file_names{imus}]);
    fprintf('(calib: %s), %s (slce = %0.3f, sfmax = %0.3f); RMSE = %0.3f\n',calibrations{icalib}, EMG_muscles{imus},scaler_lceopt,scaler_fmax,(rmse(activation_healthy,activation_EMG)))
else
    fprintf('(calib: %s), %s; RMSE = %0.3f\n',calibrations{icalib}, EMG_muscles{imus},(rmse(activation_healthy,activation_EMG)))

end
end

end


end

function plotGHStabilityAnglesIEEE(healthys, RClims)
alphabet = {'a','b','c','d','e','f','g','h','i','i'};

healthy = load(['Motions\',healthys{5},'\',healthys{3},'\',healthys{1},'.mat']);
RClim   = load(['Motions\',RClims{5},'\',RClims{3},'\',RClims{1},'.mat']);
Rx_h = healthy.data.reactions(:,1);
Ry_h = healthy.data.reactions(:,2);
Rz_h = healthy.data.reactions(:,3);
Rx_r = RClim.data.reactions(:,1);
Ry_r = RClim.data.reactions(:,2);
Rz_r = RClim.data.reactions(:,3);
t            = healthy.data.tout;
t_norm = (t - t(1)) / (t(end) - t(1));
R_h_rot = zeros(size(healthy.data.reactions));
R_r_rot = zeros(size(healthy.data.reactions));
for i = 1:size(Rx_h,1)
    R_h_rot_i = R_y(14*pi/180)' * R_z(-6.5*pi/180)' * [Rx_h(i);Ry_h(i);Rz_h(i);1];
    R_r_rot_i = R_y(14*pi/180)' * R_z(-6.5*pi/180)' * [Rx_r(i);Ry_r(i);Rz_r(i);1];
    R_h_rot(i,:) = R_h_rot_i(1:3)';
    R_r_rot(i,:) = R_r_rot_i(1:3)';
    GH_force_h(i) = norm(R_h_rot_i(1:3));
    GH_force_r(i) = norm(R_r_rot_i(1:3));
end

% --- Anatomical limits (degrees) ---
AP_lim = 14.84;   % anterior/posterior
SI_lim = 23.74;   % superior/inferior

% --- Compute angles (degrees) ---
theta_AP_h = -atan2d(R_h_rot(:,3), -R_h_rot(:,1));
theta_SI_h = atan2d(R_h_rot(:,2), -R_h_rot(:,1));

theta_AP_r = -atan2d(R_r_rot(:,3), -R_r_rot(:,1));
theta_SI_r = atan2d(R_r_rot(:,2), -R_r_rot(:,1));

% --- Ellipse for glenoid boundary ---
phi = linspace(0,2*pi,300);
ellipse_x = AP_lim * cos(phi);
ellipse_y = SI_lim * sin(phi);

% --- Figure ---
figure('Color','w','Units','inches','Position',[1 1 3.5 3]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

cmap = parula(256);

% ========= Panel A: Healthy =========
nexttile;
hold on; box on; axis equal;

plot(ellipse_x, ellipse_y, 'k','LineWidth',1.2);
ax = gca;
ax.Box = 'off';
scatter(theta_AP_h, theta_SI_h, 2, t, 'filled');

xlabel(['Ant-Pos angle (deg)',newline,'']);
ylabel('Sup-Inf angle (deg)');
title('Healthy');

xlim([-AP_lim-2 AP_lim+2]);
ylim([-SI_lim-2 SI_lim+2]);
text(0.5, -0.35, ['(',alphabet{1},')'], 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 8);

set(gca,'FontSize',8,'LineWidth',0.2);
colormap(cmap);

% ========= Panel B: RC-limited =========
nexttile;
hold on; box on; axis equal;

plot(ellipse_x, ellipse_y, 'k','LineWidth',1.2);
ax = gca;
ax.Box = 'off';
scatter(theta_AP_r, theta_SI_r, 2, t, 'filled');

xlabel(['Ant-Pos angle (deg)',newline,'']);
ylabel('Sup-Inf angle (deg)');
title(['RC-limited']);

xlim([-AP_lim-2 AP_lim+2]);
ylim([-SI_lim-2 SI_lim+2]);
text(0.5, -0.35, ['(',alphabet{2},')'], 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 8);

set(gca,'FontSize',8,'LineWidth',0.2);
colormap(cmap);
cb = colorbar;

% ========= Panel C: Compression force =========
nexttile([1 2]);
hold on; box off;

plot(t, GH_force_h/650*100, 'LineWidth',1.5,'Color','blue');
plot(t, GH_force_r/650*100, 'LineWidth',1.5,'Color','red');
xlim([-inf inf])
ylim([0 max(GH_force_r/650*100)+15])

xlabel(['Time (s)',newline,'']);
ylabel('GH force (%BW)');
legend({'Healthy','RC-limited'},'Position',[0.52,0.38,1,1],'Box','off');
text(0.5, -0.35, ['(',alphabet{3},')'], 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 8);

set(gca,'FontSize',8,'LineWidth',0.2);

% ========= Colorbar =========
cb.Layout.Tile = 'east';
cb.Label.String = 'Time (s)';
cb.FontSize = 8;

% exportgraphics(gcf,'IEEE_GH_stability2.png','Resolution',600);


end

function plot_IEEE_kinematics(kinematics, varargin)

% for isim = 1:length(varargin)
% healthys = varargin{isim};
% healthy_struct = load(['Motions\',healthys{5},'\',healthys{3},'\',healthys{1},'.mat']);
% % RClim_struct   = load(['Motions\',RClims{5},'\',RClims{3},'\',RClims{1},'.mat']);
% % healthy_struct
% name_healthy = healthys{6};
% % name_RClim = RClims{6};
% t            = healthy_struct.data.tout;
% 
% kin_exp      = kinematics.mot_struct.euler;
% kin_exp = interp1(kinematics.mot_struct.time,kin_exp,t,"spline");
% kin_exp = create_objective_traj_eul(kin_exp,'YZY',0);
% % t = linspace(0,length(t),length(t));
% 
% kin_healthy_loc  = quat2eul_motion(healthy_struct.data.trajectories,'YZY');
% kin_healthy = create_objective_traj_eul(kin_healthy_loc,'YZY',1);
% % kin_RC_loc     = quat2eul_motion(RClim_struct.data.trajectories,'YZY');
% % kin_RC  = create_objective_traj_eul(kin_RC_loc,'YZY',1);
% GH_healthy = kin_healthy_loc(:,8);
% % GH_RClim = kin_RC_loc(:,8);
% 
% 
% fro = [1,41];
% sca = [101,135];
% sag = [204, 236];
% [SCHR_frontal,info_frontal] = compute_SCHR(kin_healthy(fro(1):fro(2),8)*180/pi,kin_healthy(fro(1):fro(2),5)*180/pi,GH_healthy(fro(1):fro(2))*180/pi);
% [SCHR_scapular,info_scapular] = compute_SCHR(kin_healthy(sca(1):sca(2),8)*180/pi,kin_healthy(sca(1):sca(2),5)*180/pi,GH_healthy(sca(1):sca(2))*180/pi);
% [SCHR_sagittal,info_sagittal] = compute_SCHR(kin_healthy(sag(1):sag(2),8)*180/pi,kin_healthy(sag(1):sag(2),5)*180/pi,GH_healthy(sag(1):sag(2))*180/pi);
% % SCHR_frontal
% SCHR_scapular
% SCHR_sagittal
% 
% [SCHR_frontal_RClim,info_frontal_RClim] = compute_SCHR(kin_RC(fro(1):fro(2),8)*180/pi,kin_RC(fro(1):fro(2),5)*180/pi,GH_RClim(fro(1):fro(2))*180/pi);
% [SCHR_scapular_RClim,info_scapular_RClim] = compute_SCHR(kin_RC(sca(1):sca(2),8)*180/pi,kin_RC(sca(1):sca(2),5)*180/pi,GH_RClim(sca(1):sca(2))*180/pi);
% [SCHR_sagittal_RClim,info_sagittal_RClim] = compute_SCHR(kin_RC(sag(1):sag(2),8)*180/pi,kin_RC(sag(1):sag(2),5)*180/pi,GH_RClim(sag(1):sag(2))*180/pi);
% SCHR_frontal_RClim
% SCHR_scapular_RClim
% SCHR_sagittal_RClim
% print_SCHR(['Frontal / '  name_healthy], SCHR_frontal,info_frontal);
% print_SCHR(['Scapular / ' name_healthy], SCHR_scapular,info_scapular);
% print_SCHR(['Sagittal / ' name_healthy], SCHR_sagittal,info_sagittal);

% print_SCHR(['Frontal / ',name_RClim], SCHR_frontal_RClim,info_frontal_RClim);
% print_SCHR(['Scapular / ',name_RClim], SCHR_scapular_RClim,info_scapular_RClim);
% print_SCHR(['Sagittal / ',name_RClim], SCHR_sagittal_RClim,info_sagittal_RClim);

% ---------- Labels ----------
labels = { ...
    'Clavicle protraction/retraction', 'Clavicle elevation','Clavicle axial rotation', ...
    'Scapula internal/external rotation','Scapula upward/downward rotation','Scapula anterior/posterior tilting', ...
    'Humerus plane of elevation','Humerus elevation','Humerus axial rotation'};

% ---------- Figure setup ----------
figure('Color','w','Units','inches','Position',[1 1 7.16 4.3]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
% 
% % ---------- Line styles ----------
lw = 1.6;
exp_style     = {'Color',[0.5 0.5 0.5],'LineWidth',1.2};                 % IK (feasible)
healthy_style = {'-','Color','blue','LineWidth',lw};
healthy_interp_style = {'--','Color','blue','LineWidth',0.8};
RC_style      = {'-','Color','red','LineWidth',lw};
RC_interp_style = {'--','Color','red','LineWidth',0.8};
alphabet = {'a','b','c','d','e','f','g','h','i','i'};
names = {'Experimental trajectory'};
for i = 1:9
    nexttile; hold on; box on;

    for isim = 1:length(varargin)
    healthys = varargin{isim};
    healthy_struct = load(['Motions\',healthys{5},'\',healthys{3},'\',healthys{1},'.mat']);
    % RClim_struct   = load(['Motions\',RClims{5},'\',RClims{3},'\',RClims{1},'.mat']);
    % healthy_struct
    name_healthy = healthys{6};
    if i == 1
        names{end+1} = name_healthy;
    end
    % name_RClim = RClims{6};
    t            = healthy_struct.data.tout;
    
    kin_exp      = kinematics.mot_struct.euler;
    kin_exp = interp1(kinematics.mot_struct.time,kin_exp,t,"spline");
    kin_exp = create_objective_traj_eul(kin_exp,'YZY',0);
    % t = linspace(0,length(t),length(t));
    
    kin_healthy_loc  = quat2eul_motion(healthy_struct.data.trajectories,'YZY');
    kin_healthy = create_objective_traj_eul(kin_healthy_loc,'YZY',1);
    % kin_RC_loc     = quat2eul_motion(RClim_struct.data.trajectories,'YZY');
    % kin_RC  = create_objective_traj_eul(kin_RC_loc,'YZY',1);
    GH_healthy = rotyzy(healthy_struct.data.trajectories(:,9:12));
    % GH_RClim = kin_RC_loc(:,8);
    
    if i == 1
        fro = [1,41];
        sca = [101,135];
        sag = [204, 236];
        [SCHR_frontal,info_frontal] = compute_SCHR(kin_healthy(fro(1):fro(2),8)*180/pi,kin_healthy(fro(1):fro(2),5)*180/pi,GH_healthy(fro(1):fro(2))*180/pi);
        [SCHR_scapular,info_scapular] = compute_SCHR(kin_healthy(sca(1):sca(2),8)*180/pi,kin_healthy(sca(1):sca(2),5)*180/pi,GH_healthy(sca(1):sca(2))*180/pi);
        [SCHR_sagittal,info_sagittal] = compute_SCHR(kin_healthy(sag(1):sag(2),8)*180/pi,kin_healthy(sag(1):sag(2),5)*180/pi,GH_healthy(sag(1):sag(2))*180/pi);
        print_SCHR(['Frontal / '  name_healthy], SCHR_frontal,info_frontal);
        print_SCHR(['Scapular / ' name_healthy], SCHR_scapular,info_scapular);
        print_SCHR(['Sagittal / ' name_healthy], SCHR_sagittal,info_sagittal);
    end

    if i ~= 3
        plot(t, rad2deg(kin_exp(:,i)), exp_style{:}); hold on
    end

    if i == 7 || i == 9
        kin_healthy_interp = fillmissing(kin_healthy,'linear');
        plot(t, rad2deg(kin_healthy_interp(:,i))); hold on
        plot(t, rad2deg(kin_healthy(:,i))); hold on

        % kin_RC_interp = fillmissing(kin_RC,'linear');
        % plot(t, rad2deg(kin_RC_interp(:,i)), RC_interp_style{:}); hold on
        % plot(t, rad2deg(kin_RC(:,i)), RC_style{:});
    else
        plot(t, rad2deg(kin_healthy(:,i)));
        % plot(t, rad2deg(kin_RC(:,i)), RC_style{:});
    end
    

    title(labels{i});
    xlabel(['Time (s)',newline,'']);
    ylabel('Angle (deg)');

    xlim([t(1) t(end)]);
    set(gca,'FontSize',8,'LineWidth',0.2);
    text(0.5, -0.55, ['(',alphabet{i},')'], 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 8);
    end
end

% "predicted" annotation
annotation('line', [0.96 0.96], [0.41 0.97], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);
annotation('line', [0.95 0.96], [0.97 0.97], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);
annotation('line', [0.95 0.96], [0.41 0.41], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);

annotation('textbox', [0.96 0.67 0.0 0.0], ...
    'String', 'Predicted', ...
    'EdgeColor', 'none', ...
    'Rotation', 90, ...
    'FontSize', 15, ...
    'FontAngle', 'italic', ...
    'HorizontalAlignment', 'center');


% "tracked" annotation
annotation('line', [0.96 0.96], [0.10 0.35], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);
annotation('line', [0.95 0.96], [0.10 0.10], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);
annotation('line', [0.95 0.96], [0.35 0.35], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1);

annotation('textbox', [0.96 0.235 0.0 0.0], ...
    'String', 'Tracked', ...
    'EdgeColor', 'none', ...
    'Rotation', 90, ...
    'FontSize', 15, ...
    'FontAngle', 'italic', ...
    'HorizontalAlignment', 'center');
% names
% ---------- Legend ----------
lg = legend(names,'NumColumns',3); %, ...
             % 'Orientation','horizontal', ...
             % 'Location','southoutside');
lg.FontSize = 8;
lg.Box = 'off';
lg.Layout.Tile = 'south';
lg.Orientation = 'horizontal';
set(gca, 'LineWidth',0.2)

figure
plot(t,rad2deg(GH_healthy)); hold on
xline(t(sca(1))); hold on
xline(t(sag(1)))
% plot(t,rad2deg(kin_healthy_loc(:,7)));

% exportgraphics(gcf,'zoufa5_new.png','Resolution',600);
    
end

function print_SCHR(name, SCHR, info)
fprintf('\n%s\n', name);
fprintf('%-11s %4s %8s %6s %9s %8s %7s %7s %7s %7s %7s\n', ...
    'phase','n','slope','R2','endpoint','range','dTH','dTS','dGH','revTS','maxTS');
for i = 1:numel(SCHR)
    fprintf('%4.0f-%-6.0f %4d %8.2f %6.3f %9.2f %8.2f %7.1f %7.1f %7.1f %7.2f %7.2f\n', ...
        info.phase(i,1), info.phase(i,2), info.n(i), info.slope(i), ...
        info.R2(i), info.endpoint(i), info.range(i), ...
        info.netTH(i), info.netTS(i), info.netGH(i), info.fracRevTS(i), info.maxTS(i));
end
end

function elev=rotyzy(quat)
% calculates the Euler angles around the y,z, and new y axes
% from the rotation matrix R
num_data = size(quat,1);
elev = zeros(1,num_data);
for i = 1:num_data
r = quat2rotm(quat(i,:));
z1 = acos(r(2,2));
if (z1==0)
    y=acos(r(1,1));
	z=z1;
	ya=0.0;
	return;
end
sy = r(3,2)/sin(z1);
cy = -r(1,2)/sin(z1);
y1 = atan2(sy,cy);
sya = r(2,3)/sin(z1);
cya = r(2,1)/sin(z1);
ya1 = atan2(sya,cya);
z2 = -z1;
sy = r(3,2)/sin(z2);
cy = -r(1,2)/sin(z2);
y2 = atan2(sy,cy);
sya = r(2,3)/sin(z2);
cya = r(2,1)/sin(z2);
ya2 = atan2(sya,cya);
if (0 <= z1 && z1 <= pi)
   y = y1;
   z = z1;
   ya = ya1;
else
   y = y2;
   z = z2;
   ya = ya2;
end
% if z<(10*pi/180)
%     ya=y+ya;
%     y=0;
% end
elev(i) = z;

end
end

% function [SCHR] = compute_SCHR(TH, TS, GH)
% 
% TH = TH(:);
% TS = TS(:);
% 
% validIdx = ~isnan(TH) & ~isnan(TS);
% TH = TH(validIdx);
% TS = TS(validIdx);
% 
% % Define phase limits
% phases = [ ...
%     min(TH) 30; 
%     30 60; 
%     60 90;
%     min(TH) 90];
% 
% nPhases = size(phases,1);
% 
% SCHR = zeros(nPhases,1);
% phaseData = struct();
% 
% for i = 1:nPhases
% 
%     lower = phases(i,1);
%     upper = phases(i,2);
% 
%     idx = TH >= lower & TH <= upper;
% 
%     if sum(idx) < 2
%         SCHR(i) = NaN;
%         continue;
%     end
% 
%     deltaGH = max(GH(idx)) - min(GH(idx));
%     deltaTS = max(TS(idx)) - min(TS(idx));
% 
%     % Avoid division by zero
%     if abs(deltaTS) < 1e-6
%         SCHR(i) = NaN;
%     else
%         SCHR(i) = round(deltaGH / deltaTS,2);
%     end
% 
% end
% 
% end
function [SCHR, info] = compute_SCHR(TH, TS, GH, phaseEdges)
%COMPUTE_SCHR  Phase-specific scapulohumeral rhythm from simulated kinematics.
%
%   SCHR = COMPUTE_SCHR(TH, TS, GH) returns a 4x1 vector of SHR values for
%   the phases [min(TH)-30, 30-60, 60-90, min(TH)-90], estimated as the
%   least-squares slope of GH against TS within each phase window.
%
%   [SCHR, INFO] = COMPUTE_SCHR(...) also returns a struct of per-phase
%   diagnostics, including alternative estimators for comparison.
%
%   SCHR = COMPUTE_SCHR(TH, TS, GH, PHASEEDGES) uses custom phase limits,
%   given as an Nx2 matrix of [lower upper] bounds on TH. NaN in the lower
%   bound is replaced by min(TH).
%
%   INPUTS (all in degrees, equal length, same simulation window)
%     TH  thoracohumeral elevation
%     TS  thoracoscapular upward rotation
%     GH  glenohumeral elevation, taken directly from the relative
%         humerus-scapula rotation (NOT TH - TS)
%
%   INFO FIELDS (one row per phase)
%     .slope        least-squares dGH/dTS                        [primary]
%     .R2           coefficient of determination of that fit
%     .endpoint     (GH(end)-GH(1)) / (TS(end)-TS(1))            [check]
%     .range        (max-min GH) / (max-min TS)                  [legacy]
%     .netTH/.netGH/.netTS   signed change over the phase
%     .spanTS       max(TS)-min(TS), the conditioning of the slope
%     .fracRevTS    fraction of samples where TS moves opposite to its net
%                   direction (flags the scapular setting phase)
%     .n            samples in the phase
%     .phase        [lower upper] bounds used
%
%   NOTES
%     The slope estimator is preferred over the range ratio because it uses
%     all samples, preserves sign, and is unaffected by non-monotonic TS
%     (e.g. transient downward rotation during scapular setting), which
%     inflates a max-min denominator and biases SHR downward.

% ---------------------------------------------------------------- inputs
TH = TH(:); TS = TS(:); GH = GH(:);

if ~isequal(numel(TH), numel(TS), numel(GH))
    error('compute_SCHR:sizeMismatch', ...
        'TH (%d), TS (%d) and GH (%d) must be the same length.', ...
        numel(TH), numel(TS), numel(GH));
end

% Filter jointly so the three signals cannot become misaligned. There
% should be no NaN if all three come from the same simulation window, so
% warn rather than fail silently.
bad = isnan(TH) | isnan(TS) | isnan(GH);
if any(bad)
    warning('compute_SCHR:nanFound', ...
        '%d of %d samples contain NaN and were removed from all signals.', ...
        nnz(bad), numel(bad));
    TH = TH(~bad); TS = TS(~bad); GH = GH(~bad);
end

if nargin < 4 || isempty(phaseEdges)
    phaseEdges = [NaN 30; 30 60; 60 90; NaN 90];
end
phaseEdges(isnan(phaseEdges(:,1)), 1) = min(TH);

nPhases   = size(phaseEdges, 1);
minSpanTS = 1e-3;   % deg; below this the slope is ill-conditioned
minSamp   = 3;      % samples needed for a meaningful regression

% --------------------------------------------------------------- outputs
SCHR = nan(nPhases, 1);
z    = nan(nPhases, 1);
info = struct('phase', phaseEdges, 'n', z, 'slope', z, 'R2', z, ...
              'endpoint', z, 'range', z, 'netTH', z, 'netTS', z, ...
              'netGH', z, 'spanTS', z, 'fracRevTS', z, 'maxTS', z);

% ----------------------------------------------------------------- loop
for i = 1:nPhases

    idx = TH >= phaseEdges(i,1) & TH <= phaseEdges(i,2);
    n   = nnz(idx);
    info.n(i) = n;

    if n < minSamp
        continue
    end

    ts = TS(idx);  gh = GH(idx);  th = TH(idx);

    info.netTH(i)  = th(end) - th(1);
    info.netTS(i)  = ts(end) - ts(1);
    info.netGH(i)  = gh(end) - gh(1);
    info.spanTS(i) = max(ts) - min(ts);
    info.maxTS(i) = ts(end);

    % How much of the phase moves against the net TS direction. A large
    % value means max-min is not measuring the net rotation.
    dts = diff(ts);
    if info.netTS(i) ~= 0 && ~isempty(dts)
        info.fracRevTS(i) = mean(sign(dts) == -sign(info.netTS(i)));
    end

    if info.spanTS(i) < minSpanTS
        continue    % scapula effectively stationary: SHR undefined
    end

    % Primary estimator: least-squares slope of GH on TS.
    A = [ts, ones(n,1)];
    p = A \ gh;
    info.slope(i) = p(1);

    ghFit    = A * p;
    SSres    = sum((gh - ghFit).^2);
    SStot    = sum((gh - mean(gh)).^2);
    info.R2(i) = 1 - SSres / max(SStot, eps);

    % Secondary estimators, for the sensitivity table.
    if abs(info.netTS(i)) > minSpanTS
        info.endpoint(i) = info.netGH(i) / info.netTS(i);
    end
    info.range(i) = (max(gh) - min(gh)) / info.spanTS(i);

    SCHR(i) = info.slope(i);
end

SCHR = round(SCHR, 2);

end

function plot_EMG_healthy_RClim_IEEE(EMG_struct, OS_model, varargin)

figure('Color','w','Units','inches','Position',[1 1 3.5 3.5]);

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% ---------- Line styles ----------
lw = 1.2;
EMG_style       = {'Color',[0.5 0.5 0.5],'LineWidth',0.8};            % black solid
results_style  = {{'-','Color','blue','LineWidth',lw},{'-','Color','red','LineWidth',lw}};
alphabet = {'a','b','c','d'};
EMG_muscles = {'IntermediateDelt','Infrasp','UpperTrap','Serrupper'};
model_names = {["delt_scap_8","delt_scap_9","delt_scap_10","delt_scap11"],["infra_3","infra_4"],["trap_clav_1"],["serr_ant_2","serr_ant_3","serr_ant_4"]};
figure_names = {'Lateral Deltoid','Infraspinatus',' Clavicular trapezius','Serratus anterior'};
num_res = length(varargin);
emg_data = load(EMG_struct);
model = load(OS_model);
muscles = model.model.muscles;
num_muscles = length(muscles);

for i = 1:num_muscles
    muscle_names{i} = muscles{i}.osim_name;
end
mus_index = {};
for igroup = 1:length(model_names)
    current_group = model_names{igroup};
    group_indeces = [];
        for imus_in_group = 1:length(current_group)
            group_indeces = [group_indeces,find(strcmp(muscle_names,current_group(imus_in_group)))];
        end
    mus_index{end+1} = int16(group_indeces);
end

for i = 1:length(EMG_muscles)
    nexttile; hold on; box on;
    legend_names = {};
    activation_rmse = [];
     for ires = 1:num_res
        iresult = varargin{ires};
        file_name = iresult{1};
        motion_name = iresult{3};
        participant = iresult{5};
        plot_name = iresult{6};
        result = load(['Motions\',participant,'\',motion_name,'\',file_name,'.mat']);
        % result = load(['Motions\',participant,'\All_motions\res_',rot_type,'_',motion_name,'.mat']);
        current_index = mus_index{i};
        activations = zeros(size(result.data.activations(:,1)));
        excitations = zeros(size(result.data.excitations(:,1)));
        for ielement = 1:length(current_index)
            activations = activations+result.data.activations(:,current_index(ielement));
            excitations = excitations+result.data.excitations(:,current_index(ielement));
        end
        activations = activations/length(current_index);
        excitations = excitations/length(current_index);
        tout = result.data.tout;
        time = linspace(0,tout(end),length(result.data.tout));
        activation_rmse = [activation_rmse, activations];

        legend_names{end+1} = plot_name;
        legend_names{end+1} = '';

        if ires == 1
            rsmpl_simulation = linspace(0,100,length(activations));
            current_emg_rsmpld = zeros(length(activations),6);
            for imot = 1:6
                try
                    current_emg = emg_data.data.(['num_',num2str(imot)]).(EMG_muscles{i});
                end
                time_emg = linspace(0,tout(end),length(current_emg));
                % current_emg_rsmpld(:,imot) = spline(time_emg,current_emg,rsmpl_simulation);
            end
        
            plot(time_emg, current_emg, EMG_style{:}); hold on

            % Get current y-limits (so shading spans full plot height)
            yl = [0 1];
            
            % Logical vector where EMG is NaN
            isNaN = isnan(current_emg);
            % size(isNaN)
            
            % Find start and end indices of NaN regions
            d = diff([false; isNaN'; false]);
            nanStart = find(d == 1);
            nanEnd   = find(d == -1) - 1;
            
            % Add shaded patches
            hold on
            for k = 1:length(nanStart)
                xPatch = [time_emg(nanStart(k)) time_emg(nanEnd(k)) time_emg(nanEnd(k)) time_emg(nanStart(k))];
                yPatch = [yl(1) yl(1) yl(2) yl(2)];

                xregion(xPatch(1),xPatch(2), 'FaceColor', [0.75 0.75 0.75], ...
                    'FaceAlpha', 0.15, ...
                    'EdgeColor', 'none');
            end
            % hold off
            hold on

        end
        cur_res_style = results_style{ires};
        % ylim([0 max([max(activations),max(current_emg_rsmpld)])])
        plot(time,excitations,cur_res_style{:}); hold on %,'Linestyle',line_styles{ires}
        

     end

     
    title(figure_names{i});
    xlabel(['Time (s)',newline,'']);
    ylabel('Excitation (s)');
    text(0.5, -0.35, ['(',alphabet{i},')'], 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 8);

    % ylim([0 0.5]);
    % xlim([t(1) t(end)]);

    set(gca,'FontSize',8,'LineWidth',0.2);

end
lg = legend({'EMG','EMG not defined','','Healthy','RC-limited'});
lg.Layout.Tile = 'south';
lg.Orientation = 'horizontal';
lg.ItemTokenSize = 10;
lg.FontSize = 8;
lg.Box = 'off';

% exportgraphics(gcf,'IEEE_emg_healthy_RClim2.png','Resolution',600);

end

function plot_EMG_optim_IEEE(EMG_struct, OS_model, varargin)

figure('Color','w','Units','inches','Position',[1 1 3.5 3.5]);
alphabet = {'a','b','c','d'};

tt=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% ---------- Line styles ----------
lw = 1.2;
EMG_style       = {'Color',[0.5 0.5 0.5],'LineWidth',0.8};            % black solid
results_style  = {{'-','Color',[245 190 40]/255,'LineWidth',lw},{'-','Color',[166 20 146]/255,'LineWidth',lw}};

EMG_muscles = {'IntermediateDelt','Infrasp','UpperTrap','Serrupper'};
model_names = {["delt_scap_8","delt_scap_9","delt_scap_10","delt_scap11"],["infra_3","infra_4","infra_5"],["trap_clav_1"],["serr_ant_2","serr_ant_3","serr_ant_4"]};
figure_names = {'Lateral Deltoid','Infraspinatus',' Clavicular trapezius','Serratus anterior'};
num_res = length(varargin);
emg_data = load(EMG_struct);
model = load(OS_model);
muscles = model.model.muscles;
num_muscles = length(muscles);

for i = 1:num_muscles
    muscle_names{i} = muscles{i}.osim_name;
end
mus_index = {};
for igroup = 1:length(model_names)
    current_group = model_names{igroup};
    group_indeces = [];
        for imus_in_group = 1:length(current_group)
            group_indeces = [group_indeces,find(strcmp(muscle_names,current_group(imus_in_group)))];
        end
    mus_index{end+1} = int16(group_indeces);
end

for i = 1:length(EMG_muscles)
    nexttile; hold on; box on;
    legend_names = {};
    activation_rmse = [];
     for ires = 1:num_res
        iresult = varargin{ires};
        file_name = iresult{1};
        motion_name = iresult{3};
        participant = iresult{5};
        plot_name = iresult{6};
        result = load(['Motions\',participant,'\',motion_name,'\',file_name,'.mat']);
        % result = load(['Motions\',participant,'\All_motions\res_',rot_type,'_',motion_name,'.mat']);
        current_index = mus_index{i};
        activations = zeros(size(result.data.activations(:,1)));
        excitations = zeros(size(result.data.excitations(:,1)));
        for ielement = 1:length(current_index)
            activations = activations+result.data.activations(:,current_index(ielement));
            excitations = excitations+result.data.excitations(:,current_index(ielement));
        end
        activations = excitations/length(current_index);
        excitations = excitations/length(current_index);
        tout = result.data.tout;
        time = linspace(0,tout(end),length(result.data.tout));
        activation_rmse = [activation_rmse, activations];

        legend_names{end+1} = plot_name;
        legend_names{end+1} = '';

        % plot(time,activations,'Color',line_colors_act(ires,:),'LineWidth',1.5) %,'Linestyle',line_styles{ires}
        % hold on
        % results_style{ires,:}
        if ires == 1
            rsmpl_simulation = linspace(0,100,length(activations));
            current_emg_rsmpld = zeros(length(activations),6);
            for imot = 1:6
                try
                    current_emg = emg_data.data.(['num_',num2str(imot)]).(EMG_muscles{i});
                end
                time_emg = linspace(0,tout(end),length(current_emg));
                % current_emg_rsmpld(:,imot) = spline(time_emg,current_emg,rsmpl_simulation);
            end
        
            plot(time_emg, current_emg, EMG_style{:}); hold on

            % Get current y-limits (so shading spans full plot height)
            yl = [0 1];
            
            % Logical vector where EMG is NaN
            isNaN = isnan(current_emg);
            % size(isNaN)
            
            % Find start and end indices of NaN regions
            d = diff([false; isNaN'; false]);
            nanStart = find(d == 1);
            nanEnd   = find(d == -1) - 1;
            
            % Add shaded patches
            hold on
            for k = 1:length(nanStart)
                xPatch = [time_emg(nanStart(k)) time_emg(nanEnd(k)) time_emg(nanEnd(k)) time_emg(nanStart(k))];
                yPatch = [yl(1) yl(1) yl(2) yl(2)];

                patch(xPatch, yPatch, [0.8 0.8 0.8], ...
                    'FaceAlpha', 0.15, ...
                    'EdgeColor', 'none');
            end
            % hold off
            hold on

        end
        xlim([0 tout(end)])

        cur_res_style = results_style{ires};
        plot(time,excitations,cur_res_style{:}); hold on %,'Linestyle',line_styles{ires}
        

     end

     
    title(figure_names{i});
    xlabel(['Time (s)',newline,'']);
    ylabel('Excitation (-)');
    text(0.5, -0.35, ['(',alphabet{i},')'], 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 8);

    ylim([0 0.35]);
    % xlim([t(1) t(end)]);

    set(gca,'FontSize',8,'LineWidth',0.2);
    % if i == 3
        
    % end

end

lg = legend({'EMG','Original parameters','Adjusted parameters'});
lg.Layout.Tile = 'south';
lg.Orientation = 'horizontal';
lg.FontSize = 8;
lg.Box = 'off';
lg.ItemTokenSize = 10;

% exportgraphics(gcf,'IEEE_emg_orig_adjusted.png','Resolution',600);

end


function computational_performance_eul_vs_quat()
motions = {'Elevation', 'Scabduction', 'Flexion'};
participants = {'par1','par2','par3'};
rot_type = {'euler','quat'};
rot_type_weights = {'50','200'};
itersE = [];
itersE_GL = [];
itersQ = [];
timeE = [];
timeE_GL = [];
timeQ = [];
for ipar = 1:length(participants)
    for imot = 1:length(motions)
        for irot = 1:2
            struct_path = ['Motions\',participants{ipar},'\',motions{imot},'\res_',rot_type{irot},'_',motions{imot},'_',rot_type_weights{irot},'.mat'];
            load(struct_path);
            num_iter = length(data.objective_value);
            time2solve = data.time2sol;
            if strcmp(rot_type{irot},'euler')
                trajectory = data.trajectories;
                time = data.tout;
                GL_pos = find_gimbal_lock(time,trajectory(:,8));
                if isempty(GL_pos)
                    itersE = [itersE;num_iter];
                    timeE = [timeE;time2solve];
                else
                    itersE_GL = [itersE_GL;num_iter];
                    timeE_GL = [timeE_GL;time2solve];
                end
            elseif strcmp(rot_type{irot},'quat')
                itersQ = [itersQ;num_iter];
                timeQ = [timeQ;time2solve];
            end
        end
    end
end
% itersE
plot_computational_performance(itersE,itersQ,itersE_GL,...
                                        timeE,timeQ,timeE_GL)
end

function plot_computational_performance(itersEo,itersQ,itersE_GL,...
                                        timeE,timeQ,timeE_GL)

figure('Color','w','Units','inches','Position',[1 1 3.5 2.5]);
tiledlayout(1,2,'TileSpacing','compact','Padding','loose')

% ===================== ITERATIONS =====================
nexttile; hold on;
itersE = [itersEo;itersE_GL];
timeE = [timeE;timeE_GL];
n = length(itersE);

% Draw connecting lines first (background)
for i = 1:n
    plot([1 2],[itersE(i) itersQ(i)],...
        'Color',[0.8 0.8 0.8],'LineWidth',1);
end

% Scatter points
scatter(ones(n,1), itersE, 50, 'k','filled');
scatter(2*ones(n,1), itersQ, 50, 'k','filled');

% Highlight GL cases in Euler
% scatter(ones(length(itersE_GL),1), itersE_GL,...
%     70,'r','x','LineWidth',1.5);
scatter(ones(length(itersE_GL),1), itersE_GL, ...
    70,'r','x','LineWidth',1.5);
ax = gca;
ax.TickLabelInterpreter = 'tex';
ax.YAxis.Exponent = 0;

set(gca,'YScale','log');
xlim([0.7 2.3])
xticks([1 2])
xticklabels({'Euler\newline angles','Quaternions'})
ylim([200,2500])
yticks([200 500 1000 2000])
yticklabels({'200','500','1000','2000'})
ylabel('Iterations')
t = title('Solver Iterations');

ax = gca;
ax.Position(2) = ax.Position(2) + 0.05;

set(gca,'FontSize',8,'Box','off','YMinorTick','off')


% ===================== TIME =====================
nexttile; hold on;

n = length(timeE);

% Connecting lines
for i = 1:n
    plot([1 2],[timeE(i) timeQ(i)],...
        'Color',[0.8 0.8 0.8],'LineWidth',1);
end

% Scatter
scatter(ones(n,1), timeE, 50, 'k','filled');
scatter(2*ones(n,1), timeQ, 50, 'k','filled');

% Highlight GL
hGL = scatter(ones(length(timeE_GL),1), timeE_GL,...
    70,'r','x','LineWidth',1.5);

set(gca,'YScale','log');

xlim([0.7 2.3])
xticks([1 2])
ylim([400 9000])
yticks([500 1000 2000 4000 8000])
yticklabels({'500','1000','2000','4000','8000'})


xticklabels({'Euler\newline angles','Quaternions'})
ylabel('Time to solve [s]')
t = title('Solver Time');

legend(hGL,['Gimbal lock',newline, 'occurrence'],'Position',[0.4 0.45 0.9 0.5],'Box','off')

set(gca,'FontSize',8,'Box','off','YMinorTick','on')

ax = gca;
ax.Position(2) = ax.Position(2) + 0.5;
sgtitle('Computational Performance','FontWeight','bold','FontSize',10)

% exportgraphics(gcf,'comp_perf.png','Resolution',600)

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

function res = dquatdt(quat,w)
    res = 1/2 * G(quat)' * w';
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

function res = create_objective_traj_eul(trajectory,GH_seq,GL_zone)
    res = zeros(size(trajectory));
    for istep = 1:size(trajectory,1)
        scapula_thorax = YZX_seq(trajectory(istep,1:3)) * YZX_seq(trajectory(istep,4:6));
        if strcmp(GH_seq,'YZY')
            humerus_thorax = scapula_thorax * YZY_seq (trajectory(istep,7:9));
        elseif strcmp(GH_seq,'YZX')
            humerus_thorax = scapula_thorax * YZX_seq (trajectory(istep,7:9));
        end
        res(istep,4:6) = rotm2eul(scapula_thorax(1:3,1:3),'YZX');
        if GL_zone == 0
            res(istep,7:9) = rotm2eul(humerus_thorax(1:3,1:3),GH_seq);
        else
            res(istep,7:9) = rotm2yzy_shoulder(humerus_thorax(1:3,1:3));
        end
    end
    res(:,[1:3,10]) = trajectory(:,[1:3,10]);
end


function res = rotm2yzy_shoulder(R)

    cos_z = min(max(R(2,2), -1), 1);

    z_mag = acosd(cos_z); 

    if z_mag < 7.5
        y1 = nan;
        y2 = nan;

        sin_z = (R(2,1) - R(1,2)) / 2;
        z = atan2d(sin_z, cos_z);

    else
        sin_z_normal = sqrt(R(2,1)^2 + R(2,3)^2); 

        z = atan2d(sin_z_normal, R(2,2));
        y1 = atan2d(R(3,2), -R(1,2));
        y2 = atan2d(R(2,3), R(2,1));
    end
    res = [y1,z,y2]*pi/180;
end

function res = create_objective_traj_quat(trajectory,OS_model)
    res = zeros(size(trajectory));
    model = load(OS_model);
    AC_offset = model.model.joints{1,5}.location;
    for istep = 1:size(trajectory,1)
        AC_pos = Qrm(trajectory(istep,1:4)) * [AC_offset(1);AC_offset(2);AC_offset(3);1];
        scapula_thorax = mulQuat(trajectory(istep,1:4),trajectory(istep,5:8));
        humerus_thorax = mulQuat(scapula_thorax,trajectory(istep,(9:12)));
        res(istep,1:3) = AC_pos(1:3);
        res(istep,5:8) = scapula_thorax;
        res(istep,9:12) = humerus_thorax;
    end
    res(:,13) = trajectory(:,13);
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

function results = computational_performance_RMS_angles2(cfg)
%COMPUTATIONAL_PERFORMANCE_RMS_ANGLES  Euler vs quaternion tracking accuracy.
%
%   Reports three complementary views of tracking error:
%
%   (1) Per-DOF Euler-component error, with mean +/- SD across trials, 95%
%       CI on the paired difference, max absolute error, bias and % of ROM.
%       Answers R1.4.
%
%   (2) The same, restricted to samples above a thoracohumeral elevation
%       threshold. Sweeping the threshold shows that the quaternion penalty
%       in the first and third YZY axes is a conversion artefact, not a
%       tracking failure. Answers R1.5.
%
%   (3) Per-segment GEODESIC orientation error: the rotation angle of
%       R_sim * R_ref'. Sequence-free and singularity-free, so it scores
%       both formulations on the same physical quantity. Euler component
%       error cannot do this: near the YZY singularity the first and third
%       angles are individually undetermined, and the Euler simulation is
%       additionally being scored on its own objective function.
%
%   Robust statistics (median, 95th percentile) and the elevation at which
%   the maximum occurs are reported alongside RMSE, because the large
%   scapular maxima are isolated spikes rather than distributed error.

% ------------------------------------------------------------------ config
if nargin < 1, cfg = struct(); end
def = struct( ...
    'motions',      {{'Elevation','Scabduction','Flexion'}}, ...
    'participants', {{'par1','par2','par3'}}, ...
    'rotType',      {{'euler','quat'}}, ...
    'rotWeights',   {{'50','200'}}, ...
    'root',         'Motions', ...
    'nSamp',        100, ...
    'seq',          'YZY', ...
    'globFlag',     0, ...
    'elevCol',      8, ...
    'elevThresh',   0, ...          % deg; 0 = use all samples
    'dofCols',      [1 2 4 5 6 7 8 9], ...
    'dofNames',     {{'Clav protraction','Clav elevation', ...
                      'Scap protraction','Scap upward rot','Scap tilt', ...
                      'TH plane of elev','TH elevation','TH axial rot'}}, ...
    'segCols',      {{1:3, 4:6, 7:9}}, ...
    'segNames',     {{'Thoracoclavicular','Thoracoscapular','Thoracohumeral'}}, ...
    'segSeq',       {{'YZY','YZY','YZY'}}, ...
    'intrinsic',    true);          % see eulToRotm and the convention check
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k}) || isempty(cfg.(fn{k})), cfg.(fn{k}) = def.(fn{k}); end
end

nPar = numel(cfg.participants);
nMot = numel(cfg.motions);
nDOF = numel(cfg.dofCols);
nSeg = numel(cfg.segCols);
nTr  = nPar * nMot;

% ------------------------------------------------------------- preallocate
Z  = nan(nTr, nDOF);
Zs = nan(nTr, nSeg);
blankD = struct('rmse',Z,'maxabs',Z,'p95',Z,'med',Z,'bias',Z,'rom',Z, ...
                'nrmse',Z,'elevAtMax',Z);
blankG = struct('rmse',Zs,'maxabs',Zs,'p95',Zs,'med',Zs,'elevAtMax',Zs);

results.trial.euler = blankD;   results.trial.quat = blankD;
results.geo.euler   = blankG;   results.geo.quat   = blankG;
results.trial.label = cell(nTr,1);
results.trial.fracKept = nan(nTr,1);

iTr = 0;

% ------------------------------------------------------------------- loop
for ipar = 1:nPar
    for imot = 1:nMot

        iTr = iTr + 1;
        results.trial.label{iTr} = sprintf('%s / %s', ...
            cfg.participants{ipar}, cfg.motions{imot});

        OS = load(fullfile(cfg.root, cfg.participants{ipar}, ...
                           cfg.motions{imot}, [cfg.motions{imot} '.mat']));
        IK_glob = create_objective_traj_eul(OS.mot_struct.euler, cfg.seq, cfg.globFlag);
        assert(size(IK_glob,1) == cfg.nSamp, ...
            'IK reference has %d rows, expected %d.', size(IK_glob,1), cfg.nSamp);

        elev = rad2deg(IK_glob(:, cfg.elevCol));
        keep = elev >= cfg.elevThresh;
        results.trial.fracKept(iTr) = mean(keep);

        refDeg = rad2deg(IK_glob(:, cfg.dofCols));

        for irot = 1:numel(cfg.rotType)

            S = load(fullfile(cfg.root, cfg.participants{ipar}, cfg.motions{imot}, ...
                 ['res_' cfg.rotType{irot} '_' cfg.motions{imot} '_' cfg.rotWeights{irot} '.mat']));

            t    = S.data.tout;
            traj = interp1(t, S.data.trajectories, linspace(0, t(end), cfg.nSamp));

            if strcmp(cfg.rotType{irot}, 'quat')
                traj = quat2eul_motion(traj, cfg.seq);
            end
            simGlob = create_objective_traj_eul(traj, cfg.seq, cfg.globFlag);
            simDeg  = rad2deg(simGlob(:, cfg.dofCols));

            f = cfg.rotType{irot};

            % --- (1)-(2) Euler-component error -----------------------------
            err = simDeg(keep,:) - refDeg(keep,:);
            results.trial.(f) = putDOF(results.trial.(f), iTr, err, ...
                                       refDeg(keep,:), elev(keep));

            % --- (3) geodesic orientation error ----------------------------
            g = nan(sum(keep), nSeg);
            for iseg = 1:nSeg
                g(:,iseg) = geodesicError( ...
                    simGlob(keep, cfg.segCols{iseg}), ...
                    IK_glob(keep, cfg.segCols{iseg}), ...
                    cfg.segSeq{iseg}, cfg.intrinsic);
            end
            results.geo.(f) = putGeo(results.geo.(f), iTr, g, elev(keep));
        end
    end
end

% -------------------------------------------------------------- aggregate
results.summary.euler = aggregate(results.trial.euler);
results.summary.quat  = aggregate(results.trial.quat);
results.geoSum.euler  = aggregate(results.geo.euler);
results.geoSum.quat   = aggregate(results.geo.quat);

results.paired.dof = aggVec(results.trial.quat.rmse - results.trial.euler.rmse);
results.paired.geo = aggVec(results.geo.quat.rmse   - results.geo.euler.rmse);

% Excess error attributable to the quaternion formulation. Equal values in
% the first and third YZY axes indicate coordinate degeneracy rather than
% tracking error, since a perturbation near the singularity splits equally
% and oppositely between them.
results.excess.dof = sqrt(max(results.summary.quat.rmse_mean.^2 - ...
                              results.summary.euler.rmse_mean.^2, 0));
results.excess.geo = sqrt(max(results.geoSum.quat.rmse_mean.^2 - ...
                              results.geoSum.euler.rmse_mean.^2, 0));

results.cfg = cfg;

% ------------------------------------------------------------------ report
hdr = 'ALL SAMPLES';
if cfg.elevThresh > 0
    hdr = sprintf('EXCLUDING TH ELEVATION < %g deg', cfg.elevThresh);
end

printGeo(hdr, cfg.segNames, results.geoSum, results.paired.geo, results.excess.geo);
printDOF(cfg.dofNames, results.summary, results.paired.dof, results.excess.dof);

fprintf('\nSamples retained: %.1f%% (range %.1f-%.1f%% across trials)\n', ...
    100*mean(results.trial.fracKept), ...
    100*min(results.trial.fracKept), 100*max(results.trial.fracKept));
fprintf('n = %d trials (%d participants x %d motions)\n', nTr, nPar, nMot);
fprintf(['\nCONVENTION CHECK: geodesic error for the EULER run should be of\n' ...
         'the same order as its component errors (~1-3 deg). A value an order\n' ...
         'of magnitude larger means eulToRotm does not match the convention\n' ...
         'used by create_objective_traj_eul; toggle cfg.intrinsic and retry.\n\n']);

end