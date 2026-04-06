addpath Matlab_functions/
%%  Plot the activations of two elements affected by gimbal lock. 
OS_model = ['Motions/',participant,'/OS_model.mat'];
participant = 'par2';
motion_name = 'Elevation_yzy';
res_q1 = {['res_euler_',motion_name,'_50'],'euler',motion_name,'YZY',participant,'Euler angles'};
res_q2 = {['res_quat_',motion_name,'_200'],'quat',motion_name,'YZY',participant,'Quaternions'};
plot_IEEE_activations(OS_model,["pect_maj_c_2","delt_scap_6"],0,res_q1,res_q2)

%%
participant = 'par2';
motion_name = 'All_motions';
res_q1 = {['res_SHR_0'],'quat',motion_name,'YZY',participant,'Healthy'};
res_q2 = {['res_SHR_5'],'quat',motion_name,'YZY',participant,'RC-limited'};

plot_GH_seq = 'YZY';
% mus_group = ["delt_scap_5","delt_scap_6","delt_scap_7","delt_scap_8","subscap_5"];
mus_group = ["delt"];
OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
OS_model = ['Motions/',participant,'/OS_model.mat'];
% plot_paper_kinematics(OS_struct,OS_model, plot_GH_seq,0,res_q1,res_q2)
% plot_activations_EMG(['Motions\',participant,'\',motion_name,'\EMG_',participant,'_',motion_name,'.mat'], OS_model,'one',1,res_q1)
% plot_paper_activations(OS_model,mus_group,1,res_q1,res_q2)
% plot_conoid_length(OS_model,res_q2)
% plot_polynomials(OS_model,mus_group,1.0,res_q1)
% plot_GH_reactions(res_q2)
% examine_fvectors(res_q2,OS_model)
% plot_EMG_optim_IEEE(['Motions\',participant,'\',motion_name,'\EMG_',participant,'_',motion_name,'.mat'], OS_model,'one',1,res_q1,res_q2)
% plot_EMG_healthy_RClim_IEEE(['Motions\',participant,'\',motion_name,'\EMG_',participant,'_',motion_name,'_fig.mat'], OS_model,'one',1,res_q1,res_q2)

plot_IEEE_kinematics(OS_struct,res_q1,res_q2)
% plotGHStabilityAnglesIEEE(res_q1,res_q2)
% plotGHStabilityAnglesIEEE_sanalysis()
% plot_IEEE_kinematics_sanalysis(OS_struct)
% plot_activation_snalaysis()
% activation_RMSE
% computational_performance_eul_vs_quat
% computational_performance_RMS_angles

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

            % if ~isempty(GL_pos_prcnt)
            %     for GH_num = 1:length(GL_pos_prcnt)
            %         if GH_num == length(GL_pos_prcnt) & imus == num_in_group
            %             legend_names{end+1} = ['GL'];
            %         else
            %             legend_names{end+1} = [''];
            %         end
            %     end
            % end
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


function plot_activation_snalaysis()
clc
participant = '3NA';
motion_name = 'All_motions';
mus_group = ["subscap"];
model_names = {"subscap"};
weights = {'198141','198142','198143','198144','198145'};

OS_model = ['Motions/',participant,'/OS_model.mat'];
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
            group_indeces = [group_indeces,find(contains(muscle_names,current_group(imus_in_group)))];
        end
    mus_index{end+1} = int16(group_indeces);
end

for i = 1:length(mus_group)
    nexttile; hold on; box on;
    legend_names = {};
     for iweight = 1:length(weights)
        result = load(['Motions\',participant,'\',motion_name,'\res_quat_',motion_name,'_',weights{iweight},'.mat']);
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

        % legend_names{end+1} = plot_name;
        % legend_names{end+1} = '';

        plot(time,excitations); hold on %,'Linestyle',line_styles{ires}

     end


    title(mus_group{i});
    xlabel(['Time (s)',newline,'']);
    ylabel('Excitation (s)');
    % text(0.5, -0.35, ['(',alphabet{i},')'], 'Units', 'normalized', ...
    % 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    % 'FontName', 'Times New Roman', 'FontSize', 8);

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


end

function plot_IEEE_kinematics_sanalysis(kinematics)
clc
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
RClim_struct   = load(['Motions\3NA\All_motions\',weights{iweight},'.mat']);

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

% set(gcf, 'Units', 'inches', 'Position', [1 1 7.16 3])
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
res   = load(['Motions\3NA\All_motions\',weights{iweight},'.mat']);

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

% --- Anatomical limits (degrees) ---
AP_lim = 14.84;   % anterior/posterior
SI_lim = 23.74;   % superior/inferior

% --- Compute angles (degrees) ---
theta_AP_r = -atan2d(R_r_rot(:,3), -R_r_rot(:,1));
theta_SI_r = atan2d(R_r_rot(:,2), -R_r_rot(:,1));

% --- Ellipse for glenoid boundary ---
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
res   = load(['Motions\3NA\All_motions\',weights{iweight},'.mat']);

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
exportgraphics(gcf,'IEEE_GH_stability_sensitivity_analysis.png','Resolution',600);


end

function computational_performance_RMS_angles()
motions = {'Elevation_yzy', 'Scabduction_yzy', 'Flexion_yzy'};
participants = {'1LU','3NA','8tata'};
rot_type = {'euler','quat'};
rot_type_weights = {'50','200'};
itersE = [];
itersE_GL = [];
itersQ = [];
timeE = [];
timeE_GL = [];
timeQ = [];
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

rmse(rad2deg(qSC_eul),rad2deg(qSC_IK))
rmse(rad2deg(qSC_quat),rad2deg(qSC_IK))

rmse(rad2deg(qAC_eul),rad2deg(qAC_IK))
rmse(rad2deg(qAC_quat),rad2deg(qAC_IK))

rmse(rad2deg(qGH_eul),rad2deg(qGH_IK))
rmse(rad2deg(qGH_quat),rad2deg(qGH_IK))
end

function activation_RMSE()
motions = {'shelf_reaching','lifting_2kg','driving','drinking'};
participant = '3NA';
OS_model = 'Motions/3NA/OS_model_prediction.mat';
EMG_muscles = {'Infrasp','UpperTrap','Serrupper','IntermediateDelt','PosteriorDelt'};
model_names = {["infra_3","infra_4","infra_5"],["trap_clav_1"],["serr_ant_2","serr_ant_3","serr_ant_4"],["delt_scap11","delt_scap10","delt_scap_9","delt_scap_8"],["delt_scap_3","delt_scap_4","delt_scap_5"]};


for imus = 1:length(EMG_muscles)

model = load(OS_model);
muscles = model.model.muscles;
num_muscles = length(muscles);
results = {'0','1'};
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

    
     for ires = 1:length(results)
        emg_data = load(['Motions/3NA/',motions{imot},'/EMG_',participant,'_',motions{imot},'.mat']);
        result = load(['Motions/3NA/',motions{imot},'/res_',motions{imot},'_',results{ires},'.mat']);
        current_index = mus_index{imus};
        activations = zeros(size(result.data.activations(:,1)));
        excitations = zeros(size(result.data.excitations(:,1)));
        for ielement = 1:length(current_index)
            activations = activations+result.data.activations(:,current_index(ielement));
            excitations = excitations+result.data.excitations(:,current_index(ielement));
        end
        activations = activations/length(current_index);
        % excitations = excitations/length(current_index);
        time = linspace(0,100,length(result.data.tout));

        if ires == 1
            rsmpl_simulation = linspace(0,100,length(activations));
            current_emg_rsmpld = zeros(length(activations),6);
            for icase = 1:6
                try
                    current_emg = emg_data.data.(['num_',num2str(icase)]).(EMG_muscles{imus});
                end
                time_emg = linspace(0,100,length(current_emg));
                current_emg_rsmpld = spline(time_emg,current_emg,rsmpl_simulation);
            end

            activation_EMG = [activation_EMG;current_emg_rsmpld'];
            activation_healthy = [activation_healthy; activations];
        
        else
            activation_RClim = [activation_RClim;activations];
        end

        

     end
end

disp(['RMSE original ',EMG_muscles{imus},' = ',num2str(rmse(activation_healthy,activation_EMG))])
disp(['RMSE adjusted ',EMG_muscles{imus},' = ',num2str(rmse(activation_RClim,activation_EMG))])
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


% set(gcf, 'Units', 'inches', ' Position', [1 1 3.5 3])
exportgraphics(gcf,'IEEE_GH_stability2.png','Resolution',600);


end

function plot_IEEE_kinematics(kinematics, healthys, RClims)
clc
motion_name = 'All_motions';
healthy_struct = load(['Motions\',healthys{5},'\',healthys{3},'\',healthys{1},'.mat']);
RClim_struct   = load(['Motions\',RClims{5},'\',RClims{3},'\',RClims{1},'.mat']);
t            = healthy_struct.data.tout;

kin_exp      = kinematics.mot_struct.euler;
kin_exp = interp1(kinematics.mot_struct.time,kin_exp,t,"spline");
kin_exp = create_objective_traj_eul(kin_exp,'YZY',0);
% t = linspace(0,length(t),length(t));

kin_healthy_loc  = quat2eul_motion(healthy_struct.data.trajectories,'YZY');
kin_healthy = create_objective_traj_eul(kin_healthy_loc,'YZY',1);
kin_RC_loc     = quat2eul_motion(RClim_struct.data.trajectories,'YZY');
kin_RC  = create_objective_traj_eul(kin_RC_loc,'YZY',1);
GH_healthy = kin_healthy(:,8) - kin_healthy(:,5);
GH_RClim = kin_RC(:,8) - kin_RC(:,5);


fro = [1,41];
sca = [101,135];
sag = [204, 236];
[SCHR_frontal,phaseData] = compute_SCHR(kin_healthy(fro(1):fro(2),8)*180/pi,kin_healthy(fro(1):fro(2),5)*180/pi,GH_healthy(fro(1):fro(2))*180/pi);
[SCHR_scapular,phaseData] = compute_SCHR(kin_healthy(sca(1):sca(2),8)*180/pi,kin_healthy(sca(1):sca(2),5)*180/pi,GH_healthy(sca(1):sca(2))*180/pi);
[SCHR_sagittal,phaseData] = compute_SCHR(kin_healthy(sag(1):sag(2),8)*180/pi,kin_healthy(sag(1):sag(2),5)*180/pi,GH_healthy(sag(1):sag(2))*180/pi);
SCHR_frontal
SCHR_scapular
SCHR_sagittal
% 
[SCHR_frontal_RClim,phaseData] = compute_SCHR(kin_RC(fro(1):fro(2),8)*180/pi,kin_RC(fro(1):fro(2),5)*180/pi,GH_RClim(fro(1):fro(2))*180/pi);
[SCHR_scapular_RClim,phaseData] = compute_SCHR(kin_RC(sca(1):sca(2),8)*180/pi,kin_RC(sca(1):sca(2),5)*180/pi,GH_RClim(sca(1):sca(2))*180/pi);
[SCHR_sagittal_RClim,phaseData] = compute_SCHR(kin_RC(sag(1):sag(2),8)*180/pi,kin_RC(sag(1):sag(2),5)*180/pi,GH_RClim(sag(1):sag(2))*180/pi);
SCHR_frontal_RClim
SCHR_scapular_RClim
SCHR_sagittal_RClim



% ---------- Labels ----------
labels = { ...
    'Clavicle protraction/retraction', 'Clavicle elevation','Clavicle axial rotation', ...
    'Scapula internal/external rotation','Scapula upward/downward rotation','Scapula anterior/posterior tilting', ...
    'Humerus plane of elevation','Humerus elevation','Humerus axial rotation'};

% ---------- Figure setup ----------
figure('Color','w','Units','inches','Position',[1 1 7.16 4.5]);
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
for i = 1:9
    nexttile; hold on; box on;
    if i ~= 3
        plot(t, rad2deg(kin_exp(:,i)), exp_style{:});
    end

    if i == 7 || i == 9
        kin_healthy_interp = fillmissing(kin_healthy,'linear');
        plot(t, rad2deg(kin_healthy_interp(:,i)), healthy_interp_style{:}); hold on
        plot(t, rad2deg(kin_healthy(:,i)), healthy_style{:});

        kin_RC_interp = fillmissing(kin_RC,'linear');
        plot(t, rad2deg(kin_RC_interp(:,i)), RC_interp_style{:}); hold on
        plot(t, rad2deg(kin_RC(:,i)), RC_style{:});
    else
        plot(t, rad2deg(kin_healthy(:,i)), healthy_style{:});
        plot(t, rad2deg(kin_RC(:,i)), RC_style{:});
    end
    


    title(labels{i});
    xlabel(['Time (s)',newline,'']);
    ylabel('Angle (deg)');

    xlim([t(1) t(end)]);
    set(gca,'FontSize',8,'LineWidth',0.2);
    text(0.5, -0.35, ['(',alphabet{i},')'], 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 8);
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


% ---------- Legend ----------
lg = legend({'Experimental data','', ...
             'Healthy','', ...
             'RC-limited'}); %, ...
             % 'Orientation','horizontal', ...
             % 'Location','southoutside');
lg.FontSize = 8;
lg.Box = 'off';
lg.Layout.Tile = 'south';
lg.Orientation = 'horizontal';
set(gca, 'LineWidth',0.2)

% set(gcf, 'Units', 'inches', 'Position', [1 1 7.16 3])
% exportgraphics(gcf,'IEEE_kinematics_healthy_RClim2.png','Resolution',600);
    
end

function [SCHR, phaseData] = compute_SCHR(TH, TS, GH)

% Ensure column vectors
TH = TH(:);
TS = TS(:);

% Remove NaNs (if present)
validIdx = ~isnan(TH) & ~isnan(TS);
TH = TH(validIdx);
TS = TS(validIdx);

% Optional: smooth data slightly (recommended if noisy)
% TH = smoothdata(TH,'movmean',5);
% TS = smoothdata(TS,'movmean',5);

% Define phase limits
phases = [ ...
    min(TH) 30; 
    30 60; 
    60 90;
    min(TH) 90];

nPhases = size(phases,1);

SCHR = zeros(nPhases,1);
phaseData = struct();

for i = 1:nPhases
    
    lower = phases(i,1);
    upper = phases(i,2);
    
    % Find indices within this TH range
    idx = TH >= lower & TH <= upper;
    
    if sum(idx) < 2
        SCHR(i) = NaN;
        continue;
    end
    
    % Compute elevation change
    deltaGH = max(GH(idx)) - min(GH(idx));
    deltaTS = max(TS(idx)) - min(TS(idx));
    
    % Avoid division by zero
    if abs(deltaTS) < 1e-6
        SCHR(i) = NaN;
    else
        SCHR(i) = round(deltaGH / deltaTS,2);
    end
    
    % Store additional info
    phaseData(i).rangeTH = [lower upper];
    phaseData(i).deltaTH = deltaGH;
    phaseData(i).deltaTS = deltaTS;
    
end

end


function plot_EMG_healthy_RClim_IEEE(EMG_struct, OS_model,all_or_one,GH_plot, varargin)

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
        
            % plot(rsmpl_simulation,M,'k')
            % hold on
            % plot(rsmpl_simulation,upper_bound,'k',rsmpl_simulation,lower_bound,'k')
            % hold on
            % plot(rsmpl_simulation,simulation.data.inputs(1:end-ending_val,mus_index),'r','LineWidth',1.5)
            % patch([time' fliplr(time')], [lower_bound fliplr(upper_bound)], 'g')
            % fill([rsmpl_simulation'; flip(rsmpl_simulation')],[lower_bound; flip(upper_bound)], 'b', 'edgecolor', 'none', 'facealpha', 0.1)
            % time_emg
            % current_emg
            % figure
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

exportgraphics(gcf,'IEEE_emg_healthy_RClim2.png','Resolution',600);

end

function plot_EMG_optim_IEEE(EMG_struct, OS_model,all_or_one,GH_plot, varargin)

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
        
            % plot(rsmpl_simulation,M,'k')
            % hold on
            % plot(rsmpl_simulation,upper_bound,'k',rsmpl_simulation,lower_bound,'k')
            % hold on
            % plot(rsmpl_simulation,simulation.data.inputs(1:end-ending_val,mus_index),'r','LineWidth',1.5)
            % patch([time' fliplr(time')], [lower_bound fliplr(upper_bound)], 'g')
            % fill([rsmpl_simulation'; flip(rsmpl_simulation')],[lower_bound; flip(upper_bound)], 'b', 'edgecolor', 'none', 'facealpha', 0.1)
            % time_emg
            % current_emg
            % figure
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
% lg = legend({'EMG','Unscaled simulation','Scaled simulation'}, ...
%             'Orientation','horizontal','Location','southoutside');
%         lg.NumColumns = 2;

lg = legend({'EMG','Original parameters','Adjusted parameters'});
lg.Layout.Tile = 'south';
lg.Orientation = 'horizontal';
lg.FontSize = 8;
lg.Box = 'off';
lg.ItemTokenSize = 10;

exportgraphics(gcf,'IEEE_emg_orig_adjusted.png','Resolution',600);


end

function plot_paper_kinematics_ISG(OS_struct,OS_model,plot_GH_seq,plot_quaternion,varargin)
    num_res = length(varargin);
    legend_names = {};
    line_styles = {'--','-.','-'};
    line_colors = [0, 0, 1, 1; 1, 0, 0, 1; 0.4660, 0.6740, 0.1880, 1];
    figure
    for ires = 1:num_res
        iresult = varargin{ires};
        file_name = iresult{1};
        rot_type = iresult{2};
        motion_name = iresult{3};
        GH_seq = iresult{4};
        participant = iresult{5};
        plot_name = iresult{6};
        result = load(['Motions\',participant,'\',motion_name,'\',file_name,'.mat']);
        % result = load(['Motions\',participant,'/',motion_name,'\res_',rot_type,'_',motion_name,'_',weight,'.mat']);
        time = result.data.tout;
        Nfr = length(time);
        startG = floor(0.38*Nfr);
        endG = floor(0.62*Nfr);
        time = time(startG:endG);
        percent_of_motion = linspace(0,100,length(time));
        traj = result.data.trajectories;
        traj = traj(startG:endG,:);
        dofs_names = {'R_{Y}^{clav} in W_{F}','R_{Z}^{clav} in W_{F}','R_{X}^{clav} in W_{F}','R_{Y}^{scap} in W_{F}','R_{Z}^{scap} in W_{F}','R_{X}^{scap} in W_{F}','R_{Y}^{hum} in W_{F}','R_{Z}^{hum} in W_{F}',['R_{',GH_seq(end),'Y}^{hum} in W_{F}'],'Elbow'};
    
        if strcmp(rot_type,'euler')
            traj_obj_eul = rad2deg(create_objective_traj_eul(traj,GH_seq));
            traj_obj = traj_obj_eul;
            traj_gh = traj(:,8);
        elseif strcmp(rot_type,'quat')
            traj_in_eul = quat2eul_motion(traj,GH_seq);
            traj_obj_quat = rad2deg(create_objective_traj_eul(traj_in_eul,GH_seq));
            traj_obj = traj_obj_quat;
            for itime = 1:length(traj(:,1))
                nSC(itime) = norm(traj(itime,1:4));
                nAC(itime) = norm(traj(itime,5:8));
                nGH(itime) = norm(traj(itime,9:12));
            end
            traj_gh = traj_in_eul(:,8);
            % sum(nSC-1)
            % sum(nAC-1)
            % sum(nGH-1)
        end
        
        plot(percent_of_motion,traj_obj(:,5),'Color',line_colors(ires,:),'Linestyle',line_styles{ires},'LineWidth',3)
        % plot(traj_obj(:,8),traj_obj(:,5),'Color',line_colors(ires,:),'Linestyle',line_styles{ires},'LineWidth',3)

        legend('Healthy','RCI')
        hold on
        xlabel('% of motion','FontSize',15,'FontWeight','bold')
        ylabel('Scapula - upward rotation [deg]','FontSize',15,'FontWeight','bold')

        end


    OS_interp = interp1(OS_struct.mot_struct.time,OS_struct.mot_struct.euler,time',"spline");
    OS_interp_obj = rad2deg(create_objective_traj_eul(OS_interp,plot_GH_seq));
    % legend_names{end+1} = 'Inverse kinematics';
    plot(percent_of_motion, OS_interp_obj(:,5),'g','LineWidth',3); hold on



    hold off
    fig = gca;
    % ax = axes('Parent',fig);
    fig.FontSize = 15;
    Lgnd = legend({'RC intact','RC impairment','Healthy participant (IK)'},'FontSize',16,'Location','south'); %

    title(['Upward rotation of scapula',newline,'during humeral elevation'],'FontSize',25)

    % Lgnd.Position(1) = 0.62;
    % Lgnd.Position(2) = 0.075;
    % Lgnd.Direction = "reverse";
    % fig.Position(3) = fig.Position(3) - 1;
    fig.Position(3) = fig.Position(3);
    fig.Position(4) = fig.Position(4);
    % Lgnd.Position = 'northwest';
    % axes('FontSize',25)
    % exportgraphics(fig,['Motions\',participant,'\',motion_name,'\Results\kinematics_all_',participant,'_',motion_name,'.png'],'Resolution',600);
    % print(gcf,'-vector','-dsvg',['kinematics.svg'])
    exportgraphics(fig,'scapula_ISG.png','Resolution',600);
    
end


function [force,y_comp,z_comp] = plot_GH_reactions(result)
file_name = result{1};
rot_type = result{2};
motion_name = result{3};
GH_seq = result{4};
participant = result{5};
result = load(['Motions\',participant,'\',motion_name,'\',file_name,'.mat']);
reactions = result.data.reactions;
traj = result.data.trajectories;

for i = 1:size(reactions,1)
    reactions_rotated =  R_y(16.5*pi/180)' * R_z(-6*pi/180)' * [reactions(i,:)';1];
    reactions_in_hum(i,:)= R_y(45*pi/180)' * Qrm(traj(i,9:12))' * [-reactions(i,:)';1];
    react_rot_all(i,:) = reactions_rotated;
    fvec_norm = norm(reactions_rotated);
    force(i) = fvec_norm;
    % y_comp(i) = atan(reactions(i,2)/fvec_norm)*180/pi;
    % z_comp(i) = atan(reactions(i,3)/fvec_norm)*180/pi;
    y_comp(i) = -atan2(-reactions_rotated(2),-reactions_rotated(1))*180/pi;
    z_comp(i) = atan2(-reactions_rotated(3),-reactions_rotated(1))*180/pi;
end
end

function computational_performance_eul_vs_quat()
motions = {'Elevation_yzy', 'Scabduction_yzy', 'Flexion_yzy'};
participants = {'1LU','3NA','8tata'};
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
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

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
ylabel('Iterations')
t = title('Solver Iterations');
% t.Position(2) = t.Position(2)*1.02;

set(gca,'FontSize',8,'Box','off','YMinorTick','on')


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
yticks([500 1000 1500 2000 3000 4000])
yticklabels({'500','1000','1500','2000','3000','4000'})


xticklabels({'Euler\newline angles','Quaternions'})
ylabel('Time to solve [s]')
t = title('Solver Time');

legend(hGL,['Gimbal lock',newline, 'occurrence'],'Position',[0.4 0.45 0.9 0.5],'Box','off')

set(gca,'FontSize',8,'Box','off','YMinorTick','on')

sgtitle('Computational Performance','FontWeight','bold','FontSize',10)

% exportgraphics(gcf,'comp_perf.png','Resolution',600)

end

function plot_conoid_length(OS_model,result)

file_name = result{1};
rot_type = result{2};
motion_name = result{3};
GH_seq = result{4};
participant = result{5};
% data = load(['Motions\',motion_name,'\res_',rot_type,'_',motion_name,'_',weight,'.mat']);
struct_path = ['Motions\',participant,'\',motion_name,'\',file_name,'.mat'];
% struct_path = ['Motions\',participant,'\All_motions\res_',rot_type,'_',motion_name,'.mat'];

load(struct_path);
trajectory = data.trajectories;
if strcmp(rot_type,'quat')
    trajectory = quat2eul_motion(trajectory,GH_seq);
end
numdata = size(trajectory,1);
OS_model = load(OS_model);
for i = 1:numdata
    actual_conoid_length(i) = get_conoid_length(trajectory(i,4:6),OS_model.model);
    actual_conoid_force(i) = get_conoid_force(actual_conoid_length(i),OS_model.model);
    expected_traj(i,:) = min_conoid_length(trajectory(i,:),OS_model.model);
    expected_length(i) = get_conoid_length(expected_traj(i,4:6),OS_model.model);
    expected_force(i) = get_conoid_force(expected_length(i));
end

if strcmp(rot_type,'euler')
    expected_traj(1,1:3);
    format long
    disp(expected_traj(1,1:3))
    data.q_clavicula_init = expected_traj(1,1:3);
    save(struct_path,'data')
end

if strcmp(rot_type,'quat')
    expected_traj_init = eul2quat_motion(expected_traj,GH_seq);
    format long
    disp(expected_traj_init(1,1:4));
    data.q_clavicula_init = expected_traj_init(1,1:8);
    save(struct_path,'data')
end



figure
subplot(3,1,1)
plot(rad2deg(trajectory(:,3)))
title('actual angle')
subplot(3,1,2)
plot(actual_conoid_length)
title('actual length')
subplot(3,1,3)
plot(actual_conoid_force)
title('actual force')

figure
subplot(3,1,1)
plot(rad2deg(expected_traj(:,3)))
title('exptected angle')
subplot(3,1,2)
plot(expected_length)
title('expected length')
subplot(3,1,3)
plot(expected_force)
title('expected force')
end

function plot_objective_values(varargin)
    num_res = length(varargin);
    legend_names = {}
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
        plot(obj_values,'LineWidth',1.5); hold on
        if strcmp(rot_type,'euler')
            legend_names{end+1} = [rot_type,' ',GH_seq]
        else
            legend_names{end+1} = [rot_type]
        end
    end
    legend(legend_names)
    hold off

    
end

function plot_GH_near_GL(varargin)
    num_res = length(varargin);
    iplot = 1;
    colors_fig = {[0.4660, 0.6740, 0.1880, 1],[0, 0, 1, 1]};
    figure
    for ires = 1:num_res
        iresult = varargin{ires};
        motion_name = iresult{1};
        rot_type = iresult{2};
        weight = iresult{3};
        GH_seq = iresult{4};
        result = load(['Motions\',motion_name,'\res_',rot_type,'_',motion_name,'_',weight,'.mat']);
        time = result.data.tout;
        percent_of_motion = linspace(0,100,length(time));
        traj = rad2deg(result.data.trajectories);
        if strcmp(GH_seq,'YZY')
            dofs_names = {'GH_{y}','GH_{z}','GH_{yy}'};
        elseif strcmp(GH_seq,'YZX')
            dofs_names = {'GH_{y}','GH_{z}','GH_{x}'};
        end
        
        for i = 1:length(dofs_names)
        subplot(num_res,3,iplot)
        plot(percent_of_motion,traj(:,i+6),'Color',colors_fig{ires},'LineWidth',1.5)
        if i ==2 && strcmp(GH_seq,'YZY')
            yline(0,'-.','Color','Black','LineWidth',1.5)
            text(20,7,'\downarrow Gimbal lock')
            axis([-inf, inf, -5, inf])
        else
            axis([-inf, inf, -inf, inf])
        end
        xlabel(['% of humeral elevation'])
        ylabel('Angle [deg]')
        title(dofs_names{i})
        iplot = iplot+1;
        end
        
    
    end
    fig = gcf;
    fig.Position(3) = fig.Position(3) + 300;
    fig.Position(4) = fig.Position(4) - 150;
    exportgraphics(fig,'plot_near_GH.png','Resolution',600);
    % Lgnd.Position(1) = 0.7;
    % Lgnd.Position(2) = 0.14;
end

function plot_paper_activations(OS_model,mus_group,plot_excitation,varargin)
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
    if rem(num_in_group,2) == 0
        plot_rows = plot_rows + 1;
    end
    current_names = muscle_names(mask);
    legend_names = {};
    num_res = length(varargin);
    figure
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
        for imus = 1:num_in_group
            subplot(plot_rows,3,imus)
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
            ylabel('act [-]')
        end


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

    end % end num_res

    for imus = 1:num_in_group
        if ~isempty(GL_pos_prcnt)
            % subplot(plot_rows,2,imus)
            xline(GL_pos_prcnt,'LineWidth',1.5);
        end
    end

    if ~isempty(GL_pos_prcnt)
        for GH_num = 1:length(GL_pos_prcnt)
            if GH_num == 1
                legend_names{end+1} = ['Gimbal lock'];
            else
                legend_names{end+1} = [''];
            end
        end
    end


    fig = gcf;
    
    fig.Position(3) = fig.Position(3)+200; %200 normal  *0.85 GL fig
    fig.Position(4) = fig.Position(4)+70; %70 normal  *0.8 GL fig
    Lgnd = legend(legend_names,'FontSize',10); %
    % Lgnd.Position(1) = 0.37;
    % Lgnd.Position(2) = 0.26;
    Lgnd.IconColumnWidth = 10;
    % annotation('textarrow',[0.4,0.3],[0.28,0.4],'String','  ','FontSize',13,'Linewidth',2)
    % annotation('textarrow',[0.65,0.82],[0.3,0.35],'String','  ','FontSize',13,'Linewidth',2)
    % annotation('textbox',[0.4 0.38 0.1 0.1],'String',['Gimbal lock for',newline,'Euler angles',newline, 'based-model'],'FontSize',13,'Linewidth',2,'EdgeColor','None','FontWeight','bold')

    % sgtitle('Activations', 'Interpreter', 'none')
    
    hold off
    % print(gcf,'-vector','-dsvg',['act_GL.svg'])
    % movegui(fig,'center');
    % exportgraphics(fig,['Motions\',participant,'\',motion_name,'\Results\activations_only_',participant,'_',motion_name,'.png'],'Resolution',600);
    % exportgraphics(gca,'act_GL.jpg',...   % since R2020a
    % print(gcf,'activation_near_GL.png','-dpng','-r300') 
    % 'BackgroundColor','none')
end

function plot_polynomials(OS_model,mus_group,lceopt_scaler,varargin)
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
        plot_rows = plot_rows + 2;
    end
    current_names = muscle_names(mask);
    legend_names = {};
    num_res = length(varargin);
    indexes = find(mask);
    ilength_plot = 1;

    for imus = 1:length(indexes)
        figure
        for ires = 1:num_res
            iresult = varargin{ires};
            file_name = iresult{1};
            rot_type = iresult{2};
            motion_name = iresult{3};
            GH_seq = iresult{4};
            participant = iresult{5};
            result = load(['Motions\',participant,'\',motion_name,'\',file_name,'.mat']);
            % result = load(['Motions\',participant,'\All_motions\res_',rot_type,'_',motion_name,'.mat']);

            model_res = load(['Motions\',participant,'\OS_model.mat']);
            muscles_res = model_res.model.muscles;
            time = result.data.tout;
            traj = result.data.trajectories;
            speeds = result.data.speeds;
            numdata = size(traj,1);
            current_mus = muscles_res{indexes(imus)};
            current_name = current_mus.osim_name;
            dof_names = current_mus.dof_names;
            dof_indeces = current_mus.dof_indeces-3;
    
            if strcmp(rot_type,'euler')
                [lengths,moment_arms] = momarms(current_mus.Euler, dof_indeces, traj);
                motion = traj;
                
            elseif strcmp(rot_type,'quat')
                motion_quat_WO_real = traj(:,[2:4,6:8,10:12,13]);
                size(motion_quat_WO_real)
                [lengths,jacobian] = momarms(current_mus.Quaternion, dof_indeces, motion_quat_WO_real);
                mus_velocity = zeros(size(lengths));
                for iframe = 1:length(lengths)
                    try
                        dqdtSC = dquatdt(traj(iframe,1:4),speeds(iframe,1:3));
                        dqdtAC = dquatdt(traj(iframe,5:8),speeds(iframe,4:6));
                        dqdtGH = dquatdt(traj(iframe,9:12),speeds(iframe,7:9));
                        dqEL = speeds(iframe,10);
                        dqdt_full = [dqdtSC(2:4);dqdtAC(2:4);dqdtGH(2:4);dqEL];
                        mus_velocity(iframe) = jacobian(iframe,:) * dqdt_full(dof_indeces);
                    catch
                        mus_velocity(iframe) = 0;
                    end
                end
                try
                    fvecs_poly = fvectors(current_mus, dof_indeces, motion_quat_WO_real);
                catch
                    fvecs_poly = zeros(numdata,3);
                end
                moment_arms_full = zeros(numdata,11);
                jacobian_quat = zeros(numdata,11);
                jacobian_quat(:,dof_indeces) = jacobian;
                alljoints = {'YZX','YZX',GH_seq};
                motion = quat2eul_motion(traj, GH_seq);
    
                for iframe=1:numdata
                    for j = 1:3
                        ind3 = ((j-1)*3+1:(j-1)*3+3);
                        ind4 = ((j-1)*4+1:(j-1)*4+4);
                        JQuatInSpat = invJtrans(traj(iframe,ind4)) * jacobian_quat(iframe,ind3)';
                        moment_arms_full(iframe,ind3) = GeomJ(motion(iframe,ind3),alljoints{j})*(JQuatInSpat);
                    end
                     moment_arms_full(iframe,10) = jacobian_quat(iframe,10);
                     moment_arms_full(iframe,11) = jacobian_quat(iframe,11);
                end
                moment_arms = moment_arms_full(:,dof_indeces);
    
            end
        
            [lengths_OS, moment_arms_OS,fvecs_OS] = opensim_get_polyvalues(motion, indexes(imus), current_mus.dof_indeces,GH_seq);
            max_momarm = max([max(abs(moment_arms),[],'all'),max(abs(moment_arms_OS),[],'all')]) + 0.01;
    
            for j=1:length(current_mus.dof_indeces)
                subplot(length(current_mus.dof_indeces)+1,num_res,j*num_res+ires-num_res)
                plot(time,moment_arms_OS(:,j),time,-moment_arms(:,j)); hold on
                axis([-inf inf -max_momarm max_momarm])
            end
            subplot(length(current_mus.dof_indeces)+1,num_res,length(current_mus.dof_indeces)+1)
            % plot(time,(lengths_OS-current_mus.lslack)/current_mus.lceopt/lceopt_scaler,time,(lengths-current_mus.lslack)/current_mus.lceopt/lceopt_scaler); hold on
            kpe = 5;
            epsm0 = 0.6;
            lm = lengths-current_mus.lslack;
            fpe = (exp(kpe*(lm / current_mus.lceopt - 1)/epsm0)-1)/(exp(kpe)-1);
            plot(time,fpe*current_mus.fmax)
    
            % subplot(length(current_mus.dof_indeces)+1,num_res,num_res*length(current_mus.dof_indeces)+ilength_plot)
            % plot(time, lengths_OS, time, lengths)
            % title([rot_type,' ',GH_seq])
        
            
        legend('OpenSim momarm','Polynomial')
        sgtitle(current_name)

        figure
        subplot(3,1,1)
        plot(fvecs_poly(:,1)); hold on
        plot(fvecs_OS(:,1))
        axis([-inf inf min(fvecs_poly,[],'all')-0.1 max(fvecs_poly,[],'all')+0.1])
        subplot(3,1,2)
        plot(fvecs_poly(:,2)); hold on
        plot(fvecs_OS(:,2))
        axis([-inf inf min(fvecs_poly,[],'all')-0.1 max(fvecs_poly,[],'all')+0.1])
        subplot(3,1,3)
        plot(fvecs_poly(:,3)); hold on
        plot(fvecs_OS(:,3))
        axis([-inf inf min(fvecs_poly,[],'all')-0.1 max(fvecs_poly,[],'all')+0.1])
        sgtitle(current_name)
        legend('Polynomial','OpenSim fvec')
        end
        ilength_plot = ilength_plot + 1;
    end
    % d1 = -0.318;
    % d2 = -8.149;
    % d3 = -0.374;
    % d4 = 0.886;
    vmax_norm = 10 * current_mus.lceopt;
    vnorm = mus_velocity/vmax_norm;
    % 
    % fvce = d1 * log(d2 * vnorm + d3 + sqrt((d2 * vnorm + d3).^2 + 1)) + d4;
    figure
    velocityx = linspace(0,100,length(time));
    plot(velocityx,mus_velocity)

end

function plot_paper_kinematics(OS_struct,OS_model,plot_GH_seq,plot_quaternion,varargin)
    num_res = length(varargin);
    legend_names = {};
    line_styles = {'--','-.','-'};
    line_colors = [0, 0, 1, 1; 1, 0, 0, 1; 0.4660, 0.6740, 0.1880, 1];
    figure
    for ires = 1:num_res
        iresult = varargin{ires};
        file_name = iresult{1};
        rot_type = iresult{2};
        motion_name = iresult{3};
        GH_seq = iresult{4};
        participant = iresult{5};
        plot_name = iresult{6};
        result = load(['Motions\',participant,'\',motion_name,'\',file_name,'.mat']);
        % result = load(['Motions\',participant,'/',motion_name,'\res_',rot_type,'_',motion_name,'_',weight,'.mat']);
        time = result.data.tout;
        percent_of_motion = linspace(0,100,length(time));
        traj = result.data.trajectories;
        dofs_names = {'R_{Y}^{clav} in W_{F}','R_{Z}^{clav} in W_{F}','R_{X}^{clav} in W_{F}','R_{Y}^{scap} in W_{F}','R_{Z}^{scap} in W_{F}','R_{X}^{scap} in W_{F}','R_{Y}^{hum} in W_{F}','R_{Z}^{hum} in W_{F}',['R_{',GH_seq(end),'Y}^{hum} in W_{F}'],'Elbow'};
    
        if strcmp(rot_type,'euler')
            traj_obj_eul = rad2deg(create_objective_traj_eul(traj,GH_seq,1));
            traj_obj = traj_obj_eul;
            traj_gh = traj(:,8);
        elseif strcmp(rot_type,'quat')
            traj_in_eul = quat2eul_motion(traj,GH_seq);
            traj_obj_quat = rad2deg(create_objective_traj_eul(traj_in_eul,GH_seq,1));
            traj_obj = traj_obj_quat;
            for itime = 1:length(traj(:,1))
                nSC(itime) = norm(traj(itime,1:4));
                nAC(itime) = norm(traj(itime,5:8));
                nGH(itime) = norm(traj(itime,9:12));
            end
            traj_gh = traj_in_eul(:,8);
            rmse(nSC,ones(1,length(nSC)))
            rmse(nAC,ones(1,length(nAC)))
            rmse(nGH,ones(1,length(nGH)))
        end
        

        for i = 1:length(dofs_names)
        subplot(4,3,i)
        plot(percent_of_motion,traj_obj(:,i),'Color',line_colors(ires,:),'Linestyle',line_styles{ires},'LineWidth',1.5)
        hold on
        title(dofs_names{i})
        xlabel('% of motion')
        if i == 1 || i == 4 || i == 7 || i == 10
            ylabel('Angle [deg]')
        end
        


        end
     if strcmp(rot_type,'euler')
        legend_names{end+1} = plot_name;
    elseif strcmp(rot_type,'quat')
        legend_names{end+1} = plot_name;
     end

    subplot(4,3,11)
    plot(percent_of_motion,traj_gh,'Color',line_colors(ires,:),'Linestyle',line_styles{ires},'LineWidth',1.5); hold on
    
    end

    OS_interp = interp1(OS_struct.mot_struct.time,OS_struct.mot_struct.euler,time',"spline");
    OS_interp_obj = rad2deg(create_objective_traj_eul(OS_interp,plot_GH_seq,0));
    legend_names{end+1} = 'Inverse kinematics';

    try
        OS_interp_obj;
        traj_obj;
        rms_eul = rmse(OS_interp_obj,traj_obj_eul,"all")
        rms_quat = rmse(OS_interp_obj,traj_obj_quat,'all')
    end

    for i = 1:length(dofs_names)
        subplot(4,3,i)
        if i == 3
            continue
        else
            plot(percent_of_motion, OS_interp_obj(:,i),'g','LineWidth',1.5); hold on
            % axis([-inf inf -inf inf])  

    
            
            

        end
        % OS_interp_quat = eul2quat_motion(OS_interp,GH_seq);
        % OS_interp_quat_obj = create_objective_traj_quat(OS_interp_quat);
    end



    hold off
    fig = gcf;
    Lgnd = legend(legend_names,'FontSize',15); %
    Lgnd.Position(1) = 0.62;
    Lgnd.Position(2) = 0.075;
    % Lgnd.Direction = "reverse";
    % fig.Position(3) = fig.Position(3) - 1;
    fig.Position(3) = fig.Position(3)*1.2-50;
    fig.Position(4) = fig.Position(4)*1.2;

    % exportgraphics(fig,['Motions\',participant,'\',motion_name,'\Results\kinematics_all_',participant,'_',motion_name,'.png'],'Resolution',600);
    % print(gcf,'-vector','-dsvg',['kinematics.svg'])
    % exportgraphics(fig,'kinematics_DC.png','Resolution',600);

    if plot_quaternion == 1
    figure
        for ires = 1:num_res
        iresult = varargin{ires};
        file_name = iresult{1};
        rot_type = iresult{2};
        motion_name = iresult{3};
        GH_seq = iresult{4};
        participant = iresult{5};
        result = load(['Motions\',participant,'\',motion_name,'\',file_name,'.mat']);
        time = result.data.tout;
        percent_of_motion = linspace(0,100,length(time));
        traj = result.data.trajectories;
    
        if strcmp(rot_type,'quat')
            result_quat_obj = create_objective_traj_quat(traj,OS_model);
        else
            continue
        end

        OS_interp = interp1(OS_struct.mot_struct.time,OS_struct.mot_struct.euler,time',"spline");
        OS_interp_quat = eul2quat_motion(OS_interp,plot_GH_seq);
        OS_interp_quat_obj = create_objective_traj_quat(OS_interp_quat,OS_model);
        
        subplot(3,4,1)
        plot(percent_of_motion,result_quat_obj(:,1),'Color',[1,0,0,1],'LineStyle','-.',"LineWidth",1.5); hold on
        plot(percent_of_motion, OS_interp_quat_obj(:,1),'g','LineWidth',1.5)
        title('AC_{pos} in W_{F}^X')
        ylabel('position [m]')
        xlabel('% of motion')
        subplot(3,4,2)
        plot(percent_of_motion,result_quat_obj(:,2),'Color',[1,0,0,1],'LineStyle','-.',"LineWidth",1.5); hold on
        plot(percent_of_motion, OS_interp_quat_obj(:,2),'g','LineWidth',1.5)
        title('AC_{pos} in W_{F}^Y')
        xlabel('% of motion')
        subplot(3,4,3)
        plot(percent_of_motion,result_quat_obj(:,3),'Color',[1,0,0,1],'LineStyle','-.',"LineWidth",1.5); hold on
        plot(percent_of_motion, OS_interp_quat_obj(:,3),'g','LineWidth',1.5)
        title('AC_{pos} in W_{F}^Z')
        xlabel('% of motion')

        quat_plot_names = {'Q_0^{scap} in W_{F}','Q_1^{scap} in W_{F}','Q_2^{scap} in W_{F}','Q_3^{scap} in W_{F}','Q_0^{hum} in W_{F}','Q_1^{hum} in WF','Q_2^{hum} in W_{F}','Q_3^{hum} in W_{F}'};
        for iquat = 1:8
            subplot(3,4,iquat+4)
            plot(percent_of_motion,result_quat_obj(:,iquat+4),'Color',[1,0,0,1],'LineStyle','-.',"LineWidth",1.5); hold on
            plot(percent_of_motion, OS_interp_quat_obj(:,iquat+4),'g','LineWidth',1.5);
            title(quat_plot_names{iquat})
            if iquat == 1 || iquat == 5
                ylabel('[-]')
            end
            % if iquat == 1
            xlabel(['% of motion',newline])
            % else
            %     xlabel(['% of motion'])
            % end

        end

        end

        % xlabel('% of motion','Position',[-10,-0.1])

    fig = gcf;
    fig.Position(3) = fig.Position(3) + 250;
    fig.Position(4) = fig.Position(4) + 0;
    Lgnd = legend({['Quaternion',newline,'model'],['Inverse',newline,'kinematics']}); %,'FontSize',15
    Lgnd.Position(1) = 0.75;
    Lgnd.Position(2) = 0.7;
    exportgraphics(fig,'kinematics_quat_only.png','Resolution',600);
    % exportgraphics(fig,['Motions\',participant,'\',motion_name,'\Results\kinematics_quat_',participant,'_',motion_name,'.png'],'Resolution',600);

    end
    
end

function plot_elipsoid_eq(OS_model,varargin)
    ires = 1;
    iresult = varargin{ires};
    motion_name = iresult{1};
    rot_type = iresult{2};
    weight = iresult{3};
    GH_seq = iresult{4};
    participant = iresult{5};
    result = load(['Motions\',participant,'/',motion_name,'\res_',rot_type,'_',motion_name,'_',weight,'.mat']);
    time = result.data.tout;
    traj = result.data.trajectories;
    model = load(OS_model);
    numdata = length(time);
    Epos = model.model.thoracic_wall.translation;
    Edim = model.model.thoracic_wall.dimensions;
    for i = 1:numdata
        [iTS,iAI] = contact_points_position(traj(i,:),model.model);
        eq_TS(i) = elips_eq(iTS,Epos,Edim);
        eq_AI(i) = elips_eq(iAI,Epos,Edim);
    end

    figure
    plot(eq_TS); hold on
    plot(eq_AI)
    legend('TS','AI')

    figure
    plot(sqrt(eq_TS.^2+eq_AI.^2))

end

function plot_SCx(folders, OS_model)
    GH_names = {'GHy','GHz','GHyy'};
    num_coords = 10;
    model = load(OS_model);
    euler_struct = load(folders{1});
    time = euler_struct.data.tout;
    traj_eul = euler_struct.data.trajectories;
    quat_struct = load(folders{2});
    traj_quat_orig = quat_struct.data.trajectories;
    % traj_quat_eul = quat2eul_motion(traj_quat_orig);

    for i = 1:length(time)
        newSCx_eul(i,:) = min_conoid_length(traj_eul(i,:));
        length_eul(i) = get_conoid_length(newSCx_eul(i,4:6),model.model);
        force(i) = get_conoid_force(length_eul(i));
        % newSCx_quat(i,:) = min_conoid_length(traj_quat_eul(i,:));
    end

    figure
    subplot(2,1,1)
    plot(time,rad2deg(newSCx_eul(:,3)),time,rad2deg(traj_eul(:,3)),'*')
    legend('Min length','simulation')
    subplot(2,1,2)
    plot(time,length_eul)
    title('eul')
    figure
    plot(get_conoid_force)

    % figure
    % plot(time,rad2deg(newSCx_quat(:,3)),time,rad2deg(traj_quat_eul(:,3)),'*')
    % legend('Min length','simulation')
    % title('quat')
    

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

function plot_activations_EMG(EMG_struct, OS_model,all_or_one,GH_plot, varargin)
line_styles = {'--','-.','-'};
line_colors_act = [0, 0, 1, 1; 1, 0, 0, 1; 0.4660, 0.6740, 0.1880, 1];
line_colors_exc = [0, 0, 1, 0.5; 1, 0, 0, 0.5; 0.4660, 0.6740, 0.1880, 1];


EMG_muscles = {'AnteriorDelt','IntermediateDelt','PosteriorDelt','Infrasp','Suprasp','MiddleTrap','UpperTrap','Serrupper'};
% EMG_muscles = {'AnteriorDelt','IntermediateDelt','PosteriorDelt','Suprasp','MiddleTrap','Serrupper'};
if contains(EMG_struct,'flexion')
    ignore = {[],[2],[9],[9],[9],[],[],[6]}; %flexion
    motion_name_figure = 'flexion';
elseif contains(EMG_struct, 'scabduction')
    ignore =      {[],[2,8],[],[9],[4],[3],[],[6]}; %scabduction
    motion_name_figure = 'scapular plane elevation';
elseif contains(EMG_struct,'elevation')
    ignore = {[],[],[9],[9],[4,1],[],[],[6]}; %elevation
    motion_name_figure = 'abduction';
end %,"delt_clav_1","delt_clav_2"
model_names = {["delt_clav_1","delt_clav_2"],["delt_scap_9","delt_scap_10","delt_scap11"],["delt_scap_3","delt_scap_4","delt_scap_5"],["infra_3","infra_4","infra_5"],["trap_scap_1","trap_scap_2"],["trap_scap_6","trap_scap_7","trap_scap_8"],["trap_clav_1"],["serr_ant_2","serr_ant_3","serr_ant_4"]};
figure_names = {'Anterior deltoid', 'Lateral deltoid', 'Posterior deltoid', 'Infraspinatus',['Scapular trapezius',newline,'- descending'], ['Scapular trapezius',newline,'- transverse'],' Clavicular trapezius','Serratus anterior'};
% model_names = {["delt_scap11"],["delt_scap_9","delt_scap10"],["delt_scap_3","delt_scap_4"],"trap_scap_2","trap_scap_6",["serr_ant_3","serr_ant_4"]};
% figure_names = {'Anterior deltoid', 'Lateral deltoid', 'Posterior deltoid',['Scapular trapezius',newline,'- descending'], ['Scapular trapezius',newline,'- transverse'],'Serratus anterior'};
num_res = length(varargin);
% simulation = load(folders{2});
% ending_val = 0;


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
figure
tiledlayout(3,3);


for i = 1:length(EMG_muscles)
    nexttile
    legend_names = {};
    activation_rmse = [];
     for ires = 1:num_res
        iresult = varargin{ires};
        [fvec(ires,:),y_comp(ires,:),z_comp(ires,:)] = plot_GH_reactions(iresult);
        iresult = varargin{ires};
        file_name = iresult{1};
        rot_type = iresult{2};
        motion_name = iresult{3};
        GH_seq = iresult{4};
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
        time = linspace(0,100,length(result.data.tout));
        activation_rmse = [activation_rmse, activations];

        legend_names{end+1} = plot_name;
        legend_names{end+1} = '';

        plot(time,activations,'Color',line_colors_act(ires,:),'LineWidth',1.5) %,'Linestyle',line_styles{ires}
        hold on
        plot(time,excitations,'Color',line_colors_exc(ires,:),'LineWidth',1.5) %,'Linestyle',line_styles{ires}
        hold on

     end

    if strcmp(all_or_one,'all') 
    current_emg_rsmpld = zeros(length(rsmpl_simulation),14);
    for ipar = 1:14
        current_emg = emg_data.data.(['par_',num2str(ipar)]).(EMG_muscles{i});
        time_emg = linspace(0,100,length(current_emg));
        current_emg_rsmpld(:,ipar) = spline(time_emg,current_emg,rsmpl_simulation);
        % plot(rsmpl_simulation,current_emg_rsmpld(:,ipar))
        % hold on
    end
    current_emg_rsmpld(:,ignore{i}) = [];
    [S,M] = std(current_emg_rsmpld,0,2);
    upper_bound = M+S;
    lower_bound = M-S;
    % plot(rsmpl_simulation,M,'k')
    % hold on
    % plot(rsmpl_simulation,upper_bound,'k',rsmpl_simulation,lower_bound,'k')
    % hold on
    % plot(rsmpl_simulation,simulation.data.inputs(1:end-ending_val,mus_index),'r','LineWidth',1.5)
    % patch([time' fliplr(time')], [lower_bound fliplr(upper_bound)], 'g')
    fill([rsmpl_simulation'; flip(rsmpl_simulation')],[lower_bound; flip(upper_bound)], 'b', 'edgecolor', 'none', 'facealpha', 0.1);current_emg_rsmpld(:,ignore{i}) = [];
    [S,M] = std(current_emg_rsmpld,0,2);
    upper_bound = M+S;
    lower_bound = M-S;
    % plot(rsmpl_simulation,M,'k')
    % hold on
    % plot(rsmpl_simulation,upper_bound,'k',rsmpl_simulation,lower_bound,'k')
    % hold on
    % plot(rsmpl_simulation,simulation.data.inputs(1:end-ending_val,mus_index),'r','LineWidth',1.5)
    % patch([time' fliplr(time')], [lower_bound fliplr(upper_bound)], 'g')
    fill([rsmpl_simulation'; flip(rsmpl_simulation')],[lower_bound; flip(upper_bound)], 'b', 'edgecolor', 'none', 'facealpha', 0.1);
    
    axis([-inf inf 0 inf])
    xlabel(['% of motion'])
    ylabel('act [-]')
    legend_names{end+1} = 'filtered EMG';

    

    elseif strcmp(all_or_one,'one')
    rsmpl_simulation = linspace(0,100,length(activations));
    current_emg_rsmpld = zeros(length(activations),6);
    for imot = 1:6
        try
            current_emg = emg_data.data.(['num_',num2str(imot)]).(EMG_muscles{i});
        end
        time_emg = linspace(0,100,length(current_emg));
        current_emg_rsmpld(:,imot) = spline(time_emg,current_emg,rsmpl_simulation);
        % plot(rsmpl_simulation,current_emg_rsmpld(:,ipar))
        % hold on
    end
    % current_emg_rsmpld(:,ignore{i}) = [];
    [S,M] = std(current_emg_rsmpld,0,2);
    upper_bound = M+1*S;
    lower_bound = M-1*S;
    % plot(rsmpl_simulation,M,'k')
    % hold on
    % plot(rsmpl_simulation,upper_bound,'k',rsmpl_simulation,lower_bound,'k')
    % hold on
    % plot(rsmpl_simulation,simulation.data.inputs(1:end-ending_val,mus_index),'r','LineWidth',1.5)
    % patch([time' fliplr(time')], [lower_bound fliplr(upper_bound)], 'g')
    % fill([rsmpl_simulation'; flip(rsmpl_simulation')],[lower_bound; flip(upper_bound)], 'b', 'edgecolor', 'none', 'facealpha', 0.1)
    
    % length(current_emg_rsmpld(:,1))
    plot(rsmpl_simulation, current_emg_rsmpld(:,1))
    title(figure_names{i},'Interpreter','None','FontSize',10)
    axis([-inf inf 0 inf])
    xlabel(['% of motion'])
    ylabel('act [-]')
    legend_names{end+1} = 'filtered EMG';
        % current_emg = emg_data.data.(['par_',participant_number]).(EMG_muscles{i});
        % time_emg = linspace(0,100,length(current_emg));
        % plot(time_emg,current_emg,'Color',line_colors(3,:))
        % title(figure_names{i},'Interpreter','None')
        % axis([-inf inf 0 inf])
        % xlabel(['% of motion'])
        % ylabel('act [-]')
        % legend_names{end+1} = 'filtered EMG';
        
    end
    % RMS = sqrt(mean(( - ).^2));

    for ires = 1:num_res
        RMS(ires) = rmse(current_emg_rsmpld(:,1),activation_rmse(:,ires));
        MAE(ires) = sum(abs(current_emg_rsmpld(:,1)-activation_rmse(:,ires)))/length(activations);
        [r(ires), p] = corr(current_emg_rsmpld(:,1), activation_rmse(:,ires));
    end
    
    % RMS = rmse([1,1,1,1],[0,0,0,0])
    % [r, p] = corr(current_emg_rsmpld(:,1), activations); 
    % fprintf(muscle_names)
    % fprintf('RMS = %.3f -> %.3f \n', RMS(1), RMS(2))
    % fprintf('MAE = %.3f \n', MAE)
    % fprintf('Pearson correlation r = %.3f (p = %.3f)\n', r, p);
    % nrmse = rmse / (max(emg) - min(emg));  % normalized RMSE
    % fprintf('RMSE = %.3f', rmse);
    try
        title([figure_names{i}, newline, '|RMSE| ',num2str(RMS(1),'%.3f'),'->',num2str(RMS(2),'%.3f'), newline, '|MAE| ',num2str(MAE(1),'%.3f'),'->',num2str(MAE(2),'%.3f'),newline, '|r| ',num2str(r(1),'%.3f'),'->',num2str(r(2),'%.3f')],'Interpreter','None')
    catch
        title([figure_names{i}],'Interpreter','None')
    end

end




hold off

fig = gcf;
Lgnd = legend(legend_names,'FontSize',10); %
Lgnd.Position(1) = 0.02;
Lgnd.Position(2) = 0.65;
Lgnd.IconColumnWidth = 10;
fig.Position(3) = fig.Position(3) + 200;
fig.Position(4) = fig.Position(4) + 100;

ires = GH_plot;
axes('Position',[.7 .075 .25 .25])
box on
force = fvec(ires,:)/650*100;
numdata = size(force,2);
force1 = force(1:floor(numdata/3));
force2 = force(floor(numdata/3):floor(numdata*2/3));
force3 = force(floor(numdata*2/3):end);
time1 = time(1:floor(numdata/3));
time2 = time(floor(numdata/3):floor(numdata*2/3));
time3 = time(floor(numdata*2/3):end);
plot(time1,force1,time2,force2,time3,force3,'LineWidth',1.5); hold on
% legend('1','2','3','Location','northeast')
title('GH force/stability')
ylabel('Reaction force [%BW]')
xlabel('% of motion')

axes('Position',[.825 .2175 .1 .1])
box on
y_cur = y_comp(ires,:);
z_cur = z_comp(ires,:);
y_cur1 = y_cur(1:floor(numdata/3));
y_cur2 = y_cur(floor(numdata/3):floor(numdata*2/3));
y_cur3 = y_cur(floor(numdata*2/3):end);
z_cur1 = z_cur(1:floor(numdata/3));
z_cur2 = z_cur(floor(numdata/3):floor(numdata*2/3));
z_cur3 = z_cur(floor(numdata*2/3):end);
plot(z_cur1,y_cur1,z_cur2,y_cur2,z_cur3,y_cur3,'LineWidth',1.5); hold on
a=15; % horizontal radius
b=20; % vertical radius
t=-pi:0.01:pi;
x=a*cos(t);
y=b*sin(t);
plot(x,y,'Color','k'); hold on
set(gca,'XColor', 'none','YColor','none','color','none')
% title('stability')
axis equal
    
% exportgraphics(fig,'EMG_stab.png','Resolution',600);
% print(gcf,'-vector','-dsvg',['emg_act.svg']) % svg % pdf
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

    % if strcmp(current_name,'serr_ant_7')
    if 1==1

        figure
        tiledlayout(3,3);
        nexttile([1,3])
        plot(time,lengths_eul,'--','LineWidth',1.5);hold on
        plot(time,lengths_quat,'-.','LineWidth',1.5);hold on
        plot(time,lengths_OS_euler,'g','LineWidth',1.5);hold off
        title('Muscle length and moment arms approximation')
        ylabel('length [m]')
        xlabel('time [s]')
        axis([-inf inf 0 0.25])

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
        Lgnd = legend(['Rotation sequences',newline,'model'],'Quaternion model','OpenSim');
        Lgnd.FontSize = 9;
        Lgnd.Position(1) = 0.65;
        Lgnd.Position(2) = 0.65;
        
        length_numdata = numel(lengths_OS_euler)
        MA_numdata = numel(dLdq_euler)
        LRMS_eul = sqrt(sum((lengths_OS_euler-lengths_eul').^2,'all')/length_numdata)*1000
        LRMS_quat = sqrt(sum((lengths_OS_quat-lengths_quat').^2,'all')/length_numdata)*1000
        max_length_error_eul = max(abs(lengths_OS_euler-lengths_eul'))*1000
        max_length_error_quat = max(abs(lengths_OS_quat-lengths_quat'))*1000
        RRMS_eul = sqrt(sum((dLdq_euler+jacobian_eul(:,1:6)).^2,'all')/MA_numdata)*1000
        RRMS_quat = sqrt(sum((dLdq_quat+JQuatInJEul(:,1:6)).^2,'all')/MA_numdata)*1000
        max_ma_eul = max(max(abs(dLdq_euler+jacobian_eul(:,1:6))))*1000
        max_ma_quat = max(max(abs(dLdq_quat+JQuatInJEul(:,1:6))))*1000

        exportgraphics(fig,'Serr_ant_3_approx.png','Resolution',600);
        
        

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

% function res = rotm2yzy_shoulder(R)
%     cos_z = min(max(R(2,2), -1), 1);
%     z_mag = acosd(cos_z);
% 
%     % Always compute the normal extraction
%     sin_z_normal = sqrt(R(2,1)^2 + R(2,3)^2);
%     z_normal = atan2d(sin_z_normal, cos_z);
%     y1_normal = atan2d(R(3,2), -R(1,2));
%     y2_normal = atan2d(R(2,3), R(2,1));
% 
%     % Always compute the gimbal lock extraction
%     sin_z_gl = (R(2,1) - R(1,2)) / 2;
%     z_gl = atan2d(sin_z_gl, cos_z);
%     y1_gl = 0;
%     y2_gl = 0;
% 
%     % Smooth blending using a sigmoid weight
%     % Below ~5 deg: mostly gimbal lock solution
%     % Above ~10 deg: mostly normal solution
%     threshold = 8;  % center of transition
%     width = 0.01;      % controls sharpness (smaller = sharper)
%     w = 1 / (1 + exp(-(z_mag - threshold) / width));
% 
%     y1 = w * y1_normal + (1 - w) * y1_gl;
%     z  = w * z_normal  + (1 - w) * z_gl;
%     y2 = w * y2_normal + (1 - w) * y2_gl;
% 
%     if z_mag < 7
%         y1 = nan;
%         y2 = nan;
%     end
% 
%     res = [y1, z, y2] * pi/180;
% end

function res = rotm2yzy_shoulder(R)
    % ROTM2YZY_SHOULDER extracts YZY Euler angles from a rotation matrix R.
    % y1 = Plane of Elevation, z = Elevation, y2 = Axial Rotation.
    % All outputs are in degrees.

    % 1. Check the magnitude of the elevation angle to identify gimbal lock
    % In a YZY sequence, R(2,2) represents cos(z). 
    % We clamp it between -1 and 1 to prevent complex numbers from floating-point errors.
    cos_z = min(max(R(2,2), -1), 1);

    % acosd gives the absolute magnitude of the elevation angle [0, 180]
    z_mag = acosd(cos_z); 

    if z_mag < 7.5
        % --- GIMBAL LOCK ZONE (< 10 degrees) ---
        % As requested, ignore plane of elevation and axial rotation
        y1 = nan;
        y2 = nan;

        % When y1 = 0 and y2 = 0, the matrix collapses to a pure Z rotation.
        % R(2,1) becomes sin(z) and R(1,2) becomes -sin(z).
        % We subtract them and divide by 2 for numerical robustness.
        sin_z = (R(2,1) - R(1,2)) / 2;

        % atan2d smoothly calculates the signed angle, allowing negative elevation
        z = atan2d(sin_z, cos_z);

    else
        % --- NORMAL EXTRACTION (>= 10 degrees) ---
        % Standard YZY assumes positive elevation
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


function [lengths, minusdLdq, GHfvecs] = opensim_get_polyvalues(angles, iMus, Dofs, GH_seq)
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
% osimfile = ['Motions','das3_',GH_seq,'_scaled.osim'];
osimfile = ['Motions\3NA\scaled_3NA_GH.osim'];
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

% if nargout>3
GHfvecs = zeros(size(angles,1),3);
% end

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
    if nargout>2
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
        try
            scap_pt = fdarray.get(scap_pt_index);
            fvec = scap_pt.direction();
            
            % transform to the scapular coordinate frame
            SimEn.transform(state,groundbody,fvec,scapulabody,fvec2);
            GHfvecs(istep,1)=fvec2.get(0);
            GHfvecs(istep,2)=fvec2.get(1);
            GHfvecs(istep,3)=fvec2.get(2);
        end
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

function pfvectors = fvectors(musmodel,dof_indeces,angles)
% plot momentarm-angle data
angles = [angles,ones(size(angles,1),1)*120*pi/180];
% choose a subset of "angles" that contains only 100 points
indeces = 1:size(angles,1);
sangles = angles(:,dof_indeces);

% ...or use all angles
%indeces = 1:size(angles,1);
%sangles = angles;

% calculate fvectors from polynomial
pfvectors = zeros(length(indeces),3);
for iframe = 1:length(indeces)
    for i=1:musmodel.xparam_count
        % add this term's contribution to the x vector
        term = musmodel.xcoefs(i);
        for j=1:size(musmodel.xparams,2)
            for k=1:musmodel.xparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end
        pfvectors(iframe,1) = pfvectors(iframe,1) + term;
    end
    for i=1:musmodel.yparam_count
        % add this term's contribution to the y vector
        term = musmodel.ycoefs(i);
        for j=1:size(musmodel.yparams,2)
            for k=1:musmodel.yparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end
        pfvectors(iframe,2) = pfvectors(iframe,2) + term;
    end
    for i=1:musmodel.zparam_count
        % add this term's contribution to the z vector
        term = musmodel.zcoefs(i);
        for j=1:size(musmodel.zparams,2)
            for k=1:musmodel.zparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end
        pfvectors(iframe,3) = pfvectors(iframe,3) + term;
    end
end



end

function [L,pmoment_arms] = momarms(musmodel, dof_indeces, angles)
% ...or use all angles
angles = [angles,ones(size(angles,1),1)*120*pi/180];
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