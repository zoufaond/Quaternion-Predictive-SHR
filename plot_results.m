%% Run this first
addpath Matlab_functions/


%% Plot the computational performance and display the RMS errors - Fig. 2
computational_performance_eul_vs_quat

%% Print the tracking statistics, EMG validations and validation tracking
res = tracking_statistics();
calibration_activation_RMSE
res = RMSE_calibration;

%%  Plot the activations of two elements affected by gimbal lock - Fig. 3

participant = 'par2';
OS_model = ['Motions/',participant,'/OS_model.mat'];
motion_name = 'Elevation';
res_q1 = {['res_euler_',motion_name,'_50'],'euler',motion_name,'YZY',participant,'Euler angles'};
res_q2 = {['res_quat_',motion_name,'_200'],'quat',motion_name,'YZY',participant,'Quaternions'};
plot_IEEE_activations(OS_model,["pect_maj_c_2","delt_scap_6"],0,res_q1,res_q2)

%% PLOT EMG validation - Fig. 4
participant = 'par2';
plot_GH_seq = 'YZY';
motion_name = 'lifting_5kg';
OS_model = ['Motions/',participant,'/OS_model.mat'];
OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
res_q1 = {['res_',motion_name,'_0'],'quat',motion_name,'YZY',participant,'Generic',[245 190 40]/255,1.2};
res_q2 = {['res_',motion_name,'_All_elevations'],'quat',motion_name,'YZY',participant,'Fro+Sca+Sag',[166 20 146]/255,1.2};
res_q3 = {['res_',motion_name,'_drinking_shelf_reaching'],'quat',motion_name,'YZY',participant,'ADL',[0 130 130]/255,1.2};
plot_EMG_optim_IEEE2(['Motions\',participant,'\',motion_name,'\EMG_',participant,'_',motion_name,'.mat'], OS_model,0,'Figures/EMG_validation',res_q1,res_q2,res_q3)

%% MAIN - GENERIC VS CALIBRATION (FRO+SCA+SAG calibration) - Fig. 5
participant = 'par2';
OS_model = ['Motions/',participant,'/OS_model.mat'];
motion_name = 'All_elevations';
OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
res_q1 = {['res_SHR_0'],'quat',motion_name,'YZY',participant,'Calibrated model (Fro+Sca+Sag)',[0 0 1],1.6};
res_q2 = {['res_SHR_6'],'quat',motion_name,'YZY',participant,'Generic model',[128 51 26]/255,1.6};
plot_IEEE_kinematics(OS_struct,0,'Figures/Healthy_gen_vs_calib',res_q1,res_q2);

%% Graded cuff limitation sensitivity analysis - Fig. 6
participant = 'par2';
OS_model = ['Motions/',participant,'/OS_model.mat'];
motion_name = 'All_elevations';
OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
plot_IEEE_kinematics_2DoFs_RCSA(OS_struct,0,'Figures/RC_limitations_SA')

%% HEALTHY vs. RC 0% activation during predictive simulation and glenoid reaction force - Fig. 7, 8
participant = 'par2';
OS_model = ['Motions/',participant,'/OS_model.mat'];
motion_name = 'All_elevations';
OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
res_q1 = {['res_SHR_0'],'quat',motion_name,'YZY',participant,'Healthy','blue',1};
res_q2 = {['res_SHR_RClim0_All_elevations_wGH10'],'quat',motion_name,'YZY',participant,'Supra+infra 0%','red',1};
plot_EMG_optim_IEEE2(['Motions\',participant,'\',motion_name,'\EMG_',participant,'_',motion_name,'.mat'], OS_model,0,'Figures/Act_healthy_vs_RClim',res_q1,res_q2);
plotGHStabilityAnglesIEEE(res_q1,res_q2)


%% Supplementary figures - Fig. S4, S5
cap_colors_RClim = [
    0.00 0.00 1.00;
    0.35 0.00 0.85;
    0.60 0.00 0.65;
    0.82 0.00 0.40;
    1.00 0.00 0.00];

participant = 'par2';
OS_model = ['Motions/',participant,'/OS_model.mat'];
motion_name = 'All_elevations';
OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
res_q1 = {['res_SHR_6'],'quat',motion_name,'YZY',participant,'Healthy Generic',[128 51 26]/255,1.2};
res_q2 = {['res_SHR_0'],'quat',motion_name,'YZY',participant,'Healthy Fro+Sca+Sag',[0 0 1],1.2};
res_q3 = {['res_SHR_Healthy_Elevation_Scabduction_wGH2'],'quat',motion_name,'YZY',participant,'Healthy Fro+Sca',[1.00 0.00 0.00],0.6};
res_q4 = {['res_SHR_Healthy_Elevation_Flexion_wGH2'],'quat',motion_name,'YZY',participant,'Healthy Fro+Sag',        [0.93 0.69 0.13],0.6};
res_q5 = {['res_SHR_Healthy_Scabduction_Flexion_wGH2'],'quat',motion_name,'YZY',participant,' Healthy Sca+Sag',   [0.49 0.18 0.56],0.6};
res_q6 = {['res_SHR_Healthy_drinking_shelf_reaching_wGH2'],'quat',motion_name,'YZY',participant,'Healthy ADLs',               [0.00 0.50 0.50],0.6};
plot_IEEE_kinematics(OS_struct,0,'Figures/Healthy_calibrations_SA',res_q1,res_q2,res_q3,res_q4,res_q5,res_q6);

participant = 'par2';
OS_model = ['Motions/',participant,'/OS_model.mat'];
motion_name = 'All_elevations';
OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
res_q1 = {['res_SHR_RClim0_Generic_wGH10'],'quat',motion_name,'YZY',participant,'Sup+Inf 0% Generic',[128 51 26]/255,1.2};
res_q2 = {['res_SHR_RClim0_All_elevations_wGH10'],'quat',motion_name,'YZY',participant,'Sup+Inf 0% Fro+Sca+Sag',[0 0 1],1.2};
res_q3 = {['res_SHR_RClim0_Elevation_Scabduction_wGH10'],'quat',motion_name,'YZY',participant,'Sup+Inf 0% Fro+Sca',[1.00 0.00 0.00],0.6};
res_q4 = {['res_SHR_RClim0_Elevation_Flexion_wGH10'],'quat',motion_name,'YZY',participant,'Sup+Inf 0% Fro+Sag',        [0.93 0.69 0.13],0.6};
res_q5 = {['res_SHR_RClim0_Scabduction_Flexion_wGH10'],'quat',motion_name,'YZY',participant,' Sup+Inf 0% Sca+Sag',   [0.49 0.18 0.56],0.6};
res_q6 = {['res_SHR_RClim0_drinking_shelf_reaching_wGH10'],'quat',motion_name,'YZY',participant,'Sup+Inf 0% ADLs',               [0.00 0.50 0.50],0.6};
plot_IEEE_kinematics(OS_struct,0,'Figures/RClim0_calibrations_SA',res_q1,res_q2,res_q3,res_q4,res_q5,res_q6);


%% Full figures of graded reduction - not present in the paper

% participant = 'par2';
% OS_model = ['Motions/',participant,'/OS_model.mat'];
% motion_name = 'All_elevations';
% OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
% res_q1 = {['res_SHR_0'],'quat',motion_name,'YZY',participant,'Full capacity wGH=2',cap_colors_RClim(1,:),1};
% res_q2 = {['res_SHR_supra75_All_elevations_wGH2'],'quat',motion_name,'YZY',participant,'Supra=75% wGH=2',cap_colors_RClim(2,:),1};
% res_q3 = {['res_SHR_supra50_All_elevations_wGH2'],'quat',motion_name,'YZY',participant,'Supra=50% wGH=2',cap_colors_RClim(3,:),1};
% res_q4 = {['res_SHR_supra25_All_elevations_wGH2'],'quat',motion_name,'YZY',participant,'Supra=25% wGH=2',cap_colors_RClim(4,:),1};
% res_q5 = {['res_SHR_supra0_All_elevations_wGH4'],'quat',motion_name,'YZY',participant,'Supra=0% wGH=4',cap_colors_RClim(5,:),1};
% plot_IEEE_kinematics(OS_struct,0,'',res_q1,res_q2,res_q3,res_q4,res_q5);

% participant = 'par2';
% OS_model = ['Motions/',participant,'/OS_model.mat'];
% motion_name = 'All_elevations';
% OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
% res_q1 = {['res_SHR_0'],'quat',motion_name,'YZY',participant,'Full capacity wGH=2',cap_colors_RClim(1,:),1};
% res_q2 = {['res_SHR_infra75_All_elevations_wGH2'],'quat',motion_name,'YZY',participant,'Infra=75% wGH=2',cap_colors_RClim(2,:),1};
% res_q3 = {['res_SHR_infra50_All_elevations_wGH2'],'quat',motion_name,'YZY',participant,'Infra=50% wGH=2',cap_colors_RClim(3,:),1};
% res_q4 = {['res_SHR_infra25_All_elevations_wGH4'],'quat',motion_name,'YZY',participant,'Infra=25% wGH=2',cap_colors_RClim(4,:),1};
% res_q5 = {['res_SHR_infra0_All_elevations_wGH6'],'quat',motion_name,'YZY',participant, 'Infra=0% wGH=4',cap_colors_RClim(5,:),1};
% plot_IEEE_kinematics(OS_struct,0,'',res_q1,res_q2,res_q3,res_q4,res_q5);

% participant = 'par2';
% OS_model = ['Motions/',participant,'/OS_model.mat'];
% motion_name = 'All_elevations';
% OS_struct = load(['Motions/',participant,'/',motion_name,'/',motion_name,'.mat']);
% res_q1 = {['res_SHR_0'],'quat',motion_name,'YZY',participant,'Full capacity wGH=2',cap_colors_RClim(1,:),1};
% res_q2 = {['res_SHR_RClim75_All_elevations_wGH2'],'quat',motion_name,'YZY',participant,'Supra+Infra=75% wGH=2',cap_colors_RClim(2,:),1};
% res_q3 = {['res_SHR_RClim50_All_elevations_wGH4'],'quat',motion_name,'YZY',participant,'Supra+Infra=50% wGH=4',cap_colors_RClim(3,:),1};
% res_q4 = {['res_SHR_RClim25_All_elevations_wGH6'],'quat',motion_name,'YZY',participant,'Supra+Infra=25% wGH=6',cap_colors_RClim(4,:),1};
% res_q5 = {['res_SHR_RClim0_All_elevations_wGH10'],'quat',motion_name,'YZY',participant,' Supra+Infra=0% wGH=10',cap_colors_RClim(5,:),1};
% plot_IEEE_kinematics(OS_struct,0,'',res_q1,res_q2,res_q3,res_q4,res_q5);

%% Print joint torque capacity - Supplementary Section SVIII
res = plot_torque_capacity();

function results = plot_torque_capacity(cfg)

if nargin < 1, cfg = struct(); end
def = struct( ...
    'participant', 'par2', ...
    'motion',      'All_elevations', ...
    'root',        'Motions', ...
    'resultFile',  'All_elevations_params_calibration', ...
    'refSet',      2, ...   % index into calibDirs, the set the solution used
    'calibDirs',   {{'', 'All_elevations', 'Elevation_Scabduction', ...
                     'Elevation_Flexion', 'Scabduction_Flexion', ...
                     'drinking_shelf_reaching'}}, ...
    'calibNames',  {{'Uncalibrated','All three planes','Frontal + scapular', ...
                     'Frontal + sagittal','Scapular + sagittal', ...
                     'Drinking + shelf'}}, ...
    'plotDofs',    [1 2 4 5 6 7 8 9], ...
    'dofNames',    {{'SC_y','SC_z', ...
                     'AC_y','AC_z','AC_x' ...
                     'GH_y','GH_z','GH_{yy}'}}, ...
    'colors',      [0.47 0.67 0.19;
                    0.00 0.00 1.00;
                    1.00 0.00 0.00;
                    0.93 0.69 0.13;
                    0.49 0.18 0.56;
                    0.00 0.50 0.50]);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k}) || isempty(cfg.(fn{k})), cfg.(fn{k}) = def.(fn{k}); end
end

nCal = numel(cfg.calibDirs);
nDof = numel(cfg.plotDofs);

% ------------------------------------------------------------ load once
res   = load(fullfile(cfg.root, cfg.participant, cfg.motion, [cfg.resultFile '.mat']));
model = load(fullfile(cfg.root, cfg.participant, 'OS_model.mat'));

time    = res.data.tout;
numdata = numel(time);
traj    = res.data.trajectories;
acts    = res.data.activations;
speeds  = res.data.speeds;

motion_eul = quat2eul_motion(traj, 'YZY');

% ------------------------------------------------------------ preallocate
results.time     = time;
results.names    = cfg.calibNames;
results.dofName  = cfg.dofNames;
results.required = nan(numdata, 11);
results.capPos   = nan(numdata, 11, nCal);
results.capNeg   = nan(numdata, 11, nCal);

for ical = 1:nCal
    muscles = applyScalers(model.model.muscles, cfg, ical);
    [req, cpos, cneg] = momentsForSet(muscles, traj, acts, speeds, ...
                                      motion_eul, numdata);
    if ical == cfg.refSet
        results.required = req;
    end
    results.capPos(:,:,ical) = cpos;
    results.capNeg(:,:,ical) = cneg;
end

results.margin = min(results.capPos - results.required, ...
                     results.required - results.capNeg);

% ------------------------------------------------------------------ plot
fig = figure('Color','w','Units','inches','Position',[1 1 7.16 4.0]);
tiledlayout(2, nDof, 'TileSpacing','compact','Padding','compact');

hLeg = gobjects(nCal,1);
for id = 1:nDof
    d = cfg.plotDofs(id);
    nexttile(id); hold on; box on;

    for ical = 1:nCal
        c = cfg.colors(ical,:);
        h = plot(time, results.capPos(:,d,ical), '-', 'Color', c, 'LineWidth', 0.8);
        plot(time, results.capNeg(:,d,ical), '-', 'Color', c, 'LineWidth', 0.8);
        if id == 1, hLeg(ical) = h; end
    end
    plot(time, results.required(:,d), 'k-', 'LineWidth', 1.6);

    title(cfg.dofNames{id}, 'FontSize', 8, 'FontWeight','normal');
    if id == 1, ylabel('Moment (N m)', 'FontSize', 8); end
    set(gca,'FontSize',8,'XTickLabel',[]);
    xlim([time(1) time(end)]);
end

for id = 1:nDof
    d = cfg.plotDofs(id);
    nexttile(nDof + id); hold on; box on;

    for ical = 1:nCal
        plot(time, results.margin(:,d,ical), '-', ...
             'Color', cfg.colors(ical,:), 'LineWidth', 1.0);
    end
    yline(0, 'k:', 'LineWidth', 0.8);

    xlabel('Time (s)', 'FontSize', 8);
    if id == 1, ylabel('Margin (N m)', 'FontSize', 8); end
    set(gca,'FontSize',8);
    xlim([time(1) time(end)]);
end

lg = legend(hLeg, cfg.calibNames, 'NumColumns', 3);
lg.FontSize = 8;  lg.Box = 'off';  lg.Layout.Tile = 'south';

% ----------------------------------------------------------- print table
fprintf('\nREQUIRED MOMENT, CAPACITY BOUNDS AND FEASIBILITY MARGIN [N m]\n');
fprintf('%s\n', repmat('-', 1, 108));
fprintf('%-22s %-22s %9s %9s %9s %10s\n', ...
    'Coordinate','Parameter set','cap pos','cap neg','margin','infeasible');

for id = 1:nDof
    d = cfg.plotDofs(id);
    fprintf('%-22s required moment: %.2f to %.2f N m\n', ...
        cfg.dofNames{id}, min(results.required(:,d)), max(results.required(:,d)));
    for ical = 1:nCal
        m = results.margin(:,d,ical);
        fprintf('%-22s %-22s %9.2f %9.2f %9.2f %9.1f%%\n', ...
            '', cfg.calibNames{ical}, ...
            max(results.capPos(:,d,ical)), ...
            min(results.capNeg(:,d,ical)), ...
            min(m), 100*mean(m < 0));
    end
    fprintf('\n');
end

fprintf(['The infeasible column gives the fraction of the motion for which\n' ...
         'the required moment lies outside the capacity bounds.\n\n']);
exportgraphics(fig, 'Figures/torque_capacity_calibrations.pdf', 'ContentType','vector', 'BackgroundColor','white');
end

% ========================================================================
function muscles = applyScalers(muscles, cfg, ical)

if isempty(cfg.calibDirs{ical})
    return
end

C = load(fullfile(cfg.root, cfg.participant, cfg.calibDirs{ical}, ...
                  [cfg.calibDirs{ical} '_calibrated_params.mat']));
p = C.calibrated_params;

map = { ...
    'trapclav',  {'trap_clav'}; ...
    'deltclav',  {'delt_clav'}; ...
    'deltscap',  {'delt_scap'}; ...
    'infra',     {'infra'};     ...
    'trapscap',  {'trapscap'};  ...
    'serr',      {'serr_ant'}   };

for ig = 1:size(map,1)
    fF = ['fmax_scaler_'   map{ig,1}];
    fL = ['lceopt_scaler_' map{ig,1}];
    if ~isfield(p, fF) && ~isfield(p, fL), continue, end

    sF = 1; if isfield(p, fF), sF = p.(fF); end
    sL = 1; if isfield(p, fL), sL = p.(fL); end

    for k = 1:numel(muscles)
        if startsWithAny(muscles{k}.name, map{ig,2})
            muscles{k}.fmax   = muscles{k}.fmax   * sF;
            muscles{k}.lceopt = muscles{k}.lceopt * sL;
        end
    end
end
end

function tf = startsWithAny(name, prefixes)
tf = false;
for i = 1:numel(prefixes)
    if strncmp(name, prefixes{i}, numel(prefixes{i})), tf = true; return, end
end
end

%
function [required, capPos, capNeg] = momentsForSet(muscles, traj, acts, ...
                                                    speeds, motion_eul, numdata)

nmus     = numel(muscles);
required = zeros(numdata, 11);
capPos   = zeros(numdata, 11);
capNeg   = zeros(numdata, 11);

alljoints = {'YZX','YZX','YZY'};

for iframe = 1:numdata

    R       = zeros(nmus, 11);
    Fact    = zeros(nmus, 1);
    Fmaxact = zeros(nmus, 1);
    Fpas    = zeros(nmus, 1);

    for imus = 1:nmus

        m = muscles{imus};
        dof_indeces = m.dof_indeces - 3;

        motion_quat = traj(iframe, [2:4, 6:8, 10:12, 13]);
        [len, jac]  = momarms(m.Quaternion, dof_indeces, motion_quat);

        jacobian_quat = zeros(1,11);
        jacobian_quat(dof_indeces) = jac;

        for j = 1:3
            ind3 = (j-1)*3 + (1:3);
            ind4 = (j-1)*4 + (1:4);
            JQuatInSpat = invJtrans(traj(iframe,ind4)) * jacobian_quat(ind3)';
            R(imus,ind3) = GeomJ(motion_eul(iframe,ind3), alljoints{j}) * JQuatInSpat;
        end
        R(imus,10) = jacobian_quat(10);
        R(imus,11) = jacobian_quat(11);

        Fact(imus)    = muscle_force(acts(iframe,imus), len, 0, ...
                                     m.fmax, m.lceopt, m.lslack, m.vmax);
        Fmaxact(imus) = muscle_force(1, len, 0, ...
                                     m.fmax, m.lceopt, m.lslack, m.vmax);
        Fpas(imus)    = muscle_force(0, len, 0, ...
                                     m.fmax, m.lceopt, m.lslack, m.vmax);
    end
    R = -R;
    required(iframe,:) = R' * Fact;

    for d = 1:11
        pos = R(:,d) > 0;
        F = Fpas;  F(pos) = Fmaxact(pos);
        capPos(iframe,d) = R(:,d)' * F;

        neg = R(:,d) < 0;
        F = Fpas;  F(neg) = Fmaxact(neg);
        capNeg(iframe,d) = R(:,d)' * F;
    end
end
end


function force = muscle_force(act, lmt, vce, fmax, lceopt, lslack, vmax)

lm = lmt - lslack;

f_gauss = 0.25;
kpe = 5;
epsm0 = 0.6;
fpe = (exp(kpe*(lm / lceopt - 1)/epsm0)-1)/(exp(kpe)-1);
flce = (exp(-(lm / lceopt - 1)^2 / f_gauss));

d1 = -0.318;
d2 = -8.149;
d3 = -0.374;
d4 = 0.886;
vmax_norm = vmax * lceopt;
vnorm = vce/vmax_norm;

fvce = d1 * log(d2 * vnorm + d3 + sqrt((d2 * vnorm + d3)^2 + 1)) + d4;

force = (flce * act * fvce +  fpe) * fmax;
    
end

function [L,pmoment_arms] = momarms(musmodel, dof_indeces, angles)

angles = [angles, 120*pi/180];
indeces = 1:size(angles,1);
sangles = angles(:,dof_indeces);

% calculate moment arms from polynomial
pmoment_arms = zeros(length(indeces),length(dof_indeces));
% disp(musmodel.lparam_count)
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
    figure('Color','w',"Units","inches",'Position',[1 1 3.5 1.8])
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
            xlabel(['% of motion'])
            ylabel('Activation [-]')
            text(0.5, -0.35, ['(',alphabet{imus},')'], 'Units', 'normalized', ...
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

    exportgraphics(gcf,'Figures/IEEE_GL_activation.pdf', 'ContentType','vector', 'BackgroundColor','white');

end



function calibration_activation_RMSE()
motions = {'driving','lifting_5kg'};
calibrations = {'0','All_elevations','Elevation_Scabduction','Elevation_Flexion','Scabduction_Flexion','drinking_shelf_reaching'};
participant = 'par2';
OS_model = 'Motions/par2/OS_model_prediction.mat';
EMG_muscles = {'Infrasp','UpperTrap','Serrupper','IntermediateDelt','PosteriorDelt','AnteriorDelt','MiddleTrap-TS','MiddleTrap-Spin'};
calib_file_names = {'infra','trapclav','serr','deltscap','deltscap','deltclav','trapscap','trapscap'};
model_names = {["infra_2","infra_3","infra_4"],["trap_clav_1"],["serr_ant_2","serr_ant_3","serr_ant_4"],["delt_scap11","delt_scap10","delt_scap_9","delt_scap_8"],["delt_scap_3","delt_scap_4","delt_scap_5"],["delt_clav_1","delt_clav_2"],["trap_scap_5","trap_scap_6"],["trap_scap_1","trap_scap_2"]};


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
    R_h_rot_i = R_y(13*pi/180)' * R_z(-6.5*pi/180)' * [Rx_h(i);Ry_h(i);Rz_h(i);1];
    R_r_rot_i = R_y(13*pi/180)' * R_z(-6.5*pi/180)' * [Rx_r(i);Ry_r(i);Rz_r(i);1];
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
title(['Supra+infra 0%']);

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
legend({'Healthy','Supra+infra 0%'},'Position',[0.45,0.38,1,1],'Box','off');
text(0.5, -0.35, ['(',alphabet{3},')'], 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 8);

set(gca,'FontSize',8,'LineWidth',0.2);

% ========= Colorbar =========
cb.Layout.Tile = 'east';
cb.Label.String = 'Time (s)';
cb.FontSize = 8;

exportgraphics(gcf,'Figures/IEEE_GH_stability.pdf','ContentType','vector', 'BackgroundColor','white');


end

function plot_IEEE_kinematics(kinematics,save,filename, varargin)

labels = { ...
    'Clavicular protraction/retraction', 'Clavicular elevation','Clavicular axial rotation', ...
    'Thoracoscapular protraction/retraction','Thoracoscapular upward/downward rotation','Thoracoscapular anterior/posterior tilting', ...
    'Thoracohumeral plane of elevation','Thoracohumeral elevation','Thoracohumeral axial rotation'};

alphabet = {'a','b','c','d','e','f','g','h','i'};

nsim = numel(varargin);
% lw   = 1.6;

% ---------- Figure setup ----------
fig = figure('Color','w','Units','inches','Position',[1 1 7.16 3.5]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

% Legend handles. First entry is the experimental trajectory.
hLeg  = gobjects(nsim+1,1);
names = cell(nsim+1,1);
names{1} = 'Experimental data';

for i = 1:9
    nexttile; hold on; box on;

    for isim = 1:nsim

        cfg    = varargin{isim};
        fname  = cfg{1};
        motion = cfg{3};
        par    = cfg{5};
        lname  = cfg{6};
        lcolor = cfg{7};
        lw = cfg{8};

        S = load(fullfile('Motions', par, motion, [fname '.mat']));
        t = S.data.tout;

        kin_loc = quat2eul_motion(S.data.trajectories,'YZY');
        kin_sim = create_objective_traj_eul(kin_loc,'YZY',1);

        % Experimental trajectory, plotted once per panel from the first
        % simulation time base
        if isim == 1 && i ~= 3
            kin_exp = interp1(kinematics.mot_struct.time, ...
                              kinematics.mot_struct.euler, t, 'spline');
            kin_exp = create_objective_traj_eul(kin_exp,'YZY',0);
            h = plot(t, rad2deg(kin_exp(:,i)), ...
                     'Color',[0.5 0.5 0.5],'LineWidth',1.2);
            if i == 1, hLeg(1) = h; end
        end

        if i == 2
            if isim == 1
                fprintf('IK - Max clavicular depression = %2.2f\n',rad2deg(min(kin_exp(:,i))));
            end
            fprintf('%s - Max clavicular depression = %2.2f\n',lname,rad2deg(min(kin_sim(:,i))));
        end
        if i == 5
            if isim == 1
                fprintf('IK - Max thoracoscapular downward rotation = %2.2f\n',rad2deg(min(kin_exp(:,i))));
            end
            fprintf('%s - Max thoracoscapular downward rotation = %2.2f\n',lname,rad2deg(min(kin_sim(:,i))));
        end

        % Simulated trajectory. Panels 7 and 9 hold the degrees of freedom
        % that are interpolated through the low elevation region, shown as
        % a dashed line underneath the solid trace.
        if i == 7 || i == 9
            kin_interp = fillmissing(kin_sim,'linear');
            plot(t, rad2deg(kin_interp(:,i)), ...
                 '--','Color',lcolor,'LineWidth',0.8);
            h = plot(t, rad2deg(kin_sim(:,i)), ...
                     '-','Color',lcolor,'LineWidth',lw);
        else
            h = plot(t, rad2deg(kin_sim(:,i)), ...
                     '-','Color',lcolor,'LineWidth',lw);
        end

        if i == 1
            hLeg(isim+1) = h;
            names{isim+1} = lname;
        end

        % Coordination values, printed once per simulation
        if i == 1
            GH  = rotyzy(S.data.trajectories(:,9:12));
            fro = [1,41];
            sca = [101,135];
            sag = [204,236];

            [S_fro,i_fro] = compute_SCHR(kin_sim(fro(1):fro(2),8)*180/pi, ...
                                         kin_sim(fro(1):fro(2),5)*180/pi, ...
                                         GH(fro(1):fro(2))*180/pi);
            [S_sca,i_sca] = compute_SCHR(kin_sim(sca(1):sca(2),8)*180/pi, ...
                                         kin_sim(sca(1):sca(2),5)*180/pi, ...
                                         GH(sca(1):sca(2))*180/pi);
            [S_sag,i_sag] = compute_SCHR(kin_sim(sag(1):sag(2),8)*180/pi, ...
                                         kin_sim(sag(1):sag(2),5)*180/pi, ...
                                         GH(sag(1):sag(2))*180/pi);

            print_SCHR(['Frontal / '  lname], S_fro, i_fro);
            print_SCHR(['Scapular / ' lname], S_sca, i_sca);
            print_SCHR(['Sagittal / ' lname], S_sag, i_sag);
        end
    end

    % ---------- Panel formatting ----------
    title(labels{i},'FontSize',8,'FontWeight','normal');
    xlabel(['Time (s)']);
    ylabel('Angle (deg)');
    xlim([t(1) t(end)]);
    set(gca,'FontSize',7,'LineWidth',0.2);

    text(0.01,1.04,['(' alphabet{i} ')'], ...
             'Units','normalized','FontSize',8, ...
             'VerticalAlignment','top','HorizontalAlignment','left');
end

% ---------- Predicted and tracked annotations ----------
annotation('line',[0.96 0.96],[0.41 0.97],'Color',[0.6 0.6 0.6],'LineWidth',1);
annotation('line',[0.95 0.96],[0.97 0.97],'Color',[0.6 0.6 0.6],'LineWidth',1);
annotation('line',[0.95 0.96],[0.41 0.41],'Color',[0.6 0.6 0.6],'LineWidth',1);
annotation('textbox',[0.96 0.67 0 0],'String','Predicted','EdgeColor','none', ...
           'Rotation',90,'FontSize',15,'FontAngle','italic', ...
           'HorizontalAlignment','center');

annotation('line',[0.96 0.96],[0.10 0.35],'Color',[0.6 0.6 0.6],'LineWidth',1);
annotation('line',[0.95 0.96],[0.10 0.10],'Color',[0.6 0.6 0.6],'LineWidth',1);
annotation('line',[0.95 0.96],[0.35 0.35],'Color',[0.6 0.6 0.6],'LineWidth',1);
annotation('textbox',[0.96 0.235 0 0],'String','Tracked','EdgeColor','none', ...
           'Rotation',90,'FontSize',15,'FontAngle','italic', ...
           'HorizontalAlignment','center');

% ---------- Legend ----------
valid = isgraphics(hLeg);
lg = legend(hLeg(valid), names(valid), 'numColumns', 4);
lg.FontSize    = 8;
lg.Box         = 'off';
lg.Layout.Tile = 'south';
lg.Orientation = 'horizontal';

if save == 1
    exportgraphics(fig, [filename '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
end
end


function plot_IEEE_kinematics_2DoFs_RCSA(kinematics,save,filename)

% ---------- Labels ----------
labels = { ...
    ['Clavicular',newline,'elevation (deg)'], ...
    ['Thoracoscapular',newline,'upward',newline,'rotation (deg)']};

case_titles = {'Supraspinatus', 'Infraspinatus', 'Supraspinatus + infraspinatus'};

RC_cases    = {'supra','infra','RClim'};
limitations = {'100','75','50','25','0'};
wGHs = { {'2','2','2','2','4'}, ...
         {'2','2','2','4','6'}, ...
         {'2','2','4','6','10'} };

cap_colors = [
    0.00 0.00 1.00;
    0.35 0.00 0.85;
    0.60 0.00 0.65;
    0.82 0.00 0.40;
    1.00 0.00 0.00];

DoFs     = [2, 5];
alphabet = {'a','b','c','d','e','f'};
lw       = 1.6;

% ---------- Figure setup ----------
fig = figure('Color','w','Units','inches','Position',[1 1 7.16 3.5/3*2]);
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

hLeg    = gobjects(6,1);   % one handle per legend entry
legDone = false;
ipanel  = 0;

for idof = 1:2
    for ilim = 1:3

        ipanel = ipanel + 1;
        nexttile; hold on; box on;

        wGHs_cur = wGHs{ilim};

        for isim = 1:5

            if isim == 1
                S = load('Motions\par2\All_elevations\res_SHR_0.mat');
            else
                S = load(['Motions\par2\All_elevations\res_SHR_', ...
                          RC_cases{ilim}, limitations{isim}, ...
                          '_All_elevations_wGH', wGHs_cur{isim}, '.mat']);
            end

            t = S.data.tout;

            kin_loc = quat2eul_motion(S.data.trajectories,'YZY');
            kin_sim = create_objective_traj_eul(kin_loc,'YZY',1);

            DoF = DoFs(idof);

            if isim == 1
                kin_exp = interp1(kinematics.mot_struct.time, ...
                                  kinematics.mot_struct.euler, t, 'spline');
                kin_exp = create_objective_traj_eul(kin_exp,'YZY',0);

                h = plot(t, rad2deg(kin_exp(:,DoF)), ...
                         'Color',[0.5 0.5 0.5],'LineWidth',1.2);
                if ~legDone, hLeg(1) = h; end

                h = plot(t, rad2deg(kin_sim(:,DoF)), ...
                         'Color',cap_colors(isim,:),'LineWidth',1.8);
                if ~legDone, hLeg(2) = h; end
            else
                h = plot(t, rad2deg(kin_sim(:,DoF)), ...
                         'Color',cap_colors(isim,:),'LineWidth',1.2);
                if ~legDone, hLeg(isim+1) = h; end
            end
        end

        legDone = true;

        % ---------- Panel formatting ----------
        set(gca,'FontSize',8,'LineWidth',0.5,'TickDir','out','Layer','top');
        xlim([t(1) t(end)]);

        if idof == 1
            title(case_titles{ilim},'FontSize',8,'FontWeight','normal');
            ylim([0,20])
        end
        if idof == 2
            xlabel('Time (s)','FontSize',8);
            ylim([0,50])
        else
            set(gca,'XTickLabel',[]);
        end
        if ilim == 1
            ylabel(labels{idof},'FontSize',8);
        end

        % Panel letter
        text(0.02,0.96,['(' alphabet{ipanel} ')'], ...
             'Units','normalized','FontSize',8, ...
             'VerticalAlignment','top','HorizontalAlignment','left');

        % Weight annotation, one line per panel
        text(0.98,1.02, ...
             ['w_{GH} = ' strjoin(wGHs_cur,', ')], ...
             'Units','normalized','FontSize',7,'Color',[0.35 0.35 0.35], ...
             'VerticalAlignment','top','HorizontalAlignment','right');
    end
end

ax = findall(gcf,'Type','axes');
ax = flipud(ax);
for irow = 1:2
    idx = (irow-1)*3 + (1:3);
    yl  = cell2mat(get(ax(idx),'YLim'));
    set(ax(idx),'YLim',[min(yl(:,1)) max(yl(:,2))]);
end

% ---------- Legend ----------
names = {'Experimental trajectory', ...
         'Full capacity (100%)', ...
         '75% capacity', ...
         '50% capacity', ...
         '25% capacity', ...
         'No capacity (0%)'};

lg = legend(hLeg, names, 'NumColumns', 3);
lg.FontSize    = 8;
lg.Box         = 'off';
lg.Layout.Tile = 'south';

if save == 1
    exportgraphics(fig, [filename '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
end

end


function print_SCHR(name, SCHR, info)
fprintf('\n%s\n', name);
fprintf('%-11s %4s %8s %6s %9s %8s %7s %7s %7s %7s %7s %7s\n', ...
    'phase','n','slope','R2','endpoint','range','dTH','dTS','dGH','revTS','maxTS','minTS');
for i = 1:numel(SCHR)
    fprintf('%4.0f-%-6.0f %4d %8.2f %6.3f %9.2f %8.2f %7.1f %7.1f %7.1f %7.2f %7.2f %7.2f\n', ...
        info.phase(i,1), info.phase(i,2), info.n(i), info.slope(i), ...
        info.R2(i), info.endpoint(i), info.range(i), ...
        info.netTH(i), info.netTS(i), info.netGH(i), info.fracRevTS(i), info.maxTS(i), info.minTS(i));
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


function [SCHR, info] = compute_SCHR(TH, TS, GH, phaseEdges)

TH = TH(:); TS = TS(:); GH = GH(:);

if ~isequal(numel(TH), numel(TS), numel(GH))
    error('compute_SCHR:sizeMismatch', ...
        'TH (%d), TS (%d) and GH (%d) must be the same length.', ...
        numel(TH), numel(TS), numel(GH));
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
              'netGH', z, 'spanTS', z, 'fracRevTS', z, 'maxTS', z, 'minTS',z);

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
    info.minTS(i) = ts(1);

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

function plot_EMG_optim_IEEE2(EMG_struct, OS_model, save, filename, varargin)

fig = figure('Color','w','Units','inches','Position',[1 1 3.5 2.4]);
alphabet = {'a','b','c','d'};

tt = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% ---------- Line styles ----------
lw = 1.2;
EMG_style     = {'Color',[0.5 0.5 0.5],'LineWidth',0.8};
results_style = { {'-','Color',[245 190 40]/255,'LineWidth',lw}, ...
                  {'-','Color',[166 20 146]/255,'LineWidth',lw}, ...
                  {'-','Color',[0 130 130]/255,'LineWidth',lw} };

EMG_muscles = {'IntermediateDelt','Infrasp','UpperTrap','Serrupper'};
model_names = { ["delt_scap_8","delt_scap_9","delt_scap_10","delt_scap11"], ...
                ["infra_2","infra_3","infra_4"], ...
                ["trap_clav_1"], ...
                ["serr_ant_2","serr_ant_3","serr_ant_4"] };
figure_names = {'Lateral deltoid','Infraspinatus','Clavicular trapezius','Serratus anterior'};

num_res  = numel(varargin);
emg_data = load(EMG_struct);
model    = load(OS_model);
muscles  = model.model.muscles;

muscle_names = cell(1, numel(muscles));
for i = 1:numel(muscles)
    muscle_names{i} = muscles{i}.osim_name;
end

mus_index = cell(1, numel(model_names));
for igroup = 1:numel(model_names)
    current_group = model_names{igroup};
    idx = [];
    for k = 1:numel(current_group)
        idx = [idx, find(strcmp(muscle_names, current_group(k)))]; %#ok<AGROW>
    end
    mus_index{igroup} = int16(idx);
end

ax = gobjects(4,1);

for i = 1:numel(EMG_muscles)

    ax(i) = nexttile; hold on; box on;
    legend_names = {'EMG'};

    for ires = 1:num_res

        r = varargin{ires};
        result = load(fullfile('Motions', r{5}, r{3}, [r{1} '.mat']));

        idx = mus_index{i};
        excitations = zeros(size(result.data.excitations(:,1)));
        for k = 1:numel(idx)
            excitations = excitations + result.data.excitations(:, idx(k));
        end
        excitations = excitations / numel(idx);

        tout = result.data.tout;
        time = linspace(0, tout(end), numel(tout));

        legend_names{end+1} = r{6}; %#ok<AGROW>

        % Measured EMG and the shaded pauses, drawn once per panel
        if ires == 1
            for imot = 1:6
                try
                    current_emg = emg_data.data.(['num_' num2str(imot)]).(EMG_muscles{i});
                end
            end
            time_emg = linspace(0, tout(end), numel(current_emg));
            plot(time_emg, current_emg, EMG_style{:});

            isNaN    = isnan(current_emg);
            d        = diff([false; isNaN(:); false]);
            nanStart = find(d ==  1);
            nanEnd   = find(d == -1) - 1;
            for k = 1:numel(nanStart)
                p = patch([time_emg(nanStart(k)) time_emg(nanEnd(k)) ...
                           time_emg(nanEnd(k))   time_emg(nanStart(k))], ...
                          [0 0 1 1], [0.8 0.8 0.8], ...
                          'FaceAlpha', 0.15, 'EdgeColor', 'none');
                % Keep the patch out of the legend and out of the automatic
                % y limits, which would otherwise be driven to 1 by its
                % height rather than by the signals.
                p.Annotation.LegendInformation.IconDisplayStyle = 'off';
                set(p, 'HandleVisibility', 'off');
            end
        end
        icolor = r{7};
        ilw = r{8};
        % cur_style = results_style{ires};
        plot(time, excitations,'Color',icolor,'Linewidth',ilw);
    end

    xlim([0 tout(end)]);
    set(gca, 'FontSize', 7, 'LineWidth', 0.2, 'TickDir', 'out');

    % Panel letter and muscle name inside the axes, so no title row is used
    text(0.03, 0.95, ['(' alphabet{i} ') ' figure_names{i}], ...
        'Units', 'normalized', 'FontSize', 8, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    % X tick labels only on the bottom row. Y tick labels kept on every
    % panel, because the y scales differ between muscles.
    if i <= 2, set(gca, 'XTickLabel', []); end
    if i >= 2, xlabel('Time (s)', 'FontSize',8); end
    if i == 1 || i == 3, ylabel('Excitation (-)','FontSize',8); end
end

% Y limits set per panel, from the plotted signals only. The patches are
% excluded because their height is fixed and would set every panel to 1.
for i = 1:4
    ymax = 0;
    h = findobj(ax(i), 'Type', 'line');
    for k = 1:numel(h)
        y = get(h(k), 'YData');
        ymax = max(ymax, max(y(~isnan(y))));
    end
    if ymax == 0, ymax = 1; end
    set(ax(i), 'YLim', [0 1.25*ymax]);

    % Patches are drawn with a fixed height of 1, so rescale to this axis
    hp = findobj(ax(i), 'Type', 'patch');
    for k = 1:numel(hp)
        yd = get(hp(k), 'YData');
        yd(yd > 0) = 1.05*ymax;
        set(hp(k), 'YData', yd);
    end
end

lg = legend(legend_names);
lg.Layout.Tile   = 'south';
% lg.Orientation   = 'horizontal';
lg.NumColumns    = 4;
lg.FontSize      = 8;
lg.Box           = 'off';
lg.ItemTokenSize = 8;
% drawnow
% lg.Units = 'normalized';
% p = lg.Position;
% lg.Position = [p(1), p(2) + 0.55, p(3), p(4)];

if save == 1
    exportgraphics(fig, [filename '.pdf'], ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
end

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

figure('Color','w','Units','inches','Position',[1 1 3.5 2]);
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

exportgraphics(gcf,'Figures/comp_perf.pdf', 'ContentType','vector', 'BackgroundColor','white')

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

    if z_mag < 20
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



function t = tcrit(df)
tab = [12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228 ...
        2.201 2.179 2.160 2.145 2.131 2.120 2.110 2.101 2.093 2.086];
t = nan(size(df));
for i = 1:numel(df)
    if df(i) < 1,       t(i) = NaN;
    elseif df(i) <= 20, t(i) = tab(df(i));
    else,               t(i) = 1.96;
    end
    % t(i)
end
end


function results = RMSE_calibration(cfg)

if nargin < 1, cfg = struct(); end
def = struct( ...
    'motions',      {{'driving','lifting_5kg'}}, ...
    'calibrations', {{'0','All_elevations','Elevation_Scabduction', ...
                      'Elevation_Flexion','Scabduction_Flexion', ...
                      'drinking_shelf_reaching'}}, ...
    'calibNames',   {{'Uncalibrated','All three planes','Frontal + scapular', ...
                      'Frontal + sagittal','Scapular + sagittal', ...
                      'Drinking + shelf'}}, ...
    'participant',  'par2', ...
    'root',         'Motions', ...
    'nSamp',        100, ...
    'seq',          'YZY', ...
    'globFlag',     0, ...
    'dofCols',      [1 2 4 5 6 7 8 9], ...
    'dofNames',     {{'Clav protraction','Clav elevation', ...
                      'Scap protraction','Scap upward rot','Scap tilt', ...
                      'TH plane of elev','TH elevation','TH axial rot'}}, ...
    'segCols',      {{4:6, 7:9}}, ...
    'segNames',     {{'Thoracoscapular','Thoracohumeral'}}, ...
    'segSeq',       {{'YZX','YZY'}}, ...
    'intrinsic',    true);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k}) || isempty(cfg.(fn{k})), cfg.(fn{k}) = def.(fn{k}); end
end

nMot   = numel(cfg.motions);
nCal   = numel(cfg.calibrations);
nDOF   = numel(cfg.dofCols);
nSeg   = numel(cfg.segCols);

% ------------------------------------------------------------- preallocate
% Third dimension is the calibration set, so trials can be averaged later.
D = nan(nMot, nDOF, nCal);
G = nan(nMot, nSeg, nCal);

results.trial.rmse   = D;
results.trial.maxabs = D;
results.trial.bias   = D;
results.trial.rom    = D;
results.trial.geo    = G;
results.trial.motion = cfg.motions;
results.calibNames   = cfg.calibNames;

% ------------------------------------------------------------------- loop
for imot = 1:nMot

    motion = cfg.motions{imot};

    OS = load(fullfile(cfg.root, cfg.participant, motion, [motion '.mat']));
    IK_glob = create_objective_traj_eul(OS.mot_struct.euler, cfg.seq, cfg.globFlag);
    assert(size(IK_glob,1) == cfg.nSamp, ...
        'IK reference for %s has %d rows, expected %d.', ...
        motion, size(IK_glob,1), cfg.nSamp);

    refDeg = rad2deg(IK_glob(:, cfg.dofCols));

    for ical = 1:nCal

        fpath = fullfile(cfg.root, cfg.participant, motion, ...
                         ['res_' motion '_' cfg.calibrations{ical} '.mat']);

        if ~isfile(fpath)
            warning('Missing file: %s', fpath);
            continue
        end

        S    = load(fpath);
        t    = S.data.tout;
        traj = interp1(t, S.data.trajectories, linspace(0, t(end), cfg.nSamp));
        traj = quat2eul_motion(traj, cfg.seq);

        simGlob = create_objective_traj_eul(traj, cfg.seq, cfg.globFlag);
        simDeg  = rad2deg(simGlob(:, cfg.dofCols));

        % --- Euler component error -------------------------------------
        err = simDeg - refDeg;
        results.trial.rmse(imot,:,ical)   = sqrt(mean(err.^2, 1));
        results.trial.maxabs(imot,:,ical) = max(abs(err), [], 1);
        results.trial.bias(imot,:,ical)   = mean(err, 1);
        results.trial.rom(imot,:,ical)    = max(refDeg,[],1) - min(refDeg,[],1);

        % --- geodesic orientation error ---------------------------------
        for iseg = 1:nSeg
            g = geodesicError(simGlob(:, cfg.segCols{iseg}), ...
                              IK_glob(:, cfg.segCols{iseg}), ...
                              cfg.segSeq{iseg});
            results.trial.geo(imot,iseg,ical) = sqrt(mean(g.^2));
        end
    end
end

% -------------------------------------------------------------- aggregate
% Mean across the validation trials, one value per calibration set.
results.mean.rmse   = squeeze(mean(results.trial.rmse,   1, 'omitnan')).';
results.mean.maxabs = squeeze(max( results.trial.maxabs, [], 1)).';
results.mean.bias   = squeeze(mean(results.trial.bias,   1, 'omitnan')).';
results.mean.rom    = squeeze(mean(results.trial.rom,    1, 'omitnan')).';
results.mean.geo    = squeeze(mean(results.trial.geo,    1, 'omitnan')).';

results.cfg = cfg;

% ------------------------------------------------------------------ report
printBlock('GEODESIC ORIENTATION ERROR [deg]', ...
           cfg.calibNames, cfg.segNames, results.mean.geo);

printBlock('EULER COMPONENT RMSE [deg]', ...
           cfg.calibNames, cfg.dofNames, results.mean.rmse);

printBlock('BIAS [deg]', ...
           cfg.calibNames, cfg.dofNames, results.mean.bias);

fprintf('\nAveraged over %d validation trials: %s\n', ...
        nMot, strjoin(cfg.motions, ', '));

end

% ========================================================================
function ang = geodesicError(eulSim, eulRef, seq)
n = size(eulSim,1);
ang = nan(n,1);
for i = 1:n
    Rs = eul2rotm(eulSim(i,:), seq);
    Rr = eul2rotm(eulRef(i,:), seq);
    axang  = rotm2axang(Rs * Rr.');
    ang(i) = rad2deg(axang(4));
end
end


% ========================================================================
function printBlock(hdr, rowNames, colNames, M)
% M is nRows by nCols
fprintf('\n%s\n%s\n', hdr, repmat('-', 1, 20 + 11*numel(colNames)));
fprintf('%-20s', 'Calibration');
for j = 1:numel(colNames)
    fprintf('%11s', shorten(colNames{j}));
end
fprintf('\n');
for i = 1:size(M,1)
    fprintf('%-20s', rowNames{i});
    for j = 1:size(M,2)
        fprintf('%11.2f', M(i,j));
    end
    fprintf('\n');
end
end

function s = shorten(s)
if numel(s) > 10, s = s(1:10); end
end

function results = tracking_statistics(cfg)

if nargin < 1, cfg = struct(); end
def = struct( ...
    'motions',      {{'Elevation','Scabduction','Flexion'}}, ...
    'motionNames',  {{'Frontal','Scapular','Sagittal'}}, ...
    'participants', {{'par1','par2','par3'}}, ...
    'rotType',      {{'euler','quat'}}, ...
    'rotWeights',   {{'50','200'}}, ...
    'root',         'Motions', ...
    'nSamp',        100, ...
    'seq',          'YZY', ...
    'globFlag',     0, ...
    'dofCols',      [1 2 4 5 6 7 8 9], ...
    'dofNames',     {{'Clav protraction','Clav elevation', ...
                      'Scap protraction','Scap upward rot','Scap tilt', ...
                      'TH plane of elev','TH elevation','TH axial rot'}}, ...
    'segCols',      {{4:6, 7:9}}, ...
    'segNames',     {{'Thoracoscapular','Thoracohumeral'}}, ...
    'segSeq',       {{'YZX','YZY'}}, ...
    'intrinsic',    true);
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
blank  = struct('rmse',Z,'maxabs',Z,'bias',Z,'rom',Z);
blankG = struct('rmse',Zs,'maxabs',Zs);

results.trial.euler = blank;   results.trial.quat = blank;
results.geo.euler   = blankG;  results.geo.quat   = blankG;
results.trial.motionIdx = nan(nTr,1);
results.trial.label     = cell(nTr,1);

iTr = 0;

% ------------------------------------------------------------------- loop
for ipar = 1:nPar
    for imot = 1:nMot

        iTr = iTr + 1;
        results.trial.motionIdx(iTr) = imot;
        results.trial.label{iTr} = sprintf('%s / %s', ...
            cfg.participants{ipar}, cfg.motionNames{imot});

        OS = load(fullfile(cfg.root, cfg.participants{ipar}, ...
                           cfg.motions{imot}, [cfg.motions{imot} '.mat']));
        IK_glob = create_objective_traj_eul(OS.mot_struct.euler, cfg.seq, cfg.globFlag);
        refDeg  = rad2deg(IK_glob(:, cfg.dofCols));

        for irot = 1:2

            S = load(fullfile(cfg.root, cfg.participants{ipar}, cfg.motions{imot}, ...
                 ['res_' cfg.rotType{irot} '_' cfg.motions{imot} '_' ...
                  cfg.rotWeights{irot} '.mat']));

            t    = S.data.tout;
            traj = interp1(t, S.data.trajectories, linspace(0, t(end), cfg.nSamp));
            if strcmp(cfg.rotType{irot}, 'quat')
                traj = quat2eul_motion(traj, cfg.seq);
            end
            simGlob = create_objective_traj_eul(traj, cfg.seq, cfg.globFlag);
            simDeg  = rad2deg(simGlob(:, cfg.dofCols));

            f   = cfg.rotType{irot};
            err = simDeg - refDeg;

            results.trial.(f).rmse(iTr,:)   = sqrt(mean(err.^2, 1));
            results.trial.(f).maxabs(iTr,:) = max(abs(err), [], 1);
            results.trial.(f).bias(iTr,:)   = mean(err, 1);
            results.trial.(f).rom(iTr,:)    = max(refDeg,[],1) - min(refDeg,[],1);

            for iseg = 1:nSeg
                g = geodesicError(simGlob(:, cfg.segCols{iseg}), ...
                                  IK_glob(:, cfg.segCols{iseg}), ...
                                  cfg.segSeq{iseg});
                results.geo.(f).rmse(iTr,iseg)   = sqrt(mean(g.^2));
                results.geo.(f).maxabs(iTr,iseg) = max(g);
            end
        end
    end
end

results.cfg = cfg;

% ================================================================ report

% ---- (1) per component, all nine simulations ---------------------------
fprintf('\nPER COMPONENT, ALL %d SIMULATIONS (%d participants x %d tasks)\n', ...
        nTr, nPar, nMot);
fprintf('%s\n', repmat('-',1,116));
fprintf('%-18s %8s %16s %16s %20s %8s %8s\n', ...
    'Coordinate','ROM','RMSE Euler','RMSE Quat','Difference [95% CI]', ...
    'MaxAE E','MaxAE Q');
allIdx = true(nTr,1);
printComponents(cfg.dofNames, results, allIdx);

fprintf('\n%-18s %10s %10s %10s %10s\n', ...
    'Coordinate','bias E','bias Q','sd(bias) E','sd(bias) Q');
for j = 1:nDOF
    fprintf('%-18s %10.2f %10.2f %10.2f %10.2f\n', cfg.dofNames{j}, ...
        mean(results.trial.euler.bias(:,j)), mean(results.trial.quat.bias(:,j)), ...
        std(results.trial.euler.bias(:,j)),  std(results.trial.quat.bias(:,j)));
end

% ---- (2) per component, split by task ---------------------------------
for imot = 1:nMot
    fprintf('\nPER COMPONENT, %s PLANE ELEVATION (%d participants)\n', ...
            upper(cfg.motionNames{imot}), nPar);
    fprintf('%s\n', repmat('-',1,116));
    fprintf('%-18s %8s %16s %16s %20s %8s %8s\n', ...
        'Coordinate','ROM','RMSE Euler','RMSE Quat','Difference [95% CI]', ...
        'MaxAE E','MaxAE Q');
    printComponents(cfg.dofNames, results, results.trial.motionIdx == imot);
end

% ---- (3) geodesic orientation error -----------------------------------
fprintf('\nGEODESIC ORIENTATION ERROR, ALL SIMULATIONS\n');
fprintf('%s\n', repmat('-',1,86));
fprintf('%-20s %16s %16s %20s\n', ...
    'Segment','Euler','Quaternion','Difference [95% CI]');
printGeoValidation(cfg.segNames, results, true(nTr,1));

for imot = 1:nMot
    fprintf('\nGEODESIC ORIENTATION ERROR, %s PLANE ELEVATION\n', ...
            upper(cfg.motionNames{imot}));
    fprintf('%s\n', repmat('-',1,86));
    fprintf('%-20s %16s %16s %20s\n', ...
        'Segment','Euler','Quaternion','Difference [95% CI]');
    printGeoValidation(cfg.segNames, results, results.trial.motionIdx == imot);
end

fprintf(['\nMarker RMSE and maximum marker error are a separate quantity\n' ...
         'and are taken from the OpenSim inverse kinematics output.\n\n']);

end

% ========================================================================
function printComponents(names, results, idx)
for j = 1:numel(names)
    e = results.trial.euler.rmse(idx,j);
    q = results.trial.quat.rmse(idx,j);
    [lo,hi] = ci95(q - e);
    fprintf('%-18s %8.1f %7.2f +/- %-6.2f %7.2f +/- %-6.2f %7.2f [%5.2f,%5.2f] %8.2f %8.2f\n', ...
        names{j}, mean(results.trial.euler.rom(idx,j)), ...
        mean(e), std(e), mean(q), std(q), mean(q-e), lo, hi, ...
        max(results.trial.euler.maxabs(idx,j)), ...
        max(results.trial.quat.maxabs(idx,j)));
end
end



% ========================================================================
function [lo,hi] = ci95(d)

n  = numel(d);
m  = mean(d);
se = std(d) / sqrt(n);
lo = m - tcrit(n-1)*se;
hi = m + tcrit(n-1)*se;
end

function printGeoValidation(names, results, idx)
for j = 1:numel(names)
    e = results.geo.euler.rmse(idx,j);
    q = results.geo.quat.rmse(idx,j);
    [lo,hi] = ci95(q - e);
    fprintf('%-20s %7.2f +/- %-6.2f %7.2f +/- %-6.2f %7.2f [%5.2f,%5.2f]\n', ...
        names{j}, mean(e), std(e), mean(q), std(q), mean(q-e), lo, hi);
end
end

