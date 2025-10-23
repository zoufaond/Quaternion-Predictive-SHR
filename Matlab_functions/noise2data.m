function [time_training_full,trajectory_training_full] = noise2data(trajectory_training, clavicle_xrot_vals,noise_upper,noise_lower,noise_nearGH,num_noised, OS_model,osim_file,GH_seq)

trajectory_training_full = [];
numdata_training = length(trajectory_training(:,1));

for i=1:numdata_training

    NU = noise_upper;
    NL = noise_lower;

    for inoise = 1:num_noised
    if strcmp(GH_seq,'YZY') && abs(trajectory_training(i,8)) < deg2rad(15)
        % scaleY = (1.1-trajectory_training(i,8)/max(trajectory_training(:,8)));
        NU(7) = noise_nearGH(1);
        NU(9) = noise_nearGH(1);
        NL(7) = noise_nearGH(2);
        NL(9) = noise_nearGH(2);
        noise = NL + (NU-NL).*rand(1,11);
    else
        noise = NL + (NU-NL).*rand(1,11);
    end
        trajectory_noised = trajectory_training(i,:) + noise;

        % Check if any of the insertion points of serratus and lats muscles
        % are in the wrapping object, if so, ignore this sample (OpenSim would throw NaN value of moment arm)
        if is_in_thorax_humerus(trajectory_noised,OS_model) == 1
            continue
        end

        if is_in_thorax_scap(trajectory_noised,OS_model) == 1
            continue
        end
        % 
        % Check if the contact points are close to thorax elipsoid
        [is_scap_close,~] = is_close_to_thorax(trajectory_noised,OS_model,0.3,0.5,0,0.1);
        if  is_scap_close == 0
            continue
        end
        % momarm = is_serr_momarm_nan(trajectory_noised,osim_file);
        % trajectory_training_full = [trajectory_training_full; trajectory_noised]; 
        
        for ival = 1:length(clavicle_xrot_vals)
        trajectory_rotclav = change_clavx(trajectory_noised,clavicle_xrot_vals(ival));
        trajectory_training_full = [trajectory_training_full; trajectory_rotclav]; 
        end
    end

end

time_training_full = linspace(0,1,size(trajectory_training_full,1));
% data2mot(trajectory_training_full,time_training_full,[motion_name,'_training.mot'],'euler', 'struct');

end

function [res,Eeq1,Eeq2] = is_close_to_thorax(q,OS_model,TSmean,TSrange,AImean,AIrange)
    [pos_TS,pos_AI] = contact_points_position(q,OS_model);
    Epos = OS_model.thoracic_wall.translation;
    Edim = OS_model.thoracic_wall.dimensions;
    Eeq1 = elips_eq(pos_TS,Epos,Edim);
    Eeq2 = elips_eq(pos_AI,Epos,Edim);

    % if abs(TSmean - Eeq1) < TSrange && abs(AImean - Eeq2) < AIrange
    %     res = 1;
    % else
    %     res = 0;
    % end
    val = sqrt(Eeq1^2+Eeq2^2);
    mean_val = TSmean;
    range_val = TSrange;
    if abs(mean_val - val) < range_val
        res = 1;
    else
        res = 0;
    end
end

function res = is_in_thorax_scap(q,OS_model)

    Epos = OS_model.serr_ant_wrapping.translation;
    Edim = OS_model.serr_ant_wrapping.dimensions;
    for i = 25:36
        current_mus = OS_model.muscles{i};
        IP_scap = scapula_insertion_pos(q,current_mus,OS_model);
        Eeq = elips_eq(IP_scap,Epos,Edim);
        res = 0;
        if Eeq <0
            res = 1;
            break
        end
    end
end

function res = is_in_thorax_humerus(q,OS_model)
    Epos_all = [OS_model.lat_dorsi_sup_wrapping.translation;
                OS_model.lat_dorsi_mid_wrapping.translation;
                OS_model.lat_dorsi_inf_wrapping.translation];

    Edim_all = [OS_model.lat_dorsi_sup_wrapping.dimensions;
                OS_model.lat_dorsi_mid_wrapping.dimensions;
                OS_model.lat_dorsi_inf_wrapping.dimensions];
    which_elips = [1,1,1,2,3,3];
    
    j = 1;
    for imus = 89:94
        current_mus = OS_model.muscles{imus};
        IP_hum = humerus_insertion_pos(q,current_mus,OS_model);
        Eeq = elips_eq(IP_hum,Epos_all(which_elips(j),:),Edim_all(which_elips(j),:));
        res = 0;
        if Eeq <0
            res = 1;
            break
        end
        j = j +1;
    end
end

function pos = scapula_insertion_pos(q,muscle,OS_model)
    jnts = OS_model.joints;
    offset_thorax = jnts{1,2}.location;
    offset_clavicle = jnts{1,5}.location;
    insertion = muscle.origin_position;
    
    pos = T_trans(offset_thorax) * R_y(q(1)) * R_z(q(2)) * R_x(q(3)) * T_trans(offset_clavicle) * R_y(q(4)) * R_z(q(5)) * R_x(q(6)) * position(insertion);
end

function pos = humerus_insertion_pos(q,muscle,OS_model)
    jnts = OS_model.joints;
    offset_thorax = jnts{1,2}.location;
    offset_clavicle = jnts{1,5}.location;
    offset_scapula = jnts{1,8}.location;
    insertion = muscle.insertion_position;
    pos = T_trans(offset_thorax) * R_y(q(1)) * R_z(q(2)) * R_x(q(3)) * T_trans(offset_clavicle) * R_y(q(4)) * R_z(q(5)) * R_x(q(6)) *T_trans(offset_scapula) * R_y(q(7)) * R_z(q(8)) * R_x(q(9)) * position(insertion);
end

function T = T_trans(vec)
    T = [1,0,0,vec(1);
         0,1,0,vec(2);
         0,0,1,vec(3);
         0,0,0,1];
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

function r = position(vec)
    r = [vec(1);vec(2);vec(3);1];
end

function res = R_scap_glob(jnt1,jnt2)
    res = YZX_seq(jnt1) * YZX_seq(jnt2);
end

function res = YZX_seq(phi_vec)
    res = R_y(phi_vec(1)) * R_z(phi_vec(2)) * R_x(phi_vec(3));
end

function res = is_serr_momarm_nan(angles,osimfile)
    for i = 25:36
        res = opensim_get_polyvalues(angles, i, osimfile)
    end
end

function minusdLdq = opensim_get_polyvalues(angles, iMus, osimfile)
import org.opensim.modeling.*
Mod = Model(osimfile);
SimEn = SimbodyEngine();
SimEn.connectSimbodyEngineToModel(Mod);
groundbody = Mod.getBodySet().get('thorax');
scapulabody = Mod.getBodySet().get('scapula_r');
% initialize the system to get the initial state
state = Mod.initSystem;

Mod.equilibrateMuscles(state);

% If we only want one moment arm, put it in a cell
CoordSet = Mod.getCoordinateSet();

% get the muscle
% currentMuscle = Mod.getMuscles().get(Mus);
currentMuscle = Mod.getMuscles().get(iMus-1);
% if ~strcmp(fixname(char(currentMuscle.getName())),Mus)
%     error('Current muscle name is incorrect')
% end
angles = [0 0 0 angles];
minusdLdq = zeros(1,6);
nDofs = CoordSet.getSize;
    
% set dof values for this step
for idof = 1:nDofs
    currentDof = CoordSet.get(idof-1);    
    currentDof.setValue(state,angles(idof),1);
end

% get moment arms
Dofs = [4,5,6,7,8,9];
for idof = 1:6 %only serratus

    currentDof = CoordSet.get(Dofs(idof)-1);
    minusdLdq(idof) = currentMuscle.computeMomentArm(state,currentDof);
end
end