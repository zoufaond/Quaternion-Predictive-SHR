function [time_training_full,trajectory_training_full] = noise2data(trajectory_training, clavicle_xrot_vals,noise_upper,noise_lower,noise_nearGH,num_noised, OS_model,GH_seq, thorax_dimensions)

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
        NL(7) = -noise_nearGH(1);
        NL(9) = -noise_nearGH(1);
        noise = NL + (NU-NL).*rand(1,11);
    else
        noise = NL + (NU-NL).*rand(1,11);
    end
        trajectory_noised = trajectory_training(i,:) + noise;

        % Check if any of the insertion points of serratus and lats muscles
        % are in the wrapping object, if so, ignore this sample (OpenSim would throw NaN value of moment arm)
        if is_in_thorax_humerus(trajectory_noised,OS_model) == 1 || is_in_thorax_scap(trajectory_noised,OS_model) == 1
            continue
        end
        % 
        % Check if the contact points are close to thorax elipsoid
        [is_scap_close,~] = is_close_to_thorax(trajectory_noised,OS_model,0.1,0.1,thorax_dimensions);
        if  is_scap_close == 0
            continue
        end

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

function [c,ceq] = mycon(x,traj)
    Rscap = R_scap_glob(traj(1:3),traj(4:6));
    ceq = Rscap - R_scap_glob([traj(1),traj(2),x(1)],[x(2),x(3),x(4)]);
    c = [];
end

function res = conoid_length(AC)
    O = [0.1165,-0.0041,0.0143];
    I = T_trans([0.1575,0,0]) * YZX_seq(AC) * position([-0.0536, -0.0009, -0.0266]);
    res = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);
end

function [res,Eeq1,Eeq2] = is_close_to_thorax(q,OS_model,mean_val,range_val,thorax_dimensions)
    [pos_TS,pos_AI] = contact_points_position(q,OS_model);
    Epos = [0 -0.1521 0.0621];
    Edim = thorax_dimensions;
    Eeq1 = elips_eq(pos_TS,Epos,Edim);
    Eeq2 = elips_eq(pos_AI,Epos,Edim);
    val = sqrt(Eeq1^2+Eeq2^2);
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