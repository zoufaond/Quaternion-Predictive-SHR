function res = min_conoid_length(traj, OS_model)
    fun = @(x) conoid_length(x(2:4),OS_model);
    % fun = @(x) sum((Rscap - R_scap_glob([traj(1),traj(2),x(1)],[x(2),x(3),x(4)])).^2,"all");
    nonlcon = @(x) mycon(x,traj);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = ones(1,4)*(-pi/2);
    ub = ones(1,4)*pi/2;
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

function res = conoid_length(AC,OS_model)
    O = OS_model.conoid_origin;
    clavicle_offset = OS_model.joints{5}.location;
    Ipos = OS_model.conoid_insertion;
    I = T_trans(clavicle_offset) * YZX_seq(AC) * position(Ipos);
    res = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);
end

function res = R_scap_glob(jnt1,jnt2)
    res = YZX_seq(jnt1) * YZX_seq(jnt2);
end

function res = YZX_seq(phi_vec)
    res = R_y(phi_vec(1)) * R_z(phi_vec(2)) * R_x(phi_vec(3));
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