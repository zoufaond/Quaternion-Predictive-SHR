function [pos_TS,pos_AI] = contact_points_position(q,OS_model)
    jnts = OS_model.joints;
    offset_thorax = jnts{1,2}.location;
    offset_clavicle = jnts{1,5}.location;
    TS = OS_model.TScontact.translation;
    AI = OS_model.AIcontact.translation;
    pos_TS = T_trans(offset_thorax) * R_y(q(1)) * R_z(q(2)) * R_x(q(3)) * T_trans(offset_clavicle) * R_y(q(4)) * R_z(q(5)) * R_x(q(6)) * position(TS);
    pos_AI = T_trans(offset_thorax) * R_y(q(1)) * R_z(q(2)) * R_x(q(3)) * T_trans(offset_clavicle) * R_y(q(4)) * R_z(q(5)) * R_x(q(6)) * position(AI);
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

function T = T_trans(vec)
    T = [1,0,0,vec(1);
         0,1,0,vec(2);
         0,0,1,vec(3);
         0,0,0,1];
end