function res = get_conoid_length(q_AC,OS_model)

O = OS_model.conoid_origin;
clavicle_offset = OS_model.joints{5}.location;
Ipos = OS_model.conoid_insertion;
I = T_trans(clavicle_offset) * YZX_seq(q_AC) * position(Ipos);
res = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);

end

function res = YZX_seq(phi_vec)
    res = R_y(phi_vec(1)) * R_z(phi_vec(2)) * R_x(phi_vec(3));
end

function trans_x = T_trans(vec)
    trans_x = [1,0,0,vec(1);
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