function res = get_conoid_length(q_AC,model)

O_pos = model.conoid_origin;
I_pos = model.conoid_insertion;
offset_clavicle = model.offset_clavicula;
        
O = position(O_pos(1), O_pos(2), O_pos(3));
TC_S = T_trans(offset_clavicle);
RC_S = R_y(q_AC(1)) * R_z(q_AC(2)) * R_x(q_AC(3));
I = TC_S * RC_S * position(I_pos(1), I_pos(2), I_pos(3));


res = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);

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

function r = position(x,y,z)
    r = [x;y;z;1];
end