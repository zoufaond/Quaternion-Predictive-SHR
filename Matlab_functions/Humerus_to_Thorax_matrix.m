function mat = Humerus_to_Thorax_matrix(q)
    mat = Ryzx(q(1:3))*Ryzx(q(4:6))*Ryzyy(q(7:9));
end

function yzyy = Ryzyy(q)
    yzyy = Ry(q(1))*Rz(q(2))*Ry(q(3));
end

function yzx = Ryzx(q)
    yzx = Ry(q(1))*Rz(q(2))*Rx(q(3));
end

function rot_x = Rx(phi)
    rot_x = zeros(3,3);
    rot_x(1,1) = 1;
    rot_x(2,2) = cos(phi);
    rot_x(2,3) = -sin(phi);
    rot_x(3,2) = sin(phi);
    rot_x(3,3) = cos(phi);
end


function rot_y = Ry(phi)
    rot_y = zeros(3,3);
    rot_y(2,2) = 1;
    rot_y(1,1) = cos(phi);
    rot_y(3,1) = -sin(phi);
    rot_y(1,3) = sin(phi);
    rot_y(3,3) = cos(phi);
end

function rot_z = Rz(phi)
    rot_z = zeros(3,3);
    rot_z(3,3) = 1;
    rot_z(1,1) = cos(phi);
    rot_z(2,1) = sin(phi);
    rot_z(1,2) = -sin(phi);
    rot_z(2,2) = cos(phi);
end
