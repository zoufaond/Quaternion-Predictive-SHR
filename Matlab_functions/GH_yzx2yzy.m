function motion_new = GH_yzx2yzy(motion)

motion_new = motion;
GH_rotm = eul2rotm(motion(:,7:9),'YZX');
GH_YZX = rotm2eul(GH_rotm,'YZY');

motion_new(:,7:9) = GH_YZX;

end