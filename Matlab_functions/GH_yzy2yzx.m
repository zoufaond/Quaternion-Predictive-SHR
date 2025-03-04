function motion_new = GH_yzy2yzx(motion)

motion_new = motion;
GH_rotm = eul2rotm(motion(:,7:9),'YZY');
GH_YZX = rotm2eul(GH_rotm,'YZX');

motion_new(:,7:9) = GH_YZX;

end