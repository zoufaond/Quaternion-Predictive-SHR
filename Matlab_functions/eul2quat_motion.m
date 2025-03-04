function res = eul2quat_motion(angles, GH_seq)

SC_EUL = angles(:,1:3);
AC_EUL = angles(:,4:6);
GH_EUL = angles(:,7:9);

SC_Q = eul2quat(SC_EUL,'YZX');
AC_Q = eul2quat(AC_EUL,'YZX');
GH_Q = eul2quat(GH_EUL,GH_seq);

EL_x = angles(:,10);

res = [SC_Q,AC_Q,GH_Q,EL_x];
end