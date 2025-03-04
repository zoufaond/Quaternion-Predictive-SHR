function res = quat2eul_motion(quat,GH_seq)

SC_Q = quat(:,1:4);
AC_Q = quat(:,5:8);
GH_Q = quat(:,9:12);

SC_RM = quat2rotm(SC_Q);
AC_RM = quat2rotm(AC_Q);
GH_RM = quat2rotm(GH_Q);

SC_Eul = rotm2eul(SC_RM,'YZX');
AC_Eul = rotm2eul(AC_RM,'YZX');
GH_Eul = rotm2eul(GH_RM,GH_seq);

EL_x = quat(:,13);


res = [SC_Eul,AC_Eul,GH_Eul,EL_x];


end
