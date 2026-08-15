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
% GH_Eul = unwrap_symmetric_eul(GH_Eul);

EL_x = quat(:,13);


res = [SC_Eul,AC_Eul,GH_Eul,EL_x];


end

function eul = unwrap_symmetric_eul(eul)
% For symmetric sequences such as YZY the second angle is returned in
% [0, pi]. The representation (a1 + pi, -a2, a3 + pi) is equivalent.
% This selects at each sample whichever representation is continuous
% with the previous one, so the second angle can become negative.

for i = 2:size(eul,1)
    alt = [eul(i,1) + pi, -eul(i,2), eul(i,3) + pi];

    d_std = ang_dist(eul(i,:), eul(i-1,:));
    d_alt = ang_dist(alt,      eul(i-1,:));

    if d_alt < d_std
        eul(i,:) = alt;
    end
end
end

function d = ang_dist(a, b)
diff = a - b;
diff = atan2(sin(diff), cos(diff));
d = sum(abs(diff));
end