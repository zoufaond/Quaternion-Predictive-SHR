function data2mot(angles,time_simulation,name,type, output_from,GH_seq)

time = linspace(0,time_simulation(end),size(angles,1));

if strcmp('quaternion', type)
    if strcmp('simulink',  output_from)
        angles = out.simdata_Q.signals.values(:,:)';
    elseif strcmp('struct', output_from)
        angles = angles;
    end
    SC_Q = angles(:,1:4);
    AC_Q = angles(:,5:8);
    GH_Q = angles(:,9:12);

    SC_RM = quat2rotm(SC_Q);
    AC_RM = quat2rotm(AC_Q);
    GH_RM = quat2rotm(GH_Q);

    SC_Eul = rotm2eul(SC_RM,'YZX');
    AC_Eul = rotm2eul(AC_RM,'YZX');
    GH_Eul = rotm2eul(GH_RM,GH_seq);

    EL_x = angles(:,13);
    PS_y = angles(:,14);

    angles_OS = [SC_Eul,AC_Eul,GH_Eul,EL_x,PS_y]*180/pi;

elseif strcmp('euler', type)
    if strcmp('simulink',  output_from)
        angles = out.simdata_Eul.signals.values(:,:)';
    elseif strcmp('struct', output_from)
        angles = angles;
    end
    angles_OS = angles*180/pi;
end

angles_OS = angles_OS';

fid = fopen(name,'w');
fprintf(fid, 'simulation\n');
fprintf(fid,'nRows=%i\n',size(angles_OS,2));
fprintf(fid,'nColumns=%i\n',size(angles_OS,1)+4);
fprintf(fid,'endheader\n');
if strcmp(GH_seq,'YZY')
    fprintf(fid,'time  TH_x TH_y TH_z SC_y  SC_z  SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_yy  EL_x  PS_y\n');
elseif  strcmp(GH_seq,'YZX')
    fprintf(fid,'time  TH_x TH_y TH_z SC_y  SC_z  SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_x  EL_x  PS_y\n');
end

for i=1:size(angles_OS,2)
    for j=1:size(angles_OS,1)
        if j==1
            fprintf(fid, '%4f  0.000000  0.000000  0.000000  %4f', time(i),angles_OS(j,i));
        elseif j == size(angles_OS,1)
            if size(angles_OS,1) == 10
                fprintf(fid, '  %4f  100.000000\n', angles_OS(j,i));
            elseif size(angles_OS,1) == 11
                fprintf(fid, '  %4f\n', angles_OS(j,i));
            end
        else
            fprintf(fid, '  %4f', angles_OS(j,i));
        end
    end
end
fclose(fid);

end