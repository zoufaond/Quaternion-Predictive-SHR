function data2struct(angles,time,path_name, GH_seq)

quats = eul2quat_motion(angles, GH_seq);

mot_struct.euler = angles(:,1:10); % PSy not included
mot_struct.quat = quats;
mot_struct.time = time;

save(path_name,'mot_struct');
end