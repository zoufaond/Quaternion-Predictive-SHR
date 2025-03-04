clear all;
motion_training = 'Motions\Scabduction_yzx\Scabduction_yzx_training.mot';
motion_traj = 'Motions\Scabduction_yzx\Scabduction_yzx.mot';
[motion_path,motion_name,extension] = fileparts(motion_training);
osim_file = 'das3_noserrant12_yzx.osim';
mydir = 'Polyfiles';
musclepoly_file = 'musclepoly';

% das3_polynomials(osim_file,mydir,motion_training,musclepoly_file);
% model = das3_readosim(osim_file,[motion_path,'\',mydir '\' musclepoly_file]);
% save([motion_path,'\','OS_model'],'model')
motion2struct(motion_traj)