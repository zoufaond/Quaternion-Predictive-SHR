clearvars
%initial conditions Euler
SC = [-21.784 6.303 1.000]*pi/180;
AC = [46.295 4.899 -1.180]*pi/180;
GH = [0.000 8 0.000]*pi/180;
EL_x = 5.000*pi/180;

% initial conditions quat
SC_Q = eul2quat(SC,'YZX');
AC_Q = eul2quat(AC,'YZX');
GH_Q = eul2quat(GH,'YZY');

params.initCond_Eul = [SC,AC,GH,EL_x,zeros(1,10)]';
params.initCond_Q = [SC_Q,AC_Q,GH_Q,EL_x,zeros(1,10)]';

%global constants
params.model.g = 9.81;
params.model.c = 0.5;
params.model.k = 0;

%offset translations and rotations in each joints (measured from parents)
params.model.offset_thorax = [0.025 0.0068 0.0061];
params.model.offset_clavicula = [0.1575 0 0];
params.model.offset_scapula = [0.0005 -0.0458 -0.008];
params.model.offset_humerus = [0.0058 -0.2907 0.0049];
params.model.offset_humerus_rot = [0 0 0.32318];
params.model.EL_rot_axis = [0.969 -0.247 0];
[params.model.ELx_ax,params.model.ELx_angle] = vecs2ax_angle(params.model.EL_rot_axis,[0,0,1]);
params.model.offset_ulna = [0.0067234 -0.00931306 -0.0068107];
params.model.PSY_rot_axis = [0.182 0.98227 -0.044946];
[params.model.PSy_ax,params.model.PSy_angle] = vecs2ax_angle(params.model.PSY_rot_axis,[0,0,1]);
params.model.offset_radius = [-0.0333 -0.2342 0.0092];
params.model.offset_hand = [0 0 0];

% center of mass of segments
params.model.com_clavicula = [0.0747 -0.0053 0.0009];
params.model.com_scapula = [-0.0545 -0.0457 0.0112];
params.model.com_humerus = [0.00767448 -0.0794264 -0.00270681];
params.model.com_ulna = [-0.0199622 -0.0594577 0.0101954];
params.model.com_radius = [-0.00736417 -0.118863 0.00170825];
params.model.com_hand = [0.0125107 -0.105531 0.00134913];

% segments inertias
params.model.I_clavicula = [6.4e-06 2.63e-05 2.43e-05 0 0 0];
params.model.I_scapula = [0.001 0.001 0.001 0 0 0];
params.model.I_humerus = [0.0132 0.001988 0.0132 0 0 0];
params.model.I_ulna = [0.0030585 0.00045325 0.0030585 0 0 0];
params.model.I_radius = [0.0030585 0.00045325 0.0030585 0 0 0];
params.model.I_hand = [0.0006387 0.0001904 0.0006387 0 0 0];

% segment masses
params.model.mass_clavicula = 0.156;
params.model.mass_scapula = 0.7054;
params.model.mass_humerus = 2.0519;
params.model.mass_ulna = 0.5464;
params.model.mass_radius = 0.5464;
params.model.mass_hand = 0.525;

% recalculated transformation w.r.t. segments coms
params.model.rigid0P = params.model.offset_thorax;
params.model.rigid1C = -params.model.com_clavicula;
params.model.rigid1P = -params.model.com_clavicula+params.model.offset_clavicula;
params.model.rigid2C = -params.model.com_scapula;
params.model.rigid2P = -params.model.com_scapula+params.model.offset_scapula;
params.model.rigid3C = -params.model.com_humerus;
params.model.rigid3P = -params.model.com_humerus+params.model.offset_humerus;
params.model.rigid4C = -params.model.com_ulna;
params.model.rigid4P = -params.model.com_ulna+params.model.offset_ulna;
params.model.rigid5C = -params.model.com_radius;
params.model.rigid5P = -params.model.com_radius+params.model.offset_radius;
params.model.rigid6C = -params.model.com_hand;

%% contact_data
%kontaktni body na scapule
params.model.contTS = [-0.121876 -0.0228 0.0359];
params.model.contAI = [-0.116629 -0.140236 0.0359];
% translation
params.model.elips_trans = [0 -0.1521 0.0621];
params.model.elips_dim = [0.1446 0.2117 0.0956]*1.3;%;
% contact force Chadwick
params.model.k_contact_in = 20000;
params.model.eps_in = 0.01;
params.model.k_contact_out = 2000;
params.model.eps_out = 0.001;
params.model.first_elips_scale = 1;
params.model.second_elips_scale = 1.1;

% conoid_ligamentum
params.model.conoid_eps = 1e-3;
params.model.conoid_stiffness = 80000;
params.model.conoid_length = 0.0174;
params.model.conoid_origin = [0.1365,0.0206,0.0136];
params.model.conoid_insertion = [-0.0536,-9.0e-04,-0.0266];


%%
save('..\data_model.mat','params')
%%

function [axis,angle] = vecs2ax_angle(rot_axis,base_ax)
    perp_axis = cross(base_ax,rot_axis);
    axis = perp_axis/norm(perp_axis);
    angle = acos(dot(rot_axis,base_ax));
end