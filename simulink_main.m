clearvars
addpath Geometry\
addpath equations_of_motion\quaternion\
addpath equations_of_motion\euler\
load('data_model.mat');
model = params.model;
initEul = params.InitPosOptEul.initCondEul;
initQuat = params.InitPosOptEul.initCondQuat;



