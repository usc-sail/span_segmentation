% WRAP_MAKE_TEMPLATE - run wrap_make_template to create the set of
% templates which will be used to segment the video
% 
% Tanner Sorensen
% University of Southern California
% 
% Friday, June 16, 2017

addpath(genpath('functions')) % add path (recursively) to functions

avi_file = './demo_files/ac09132015_10_46_47.avi'; % edit appropriately
template_struct = './demo_files/template_struct.mat'; % initial template 
video_struct = avi_to_struct(avi_file); % video for making template

make_template(template_struct,video_struct,[])