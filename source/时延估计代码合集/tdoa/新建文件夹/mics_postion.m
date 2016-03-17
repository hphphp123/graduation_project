clear all
clc

%% work dir
workdir = 'E:\MATLAB\R2006a\work\springer_book';

%% file
file = sprintf('%s\\%s',workdir,'mics_harmon');

%% load file
Imp = load(file);
mics_ps = Imp.mics