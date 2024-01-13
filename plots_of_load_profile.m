%plotting of load profiles

clc;  clear all;
load('PP_L.mat');   load('PP_solar.mat');   load('PP_W.mat');
t=1:1:48;

figure(1)
plot(t(1:24),P_l(1:24))
hold
plot(t(24:end),P_l(24:end))
% savefig('load profile')

figure(2)
plot(t(1:24),P_solar(1:24))
hold
plot(t(24:end),P_solar(24:end))
% savefig('solar profile')

figure(3)
plot(t(1:24),P_W(1:24))
hold
plot(t(24:end),P_W(24:end))
% savefig('wind profile')