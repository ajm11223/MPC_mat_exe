clc; clear all; close all;
omega=10;
numc=omega^2; %전달 함수 분자
denc=[1 0.1*omega omega^2]; % 전달함수 분모
[Ac,Bc,Cc,Dc]=tf2ss(numc,denc);