clc
clear all
close all
x=15:10:100;
y=10:10:100;
[xx,yy]=meshgrid(x, y); % xx��yy����25x25�ľ���  
zz=xx.*0; % ���㺯��ֵ��zzҲ��21x21�ľ���
zz(10,2)=1;
mesh(xx, yy, zz); % ����������״ͼ


