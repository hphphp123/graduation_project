clc
clear all
close all
x=15:10:100;
y=10:10:100;
[xx,yy]=meshgrid(x, y); % xx和yy都是25x25的矩阵  
zz=xx.*0; % 计算函数值，zz也是21x21的矩阵
zz(10,2)=1;
mesh(xx, yy, zz); % 画出立体网状图


