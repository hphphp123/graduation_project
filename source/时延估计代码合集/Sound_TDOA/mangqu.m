clc
clear all
close all
a=0:10:100;
b=0:10:100;
plot(a,b,'.k','MarkerSize',1);
hold on


x=[80,100,60,70,80,90,100,100,90,100,60,90,100,60,90,100,100,80,100];
y=[10,10,20,20,20,20,20,30,50,50,60,60,60,70,70,70,80,90,90];
plot(x,y,'rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',20)
grid on;
