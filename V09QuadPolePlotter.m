clear all
close all
clc

Pole0 = load('quad0.csv');
Pole1 = load('quad1.csv');
Pole2 = load('quad2.csv');
Pole3 = load('quad3.csv');

x0 = Pole0(:,1);
y0 = Pole0(:,2);
z0 = Pole0(:,3);
x1 = Pole1(:,1);
y1 = Pole1(:,2);
z1 = Pole1(:,3);
x2 = Pole2(:,1);
y2 = Pole2(:,2);
z2 = Pole2(:,3);
x3 = Pole3(:,1);
y3 = Pole3(:,2);
z3 = Pole3(:,3);

xv = linspace(min(x0), max(x0), 101);
yv = linspace(min(y0), max(y0), 101);
[X0,Y0] = meshgrid(xv, yv);
Z0 = griddata(x0,y0,z0,X0,Y0);
xv = linspace(min(x1), max(x1), 101);
yv = linspace(min(y1), max(y1), 101);
[X1,Y1] = meshgrid(xv, yv);
Z1 = griddata(x1,y1,z1,X1,Y1);
xv = linspace(min(x2), max(x2), 101);
yv = linspace(min(y2), max(y2), 101);
[X2,Y2] = meshgrid(xv, yv);
Z2 = griddata(x2,y2,z2,X2,Y2);
xv = linspace(min(x3), max(x3), 101);
yv = linspace(min(y3), max(y3), 101);
[X3,Y3] = meshgrid(xv, yv);
Z3 = griddata(x3,y3,z3,X3,Y3);

figure(1);hold on
surf(X0, Y0, Z0);
surf(X1, Y1, Z1);
surf(X2, Y2, Z2);
surf(X3, Y3, Z3);
grid on
set(gca, 'ZLim',[-1.5 1.5])
shading interp