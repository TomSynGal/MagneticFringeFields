clear all
close all
clc

Pole0 = load('dipole0.csv');
Pole1 = load('dipole1.csv');

x0 = Pole0(:,1);
y0 = Pole0(:,2);
z0 = Pole0(:,3);
x1 = Pole1(:,1);
y1 = Pole1(:,2);
z1 = Pole1(:,3);

xv = linspace(min(x0), max(x0), 101);
yv = linspace(min(y0), max(y0), 101);
[X0,Y0] = meshgrid(xv, yv);
Z0 = griddata(x0,y0,z0,X0,Y0);
xv = linspace(min(x1), max(x1), 101);
yv = linspace(min(y1), max(y1), 101);
[X1,Y1] = meshgrid(xv, yv);
Z1 = griddata(x1,y1,z1,X1,Y1);

figure(1);hold on
surf(X0, Y0, Z0);
surf(X1, Y1, Z1);
grid on
set(gca, 'ZLim',[-1.5 1.5])
shading interp