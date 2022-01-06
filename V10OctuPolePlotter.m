clear all
close all
clc

Pole0 = load('octu0.csv');
Pole1 = load('octu1.csv');
Pole2 = load('octu2.csv');
Pole3 = load('octu3.csv');
Pole4 = load('octu4.csv');
Pole5 = load('octu5.csv');
Pole6 = load('octu6.csv');
Pole7 = load('octu7.csv');

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
x4 = Pole4(:,1);
y4 = Pole4(:,2);
z4 = Pole4(:,3);
x5 = Pole5(:,1);
y5 = Pole5(:,2);
z5 = Pole5(:,3);
x6 = Pole6(:,1);
y6 = Pole6(:,2);
z6 = Pole6(:,3);
x7 = Pole7(:,1);
y7 = Pole7(:,2);
z7 = Pole7(:,3);

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
xv = linspace(min(x4), max(x4), 101);
yv = linspace(min(y4), max(y4), 101);
[X4,Y4] = meshgrid(xv, yv);
Z4 = griddata(x4,y4,z4,X4,Y4);
xv = linspace(min(x5), max(x5), 101);
yv = linspace(min(y5), max(y5), 101);
[X5,Y5] = meshgrid(xv, yv);
Z5 = griddata(x5,y5,z5,X5,Y5);
xv = linspace(min(x6), max(x6), 101);
yv = linspace(min(y6), max(y6), 101);
[X6,Y6] = meshgrid(xv, yv);
Z6 = griddata(x6,y6,z6,X6,Y6);
xv = linspace(min(x7), max(x7), 101);
yv = linspace(min(y7), max(y7), 101);
[X7,Y7] = meshgrid(xv, yv);
Z7 = griddata(x7,y7,z7,X7,Y7);

figure(1);hold on
surf(X0, Y0, Z0);
surf(X1, Y1, Z1);
surf(X2, Y2, Z2);
surf(X3, Y3, Z3);
surf(X4, Y4, Z4);
surf(X5, Y5, Z5);
surf(X6, Y6, Z6);
surf(X7, Y7, Z7);
grid on
set(gca, 'ZLim',[-1.5 1.5])
shading interp