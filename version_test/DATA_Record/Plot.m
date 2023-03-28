%--------------------------------------------------------------------------
% Name: 

% Desc: 

% Auth: Hai-Shuo Wang

% Time: 01/26/2023

% Version 2.0:
%--------------------------------------------------------------------------

clc;
clear;
close all;
format LONG;

posvel = load('pos_vel_spin.dat');
posvel_post = load('pos_vel_spin_post.dat');



figure
plot3(posvel(:,7),posvel(:,8),posvel(:,9),LineWidth=2)
hold on
plot3(posvel(:,16),posvel(:,17),posvel(:,18),LineWidth=2)
xlabel('\omega_x')
ylabel('\omega_y')
zlabel('\omega_z')
grid on
% title('Trajectory')
set(gca,'FontSize',20,'FontWeight','bold')

figure
plot3(posvel(:,1),posvel(:,2),posvel(:,3),LineWidth=2,Color='blue')
hold on
plot3(posvel(:,10),posvel(:,11),posvel(:,12),LineWidth=2,Color='magenta')
plot3(posvel_post(:,1),posvel_post(:,2),posvel_post(:,3),LineWidth=2,LineStyle='--',Color='blue')
plot3(posvel_post(:,10),posvel_post(:,11),posvel_post(:,12),LineWidth=2,LineStyle='--',Color='magenta')
xlabel('x')
ylabel('y')
zlabel('z')
grid on
title('Trajectory')
set(gca,'FontSize',20,'FontWeight','bold')

figure
hold on
plot(posvel(:,1)-posvel_post(:,1))
plot(posvel(:,2)-posvel_post(:,2))
xlabel('t')
ylabel('\delta')
grid on
title('Trajectory')
set(gca,'FontSize',20,'FontWeight','bold')

momentum = load('angular_momentum.dat');
figure
hold on
plot(momentum(:,9),LineWidth=2)
xlabel('t')
ylabel('angular momentum')
grid on
set(gca,'FontSize',20,'FontWeight','bold')

energy = load('energy.dat');
figure
hold on
plot(energy(:,3),LineWidth=2)
xlabel('t')
ylabel('energy')
grid on
set(gca,'FontSize',20,'FontWeight','bold')

for i=1:size(momentum,1)
    momV = norm(momentum(i,7:9));
    we(i) = 2*energy(i,3)/momV;
    Id(i) = momV^2/2/energy(i,3);
end
figure
scatter(we,Id)

