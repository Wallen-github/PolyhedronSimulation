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
flyby = load("flybyOrbit.dat");

figure
plot3(posvel(:,7),posvel(:,8),posvel(:,9),LineWidth=2)
hold on
plot3(posvel(:,16),posvel(:,17),posvel(:,18),LineWidth=2)
xlabel('\omega_x (1/hour)')
ylabel('\omega_y (1/hour)')
zlabel('\omega_z (1/hour)')
grid on
% title('Trajectory')
set(gca,'FontSize',20,'FontWeight','bold')

figure
plot3(posvel(:,1),posvel(:,2),posvel(:,3),LineWidth=2)
hold on
plot3(posvel(:,10),posvel(:,11),posvel(:,12),LineWidth=2)
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid on
title('Trajectory')
set(gca,'FontSize',20,'FontWeight','bold')

figure
plot3(posvel(:,1)-flyby(1,2),posvel(:,2)-flyby(1,3),posvel(:,3)-flyby(1,4),LineWidth=2)
hold on
plot3(posvel(:,10)-flyby(1,2),posvel(:,11)-flyby(1,3),posvel(:,12)-flyby(1,4),LineWidth=2)
% plot3(-flyby(:,2),-flyby(:,3),-flyby(:,4),LineWidth=2,DisplayName='trajectory')
% scatter3(-flyby(1,2),-flyby(1,3),-flyby(1,4),'filled',DisplayName='Initial Position')
scatter3(0,0,0,'filled',DisplayName='Earth')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid on
title('Flyby Orbit')
set(gca,'FontSize',20,'FontWeight','bold')

energy = load('energy.dat');
momentum = load("angular_momentum.dat");
figure
yyaxis left
title('relative conserved quantity')
plot(flyby(:,1),energy(:,1)./abs(energy(1,1)),LineWidth = 2);
xlabel('time')
ylabel('Energy')
yyaxis right
plot(flyby(:,1),momentum(:,4)./momentum(1,4),LineWidth = 2);
ylabel('Momentum^2')
grid on
set(gca,'FontSize',20,'FontWeight','bold')

figure
yyaxis left
title('absolute conserved quantity')
plot(flyby(:,1),energy(:,1),LineWidth = 2);
xlabel('time')
ylabel('Energy')
yyaxis right
plot(flyby(:,1),momentum(:,4),LineWidth = 2);
ylabel('Momentum^2')
grid on
set(gca,'FontSize',20,'FontWeight','bold')

effectQ = load("effectquantity.dat");
figure
yyaxis left
title('relative conserved quantity')
plot(flyby(:,1),effectQ(:,1),LineWidth = 2);
xlabel('time')
ylabel('\omega_e')
yyaxis right
plot(flyby(:,1),effectQ(:,2),LineWidth = 2);
ylabel('I_d')
grid on
set(gca,'FontSize',20,'FontWeight','bold')

figure
scatter(flyby(:,1),effectQ(:,3))
xlabel('time')
ylabel('I_{tilde}')
grid on
set(gca,'FontSize',20,'FontWeight','bold')


figure
plot3(flyby(:,2),flyby(:,3),flyby(:,4),LineWidth=2,DisplayName='trajectory')
hold on
scatter3(flyby(1,2),flyby(1,3),flyby(1,4),'filled',DisplayName='Initial Position')
scatter3(0,0,0,'filled',DisplayName='Earth')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid on
title('Flyby Orbit')
set(gca,'FontSize',20,'FontWeight','bold')

figure
subplot(3,2,1)
plot(flyby(:,1),flyby(:,2))
ylabel('x (m)')
grid on
subplot(3,2,3)
plot(flyby(:,1),flyby(:,3))
ylabel('y (m)')
grid on
subplot(3,2,5)
plot(flyby(:,1),flyby(:,4))
ylabel('z (m)')
grid on
subplot(3,2,2)
plot(flyby(:,1),flyby(:,5))
ylabel('v_x (m/s)')
grid on
subplot(3,2,4)
plot(flyby(:,1),flyby(:,6))
ylabel('v_y (m/s)')
grid on
subplot(3,2,6)
plot(flyby(:,1),flyby(:,7))
ylabel('v_z (m/s)')
xlabel('time (hour)')
grid on
% 

for i=1:size(flyby,1)
    pos(i) = norm(flyby(i,2:4));
    vel(i) = norm(flyby(i,5:7));
end
figure
yyaxis left
title('Norm position and velocity')
plot(flyby(:,1),pos,LineWidth = 2);
xlabel('time')
ylabel('distance')
yyaxis right
plot(flyby(:,1),vel,LineWidth = 2);
ylabel('velocity')
grid on
set(gca,'FontSize',20,'FontWeight','bold')

