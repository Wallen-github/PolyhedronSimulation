A_vert = csvread('Apophis_vert.csv');
A_vert = A_vert(:,2:4);
A_facet = csvread('Apophis_facet.csv');


%% Inertial Frame
figure(1)
set(gcf, 'Position', [1 100 1100 700])
set(gcf,'Color','white')


% Plot
A = patch('Faces',A_facet,'Vertices',A_vert,'FaceColor',[.64 .58 .57],'LineStyle','none');
material dull
camlight
% hold on
% A.FaceLighting = 'gouraud';
% A.DiffuseStrength = 0.5;
% A.AmbientStrength = 0.3;
view(3)
axis equal