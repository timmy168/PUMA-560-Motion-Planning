clc
clear all
close all

% starting point 
A = [ 0.64  0.77  0      5 ;
      0.77 -0.64  0    -55 ; 
      0     0    -1    -60 ; 
      0     0     0      1 ];  

% via point
B = [ 0.87 -0.1   0.48  50 ;
      0.29  0.9  -0.34 -40 ; 
     -0.4   0.43  0.81  40 ; 
      0     0     0      1 ];  

% end point
C = [ 0.41 -0.29  0.87  60 ;
      0.69  0.71 -0.09  15 ; 
     -0.6   0.64  0.49 -30 ; 
      0     0     0      1 ]; 

% sampling_rate
sampling_rate = 0.002;

% Cartesian parameter
[xA ,yA ,zA ,RA ,PA, YA] = find_6D_poses(A);
[xB ,yB ,zB ,RB ,PB, YB] = find_6D_poses(B);
[xC ,yC ,zC ,RC ,PC, YC] = find_6D_poses(C);

% Path planning strategy 
% A -> A' straight line
% B       curve line
% C -> C' straight line

% Straight line A -> B1
indexA = 0;   %  the index of the data in A 
x = xB - xA; y = yB - yA; z = zB - zA;
R = RB - RA; P = PB - PA; Y = YB - YA;

% record the coordinate at each smapling time and save them into list
for t=-0.5:sampling_rate:-0.2
    indexA = indexA+1;
    h = (t+0.5)/0.5;
    dx = x*h; dy = y*h; dz = z*h;
    dR = R*h; dP = P*h; dY = Y*h;
    x_A(:,indexA) = xA+dx; y_A(:,indexA) = yA+dy; z_A(:,indexA) =zA+dz; 
    R_A(:,indexA) = RA+dR; P_A(:,indexA) = PA+dP; Y_A(:,indexA) =YA+dY;
end
for i = 1:indexA
    orientation_A(:,i) = Orientation(R_A(:,i),P_A(:,i),Y_A(:,i));
end

% Straight line B2 -> C
indexC = 0;  %  the index of the data in C
x = xC - xB; y = yC - yB; z = zC - zB;
R = RC - RB; P = PC - PB; Y = YC - YB;
for t=0.2:sampling_rate:0.5
    indexC = indexC+1;
    h = (t)/0.5;
    dx = x*h; dy = y*h; dz = z*h;
    dR = R*h; dP = P*h; dY = Y*h;
    x_C(:,indexC) = xB+dx; y_C(:,indexC) = yB+dy; z_C(:,indexC) = zB+dz; 
    R_C(:,indexC) = RB+dR; P_C(:,indexC) = PB+dP; Y_C(:,indexC) = YB+dY;
end
for i = 1:indexC
    orientation_C(:,i) = Orientation(R_C(:,i),P_C(:,i),Y_C(:,i));
end

% Curve B
% The end point of Straight A is the start point of curve B
xA = x_A(:,indexA)-xB; yA = y_A(:,indexA)-yB; zA = z_A(:,indexA)-zB;
RA = R_A(:,indexA)-RB; PA = P_A(:,indexA)-RB; YA = Y_A(:,indexA)-RB;
% The start point of Straight C is the end point of curve B 
xC = x_C(:,1)-xB; yC = y_C(:,1)-yB; zC = z_C(:,1)-zB;
RC = R_C(:,1)-RB; PC = P_C(:,1)-PB; YC = Y_C(:,1)-YB;
indexB = 0; %  the index of the data in B
for t=(-0.2+sampling_rate):sampling_rate:(0.2-sampling_rate)
    indexB = indexB+1;
    h = (t+0.2)/(2*0.2);
    x_B(:,indexB) = xB + ((xC+xA)*(2-h)*h^2-2*xA)*h+xA;
    y_B(:,indexB) = yB + ((yC+yA)*(2-h)*h^2-2*yA)*h+yA;
    z_B(:,indexB) = zB + ((zC+zA)*(2-h)*h^2-2*zA)*h+zA;
    R_B(:,indexB) = RB + ((RC+RA)*(2-h)*h^2-2*RA)*h+RA;
    P_B(:,indexB) = PB + ((PC+PA)*(2-h)*h^2-2*PA)*h+PA;
    Y_B(:,indexB) = YB + ((YC+YA)*(2-h)*h^2-2*YA)*h+YA;   
end
for i = 1:indexB
    orientation_B(:,i) = Orientation(R_B(:,i),P_B(:,i),Y_B(:,i));
end

% XYZ Position Variation
X=[x_A x_B x_C]; Y=[y_A y_B y_C]; Z=[z_A z_B z_C];
X1=[x_A x_B]; Y1=[y_A y_B]; Z1=[z_A z_B];
t=0:sampling_rate:1;
figure(1);
subplot(3,1,1);plot(t,X);
xlabel('Time(s)');ylabel('Position(cm)');
title('position of x');grid;
%%
subplot(3,1,2);plot(t,Y);
xlabel('Time(s)');ylabel('Position(cm)');
title('position of y');grid;
%%
subplot(3,1,3);plot(t,Z);
xlabel('Time(s)');ylabel('Position(cm)');
title('position of z');grid;

% XYZ Velocity Variation
% First Derivative
dt=t(2:501);
dX=diff(X)/sampling_rate;
dY=diff(Y)/sampling_rate;
dZ=diff(Z)/sampling_rate;
figure(2)
subplot(3,1,1);plot(dt,dX);
xlabel('Time(s)');ylabel('Velocity(cm/s)');
title('velocity of x');grid;
%%
subplot(3,1,2);plot(dt,dY);
xlabel('Time(s)');ylabel('Velocity(cm/s)');
title('velocity of y');grid;
%%
subplot(3,1,3);plot(dt,dZ);
xlabel('Time(s)');ylabel('Velocity(cm/s)');
title('velocity of z');grid;

% XYZ Acceleration Variation
% Second Derivative
dt2=t(3:501);
dX2=diff(dX)/sampling_rate;
dY2=diff(dY)/sampling_rate;
dZ2=diff(dZ)/sampling_rate;
figure(3);
subplot(3,1,1);plot(dt2,dX2);
xlabel('Time(s)');ylabel('Acceleration(cm/s^2)');
title('acceleration of x');grid;
%%
subplot(3,1,2);plot(dt2,dY2);
xlabel('Time(s)');ylabel('Acceleration(cm/s^2)');
title('acceleration of y');grid;
%%
subplot(3,1,3);plot(dt2,dZ2);
xlabel('Time(s)');ylabel('Acceleration(cm/s^2)');
title('acceleration of z');grid;

% 3D Path Visualize
figure(4);
x_all = [x_A x_B x_C]; y_all = [y_A y_B y_C]; z_all = [z_A z_B z_C]; 
orientation_all = [orientation_A orientation_B orientation_C];
quiver3(x_all,y_all,z_all,orientation_all(1,:),orientation_all(2,:),orientation_all(3,:),'Color', [0, 1, 1]);
xlabel('x(cm)');ylabel('y(cm)');zlabel('z(cm)');

hold on;
scatter3(x_all,y_all,z_all,'b','filled', 'SizeData', 5);
scatter3(  5, -55, -60, 'r', 'filled');  % A
scatter3( 50, -40,  40, 'r', 'filled');  % B
scatter3( 60,  15, -30, 'r', 'filled');  % C
hold off;

text(  5, -55, -60,'A( 5,-55,-60)');
text( 50, -40,  40,'B(50,-40, 40)');
text( 60,  15, -30,'C(60, 15,-30)');
title('3D path of Cartesion Motion')

% -------------------------------------------------------------------------
% Function Defination

% find_6D_poses: find x,y,z,roll,pitch,yaw from a T6 matrix
function [x,y,z,phi,theta,psi] = find_6D_poses(T6)
    % Degree and Radius Transformation
    D_to_R = pi / 180;
    R_to_D = 180 / pi;
   
    nx = T6(1,1); ny = T6(2,1); nz = T6(3,1);
    ox = T6(1,2); oy = T6(2,2); oz = T6(3,2);
    ax = T6(1,3); ay = T6(2,3); az = T6(3,3);
    px = T6(1,4); py = T6(2,4); pz = T6(3,4);
    x = px; y = py; z = pz;
    phi = atan2(ay, ax) * R_to_D;
    theta = atan2(sqrt(ax^2 + ay^2), az) * R_to_D;
    psi = atan2(oz, -nz) * R_to_D;
    % p = [ x, y, z, phi,theta, psi];
end

% Orientation: find the orientation of Z axis from raw pitch yaw
function orientation = Orientation(R,P,Y)
    % Degree and Radius Transformation
    D_to_R = pi / 180;
    R_to_D = 180 / pi;
    R = R * D_to_R;
    P = P * D_to_R;
    Y = Y * D_to_R;
    orientation = [cos(R)*sin(P) sin(R)*sin(P) cos(P)];
end