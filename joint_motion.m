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

% IK solutions of the points A, B, C
theta_A = inverse_kinematic(A); thetaA = theta_A(1,:)';
theta_B = inverse_kinematic(B); thetaB = theta_B(1,:)';
theta_C = inverse_kinematic(C); thetaC = theta_C(1,:)';

% Segment A
s1 = 0;
for t=-0.5:sampling_rate:-0.2
    s1 = s1+1;
    h = (t+0.5)/0.5;
    d_thetaA(:,s1) = thetaA+(thetaB-thetaA)*h;                   
    d_omegaA(:,s1) = (thetaB-thetaA)/0.5;
    d_alphaA(:,s1) = [0; 0; 0; 0; 0; 0];
    p1 = foward_kinematic(d_thetaA(:,s1)');
    x1(s1) = p1(1,4); y1(s1) = p1(2,4); z1(s1) = p1(3,4);
    orientation_A(:,s1) = [p1(1,3) p1(2,3) p1(3,3)];
end

% Segment C
s3 = 0;
for t=0.2:sampling_rate:0.5
    s3 = s3+1;
    h = t/0.5;
    d_thetaC(:,s3) = thetaB+(thetaC-thetaB)*h;
    d_omegaC(:,s3) = (thetaC-thetaB)/0.5;
    d_alphaC(:,s3) = [0; 0; 0; 0; 0; 0];
    p3 = foward_kinematic(d_thetaC(:,s3)');
    x3(s3) = p3(1,4); y3(s3) = p3(2,4); z3(s3) = p3(3,4);
    orientation_C(:,s3) = [p3(1,3) p3(2,3) p3(3,3)];
end

% Segment B 
dB = d_thetaA(:,s1)-thetaB;
dC = d_thetaC(:,1)-thetaB;
s2 = 0;
for t=(-0.2+sampling_rate):sampling_rate:(0.2-sampling_rate)
    s2 = s2+1;
    h = (t+0.2)/0.4;
    d_thetaB(:,s2) = ((dC+dB)*(2-h)*h^2-2*dB)*h+dB+thetaB;        
    d_omegaB(:,s2) = ((dC+dB)*(1.5-h)*2*h^2-dB)/0.2;
    d_alphaB(:,s2) = (dC+dB)*(1-h)*3*h/0.2^2;
    p2 = foward_kinematic(d_thetaB(:,s2)');
    x2(s2) = p2(1,4); y2(s2) = p2(2,4); z2(s2) = p2(3,4);
    orientation_B(:,s2) = [p2(1,3) p2(2,3) p2(3,3)];
end

% Joint angle Plot
t=-0.5:sampling_rate:0.5;
figure(1)   
theta1 = [d_thetaA(1,:) d_thetaB(1,:) d_thetaC(1,:)];                     
subplot(3,2,1); plot(t,theta1);
grid; title('joint1 angle'); ylabel('Angle');

theta2 = [d_thetaA(2,:) d_thetaB(2,:) d_thetaC(2,:)];                     
subplot(3,2,2); plot(t,theta2);
grid; title('joint2 angle'); ylabel('Angle');

theta3 = [d_thetaA(3,:) d_thetaB(3,:) d_thetaC(3,:)];                  
subplot(3,2,3); plot(t,theta3);
grid; title('joint3 angle'); ylabel('Angle');

theta4 = [d_thetaA(4,:) d_thetaB(4,:) d_thetaC(4,:)];                    
subplot(3,2,4); plot(t,theta4);
grid; title('joint4 angle'); ylabel('Angle');

theta5 = [d_thetaA(5,:) d_thetaB(5,:) d_thetaC(5,:)];                  
subplot(3,2,5); plot(t,theta5);
grid; title('joint5 angle'); ylabel('Angle');

theta6 = [d_thetaA(6,:) d_thetaB(6,:) d_thetaC(6,:)];                      
subplot(3,2,6); plot(t,theta6);
grid; title('joint6 angle'); ylabel('Angle');

% Joint Velocity Plot
figure(2)   
subplot(3,2,1); plot(t,[d_omegaA(1,:) d_omegaB(1,:) d_omegaC(1,:)]);  
grid; title('joint1'); ylabel('Angular Velocity');

subplot(3,2,2); plot(t,[d_omegaA(2,:) d_omegaB(2,:) d_omegaC(2,:)]);   
grid; title('joint2'); ylabel('Angular Velocity');

subplot(3,2,3); plot(t,[d_omegaA(3,:) d_omegaB(3,:) d_omegaC(3,:)]); 
grid; title('joint3'); ylabel('Angular Velocity');

subplot(3,2,4); plot(t,[d_omegaA(4,:) d_omegaB(4,:) d_omegaC(4,:)]);    
grid; title('joint4'); ylabel('Angular Velocity');

subplot(3,2,5); plot(t,[d_omegaA(5,:) d_omegaB(5,:) d_omegaC(5,:)]);  
grid; title('joint5'); ylabel('Angular Velocity');

subplot(3,2,6); plot(t,[d_omegaA(6,:) d_omegaB(6,:) d_omegaC(6,:)]);   
grid; title('joint6'); ylabel('Angular Velocity');

% Joint Acceleration Plot
figure(3)   
subplot(3,2,1); plot(t,[d_alphaA(1,:) d_alphaB(1,:) d_alphaC(1,:)]);
grid;title('joint1'); ylabel('Angular Acceleration');

subplot(3,2,2); plot(t,[d_alphaA(2,:) d_alphaB(2,:) d_alphaC(2,:)]);
grid;title('joint2'); ylabel('Angular Acceleration');

subplot(3,2,3); plot(t,[d_alphaA(3,:) d_alphaB(3,:) d_alphaC(3,:)]);
grid;title('joint3'); ylabel('Angular Acceleration');

subplot(3,2,4); plot(t,[d_alphaA(4,:) d_alphaB(4,:) d_alphaC(4,:)]);
grid;title('joint4'); ylabel('Angular Acceleration');

subplot(3,2,5); plot(t,[d_alphaA(5,:) d_alphaB(5,:) d_alphaC(5,:)]);
grid;title('joint5'); ylabel('Angular Acceleration');

subplot(3,2,6); plot(t,[d_alphaA(6,:) d_alphaB(6,:) d_alphaC(6,:)]);
grid;title('joint6'); ylabel('Angular Acceleration');

% 3D path visulaize
figure(4)
x_all = [x1 x2 x3]; y_all = [y1 y2 y3]; z_all = [z1 z2 z3]; 
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
title('3D path of joint space planning')