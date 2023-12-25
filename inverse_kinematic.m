% -------------------------------------------------------------------------
% Function Defination

%Inverse Kinematic Funtcion
function Theta_Solution = inverse_kinematic(T6)

    % Degree and Radius Transformation
    D_to_R = pi / 180;
    R_to_D = 180 / pi;

    % DH Model - Puma 560
    % d_list = [0, 0, 0.149, 0.433, 0, 0];
    % a_list = [0, 0.432, -0.02, 0, 0, 0];
    % alpha_list = [-90 ,0 ,90 ,-90 ,90 ,0];

    nx = T6(1,1); ny = T6(2,1); nz = T6(3,1);
    ox = T6(1,2); oy = T6(2,2); oz = T6(3,2);
    ax = T6(1,3); ay = T6(2,3); az = T6(3,3);
    px = T6(1,4); py = T6(2,4); pz = T6(3,4);

    % theta1 (2 solutions)
    P = sqrt(px*px + py*py);
    phi = atan2(py,px);
    px_5 = -14.9 / P; % d3 = 14.9
    px_6 = (1 - px_5 ^ 2) ^ 0.5;
    theta1_1 = (phi + atan2(px_5, px_6));
    theta1_2 = (phi + atan2(px_5,-px_6));

    % theta3 (2 solutions)
    M = (px^2+py^2+pz^2-43.2^2-(-2)^2-(14.9)^2-(43.3)^2)/( 2 * 43.2 );
    theta3_1 = (atan2(M, ((-2)^2+(43.3)^2-M^2)^0.5) - atan2(-2,43.3));
    theta3_2 = (atan2(M,-((-2)^2+(43.3)^2-M^2)^0.5) - atan2(-2,43.3));

    % theta 2 (4 solutions)
    A23_1 = cos(theta1_1)*px+sin(theta1_1)*py;
    A23_2 = cos(theta1_2)*px+sin(theta1_2)*py;
    C23_1 = -2 + 43.2*cos(theta3_1);
    C23_2 = -2 + 43.2*cos(theta3_2); 
    B23 = -pz; 
    theta2_1 = atan2(B23,A23_1) + acos(C23_1/(A23_1^2+B23^2)^0.5) - theta3_1; 
    theta2_2 = atan2(B23,A23_2) + acos(C23_1/(A23_2^2+B23^2)^0.5) - theta3_1; 
    theta2_3 = atan2(B23,A23_1) + acos(C23_2/(A23_1^2+B23^2)^0.5) - theta3_2; 
    theta2_4 = atan2(B23,A23_2) + acos(C23_2/(A23_2^2+B23^2)^0.5) - theta3_2; 

    theta123_list_1 = [theta1_1,theta2_1,theta3_1];
    theta123_list_2 = [theta1_2,theta2_2,theta3_1];
    theta123_list_3 = [theta1_1,theta2_3,theta3_2];
    theta123_list_4 = [theta1_2,theta2_4,theta3_2];

    % theta 4 (8 solutions)
    T46_1 = inv(A1_to_3(theta123_list_1(1),theta123_list_1(2),theta123_list_1(3))) * T6 ;
    T46_2 = inv(A1_to_3(theta123_list_2(1),theta123_list_2(2),theta123_list_2(3))) * T6 ;
    T46_3 = inv(A1_to_3(theta123_list_3(1),theta123_list_3(2),theta123_list_3(3))) * T6 ;
    T46_4 = inv(A1_to_3(theta123_list_4(1),theta123_list_4(2),theta123_list_4(3))) * T6 ;
    theta4_1 =  atan2(T46_1(2,3), T46_1(1,3));
    theta4_2 = -atan2(T46_1(2,3),-T46_1(1,3));
    theta4_3 =  atan2(T46_2(2,3), T46_2(1,3));
    theta4_4 = -atan2(T46_2(2,3),-T46_2(1,3));
    theta4_5 =  atan2(T46_3(2,3), T46_3(1,3));
    theta4_6 = -atan2(T46_3(2,3),-T46_3(1,3));
    theta4_7 =  atan2(T46_4(2,3), T46_4(1,3));
    theta4_8 = -atan2(T46_4(2,3),-T46_4(1,3));

    % theta 5 (8 solutions)
    T56_1 = inv(transformation(43.3,0,-90*pi/180,theta4_1)) * T46_1;
    T56_2 = inv(transformation(43.3,0,-90*pi/180,theta4_2)) * T46_1;
    T56_3 = inv(transformation(43.3,0,-90*pi/180,theta4_3)) * T46_2;
    T56_4 = inv(transformation(43.3,0,-90*pi/180,theta4_4)) * T46_2;
    T56_5 = inv(transformation(43.3,0,-90*pi/180,theta4_5)) * T46_3;
    T56_6 = inv(transformation(43.3,0,-90*pi/180,theta4_6)) * T46_3;
    T56_7 = inv(transformation(43.3,0,-90*pi/180,theta4_7)) * T46_4;
    T56_8 = inv(transformation(43.3,0,-90*pi/180,theta4_8)) * T46_4;
    theta5_1 = -atan2(T56_1(1,3),-T56_1(2,3));
    theta5_2 = -atan2(T56_2(1,3),-T56_2(2,3));
    theta5_3 = -atan2(T56_3(1,3),-T56_3(2,3));
    theta5_4 = -atan2(T56_4(1,3),-T56_4(2,3));
    theta5_5 = -atan2(T56_5(1,3),-T56_5(2,3));
    theta5_6 = -atan2(T56_6(1,3),-T56_6(2,3));
    theta5_7 = -atan2(T56_7(1,3),-T56_7(2,3));
    theta5_8 = -atan2(T56_8(1,3),-T56_8(2,3));

    % theta 6 (8 solutions)
    if theta5_1 ~= 0
        theta6_1 = atan2(T46_1(3,2),-T46_1(3,1));
    else
        theta6_1 = atan2(-T46_1(1,2),T46_1(2,2)) - theta4_1;
    end

    if theta5_2 ~= 0
        theta6_2 = atan2(-T46_1(3,2),T46_1(3,1));
    else
        theta6_2 = atan2(T46_1(1,2),-T46_1(2,2)) - theta4_2;
    end

    if theta5_3 ~= 0
        theta6_3 = atan2(T46_2(3,2),-T46_2(3,1));
    else
        theta6_3 = atan2(-T46_2(1,2),T46_2(2,2)) - theta4_3;
    end

    if theta5_4 ~= 0
        theta6_4 = atan2(-T46_2(3,2),T46_2(3,1));
    else
        theta6_4 = atan2(T46_2(1,2),-T46_2(2,2)) - theta4_4;
    end

    if theta5_5 ~= 0
        theta6_5 = atan2(T46_3(3,2),-T46_3(3,1));
    else
        theta6_5 = atan2(-T46_3(1,2),T46_3(2,2)) - theta4_5;
    end

    if theta5_6 ~= 0
        theta6_6 = atan2(-T46_3(3,2),T46_3(3,1));
    else
        theta6_6 = atan2(T46_3(1,2),-T46_3(2,2)) - theta4_6;
    end

    if theta5_7 ~= 0
        theta6_7 = atan2(T46_4(3,2),-T46_4(3,1));
    else
        theta6_7 = atan2(-T46_4(1,2),T46_4(2,2)) - theta4_7;
    end

    if theta5_8 ~= 0
        theta6_8 = atan2(-T46_4(3,2),T46_4(3,1));
    else
        theta6_8 = atan2(T46_4(1,2),-T46_4(2,2)) - theta4_8;
    end

    % IK solutions
    IK_solution_1 = [theta1_1,theta2_1,theta3_1,theta4_1,theta5_1,theta6_1];
    IK_solution_2 = [theta1_2,theta2_2,theta3_1,theta4_3,theta5_3,theta6_3];
    IK_solution_3 = [theta1_1,theta2_3,theta3_2,theta4_5,theta5_5,theta6_5];
    IK_solution_4 = [theta1_2,theta2_4,theta3_2,theta4_7,theta5_7,theta6_7];
    IK_solution_5 = [theta1_1,theta2_1,theta3_1,theta4_2,theta5_2,theta6_2];
    IK_solution_6 = [theta1_2,theta2_2,theta3_1,theta4_4,theta5_4,theta6_4];
    IK_solution_7 = [theta1_1,theta2_3,theta3_2,theta4_6,theta5_6,theta6_6];
    IK_solution_8 = [theta1_2,theta2_4,theta3_2,theta4_8,theta5_8,theta6_8];

    Theta_Solution = fix(IK_solution_5) * R_to_D;

end

% Transformation function
function A = transformation(d, a, alpha, theta)
    A = [ cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta)
          sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta)
          0           sin(alpha)             cos(alpha)            d
          0           0                      0                     1           ];
end

% Caluating A123
function T36 = A1_to_3(theta1,theta2,theta3)
    A1 = transformation(0    ,0    ,-90 * pi/180,theta1);
    A2 = transformation(0    ,43.2 ,0           ,theta2);
    A3 = transformation(14.9 ,-2   ,90 * pi/180 ,theta3);
    T36 = A1 * A2 * A3;
end

% Fixing theta values
function theta_list = fix(theta_list)
    for i= 1:6
        theta = theta_list(i);
        if(theta < -180)
            theta = theta + 360 ;
        elseif(theta > 180)
            theta = theta - 360;
        end
        theta_list(i) = theta;
    end
end