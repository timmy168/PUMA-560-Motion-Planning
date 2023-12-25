% -------------------------------------------------------------------------
% Function Defination

% Kinematic Function
function T6 = foward_kinematic(theta_list) % theta => list of joint angle

    % Degree and Radius Transformation
    D_to_R = pi / 180;
    R_to_D = 180 / pi;
    
    % DH Model - Puma 560 in (M)
    % d_list = [0, 0, 0.149, 0.433, 0, 0];
    % a_list = [0, 0.432, -0.02, 0, 0, 0];
    % alpha_list = [-90 ,0 ,90 ,-90 ,90 ,0];

    % A = transformation(d, a, alpha, theta)
    A1 = transformation(0    , 0    , -90 * D_to_R, theta_list(1) * D_to_R);
    A2 = transformation(0    , 43.2 , 0           , theta_list(2) * D_to_R);
    A3 = transformation(14.9 , -2   ,  90 * D_to_R, theta_list(3) * D_to_R);
    A4 = transformation(43.3 , 0    , -90 * D_to_R, theta_list(4) * D_to_R);
    A5 = transformation(0    , 0    ,  90 * D_to_R, theta_list(5) * D_to_R);
    A6 = transformation(0    , 0    , 0           , theta_list(6) * D_to_R);

    % Calculate the end-effector transformation matrix
    T6 = A1 * A2 * A3 * A4 * A5 * A6;
    nx = T6(1,1); ny = T6(2,1); nz = T6(3,1);
    ox = T6(1,2); oy = T6(2,2); oz = T6(3,2);
    ax = T6(1,3); ay = T6(2,3); az = T6(3,3);
    px = T6(1,4); py = T6(2,4); pz = T6(3,4);
    x = px; y = py; z = pz;

    phi = atan2(oz, -nz) * R_to_D;
    theta = atan2(sqrt(ax^2 + ay^2), az) * R_to_D;
    psi = atan2(ay, ax) * R_to_D;
    p = [x, y, z, phi, theta, psi];

end

% Transformation function
function A = transformation(d, a, alpha, theta)
     A = [ cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta)
           sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta)
           0           sin(alpha)             cos(alpha)            d
           0           0                      0                     1           ];
end