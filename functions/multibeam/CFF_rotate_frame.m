function T = CFF_rotate_frame(roll, pitch, yaw)


% PRELIMINARY NOTES:
% - All frames in the following are following right-hand convention (XYZ).
% - All angles are measured from an axis to the vector of interest,
% following the right-hand conventions for rotations (i.e. positive X
% towards Y, Y towards Z and Z towards X). 
% - All angles are in degrees and will be stored in data in degrees. To be
% used in trigonometric equations, they will be converted on the spot in
% radians and appended with the suffix "_rad" 
%
% FRAMES OF REFERENCE:
%
% Consider the sonar frame:
% reference: center of sonar face
% X_S: sonar forward
% Y_S: sonar right
% Z_S: acoustic axis, also known as broadside, normal to the sonar face
% plane XY

% To compute a sounding position, the sonar needs to do ray bending but to
% do that it needs to figure out the launch angle to the vertical. So we
% need to change frame. 
% 
% Consider the vessel frame:
% reference: center of sonar face
% X_V / AlongTrackDistReSonar: to vessel forward
% Y_V / AcrossTrackDistReSonar: to vessel starboard
% Z_V / DepthReSonar: vertical down

% How do we get a vector with coordinates in the S frame known
% (X_S,Y_S,Z_S) into the V frame? We need rotation matrices

% Consider a vessel rolling by 10 degrees (positive port up, Y->Z) and a
% sonar setup so that the roll offset is 10 degrees (positive when left
% handside of the sonar goes up, Y -> Z). Then the total roll angle is 20
% degrees. The angle that goes from Z_V back to Z_S is -20 degrees.


% rotation matrices for frame transformation. Ie, to calculate the
% coordiantes (x',y',z') of a point M in a frame F' when its coordinates
% (x,y,z) in a frame F are known.

% get the vessel motion here
roll_mov = 45;
pitch_mov = ;
yaw_mov = 0;

% and the sonar setup here
roll_setup = 
pitch_setup = ;
yaw_setup = ;

alpha_rad = (-roll_mov - roll_setup).*pi./180;
beta_rad  = (-pitch_mov - pitch_setup).*pi./180;
gamma_rad = (-yaw_mov - yaw_setup).*pi./180;

Rx = [ [ 1         0             0      ]; ...
       [ 0  cos(alpha_rad) sin(alpha_rad) ]; ...
       [ 0 -sin(alpha_rad) cos(alpha_rad) ] ];     % rotation about the X axis (roll)

Ry = [ [ cos(beta_rad) 0 -sin(beta_rad) ]; ...
       [        0       1        0        ]; ...
       [ sin(beta_rad) 0  cos(beta_rad) ] ]; % rotation about the Y axis (pitch)

Rz = [ [  cos(gamma_rad) sin(gamma_rad) 0 ]; ...
       [ -sin(gamma_rad) cos(gamma_rad) 0 ]; ...
       [       0           0        1 ] ];     % rotation about the Z axis (yaw)
   
% frame transformation matrix is given by:
T = Rz*Ry*Rx;   

% so imagine a point that was measured 10m in the acoustic axis of the sonar:
M=[0;0;10];

% its position in the vessel frame is actually:
T*M
% -7m on the acrosstrack axis a 7m on the vertical axis

% imagine a vessel pitching (nose up) by 45 deg, the result is 
% 7m on the along track axis, 7m on the vertical
