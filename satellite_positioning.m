% This code was written by Li Tian
% Department of Mechanical Engineering, The Hong Kong Polytechnic University, Hong Kong SAR, China
% PhD student

clear all
clc

% Input parameters 输入参数
I_filename_1 = 'eph.dat';
I_filename_2 = 'rcvr.dat';
c = 299792458; % speed of light (m/s)
We = 7.2921151467e-5; % earth's rotation rate (rad/s)
GM = 3.986005e+14; % earth's universal gravitation constant
F = -4.442807633e-10; % relativistic correction term constant
pi = 3.1415926535898;
user_initial(1,1) = -2694685.473;
user_initial(1,2) = -4293642.366;
user_initial(1,3) = 3857878.924;

% Data processing 数据处理
data_1 = load(I_filename_1);
data_2 = load(I_filename_2);
num_satellite = size(data_1,1);

   % Satellite clock error determination and position calculation
   for i = 1:num_satellite
   no = sqrt(GM)/(data_1(i,10))^3;
   n = no + data_1(i,11);  % mean motion
   
   delta_t = data_1(i,5) + data_1(i,6)*(data_1(i,1) - data_1(i,3)) + data_1(i,7) * (data_1(i,1) - data_1(i,3))^2;
   delta_t_later(i,1) = delta_t;
   t = data_1(i,1) - delta_t;
   tk = t - data_1(i,4); % time correction
   
   Mk = data_1(i,12) + n * tk;
   tran(1,1) = Mk;
   for j = 1:4
       tran(j+1,1) = Mk + data_1(i,9) * sin(tran(j,1));
   end
   Ek = tran(5,1); % eccentric anomaly
   
   vk = atan(sqrt( 1-( data_1(i,9) )^2 )*sin(Ek)/( cos(Ek)-data_1(i,9) )); % true anomaly
   
   u = vk + data_1(i,13); % argument of latitude 
   
   delta_u = data_1(i,19) * cos(2*u) + data_1(i,18) * sin(2*u); 
   delta_r = data_1(i,23) * cos(2*u) + data_1(i,22) * sin(2*u);
   delta_i = data_1(i,21) * cos(2*u) + data_1(i,20) * sin(2*u); % argement of latitude, radius and inclination correction
   
   uk = u + delta_u;
   rk = ( data_1(i,10) )^2 * ( 1- data_1(i,9) * cos(Ek) ) + delta_r;
   ik = data_1(i,15) + data_1(i,17) * tk + delta_i; % corrected argument of latitude, radius and inclination
   
   x = rk*cos(uk);
   y = rk*sin(uk); % position in orbital plane
   
   L = data_1(i,14) + ( data_1(i,16) - We ) * tk - We * data_1(i,4); % corrected longitude of ascending node
   
   satellite_position(i,1) = -(x*cos(L)-y*cos(ik)*sin(L));
   satellite_position(i,2) = -(x*sin(L)+y*cos(ik)*cos(L));
   satellite_position(i,3) = -y*sin(ik);
   
   end
   
   % User position estimation
   delta(1:4,1) = 1;
   while max(abs(delta(1:3,1))) > 10^-4
   for i = 1:num_satellite
       H(i,1) = ( user_initial(1,1) - satellite_position(i,1) ) / sqrt( (satellite_position(i,1) - user_initial(1,1))^2 + (satellite_position(i,2) - user_initial(1,2))^2 + (satellite_position(i,3) - user_initial(1,3))^2 );
       H(i,2) = ( user_initial(1,2) - satellite_position(i,2) ) / sqrt( (satellite_position(i,1) - user_initial(1,1))^2 + (satellite_position(i,2) - user_initial(1,2))^2 + (satellite_position(i,3) - user_initial(1,3))^2 );
       H(i,3) = ( user_initial(1,3) - satellite_position(i,3) ) / sqrt( (satellite_position(i,1) - user_initial(1,1))^2 + (satellite_position(i,2) - user_initial(1,2))^2 + (satellite_position(i,3) - user_initial(1,3))^2 );
       H(i,4) = -c;
       p_appox(i,1) = sqrt( (satellite_position(i,1) - user_initial(1,1))^2 + (satellite_position(i,2) - user_initial(1,2))^2 + (satellite_position(i,3) - user_initial(1,3))^2 ) - delta_t_later(i,1) * c;
   end  % calculate the H matrix
   p(1,1) = data_2(8,3);
   p(2,1) = data_2(1,3);
   p(3,1) = data_2(2,3);
   p(4,1) = data_2(3,3);
   delta = inv(H(1:4,:)) * ( p - p_appox(1:4,1) );
   user_initial(1,1) = user_initial(1,1) + delta(1,1);
   user_initial(1,2) = user_initial(1,2) + delta(2,1);
   user_initial(1,3) = user_initial(1,3) + delta(3,1);
   end
   
   fprintf('All done!\n'); 
   fprintf('The user position are determined to be %f,%f,%f',user_initial);
   