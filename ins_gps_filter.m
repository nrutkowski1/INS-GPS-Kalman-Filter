

function ins_gps_filter()


    load gps_ins_synthetic_data.mat...
        accel_readings accel_readings_noisy...
        euler_angles_readings euler_angles_readings_noisy...
        gps_readings gps_readings_noisy gps_time_stamp...
        gyro_readings gyro_readings_noisy...
        std_dev_accel...
        std_dev_gps_posn...
        std_dev_gps_vel...
        std_dev_gyro...
        std_dev_mag...
        time_stamp

%      load sensorlog_20191010_171321.mat...
%          Acceleration...
%          AngularVelocity...
%          Orientation...
%          Position
%      
%      TT = synchronize(Orientation, AngularVelocity, Acceleration, Position, 'regular', 'spline','TimeStep', seconds(0.01));
%     
%     euler_angles_readings_noisy = [TT.Z_Orientation*pi/180 TT.Y_Orientation*pi/180 TT.X_Orientation*pi/180]';
%     gyro_readings_noisy = [TT.Z_AngularVelocity TT.Y_AngularVelocity TT.X_AngularVelocity]';
%     accel_readings_noisy = [TT.X_Acceleration TT.Y_Acceleration TT.Z_Acceleration]';
%     
%     time_stamp = TT.Timestamp;
%     
%     
%     psi = euler_angles_readings_noisy(1, 1);
%     theta = euler_angles_readings_noisy(2, 1);
%     phi = euler_angles_readings_noisy(3, 1);
% 
%     Rei = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1]*...
%           [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)]*...
%           [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];

%     heading = TT.course(1)*pi/180;
% 
%     Reh = [cos(heading) sin(heading) 0; -sin(heading) cos(heading) 0; 0 0 1];
% 
%     ve = Rei'*[TT.speed(1); 0; 0];
%     vh = Reh*ve;
%     
%     wgs84 = referenceEllipsoid('wgs84');
%     
%     [x, y, z] = geodetic2ned(TT.latitude(1), TT.longitude(1), TT.altitude(1), TT.latitude(1), TT.longitude(1), TT.altitude(1), wgs84);
%     
%     gps_readings_noisy = [x; y; z; vh];
    
     n_states = 9;
     C = eye(n_states);
     time_pts = time_stamp;
     
     % dt needs to be in for loop especially for our collected data
     
     var_gps  = std_dev_gps_posn^2;
     var_gpv = std_dev_gps_vel^2;
%      
%      std_dev_posv = 3.872/2;
%      std_dev_posh = 0.9455;
%      
%      std_dev_gpv = 0.003;
%      std_dev_mag = 2.385*pi/180;
%      std_dev_gyro = 0.2932;
%      std_dev_accel = .000349;
%      
%      var_gpsv = std_dev_posv^2;
%      var_gpsh = std_dev_posh^2;
%      var_gpv = std_dev_gpv^2;
     
     R = [var_gps var_gps var_gps var_gpv var_gpv var_gpv std_dev_mag^2 std_dev_mag^2 std_dev_mag^2].*eye(n_states);
     G = [[std_dev_accel^2; std_dev_accel^2; std_dev_accel^2; std_dev_gyro^2; std_dev_gyro^2; std_dev_gyro^2].*eye(6) zeros(6); zeros(6) zeros(6)];
     
     % initialize 
     x_hat = [gps_readings_noisy(:, 1); euler_angles_readings_noisy(:, 1)];
     P = R;
     
     omega_it = [0 -2*pi/86400 0; 2*pi/86400 0 0; 0 0 0];
     omega_itt = [0; 0; 2*pi/86400];
     gt = [0; 0; 9.81];
     
     x_hat_rec = zeros(9, numel(time_pts));
     x_hat_rec(:, 1) = x_hat;
    
%      dt = 0.01;
    
     for n = 2:numel(time_pts)
         
%         psi = euler_angles_readings_noisy(1, n);
%         theta = euler_angles_readings_noisy(2, n);
%         phi = euler_angles_readings_noisy(3, n);
% 
%         Rei = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1]*...
%               [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)]*...
%               [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
% 
%         heading = TT.course(n)*pi/180;
% 
%         Reh = [cos(heading) sin(heading) 0; -sin(heading) cos(heading) 0; 0 0 1];
%         
%         ve = Rei'*[TT.speed(n); 0; 0];
%         vh = Reh*ve;
% 
%         gps_readings_noisy(4:6) = vh;
         
         
%         gps_readings_noisy(1:3) = geodetic2ned(TT.latitude(n), TT.longitude(n), TT.altitude(n), TT.latitude(1), TT.longitude(1), TT.altitude(1), wgs84);

         
         psi_m = euler_angles_readings_noisy(1, n);
         theta_m = euler_angles_readings_noisy(2, n);
         phi_m = euler_angles_readings_noisy(3, n);
         
         Rte_f = [1 0 0; 0 cos(phi_m) sin(phi_m); 0 -sin(phi_m) cos(phi_m);];
         Rte_s = [cos(theta_m) 0 -sin(theta_m); 0 1 0; sin(theta_m) 0 cos(theta_m)];
         Rte_t = [cos(psi_m) sin(psi_m) 0; -sin(psi_m) cos(psi_m) 0; 0 0 1];
             
         Rte = Rte_f*Rte_s*Rte_t;
        
         Ret = Rte';
         
         H321 = [-sin(theta_m) 0 1;
                 sin(phi_m)*cos(theta_m) cos(phi_m) 0;
                 cos(phi_m)*cos(theta_m) -sin(phi_m) 0];
             
         dt = time_pts(n) - time_pts(n-1);
         
         A = [eye(3) dt*eye(3) zeros(3);
              zeros(3) (eye(3) - 2*omega_it*dt) zeros(3);
              zeros(3) zeros(3) eye(3)];
         B = [zeros(3) zeros(3) zeros(3) zeros(3);
              Ret zeros(3) eye(3) zeros(3);
              zeros(3) inv(H321) zeros(3) -H321\Rte];
          
         Q = B*G*B';
         
         f_e = accel_readings_noisy(:, n - 1);
         omega_iee = gyro_readings_noisy(:, n - 1);
         u = [f_e; omega_iee; gt; omega_itt];
         x_minus = A*x_hat + B*u;
         
         P_minus = A*P*A' + Q;
        
         % Kalman gain
         K = P_minus*C' / (C*P_minus*C' + R);

         % measurement update
         z = [gps_readings_noisy(:, n); euler_angles_readings_noisy(:, n)];
         x_hat = x_minus + K*(z - C*x_minus);
         P = (eye(n_states) - K*C)*P_minus; 

         % sum of diagonal elements 
         trace(P);

         x_hat_rec(:, n) = x_hat;
        
     end

     figure(1);
     plot3(x_hat_rec(1, :), x_hat_rec(2, :), x_hat_rec(3, :), 'b-');
 
     figure(2);
     subplot(311);
     plot(time_stamp, x_hat_rec(1, :), 'b-', time_stamp, gps_readings(1, :), 'r-');
     ylabel('X - Position (m)');
     subplot(312);
     plot(time_stamp, x_hat_rec(2, :), 'b-', time_stamp, gps_readings(2, :), 'r-');
     ylabel('Y - Position (m)');
     subplot(313);
     plot(time_stamp, x_hat_rec(3, :), 'b-', time_stamp, gps_readings(3, :), 'r-');
     ylabel('Z - Position (m)');
     xlabel('Time (s)');
     
     figure(3);
     subplot(311);
     plot(time_stamp, x_hat_rec(4, :), 'b-', time_stamp, gps_readings(4, :), 'r-');
     ylabel('X - Velocity (m/s)');
     subplot(312);
     plot(time_stamp, x_hat_rec(5, :), 'b-', time_stamp, gps_readings(5, :), 'r-');
     ylabel('Y - Velocity (m/s)');
     subplot(313);
     plot(time_stamp, x_hat_rec(6, :), 'b-', time_stamp, gps_readings(6, :), 'r-');
     ylabel('Z - Velocity (m/s)');
     xlabel('Time (s)');
     
     figure(4);
     subplot(311);
     plot(time_stamp, x_hat_rec(7, :), 'b-', time_stamp, euler_angles_readings(1, :), 'r-');
     ylabel('\psi (rad)');
     subplot(312);
     plot(time_stamp, x_hat_rec(8, :), 'b-', time_stamp, euler_angles_readings(2, :), 'r-');
     ylabel('\theta (rad)');
     subplot(313);
     plot(time_stamp, x_hat_rec(9, :), 'b-', time_stamp, euler_angles_readings(3, :), 'r-');
     ylabel('\phi (rad)');
     xlabel('Time (s)');   

end


