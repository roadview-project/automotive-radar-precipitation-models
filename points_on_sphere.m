%%% Author: Henrik Toss - henrik.toss@ri.se
%%% Project: ROADVIEW - https://roadview-project.eu/
%%% Code available at: https://github.com/roadview-project
%%% Description:  Code to estimate/predict bin volume dependency on azimuth and elevation
%%% angle for planar array. Also includes code for deriving probability
%%% distribution over velocity bins


%% angle pdf
SF = 1e2; % upsamplingfactor for better estimate at edge angles and better resemblance to continous function in case e.g. MUSIC or ESPRIT is used instead of FFT digital beam forming
Phi_y = linspace(-1,1,16*SF+1); % Detection bins for phase angle along y-axis.
Phi_z = linspace(-1,1,8*SF+1); % Detection bins for phase angle along z-axis.
P_angles = ang_pdf(Phi_y,Phi_z);
[PHI_y, PHI_z] = meshgrid(Phi_y(2:end)-mean(diff(Phi_y))/2,Phi_z(2:end)-mean(diff(Phi_z))/2);
% [PHI_y, PHI_z] = meshgrid(Phi_y(1:end-1),Phi_z(1:end-1));
figure; contourf(PHI_y, PHI_z,P_angles', 'LineStyle','none')
colorbar
% figure; surf(PHI_y, PHI_z,P_angles', 'LineStyle','none')
figure; contourf(PHI_y, PHI_z,10*log10(P_angles'), 'LineStyle','none')
colorbar

TH =asin(PHI_y);
PHI = asin(PHI_z./cos(TH));
PHI(abs(imag(PHI))>0) = NaN;
% PHI_x = sqrt(1-PHI_y.^2-PHI_z.^2);
% PHI_y(abs(imag(PHI_x))>0) = NaN;
% PHI_z(abs(imag(PHI_x))>0) = NaN;
% PHI_x(abs(imag(PHI_x))>0) = NaN;
% [TH,PHI,R] = cart2sph(sqrt(1-PHI_y.^2-PHI_z.^2),PHI_y,PHI_z);

figure; contourf(TH, PHI,P_angles', 'LineStyle','none')
colorbar

figure; contourf(TH, PHI,10*log10(P_angles'*(SF)^2), 'LineStyle','none')
colorbar
xlabel('Azimuth angle [rad]')
ylabel('Elevation angle [rad]')
title('Estimated scaling factor [dB]')
axis equal
xlim([-80*pi/180,80*pi/180])
ylim([-60*pi/180,60*pi/180])
clim([-24,-16])

%% velocity pdf
v_radar = 2*[7 0.2 0];
v_drop = 1*[0,0,-2];
v = v_drop-v_radar;
V_edges = linspace(-26.5,26.5,512+1);

PV = vel_pdf(v,V_edges);

V_centers = V_edges-mean(diff(V_edges))/2;
V_centers = V_centers(1:end-1);
figure; plot(V_centers,PV)
xlabel('velocity bin [m/s]')
ylabel('probability of finding drop in bin')
title('PDF of single drop')


function PV = vel_pdf(v,V_edges)
    N = 1e7;
    [x, y, z, ~, ~] = fibonacci_sphere(N);
    y = y(x>=0);
    z = z(x>=0);
    x = x(x>=0);
    
    v_d = dot(repmat(v,length(x),1)',[x; y; z])';
    [N_v, ~] = histcounts(v_d,V_edges);
    PV = N_v/sum(N_v);
end

function P_angles = ang_pdf(Phi_y_edges,Phi_z_edges)
    N = 1e7;
    [x, y, z, ~, ~] = fibonacci_sphere(N);
    y = y(x>=0);
    z = z(x>=0);
%     x = x(x>=0);
    
    [N,~,~] = histcounts2(y,z,Phi_y_edges,Phi_z_edges);
    P_angles = N/sum(sum(N));
end

function [x, y, z, theta, phi] = fibonacci_sphere(num_points)
    golden_angle = pi * (3 - sqrt(5)); % Golden angle for uniform distribution
    indices = 1:num_points;
    theta = golden_angle * indices; % Azimuthal angle
    phi = real(acos(1 - 2 * (indices + 0.5) / num_points) -pi/2); % Elevation angle
        
    % Convert to Cartesian coordinates
    x = cos(phi).*cos(theta);
    y = cos(phi).*sin(theta);
    z = sin(phi);
end