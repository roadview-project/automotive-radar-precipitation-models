%%% Author: Henrik Toss - henrik.toss@ri.se
%%% Project: ROADVIEW - https://roadview-project.eu/
%%% Code available at: https://github.com/roadview-project
%%% Description:Code to estimate/predict expected radar backscatter from rain in the
%%% context of automotive radar, i.e. including close range and assuming
%%% ability to resolve signal source in the velocity dimension. The code
%%% also incudes simulations of rain and calculation of the coherent signal
%%% in different range and velocity bins from the simulated rain drops. For
%%% range bins, the simulation is repeated to show what statistical spread
%%% that should be expected. The signal-to-noise ratio needed for detection
%%% is not considered here.

clear
%% Physical constants, radar specific parameters, situation specific parameters
v_radar = 0*[14 0 0]; % radar velocity in inertial coordinate system [m/s]. x i is straight ahead and y is to the left i.e. y-coordinate>0 means radar velocity slightly turning left.

c = 3e8; % speed of light [m/s]
lambda = c/77e9;
k = 2*pi/lambda; % Wavenumber

rainRate = 11; % [mm/h]  %Rain rate in mm/hour. >7.6 mm up/h to 50 mm/h constitutes heavy rain fall according to https://en.wikipedia.org/wiki/Rain#Measurement

rangeRes = 0.2; % Radar range resolution.
V_bins = linspace(-26.5,26.5,512); % Set radar velocity range and resolution

eps = 12-1i*25; % approx water @77 GHz 20 deg https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/asl.77 (water/ice)
% eps = 3.15-1i*0.005; % approx ice @77 GHz -15 deg https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/asl.77 (water/ice)
Refl = abs(((sqrt(eps)-1)/(sqrt(eps)+1))^2); % reflectivity of water sphere
% Refl = 1; % reflectivity of PEC sphere
% complex refractive index for Mie calculations
eps_abs = sqrt(real(eps).^2+imag(eps).^2);
kappa = sqrt((eps_abs-real(eps))/2);
n = sqrt((eps_abs+real(eps))/2);
refrel = n+1i*kappa;
% Marshall-Palmer distribution
N_0 = 8000; % [m^-3 mm^-1]
Lam = 4.1*rainRate^-0.21; % [mm^-1]
maxDropSize = 10; % maximum drop diameter in mm

%change to SI units
N_0 = N_0 *10^3; % [m^-3 m^-1]
Lam = Lam*10^3; % [m^-1]
D_max = maxDropSize*10^-3; % [m]
N_drops_theory = N_0/Lam*(1-exp(-Lam*D_max)); % Total number of drops in Unit Volume [m-3]

% analytical approximation of RCS of unit volume rain - probably the
% fastest way to reach an approximation
D_ray = lambda/pi; % [m] % diameter cutoff for rayleigh approximation
C_ray = 9*N_0*Refl*pi^5/(4*(lambda).^4);
RCS_ray = C_ray*integrateExp(6,[0 D_ray],Lam);
C_optic = N_0*Refl*pi/4;
RCS_optic = C_optic*integrateExp(2,[D_ray, D_max],Lam); 
RCS_analytical = RCS_ray+RCS_optic; % analytically derived RCS per unit volume

%% calculate approximate mie scattering curve for interpolation insted of calculating for each drop (faster)
r = logspace(-12,-1,1000);
sigma = zeros(size(r)); %initiate RCS vector
for jj =1:length(r)
    % radius = 2.5e-3;
    radius = r(jj);
    x=2*pi*radius/lambda;
    [s1,s2,qext,qsca,qback,gsca]=bhmie(x,refrel,1); % only care about backscatter so one angle is enough
    sigma(jj) = qback*pi*radius^2; % RCS = backscatter efficiency* cross section area 
end
figure; loglog(k*r,sigma./(pi*r.^2))
xlim([1e-9, 1e3])
ylim([1e-35, 1e5])
grid on
xlabel('kD/2 [-]')
ylabel('\sigma/A_{drop} [-]')

% Calculation of RCS of unit volume rain based on interpolation of Mie
% backscattering
D = linspace(0, 1e-2, 1e4); % drop diameters
RCS_Mie = zeros(size(D));
for ii = 1:length(D)
            RCS_Mie(ii) = interp1(r,sigma, D(ii)/2,'linear','extrap'); %use interpolation to speed up calculations
end
RCS_integrated = sum(N_0*exp(-Lam*D).*RCS_Mie)*mean(diff(D));%/N_drops_theory;

% plot simplest form of range-RCS relation (antenna near field effects not included)
ranges = linspace(1,100)+rangeRes*0.05;
Vol_new = 2/3*pi*((ranges+rangeRes/2).^3-(ranges-rangeRes/2).^3);% Estimate volume of range bins over entire half sphere
figure("Name","RCS - range"); 
plot(ranges,10*log10(RCS_integrated*Vol_new.*(ranges.^2./(ranges.*(ranges-rangeRes)))*pi/2/2), "LineWidth",2)
xlim([0,100])
xlabel('range bin [m]')
xlabel('range bin [m]')
ylabel('<total range bin RCS> [dBm^2]')
grid on
legend(sprintfc('R = %d mm/h', rainRate),'Location','southeast')

%% calculate theoretical probability of finding a drop in a specific velocity bin
%% assume relative velocity can be separated into two components - one parallel with the radar direction and one perpendicular to it (drops fall down)
D = linspace(0, 1e-2, 10000); % drop diameters
N = N_0*(exp(-Lam*D(1:end-1))-exp(-Lam*D(2:end)))/Lam/N_drops_theory; % drop size probability distribution divided into diameter bins
N = [N N(end)]; % make sure vectors are the same size
v_bins_N = dropVelocity(D); % approximate average velocity of each drop size bin
PA = zeros(length(D),length(V_bins)); % probability distribution of finding drop in specific velocity bin
RCS_Vdrop = interp1(r,sigma, D/2,'linear','extrap');

V_edges = linspace(-26.5,26.5,length(V_bins)+1);

% use quasi Monte Carlo integration for fast generation of velocity pdf
for i = 1:length(v_bins_N)
    PV = vel_pdf(-v_radar+[0,0,-v_bins_N(i)],V_edges,1e5); % pdf over velocity bins
    PA_tmp = PV;
    if sum(PA_tmp)>0
        PA_tmp = PA_tmp/sum(PA_tmp); % make sure total probability is 1 (or 0)
    end
    PA_tmp = PA_tmp*N(i); % scale probability distribution by a factor corresponding to fraction of drops it applies to 
    PA(i,:) = PA_tmp;
end
RCS_avg = sum(PA.*RCS_Vdrop'./sum(PA,1),1);
PA = sum(PA); % sum probability of finding drop of any size in certain velocity bin
PA = PA/sum(PA); % Normalize probability to 1.

%% Monte Carlo simulation of selected ranges
ranges = rangeRes*1.5:rangeRes*4:5+rangeRes*0.5; % test model for different volume sizes and see how distribution of RCS scales with number of drops. Avoid range (0,rangeRes) to avoid division by zero
Vol_bin = 2/3*pi*((ranges+rangeRes/2).^3-(ranges-rangeRes/2).^3);% Estimate volume of range bins over entire half sphere
N_Iter = 500;
signal_ranges = zeros(length(ranges),N_Iter);
RCS_ranges = zeros(length(ranges),N_Iter);
RCS_mean_Mie = zeros(length(ranges),N_Iter);
legendCell = cell(length(ranges)*2,1);

for kk = 1:length(ranges)
    legendCell{kk*2-1} = num2str(round(ranges(kk),1),'%.2f m - Monte Carlo');
    legendCell{kk*2} = 'theoretical distribution';
    NDrops = round(((2*ranges(kk)+rangeRes).^3)*N_drops_theory); % Theoretical smallest number of drops needed to be simulated in the (cubic) volume containing sphere of radius range+rangeRes
    RCS_rain_V_MC = 0;
    for nn = 1:N_Iter
        % Generate random drop sizes from probability function of D
        U = random('Uniform',0,1,NDrops,1);
        D = -1/Lam*log(1-U*(1-exp(-Lam*D_max)));
        
        X = random('Uniform',-ranges(kk),ranges(kk), size(D));
        Y = random('Uniform',-ranges(kk),ranges(kk), size(D));
        Z = random('Uniform',-ranges(kk),ranges(kk), size(D));
        
        % get cut out of range bin for typical number of drops in range bin
        D(X<=0) = [];
        Y(X<=0) = [];
        Z(X<=0) = [];
        X(X<=0) = [];

        R = sqrt(X.^2+Y.^2+Z.^2);
        D = D(R<ranges(kk)+rangeRes/2 & R>=ranges(kk)-rangeRes/2);
        X = X(R<ranges(kk)+rangeRes/2 & R>=ranges(kk)-rangeRes/2);
        Y = Y(R<ranges(kk)+rangeRes/2 & R>=ranges(kk)-rangeRes/2);
        Z = Z(R<ranges(kk)+rangeRes/2 & R>=ranges(kk)-rangeRes/2);
        R = R(R<ranges(kk)+rangeRes/2 & R>=ranges(kk)-rangeRes/2);
        
         % Calculate size parameter
        x = k * D / 2;
        RCS_Mie = zeros(size(D));
        for ii = 1:length(x)
            RCS_Mie(ii) = interp1(r,sigma, D(ii)/2,'linear','extrap'); %use interpolation to speed up calculations
        end
        RCS_mean_Mie(kk,nn) = mean(RCS_Mie,'omitnan');
        
        RCS_P = RCS_Mie;
        % Alternative, less accurate but faster approximation of drop RCS
%         RCS_P = Refl*pi*(D/2).^2; % RCS of rain drop. Optical region assumed. if drop assumed to be good conductor set Refl = 1. Likely not good estimate for smaller drops. ToDo - apply Mie scattering theory
%         RCS_P(D*pi/lambda<1) = RCS_P(D*pi/lambda<1).*9.*(pi*D(D*pi/lambda<1)/lambda).^4;
        
        S_P = dropVelocity(D); % terminal velocity of each drop.
               
        dir_V_P = [0 0 -1]; % rain fall direction. Assume mostly down (z = -1), possible to add some x and y components for wind etc.
        dir_V_P = dir_V_P/norm(dir_V_P);
        V_P = S_P*dir_V_P;
        
        V_P_prim = V_P-v_radar; % rain velocity in radar coordinate system
        R_hat = [X Y Z]./R;
        V_rho = dot(V_P_prim',R_hat')'; % relative radial velocity. This is the velocity component detectable/reported by the radar
        
        P_Tx = 1;% arbitrary transmit power - will be cancelled later as it is related to the receive power. For now only added for clarity
        [TH,PHI,R] = cart2sph(X,Y,Z);
%         D_Tx = radDiag_patch(TH, PHI,'th');
%         D_Tx = D_Tx/sum(D_Tx)*length(D_Tx); % make sure the energy of the antenna is preserved
        D_Tx = 1; % assume isotropic antenna
        D_Rx = D_Tx; % assume antenna pattern is the same for Tx and Rx elements
        Aeff_Rx = lambda'.^2/(4*pi).*D_Rx; % effective antenna area

        W_i = P_Tx*D_Tx./(4*pi*R.^2);
        W_s = W_i.*RCS_P./(4*pi*R.^2);
        P_Rx = Aeff_Rx.*W_s;
        amp_ADC = sqrt(P_Rx./P_Tx); % amplitude registered by the ADC
        phi_ADC = angle(exp(1i*2*pi*2*R./lambda));
        signal_Rx = sum(amp_ADC.*exp(1i*2*pi*2*R./lambda));
        
        RCS_lin = abs(signal_Rx).^2*(4*pi.*ranges(kk).^2).^2; 
        
        signal_ranges(kk,nn) = signal_Rx;
        RCS_rain = 10*log10(abs(signal_Rx*4*pi./lambda).^2.*ranges(kk).^4*4*pi); % Signal as RCS compensated by the range of detection range bin
        RCS_ranges(kk,nn) = RCS_rain;
        
        if kk==length(ranges)
            %% Separate detections into velocity bins (now only done for last range bin)
            [~, I_vel] = min(abs(V_rho'-V_bins'));
            
            signal_Rx_V = zeros(size(V_bins)); %calculate how signal is separated over velocity bins
            RCS_v_mean = zeros(size(V_bins)); %investigate how drop velocity plays into how the velocities are distributed
            for ii = 1:length(V_bins)
                if isnan(mean(RCS_P(I_vel==ii)))
                    RCS_v_mean(ii) = 0;
                else
                    RCS_v_mean(ii) = mean(RCS_P(I_vel==ii));
                end
            %     amp_temp = sqrt(RCS_P(I_vel==ii)./R(I_vel==ii).^4).*radDiag_patch(TH(I_vel==ii), PHI(I_vel==ii),'th');
                amp_temp = sqrt(RCS_P(I_vel==ii)./R(I_vel==ii).^4);
                phi_temp = angle(exp(1i*2*pi*2*R(I_vel==ii)./lambda));
                signal_Rx_V(ii) = sum(amp_temp.*exp(1i*phi_temp)); 
            end
            RCS_rain_V = abs(signal_Rx_V).^2.*ranges(end).^4; % Signal as RCS compensated by the range of detection range bin calculated from theoretically recieved signal of entire simulated velocity bin
            RCS_rain_V_MC = RCS_rain_V_MC + RCS_rain_V; % Signal as RCS compensated by the range of detection range bin calculated from theoretically recieved signal of entire simulated velocity bin
        end
    end
end
RCS_rain_V_MC=RCS_rain_V_MC/N_Iter;
%% plotting of results from Monte Carlo simulation
RCS_Mie_theory = mean(RCS_mean_Mie(end,:)); % largest number of drops gives best estimate for mean RCS of single drop. RCS_Mie_theory*N_drops_theory can be compared to RCS_analytical and RCS_integrated to understand the error of those respective estimation methods and the number of points to chose for simulation

P_chi2 = chi2rnd(2,10000,1).*RCS_integrated*Vol_bin(1)*(ranges(1)^2/(ranges(1)^2-rangeRes^2/4))*pi/2/2; % mean of chi-squared should be k = 2
c_test= jet(length(ranges));
figure("Name","cdf:s per range"); h = cdfplot(RCS_ranges(1,:));
h.Color = c_test(1,:);
hold on
h = cdfplot(10*log10(P_chi2));
h.LineWidth = 2;
h.LineStyle = '--';
h.Color = c_test(1,:);
for jj = 2:length(ranges)
    P_chi2 = chi2rnd(2,10000,1).*RCS_integrated*Vol_bin(jj)*((ranges(jj)^2/(ranges(jj)^2-rangeRes^2/4)))*pi/2/2; % mean of chi-squared should be k = 2
    h = cdfplot(RCS_ranges(jj,:)); 
    h.Color = c_test(jj,:);
    h = cdfplot(10*log10(P_chi2));
    h.LineWidth = 2;
    h.LineStyle = '--';
    h.Color = c_test(jj,:);
end
% hleg = legend(legendCell, 'Location','EastOutside');
hleg = legend(legendCell, 'Location','Best');
title(hleg,'range bin')
xlabel('x = range bin RCS from rain [dBm^2]')

%% Separate detections into velocity bins (here only done for last range bin and last iteration of previous Monte Carlo simulation)
% [M, I_vel] = min(abs(V_rho'-V_bins'));
% V_det = V_bins(I_vel);
% 
% signal_Rx_V = zeros(size(V_bins)); %calculate how signal is separated over velocity bins
% RCS_v_mean = zeros(size(V_bins)); %investigate how drop velocity plays into how the velocities are distributed
% for ii = 1:length(V_bins)
%     if isnan(mean(RCS_P(I_vel==ii)))
%         RCS_v_mean(ii) = 0;
%     else
%         RCS_v_mean(ii) = mean(RCS_P(I_vel==ii));
%     end
% %     amp_temp = sqrt(RCS_P(I_vel==ii)./R(I_vel==ii).^4).*radDiag_patch(TH(I_vel==ii), PHI(I_vel==ii),'th');
%     amp_temp = sqrt(RCS_P(I_vel==ii)./R(I_vel==ii).^4);
%     phi_temp = angle(exp(1i*2*pi*2*R(I_vel==ii)./lambda));
%     signal_Rx_V(ii) = sum(amp_temp.*exp(1i*phi_temp)); 
% end
% RCS_rain_V = abs(signal_Rx_V).^2.*ranges(end).^4; % Signal as RCS compensated by the range of detection range bin calculated from theoretically recieved signal of entire simulated velocity bin

%% plot results from one simulation together with analytically predicted distribution over velocity 
figure('Name','RCS_lin-velocity'); 
scatter(V_bins,RCS_rain_V,18,'filled')
hold on
plot(V_bins,RCS_avg.*PA*N_drops_theory*Vol_bin(end),'--','LineWidth',2)
plot(V_bins,RCS_rain_V_MC,':','LineWidth',2)
legend('calculated at receiver','analytical prediction','Monte Carlo mean',Location='best')
ylabel('RCS [m^2]')
xlabel('velocity bin')
grid on

figure('Name','RCS_log-velocity'); 
scatter(V_bins,10*log10(RCS_rain_V),18,'filled')
hold on
% plot(V_bins,10*log10(RCS_avg.*PA*N_drops_theory*Vol_bin(end)),'--','LineWidth',2)
plot(V_bins,10*log10(RCS_avg.*PA*N_drops_theory*Vol_bin(end)),'LineWidth',2)
% plot(V_bins,10*log10(RCS_rain_V_MC),':','LineWidth',2)
scatter(V_bins,10*log10(RCS_rain_V_MC),18)
legend('calculated at receiver','analytical prediction','Monte Carlo mean',Location='south')
ylabel('RCS [dBm^2]')
xlabel('velocity bin [m/s]')
grid on
ylim([-60,-25])

figure('Name','drop RCS_lin-velocity'); plot(V_bins,RCS_v_mean ,V_bins,RCS_avg,'--')
legend('RCS from signal','analytical prediction',Location='north')
ylabel('average drop RCS [m^2]')
xlabel('velocity bin [m/s]')

grid on

figure('Name','drop RCS_log-velocity'); plot(V_bins,10*log10(RCS_v_mean), V_bins,10*log10(RCS_avg),'--')
legend('RCS from signal', 'analytical prediction',Location='north')
ylabel('average drop RCS [dBm^2]')
xlabel('velocity bin [m/s]')
grid on

function S_P = dropVelocity(D)
%%% drop velocity according to [D. Koutsoyiannis, A Langousis, Treatise on Water Science, 2011]
    a1 =2115; % [cm^(1-b)s^-1]
    b1 = 0.8;
    S_P = a1*(0.1*D).^b1; % terminal velocity of each drop.
%     a2 = 1767; % alternative parameters suggested by Atlas and Ulbrich 1977 for drops with diameter 0.5<D<5 mm
%     b2= 0.67;
end
function PV = vel_pdf(v,V_edges,N)
%%% function to generate N uniformly distributed points on sphere, with the
%%% same given velocity vector v and separate their projections onto their
%%% respective position vector into bins with bin edges V_edges
    [x, y, z, ~, ~] = fibonacci_sphere(N);
    y = y(x>=0);
    z = z(x>=0);
    x = x(x>=0);
    
    v_d = dot(repmat(v,length(x),1)',[x; y; z])';
    [N_v, ~] = histcounts(v_d,V_edges);
    PV = N_v/sum(N_v);
end

function [x, y, z, theta, phi] = fibonacci_sphere(num_points)
%%% Function to generate num_points number of uniformly distributed points
%%% on unit sphere.
    golden_angle = pi * (3 - sqrt(5)); % Golden angle for uniform distribution
    indices = 1:num_points;
    theta = golden_angle * indices; % Azimuthal angle
    phi = real(acos(1 - 2 * (indices + 0.5) / num_points) -pi/2); % Elevation angle
        
    % Convert to Cartesian coordinates
    x = cos(phi).*cos(theta);
    y = cos(phi).*sin(theta);
    z = sin(phi);
end

function I = integrateExp(N,X,Lam)
%%% anlytical integral of exponential function of Lam^N over interval X
    x1 = X(1);
    x2 = X(2);
    n = 0:N;
    u1 = -Lam*x1;
    u2 = -Lam*x2;
    I1 = 1/(-Lam)^(N+1)*exp(u1)*sum(u1.^n*factorial(N)./factorial(n).*(-1).^(N-n));
    I2 = 1/(-Lam)^(N+1)*exp(u2)*sum(u2.^n*factorial(N)./factorial(n).*(-1).^(N-n));
    I = I2-I1;
end


%%% Bohren and Huffman model and code taken from
%%% http://scatterlib.wikidot.com/mie
%%% adapted to only our point of interest

function [s1,s2,qext,qsca,qback,gsca]=bhmie(x,refrel,nang)
% Calculated based on Mie scattering theory  
% input:
%      x - size parameter =2pi*lambda/radius
%      refrel - refraction index in complex form for example:  1.5+0.02*i;
%      nang - number of angles for S1 and S2 function in range from 0 to pi/2
% output:
%        S1, S2 - funtion which coresponted to phase function
%        Qext - extinction efficiency
%        Qsca - scattering efficiency 
%        Qback -backscatter efficiency
%        gsca- asymmetry parameter

mxnang=1000;
nmxx=150000;

if (nang < 2)
   nang = 2;
end

s1=zeros(1,2*nang-1);
s2=zeros(1,2*nang-1);
s1_1=zeros(1,nang-1);    
s2_1=zeros(1,nang-1);
s1_2=zeros(1,nang-1);  
s2_2=zeros(1,nang-1);
d=zeros(1,nmxx);
pi0=zeros(1,nang);
pi1=ones(1,nang);

% if (nang > mxnang)
%     disp('error: nang > mxnang in bhmie')
% return
% end

dx = x;

drefrl = refrel;
y = x*drefrl;
ymod = abs(y);

%    Series expansion terminated after NSTOP terms
%    Logarithmic derivatives calculated from NMX on down

xstop = x + 4.*x^0.3333 + 2.;
nmx = max(xstop,ymod) + 15;
nmx=fix(nmx);
 
% BTD experiment 91/1/15: add one more term to series and compare resu<s
%      NMX=AMAX1(XSTOP,YMOD)+16
% test: compute 7001 wavelen>hs between .0001 and 1000 micron
% for a=1.0micron SiC grain.  When NMX increased by 1, only a single
% computed number changed (out of 4*7001) and it only changed by 1/8387
% conclusion: we are indeed retaining enough terms in series!
nstop = xstop;
%
% if (nmx > nmxx)
%     'error: nmx > nmxx=', nmxx, ' for |m|x=', ymod
%     return
% end
% Require NANG.GE.1 in order to calculate scattering intensities
% dang = 0.;
% if (nang > 1)
%     dang = .5*pi/ (nang-1);
% end
dang = .5*pi/ (nang-1);
amu = (0:nang-1)*dang;
amu = cos(amu);

% Logarithmic derivative D(J) calculated by downward recurrence
% beginning with initial value (0.,0.) at J=NMX

%*** Riccati-Bessel functions with real argument X
%    calculated by upward recurrence
%
psi0 = cos(dx);
psi1 = sin(dx);
chi0 = -sin(dx);
chi1 = cos(dx);
xi1 = psi1-chi1*1i;
qsca = 0.;
gsca = 0.;
p = -1;

for n=nmx:-1:2
    d(n-1) = ((n)/y) - (1./ (d(n)+(n)/y));
end

for n=1: nstop
    fn = (2.*n+1.)/ (n* (n+1.));
% for given N, PSI  = psi_n        CHI  = chi_n
%              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
%              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
% Calculate psi_n and chi_n
    psi = (2.*n-1.)*psi1/dx - psi0;
    chi = (2.*n-1.)*chi1/dx - chi0;
    xi = psi-chi*1i;
%
%*** Store previous values of AN and BN for use
%    in computation of g=<cos(theta)>
    if (n > 1) 
        an1 = an;
        bn1 = bn;
    end
%
%*** Compute AN and BN:
    an = (d(n)/drefrl+n/dx)*psi - psi1;
    an = an/ ((d(n)/drefrl+n/dx)*xi-xi1);
    bn = (drefrl*d(n)+n/dx)*psi - psi1;
    bn = bn/ ((drefrl*d(n)+n/dx)*xi-xi1);
%
%*** Augment sums for Qsca and g=<cos(theta)>
    qsca = qsca + (2*n+1)* (abs(an)^2+abs(bn)^2);
    gsca = gsca + ((2*n+1)/ (n*(n+1)))*(real(an)*real(bn) + imag(an)*imag(bn));

    if (n > 1)
        gsca = gsca + ((n^2-1)/n)*(real(an1)*real(an) + imag(an1)*imag(an) + real(bn1)*real(bn) + imag(bn1)*imag(bn));
    end 
%
%*** Now calculate scattering intensity pattern
%    First do angles from 0 to 90
    pi_l = pi1;
    tau = n*amu.*pi_l - (n+1)*pi0;
    s1_1 = s1_1 + fn*(an.*pi_l+bn.*tau);
    s2_1 = s2_1 + fn*(an.*tau+bn.*pi_l);
%
%*** Now do angles greater than 90 using PI and TAU from
%    angles less than 90.
%    P=1 for N=1,3,...% P=-1 for N=2,4,...
    p = -p;
    s1_2 = s1_2 + fn*p*(an.*pi_l-bn.*tau);
    s2_2 = s2_2 + fn*p*(bn.*pi_l-an.*tau);

    s1 = [s1_1 s1_2(end-1:-1:1)];
    s2 = [s2_1 s2_2(end-1:-1:1)];

    psi0 = psi1;
    psi1 = psi;
    chi0 = chi1;
    chi1 = chi;
    xi1 = psi1-chi1*1i;
%
%*** Compute pi_n for next value of n
%    For each angle J, compute pi_n+1
%    from PI = pi_n , PI0 = pi_n-1
    pi1 = ((2*n+1)*amu.*pi_l - (n+1)*pi0)/n;
    pi0 = pi_l;
end
%
%*** Have summed sufficient terms.
%    Now compute QSCA,QEXT,QBACK,and GSCA
gsca = 2*gsca/qsca;
qsca = (2/dx.^2)*qsca;
qext = (4/dx.^2).*real(s1(1));
qback = 4*(abs(s1(2*nang-1))/dx)^2;

% ss1=s1;
% ss2=s2;
% clear s1 s2
% a=find(ss1~=0);
% n=max(a);
% 
% s1=ss1(1:n);
% s2=ss2(1:n);
end

function E = radDiag_patch(th, phi, pol) % simple model for antenna diagram square patch antenna. Returns polarization of choice pol.
    % function [E_TH, E_PHI] = radDiag_patch(th, phi) 
    %%% Radiation pattern square patch antenna according to https://www.antenna-theory.com/antennas/patches/antenna.php
    lambda = 3e8/77e9;
    
    k = 2*pi/lambda;
    L = lambda/2;
    W = L;
    
    
    E_TH = sin(k*W*sin(th).*sin(phi)/2)./(k*W*sin(th).*sin(phi)/2).*cos(k*L*sin(th).*cos(phi)/2).*cos(phi); % use this polarization
    E_PHI = -sin(k*W*sin(th).*sin(phi)/2)./(k*W*sin(th).*sin(phi)/2).*cos(k*L*sin(th).*cos(phi)/2).*cos(th).*sin(phi);
    
    % handle cases of division by zero
    E_TH(th == 0) = cos(k*L*sin(th(th == 0)).*cos(phi(th == 0))/2).*cos(phi(th == 0)); 
    E_TH(phi == 0) = cos(k*L*sin(th(phi == 0)).*cos(phi(phi == 0))/2).*cos(phi(phi == 0));

    E_PHI(th == 0) = -cos(k*L*sin(th(th == 0)).*cos(phi(th == 0))/2).*cos(th(th == 0)).*sin(phi(th == 0));
    E_PHI(phi == 0) = -cos(k*L*sin(th(phi == 0)).*cos(phi(phi == 0))/2).*cos(th(phi == 0)).*sin(phi(phi == 0));
    
    if strcmp(pol,'th')
        E = E_TH;
    elseif strcmp(pol,'phi')
        E = E_PHI;
    else
        print(pol + ' is not a valid choice for polarization. Returning isotropic antenna. If other is desired chose pol = th, or pol = phi')
        E = ones(size(E_TH)); 
    end
end