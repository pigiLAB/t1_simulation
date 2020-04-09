function acetateT1(  )
load('Traj.mat');   % 'trajectory' variable as from Gromacs
% CSA tensor from trajectory, rotated to standard molecular basis
% and averaged over whole trajectory
standardsigma = [0.0056855360   0.002689758   0.08209500;
                 0.0485773200   0.012129073  -0.03207105;
                 0.0343148589  -0.054089100   0.42825150];
   traj = trajectory / 1E+10; % convert to m
   trajlen = size(trajectory, 1) - 1;
   MAXDEL = 200; % number of points that correlation function is
   % non-negligible (in this case, 200 points = 200ps).  It's a guess but the
   % plots will show it's a good guess  
   dt = 1E-12; % step size = 1ps
   B0 = 9.4; % T
   hbar = 6.626E-34 / (2 * pi); % hbar in Js
   gammaC = 10.705E+6 * 2 * pi; % gamma, in rad/s/T
   gammaH = 42.577E+6 * 2 * pi; % gamma, in rad/s/T
   muH = 1.4106E-26; % JA/T
   muC = muH * gammaC / gammaH;
   wC = gammaC * B0; % angular Larmor frequencies
   wH = gammaH * B0;
   mu0 = 4E-7 * pi;
   time = (1:(MAXDEL + 1)) * dt;
   moleculerange = 1:7;                     % nuclei at indices 1-7 are acetate
   waterrange = 8:(size(traj, 2)-1);        % nuclei at range 10-(end-1) are water
   counterionrange = size(traj, 2);         % last nucleus is sodium
   gamma(moleculerange) = [muC 0 0 muC muH muH muH];
   gamma(waterrange) = repmat([0 muH muH], 1, length(waterrange)/3);
   gamma(counterionrange) = 0;  % it's not really zero but who cares
   % loop over all carbons (n) and protons (np) on the pyruvate
   for n = moleculerange
       if(gamma(n) ~= muC)
           continue;
       end
       gammaII(n) = 0;
   for np = moleculerange
       if(gamma(np) ~= muH)
           continue;
       end
   sprintf('DOING CARBON %d, PROTON %d', n, np)
   G0 = zeros(1, MAXDEL + 1); % G's are correlation functions for F(0)-F(2)
   G1 = G0;
   G2 = G0;
   % calculate average value of F(t) F*(t')
   normF02 = 0;
   for j = 1:(trajlen - MAXDEL)
       [F00 F10 F20] = calcF(squeeze(traj(j, n, :)), squeeze(traj(j, np, :)));
       F00J(j) = F00*conj(F00);
       F02J(j) = F20*conj(F20);
       F01J(j) = F10*conj(F10);
       normF02 = normF02 + F00 * conj(F00) / (trajlen - MAXDEL);
       for del = 0:MAXDEL
           [F0t F1t F2t] = calcF(squeeze(traj(j + del, n, :)), squeeze(traj(j + del, np, :)));
           G0(1+del) = G0(1+del) + F00*conj(F0t) / (trajlen - MAXDEL);
           G1(1+del) = G1(1+del) + F10*conj(F1t) / (trajlen - MAXDEL);
           G2(1+del) = G2(1+del) + F20*conj(F2t) / (trajlen - MAXDEL);
       end
   end
   Gwiggle = G0 / normF02;
   % do some plots to show what the G's and their FT's look like.  The J's
   % will basically be the value of the FT at a frequency of 0
   figure;
   subplot(3,2,1); plot(time(1:(MAXDEL + 1)), real(G0), time(1:(MAXDEL + 1)), imag(G0));
   subplot(3,2,2); plot(time(1:(MAXDEL + 1)), real(G1), time(1:(MAXDEL + 1)), imag(G1));
   subplot(3,2,3); plot(time(1:(MAXDEL + 1)), real(G2), time(1:(MAXDEL + 1)), imag(G2));
   subplot(3,2,4); plot(((1:(MAXDEL + 1)) - MAXDEL/2 - 1) * (1/dt) / (MAXDEL + 1), abs(fftshift(fft(G0))));
   subplot(3,2,5); plot(((1:(MAXDEL + 1)) - MAXDEL/2 - 1) * (1/dt) / (MAXDEL + 1), abs(fftshift(fft(G1))));
   subplot(3,2,6); plot(((1:(MAXDEL + 1)) - MAXDEL/2 - 1) * (1/dt) / (MAXDEL + 1), abs(fftshift(fft(G2))));
   drawnow;
   % Do old school FT
   J0 = 0; J1 = 0; J2 = 0; Jwiggle = 0;
   for del = 1:MAXDEL
       J0 = J0 + G0(del) * exp(-1i * (wC - wH) * del * dt) * dt;
       J1 = J1 + G1(del) * exp(-1i * wC * del * dt) * dt;
       J2 = J2 + G2(del) * exp(-1i * (wC + wH) * del * dt) * dt;
       Jwiggle = Jwiggle + Gwiggle(del) * exp(-1i * wC * dt) * dt;
       if(del > 1)
           % the integral is from -inf to +inf, so add in terms with exp(+i w t) to
           % account for the negative part of the integral.  G's should be
           % symmetric
            J0 = J0 + G0(del) * exp(1i * (wC - wH) * del * dt) * dt;
            J1 = J1 + G1(del) * exp(1i * wC * del * dt) * dt;
            J2 = J2 + G2(del) * exp(1i * (wC + wH) * del * dt) * dt;
       end
   end
   % Abragam ch. VIII eq 88.  Assume the protons are all in thermal equilibrium
   gammaII(n) = gammaII(n) + (mu0/4/pi * gammaC * gammaH * hbar)^2 * 1/10 * real(5/4 * J0 + 9/4 * J1 + 3/2 * J2)
   
   deltaiso = 1/3 * (standardsigma(1,1) + standardsigma(2,2) + standardsigma(3,3));
   delta = standardsigma(3,3) - deltaiso;
   nu = (standardsigma(2,2) - standardsigma(1,1)) / delta;
   % Abragam ch. VIII eq 141.  I believe that what they call delta_z' (the
   % traceless part of the CSA tensor) is what is defined as delta above:
   % as per Spencer at al. (Int J Engng Sci 8:475-481), b_zz =
   % 1/2(a_zz+a_zz) - 1/3(a_zz) = delta.  The mu0's cancel
   % out, two from gammaC^2 and two from B0^2.
   gammaCSA(n) = 6/40 * gammaC^2 * B0^2 * delta^2 * (1 + nu^2 / 3) * Jwiggle;
   end
   end
   'gammaII (dipole part) = '
   gammaII
   'gammaCSA (CSA part) = '
   gammaCSA
   'final T1 ='
   1 ./ (gammaII + gammaCSA)
   
function [F0 F1 F2] =  calcF(r0, r1)
    % calcualte F(0)-F(2) as per Abragam Ch VIII eq 68
    delr = r1 - r0;
    absdelr = norm(delr);
    absdelr3 = absdelr * absdelr * absdelr;
    costh = delr(3) / absdelr;
    sinth = sqrt(1.0 - costh * costh);
    if(sinth < 1E-9)
        cosph = 1.0;
        sinph = 0.0;
    else
        cosph = delr(1) / (absdelr * sinth);          
        sinph = delr(2) / (absdelr * sinth);
    end
    F0 = (1-3*costh*costh)/absdelr3;
    F1 = sinth*costh*(cosph-1i*sinph)/absdelr3;
    F2 = sinth*sinth*(cosph*cosph-sinph*sinph-2*1i*sinph*cosph)/absdelr3;
