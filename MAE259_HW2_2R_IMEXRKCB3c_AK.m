clear
clc
%%%%%%%%%%%%%%%%%%%% Initialize the simulation parameters (user input) %%%%%%%%%%%%%%%%%%%%
L = 50; 
Tmax = 100; 
%N=2^15 gridpoints, one with N=2^13 gridpoints, and one with N=2^11 gridpoints.
Gridpoints = [2^11 2^13 2^15];
N = Gridpoints(2); 
dt = 0.04; 
PlotInt = 10; 
dx = L/N; 
x = (0:N-1)'*dx; 

% Coefficients for the IMEXRKCB3c scheme 
b = [0, 673488652607/2334033219546, 493801219040/853653026979,  184814777513/1389668723319];
b_im = b;
b_ex = b;

a_im = [0 0 0 0; 0 (3375509829940/4525919076317) 0 0; b(1) (-11712383888607531889907/32694570495602105556248) (566138307881/912153721139) 0; b(1), b(2), b(3), b(4)];;
a_ex = [0 0 0 0; 3375509829940/4525919076317 0 0 0; b(1) 272778623835/1039454778728 0 0; b(1) b(2) 1660544566939/2334033219546 0];

c = [0, 3375509829940/4525919076317,  272778623835/1039454778728, 1];
c_im = c;
c_ex = c;
bhat_im = [0, 366319659506/1093160237145, 270096253287/480244073137, 104228367309/1017021570740];
bhat_ex = [449556814708/1155810555193, 0, 210901428686/1400818478499,480175564215/1042748212601];


u = 0.15*randn(N,1); 
uhat = RC_RFFT(u,N); % Assuming RC_RFFT is a user-defined FFT function

s = 4;

kx = (2*pi/L)*[0:N/2-1]'; 
Aop = kx.^2 - kx.^4; % Assuming Aop is a diagonal matrix for KS equation
f = Aop;
% ... rest of the coefficients initialization remains the same ...
startTime = tic;
for k = 1:Tmax/dt
  %%%% ALL 3 RK SUBSTEPS %%%%
  for rk = 1:s
    uhat(fix(N/3)+1:end)=0;  % Dealias (see Section 5.7).
    r = RC_RFFTinv(uhat,N); % Assuming RC_RFFTinv is a user-defined inverse FFT function 
    r=-0.5*r.^2; 
    rhat=i*kx.*RC_RFFT(r,N);
    g(:,rk) = rhat;

    if rk == 1
        y = uhat;
    else   
        y = uhat + ((a_im(rk,rk-1) - b(rk-1))* dt * Aop .* y) + ((a_ex(rk,rk-1) - b(rk-1)) * dt * g(:,rk));
    end 
    y = y ./ (1 - a_im(rk,rk) * dt * Aop);
    uhat = uhat + b(rk) * dt * Aop .* y + b(rk) * dt .* g(:,rk);
 
    rs(k,:)=RC_RFFTinv(uhat,N)'; ts(k)=k*dt; % These variables are just used for plotting...
%  if (mod(k,PlotInt)==0) 
%      pause(0.001); RC_PlotXY(x,rs(k,:),k*dt,0,L,-1.5,1.5);
% %     % Uncomment the lines below to make some additional interesting plots.
% %     % figure(2); semilogy(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([0 3 1e-8 1e-1])
% %     % figure(3); loglog(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([3e-2 4 1e-8 1e-1])
% end

  end 
end 

TotalTime = toc(startTime)

% Post-processing and plotting...
hold off
figure(4); 
rs(:,N+1)=rs(:,1); 
xs=[0:N]*L/N;
contour(xs,ts,rs,[.25 .75 1.25],'r-'); 
hold on; 
contour(xs,ts,rs,[-.25 -.75 -1.25],'b-.')
xlabel('x');
ylabel('u');
title('2R IMEXRKCB3c scheme with 2^{13} gridpoints');