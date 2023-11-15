clear
%clc
%%%%%%%%%%%%%%%%%%%% Initialize the simulation parameters (user input) %%%%%%%%%%%%%%%%%%%%
L = 50; 
Tmax = 100; 
%N=2^15 gridpoints, one with N=2^13 gridpoints, and one with N=2^11 gridpoints.
Gridpoints = [2^11 2^13 2^15];
N = 128; %Gridpoints(3); 
dt = 0.04; 
PlotInt = 10; 
dx = L/N; 
x = (0:N-1)'*dx; 

% Coefficients for the IMEXRKCB3c scheme 
a_im = [0 0 0 0; 0 (3375509829940/4525919076317) 0 0; 0 (-11712383888607531889907/32694570495602105556248) (566138307881/912153721139) 0; 0 (673488652607/2334033219546) (493801219040/853653026979) (184814777513/1389668723319)];
a_ex = [0 0 0 0; (3375509829940/4525919076317) 0 0 0; 0 (673488652607/2334033219546) (272778623835/1039454778728) 0 ; 0 (673488652607/2334033219546) (493801219040/853653026979) (1660544566939/ 2334033219546) ];
b = [ 0; (673488652607/2334033219546); (493801219040/853653026979); (184814777513/1389668723319)];
b_im = b;
b_ex = b;
c = [0; (3375509829940/4525919076317);(272778623835/1039454778728); 1];
bhat_im = [0; (366319659506/1093160237145); (270096253287/480244073137); (104228367309/1017021570740) ];
bhat_ex = [(449556814708/1155810555193); 0; (210901428686/1400818478499); (480175564215/1042748212601)];

u = 0.15*randn(N,1); 
uhat = RC_RFFT(u,N); % User-defined FFT function
ub=0.15*randn(N,1); 
ubhat=RC_RFFT(ub,N);
s = 4;
kx = (2*pi/L)*[0:N/2-1]'; 
Aop = kx.^2 - kx.^4; % Assuming Aop is a diagonal matrix for KS equation
f = Aop;
% ... rest of the coefficients initialization remains the same ...
startTime = tic;
for k = 1:Tmax/dt
    for rk = 1:s
  %%%% ALL 3 RK SUBSTEPS %%%%
  if rk == 1
         y = uhat;
    else 
        y = uhat + (a_im(rk,rk-1)-b_im(rk-1)).*dt.* Z + (a_ex(rk,rk-1)-b_ex(rk-1)).*dt.*y;
    end
    Z = (Aop.*y)./(1-(a_im(rk,rk).*dt.*Aop));

    r=RC_RFFTinv(y + a_im(rk,rk).*dt.*Z,N);
    r=-0.5*r.*r;
    rhat = i*kx.*RC_RFFT(r,N); 
    uhat = uhat + b_im(rk).*dt.*Z + b_ex(rk).*dt.*rhat;
    ubhat = ubhat + bhat_im(rk).*dt.*Z + bhat_ex(rk).*dt.*rhat;
  
    rs(k,:)= RC_RFFTinv(uhat,N)';
    ts(k)=k*dt;
    rshat = RC_RFFTinv(ubhat,N)'; 
    
    % These variables are just used for plotting...
    % if (mod(k,PlotInt)==0)
    %     pause(0.001); RC_PlotXY(x,rs(k,:),k*dt,0,L,-1.5,1.5);
    %     %     % Uncomment the lines below to make some additional interesting plots.
    %     %     % figure(2); semilogy(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([0 3 1e-8 1e-1])
    %     %     % figure(3); loglog(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([3e-2 4 1e-8 1e-1])
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
title('[3R] IMEXRKCB3c scheme with 2^{13} gridpoints');
saveas(gcf,'3R_2_15.png','png') 
