clear
close
%%%%%%%%%%%%%%%%%%%% Initialize the simulation parameters (user input) %%%%%%%%%%%%%%%%%%%%
L = 50; 
Tmax = 100; 
%N=2^15 gridpoints, one with N=2^13 gridpoints, and one with N=2^11 gridpoints.
Gridpoints = [2^11 2^13 2^15];
N = 128; %Gridpoints(2); 
dt = 0.04; 
PlotInt = 10; 
dx = L/N; 
x = (0:N-1)'*dx; 

% Coefficients for the IMEXRKCB3c scheme 
a_im = [0 0 0 0; 0 (3375509829940/4525919076317) 0 0; 0 (-11712383888607531889907/32694570495602105556248) (566138307881/912153721139) 0; 0 (673488652607/2334033219546) (493801219040/853653026979) (184814777513/1389668723319)];
a_ex = [0 0 0 0; (3375509829940/4525919076317) 0 0 0; 0 (673488652607/2334033219546) (272778623835/1039454778728) 0 ; 0 (673488652607/2334033219546) (493801219040/853653026979) (1660544566939/ 2334033219546) ];
b = [ 0 (673488652607/2334033219546) (493801219040/853653026979) (184814777513/1389668723319)];
b_im = b;
b_ex = b;
c = [0; (3375509829940/4525919076317);(272778623835/1039454778728); 1];

u = 0.15*randn(N,1); 
uhat = RC_RFFT(u,N); % Assuming RC_RFFT is a user-defined FFT function

s = 4;

kx = (2*pi/L)*[0:N/2-1]'; 
Aop = kx.^2 - kx.^4; % Assuming Aop is a diagonal matrix for KS equation

% ... rest of the coefficients initialization remains the same ...
startTime = tic;
for k = 1:Tmax/dt
  %%%% ALL 3 RK SUBSTEPS %%%%
  for rk = 1:s
    uhat(fix(N/3)+1:end)=0;  % Dealias (see Section 5.7).
    if rk == 1
        y = uhat;
    else          
        y = uhat +  sum(dt * a_im(rk,1:rk-1) .* f(:,1:rk-1),2) + sum(dt * a_ex(rk,1:rk-1) .* g(:,1:rk-1),2); 
    end 
    f(:,rk) = (Aop.* y)./(1 - (a_im(rk,rk)*dt*Aop)); % IMEX update for f
    
    rhat= y+ dt * a_im(rk,rk)*f(:,rk);
    rhat(fix(N/3)+1:end)=0;
    r = RC_RFFTinv(rhat,N); % Assuming RC_RFFTinv is a user-defined inverse FFT function 
    r=-0.5*r.^2; 
    rhat=i*kx.*RC_RFFT(r,N); 
    g(:,rk) = rhat; % Explicit update for g (nonlinear term)
  end
  %%%% END OF RK LOOP %%%%
 
  %sum_b = sum(b*dt.*f(:,1:s),2) + sum(b*dt.*g(:,1:s),2); % Update uhat directly
  %sum(b(1:rk).*f(:,1:rk), 2) + sum(b(1:rk).*g(:,1:rk), 2)

  %uhat = uhat + sum_b; %bhat_im(1)*dt*f(:,1) + bhat_ex(1)*dt*g(:,1)+bhat_im(s)*dt*f(:,2) + bhat_ex(2)*dt*g(:,2)+bhat_im(s)*dt*f(:,) + bhat_ex(s)*dt*g(:,s);
  uhat = uhat + sum(b_im(1:s)'.*f(:,1:s), 2) + sum(b(1:s)'.*g(:,1:s), 2); %bhat_im(1)*dt*f(:,1) + bhat_ex(1)*dt*g(:,1)+bhat_im(s)*dt*f(:,2) + bhat_ex(2)*dt*g(:,2)+bhat_im(s)*dt*f(:,) + bhat_ex(s)*dt*g(:,s);

  %%Plotting
rs(k,:)=RC_RFFTinv(uhat,N)'; ts(k)=k*dt; % These variables are just used for plotting...
% if (mod(k,PlotInt)==0) 
%     pause(0.001); RC_PlotXY(x,rs(k,:),k*dt,0,L,-1.5,1.5);
%     % Uncomment the lines below to make some additional interesting plots.
%     % figure(2); semilogy(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([0 3 1e-8 1e-1])
%     % figure(3); loglog(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([3e-2 4 1e-8 1e-1])
% end

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
