function test


close all

% syms z SIGMA  x
% check it sums to one.
deltat = 0.005;
UB = 15;
% SIGMA = 0.555;
res = [];
% SIGMAArray = [0.1 : 0.04 : 1];
SIGMAArray = 1;
for SIGMA = SIGMAArray
	
% f = @(x,z) x.^2 .* (  2 * exp(-2*z-(x-z+deltat).^2/2/deltat) + ...
% 											2 * exp(-2*z-(x+z-deltat).^2/2/deltat)  ) / ...
% 													sqrt(2*pi*deltat);
	
f = @(x,z) abs(x) .* ( exp(-(x-z+deltat).^2/2/deltat) + ...
						exp(-(x+z-deltat).^2/2/deltat)  ) ...
						/ sqrt(2*pi*deltat) .* exp(-2*z);
					
% fplot(int(f, z))
% res = [res sqrt(integral2(f,-UB,UB,0,UB))];
res = [res integral2(f,-UB,UB,0,UB)];
end
res
if length(SIGMAArray) > 1.5
	plot(SIGMAArray, res);
	set(gca,'DataAspectRatio',[1,1,1])
	hold on
	plot(SIGMAArray, SIGMAArray);
% 	plot(SIGMAArray, res - SIGMAArray);
end





	
% f = @(x,z) x.^2 .* (  exp(-z.^2./(2*SIGMA^2)-(x-z+deltat).^2/2/deltat) + ...
% 											exp(-z.^2./(2*SIGMA^2)-(x+z-deltat).^2/2/deltat)  ) / ...
% 													(2*pi*SIGMA*sqrt(deltat));
% syms x
% f = x^5;
% fplot(int(f))
% 
% fun = @(x) exp(-x.^2).*log(x).^2;
%% Old test
% nSample = 1e8
% T = 10;
% 	mu = 0; sigma = sqrt(T); % Wiener process parameters (Gaussian)
% 	sample = normrnd(mu,sigma,1,nSample);	% hist(sample, 20);
% 	meanS = mean(abs(sample))
	
% 	[t,y] = ode45(@(t,y) 1/(3.141592*y),[0,T],1e-5,odeset('AbsTol',1e-6));  % given u, compute q.
% 	plot(t, y)
% 	y(end)
	
	