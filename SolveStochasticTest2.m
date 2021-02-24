function SolveStochasticTest2(nSimu, T, nTimeSteps, maxIter, ...
	Intervals, x0, q0, maxU, solver, verbose)

% For printing figure: something like
% SolveStochastic(50, 10, 50, 20, [30,20,15], [0 5], 0.5, 1, 1, 2)

%% Initialisations
close all; clc;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
format('long'); 
tv = linspace(0,T, nTimeSteps + 1);
duration = tv(2) - tv(1);    % duration of one time step
a = 1/pi;   % sqrt constant
nX = length(x0);
% control(nX, size(tv, 2)*0.8 : size(tv, 2)) = 0;   % only positive u's for this test

%% Calculate average final state with discrete-time simulations 
% uStoreAll = []; qStoreAll = []; xStoreAll = []; realQStore = []; 

if verbose >= 1, fprintf('\nStart %i simulations\n', nSimu); end
mu = 0; sigma = sqrt(duration); % Wiener process parameters (Gaussian)
% normalDraw = normrnd(mu, sigma, 1, nSimu * nTimeSteps); % draw in advance

randomDraws = normrnd(mu, sigma, nSimu, nTimeSteps);
saveX = zeros(nSimu, nTimeSteps);

for iSimu = 1 : nSimu
% 	showProgress(iSimu, nSimu, verbose);
% 	
% 	if solver == 2, uMDPStore = []; end		
% 	uStore = [];
% 	xStore = [];
% 	qStore = [];
% 	simuScore = 0;
	
	
	%% Loop on time steps
	x = x0';
	for timeStep = 1 : nTimeSteps
		
		%% Draw future state
		x = x - sign(x) * maxU * duration + randomDraws(iSimu, timeStep);
		saveX(iSimu, timeStep) = x;
	end
	
end

% finalX = 
% histogram(finalX)
% figure; histogram(finalX, 20)
% toplot = mean(abs(saveX),1);
% plot(1 : nTimeSteps, toplot);waitforbuttonpress
% hold on;

ratio = 1;
% stepsToKeep = ;
xToKeep = saveX(1 : nSimu, floor(ratio * nTimeSteps) : nTimeSteps);
finalX = xToKeep(:);
h2 = histogram(abs(finalX));
[h,p] = lillietest(abs(finalX),'Distr','exp')

% h = histogram(finalX);
% myFit = fitdist(abs(finalX), 'kernel')
% index = linspace(0, max(abs(finalX)), 1000);
% plot(index, pdf(myFit, index))
% methods(myFit)
% waitforbuttonpress
% h.Values
hold on;
lambda = 2; % apparent expected value of abs(x)
tv = linspace(0, 5, 1000);
plot(tv, h2.Values(1) * exp(- lambda * tv));
figure; 
expDist = - log(rand(1,length(finalX))) / lambda;
qqplot(finalX, expDist)
set(gca,'DataAspectRatio',[1,1,1])
figure; 
% histogram(expDist);
plot(sort(abs(finalX))); hold on;
plot(sort(expDist)); legend('final x', 'exponential');
[h,p] = lillietest(expDist,'Distr','exp')
figure;
plot(sort(expDist), sort(abs(finalX)))
% histogram(sign(finalX).* sqrt(abs(finalX)))

if verbose >= 1
pm=char(177);
% bounds = 1.96 * std(abs(finalX)) / sqrt(length(finalX));
bounds = 1.96 * std(abs(finalX)) / sqrt(nSimu);
fprintf('\n Average final state x: %.6f %c %.6f \n', mean(abs(finalX)), pm, bounds);
end


end

%% sub-functions

function dndt = odeX(t,x,tv,u,q0,q1,q2,a)  % state diff. equations

u = interp1(tv,u,t);
q1 = interp1(tv,q1,t);
q2 = interp1(tv,q2,t);
dndt = a/x-u*(q0*abs(2*q1-1)+(1-q0)*abs(2*q2-1));

end	

function dndt = odeQ(t,q,tv,U2) % given u, compute q.

q1 = q(1); q2 = q(2); 
U2 = interp1(tv,U2,t);
dndt = [4 * q1 * (1-q1)^2 * U2;
    -4 * q2^2 * (1-q2) * U2];

end
