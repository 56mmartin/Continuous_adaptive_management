function SolveStochasticTest(nSimu, T, nTimeSteps, maxIter, ...
	Intervals, x0, q0, maxU, solver, verbose)

% For printing figure: something like
% SolveStochastic(50, 10, 50, 20, [30,20,15], [0 5], 0.5, 1, 1, 2)

%% Initialisations
format('shortG'); close all; clc;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

tv = linspace(0,T, nTimeSteps + 1);
duration = tv(2) - tv(1);    % duration of one time step
a = 1/pi;   % sqrt constant
nX = length(x0);
control = maxU * ones(nX, size(tv, 2));   % only positive u's for this test
% control(nX, size(tv, 2)*0.8 : size(tv, 2)) = 0;   % only positive u's for this test

% x = 0.098;
% step = 0.0001;
% t = [step:step:1];
% y = sqrt(2*t/pi).*exp(-(x-t).^2./(2*t)) + abs(x-t).*erf(abs(x-t)./sqrt(2*t));
% y1 = sqrt(2*t/pi).*exp(-(x-t).^2./(2*t));
% y2 = abs(x-t).*erf(abs(x-t)./sqrt(2*t));
% plot(t, y); hold on
% plot(t, y1);
% plot(t, y2);
% legend('y', 'y1', 'y2')
% waitforbuttonpress
%% Calculate average final state with deterministic model
U2 = sum(control.^2, 1);	
[tq,yq] = ode23s(@(t,y) odeQ(t,y,tv,U2),[0,T],[q0,q0]);  % given u, compute q.
q1 = interp1(tq,yq(:,1),tv);
q2 = interp1(tq,yq(:,2),tv);

x = zeros(nX, size(tv, 2));
% for iX = 1 : nX		% given u, compute x.
% 	[tx,yn] = ode23s(@(t,y) odeX(t,y,tv,control(iX, :),q0,q1,q2,a),...
% 		[0,T],x0(iX));  
% 	x(iX, :) = interp1(tx,yn,tv);
% 	x(iX, :)
% end
	
plot(tv, x)

%% Calculate average final state with discrete-time simulations 
uStoreAll = []; qStoreAll = []; xStoreAll = []; realQStore = []; 
score = zeros(1, nSimu);

if verbose >= 1, fprintf('\nStart %i simulations\n', nSimu); end
mu = 0; sigma = sqrt(duration); % Wiener process parameters (Gaussian)
% normalDraw = normrnd(mu, sigma, 1, nSimu * nTimeSteps); % draw in advance

for iSimu = 1 : nSimu
	showProgress(iSimu, nSimu, verbose);
	
	if solver == 2, uMDPStore = []; end		
	uStore = [];
	xStore = [];
	qStore = [];
	simuScore = 0;
	
	% draw model - q is the proba that realQ = 1, i.e. the control is "positive".
	if rand < q0, realQ = 1; else realQ = -1; end  %#ok<*SEPEX>
	realQStore = [realQStore realQ];
	
	
	%% Loop on time steps
	x = x0';
	q = q0;
	for timeStep = 1 : nTimeSteps
		ts = tv(timeStep);
		if sign(x)*sign(2*q-1) == 0
			u = abs(control(1,timeStep));
		else
			u = -sign(x)*sign(2*q-1)*abs(control(1,timeStep));
		end
		
		%% Draw future state
		if timeStep == 1  % start storing data
			uStore = [uStore u]; xStore = [xStore x]; qStore = [qStore q];
			if solver == 2, uMDPStore = [uMDPStore uMDP]; end
		end
		xPos = x + u * duration + mu * duration;        % average future value if q = 1
		xNeg = x - u * duration + mu * duration;		    % average future value if q = 0
		
		x = x + realQ * u * duration + normrnd(mu, sigma, nX, 1);
% 		x = x + realQ * u * duration + normalDraw((iSimu - 1) * nTimeSteps + timeStep);  % real future x
		
		posteriorPos = q  *  exp(-(x-xPos)'*(x-xPos)/2/sigma^2) / (sqrt(2*pi) * sigma);
		posteriorNeg = (1-q)*exp(-(x-xNeg)'*(x-xNeg)/2/sigma^2) / (sqrt(2*pi) * sigma);
		% The prob of positive scenario being true is:
		if posteriorPos ~= 0 || posteriorNeg ~= 0   % else q does not change
			q = posteriorPos / (posteriorPos + posteriorNeg);
		end
		
		simuScore = simuScore + duration * (u'*u + x'*x);		
			
		uStore = [uStore u]; xStore = [xStore x]; qStore = [qStore q];
		if solver == 2, uMDPStore = [uMDPStore uMDP]; end
		if verbose >= 3   % online display of simulations
			xAxis = tv(1:length(xStore));
			for iX = 1 : nX	
				subplot(2*nX+1, 1,iX*2-1); 
				set(pplot(iX*3-2),'XData',xAxis, 'YData',xStore(iX, :));
				subplot(2*nX+1, 1,iX*2); 
				set(pplot(iX*3-1),'XData',xAxis, 'YData',uStore(iX, :));
				if solver == 2
					set(pplot(iX*3),'XData',xAxis, 'YData',uMDPStore(iX, :)); end
			end
			subplot(2*nX+1, 1,2*nX+1); 
			set(pplot(3*nX+1),'XData',xAxis, 'YData',qStore);			
			if solver == 1, pause(0.1); else drawnow; end % slow display down 
		end
		
		
	end
	score(iSimu) = simuScore;
	finalX(iSimu) = abs(x);
	unsignedX(iSimu) = x;
% 	if iSimu <= nSimuDisplay
		uStoreAll = [uStoreAll; uStore]; xStoreAll = [xStoreAll; xStore]; qStoreAll = [qStoreAll; qStore];
% 	end
% 	plot(xStore);
% 	hold on; plot(uStore);
% 	waitforbuttonpress
end

% if solver == 0
% 	xlswrite('StoredPolicies.xlsx', storedU);
% end
hist(unsignedX)
% toplot = mean(abs(xStoreAll),1)
% hold on;
% plot(tv, toplot);
if verbose >= 1
pm=char(177);
bounds = 1.96 * std(finalX) / sqrt(nSimu);
fprintf('\n Average final state x: %.2f %c %.2f \n', mean(finalX), pm, bounds);
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
