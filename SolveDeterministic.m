function [timeStep, bestControl] = SolveDeterministic(...
	T, signed_x0, q0, omega, threshold, maxIter, maxU, verbose, firstGuess)

% Solve deterministic optimal control problem to guide decisions on the real stochastic problem.

if verbose >= 1, fprintf('\nFind initial control...\n'); end
% Check problem with SolveStochastic(200, 10, 50, 20, [0 0 0 0], 0.2, 5, 0, 3)

%% Initialisation
nDiscrete = 100;
tv = linspace(0,T, nDiscrete);
duration = tv(2) - tv(1);    % duration of one time step
maxBestScoreTime = 10;   % time after which to stop if no improvement
bestScore = 1e8;
a = 1/pi;   % sqrt constant

signed_x0(signed_x0 == 0) = 1e-3; % zero not allowed for numerical reasons.
x0 = abs(signed_x0); % only positive x's in the deterministic version.
nX = length(x0);
% LOWER OMEGA AND constant axis for states

if exist('firstGuess') %#ok<EXIST>
	u = abs(firstGuess);
else % "aggressive" first guess
	u = maxU * ones(nX, size(tv, 2));
% 	u(1:floor(nDiscrete/2)) = maxU;
end
timeStep = 0;
if verbose >= 2
	f = figure; movegui('west');
	g = figure;
	h = figure; movegui('east');
end
Values = [];
timeSinceBestScore = 0;
while true   % Loop on policies
	%% Compute states x and q over time given the current policy u and display
	U2 = sum(u.^2, 1);	
	[tq,yq] = ode23s(@(t,y) odeQ(t,y,tv,U2),[0,T],[q0,q0]);  % given u, compute q.
	q1 = interp1(tq,yq(:,1),tv);
	q2 = interp1(tq,yq(:,2),tv);
	
	x = zeros(nX, size(tv, 2));
	for iX = 1 : nX		% given u, compute x.
		[tx,yn] = ode23s(@(t,y) odeX(t,y,tv,u(iX, :),q0,q1,q2,a),...
			[0,T],x0(iX));  
		x(iX, :) = interp1(tx,yn,tv);
		x(iX, :)
	end
	
	if verbose >= 2
		figure(f);
		for iX = 1 : nX
			plot(tv, x(iX, :),'DisplayName',['x' num2str(iX)]); hold on; 
		end
		title('States');
		plot(tv, q1,'DisplayName','q1');plot(tv, q2,'DisplayName','q2');
		drawnow; legend('show'); hold off
	end
	
	value = (sum(sum(u.^2)) + sum(sum(x.^2))) * duration * 100 / 101;
	if isnan(value)
		q1
		q2
	end
	
	%% Stop algorithm if solution not improved for 'maxBestScoreTime' time steps
	
	if value < bestScore * 0.95
		timeSinceBestScore = 0;		
	else
		if timeSinceBestScore > maxBestScoreTime
			timeStep = maxIter;
			if verbose >= 1, fprintf('Value has converged --> exit\n'); end
			break;
		else
			timeSinceBestScore = timeSinceBestScore + 1;
		end
	end
	if timeStep > 0
		valueDiff = value - Values(end);
	else 
		valueDiff = 0;
	end
	Values = [Values value];
	
	
	%% Remember best score & control
	if value < bestScore
		bestScore = value;
		bestControl = u;
	end	
	
	%% Compute costates (lambdas) given the states and policy
	% given q and u, compute lambdas
	lx = zeros(nX, size(tv, 2));
	for iX = 1 : nX		% compute lx.
		[tl,yl] = ode23s(@(t,y) odeLX(t,y,tv,x(iX, :),a),[T,0],0);	
		lx(iX, :) = interp1(tl,yl,tv);
	end
	
	LU = sum(lx.*u, 1);
	[tl,yl] = ode23s(@(t,y) odeLQ(t,y,tv,q0,q1,q2,LU,U2),[T,0],[0,0]);  % compute lq.
	lq1 = interp1(tl,yl(:,1),tv);
	lq2 = interp1(tl,yl(:,2),tv);
	
	if verbose >= 2   % display
		figure(g);  
		for iX = 1 : nX
			plot(tv, lx(iX, :),'DisplayName',['lx' num2str(iX)]); hold on;
		end
		title('Lambdas');
		plot(tv, lq1,'DisplayName','lq1');plot(tv, lq2,'DisplayName','lq2');
		drawnow; legend('show'); hold off
	end
	
%%  Find policy

	uold = u;	

% %% Update policy: depends whether second derivative >0 or <0 separately	
	secondDer = 8*lq1.*q1.*(q1-1).^2 + 8*lq2.*(q2-1).*q2.^2 + 2;
	
% 1) strictly convex Hamiltonian (second derivative >0)
	convexIndices = logical(secondDer > 0);
	
	numerator = lx.*(ones(nX, 1) * (q0*abs(1-2*q1)-(q0-1)*abs(1-2*q2)));
	denominator = ones(nX, 1)*secondDer;
	unew(:, convexIndices) = numerator(:, convexIndices)./denominator(:, convexIndices);
	
	
% 2) else: 'bang-bang' policy, compare Hamiltonian for small and big U's
	unew(:, ~convexIndices) = maxU;
	
% 	Bound policy: truncation...
	unew(unew>maxU) = maxU;
	unew(unew<0) = 0;
	
%  	Update u linearly
	coeff = exp(-omega*(timeStep+1));
	u = uold + (unew - uold) * coeff;


%% Display policy
	if verbose >= 2
		for iX = 1 : nX
			figure(h); 
			plot(tv, unew(iX, :), 'DisplayName','unew');
			hold on; title('Policies');
			plot(tv, u(iX, :), '-o', 'DisplayName','u');
			plot(tv, uold(iX, :), 'DisplayName','uold');
			ylim([0 maxU]);
			drawnow; legend('show'); hold off;
		end
	end
	
	
	%% Print current solution KPIs and stop if good converged
	timeStep = timeStep + 1;
	if verbose >= 1
		policyDiff = sum(sum((u - uold).^2));
		if timeStep == 1
		disp('          step | total cost | cost change | policy change');
		end
		disp([timeStep value valueDiff policyDiff coeff]);
	end

	if sum(sum((u - uold).^2)) < threshold || timeStep > maxIter
			break
	end

end

%% Display states and policy over time
if verbose >= 2
	figure; movegui('south'); 
	clf
	subplot(1,4,1)
	plot(tv,x)
	title('x')
	subplot(1,4,2)
	plot(tv,q1)
	title('q1')
	subplot(1,4,3)
	plot(tv,q2)
	title('q2')
	subplot(1,4,4)
	plot(tv,u)
	title('u')

	figure; movegui('southeast');  plot(Values);
end

for iX = 1 : nX	
	if (q0 - 0.5) * signed_x0(iX) > 0
		bestControl(iX, :) = -bestControl(iX, :);  % adjust sign of control depending on sign of (q-0.5)*x   
	end
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
	
function dldt = odeLX(t,lx,tv,x,a) % costate diff. equations

x = interp1(tv,x,t);
dldt = a * lx / x^2 - 2 * x;
end
	
function dldt = odeLQ(t,lq,tv,q0,q1,q2,LU,U2) % costate diff. equations

lq1 = lq(1); lq2 = lq(2);

U2 = interp1(tv,U2,t);
LU = interp1(tv,LU,t);
q1 = interp1(tv,q1,t);
q2 = interp1(tv,q2,t);
eps = 1e-6;
dldt = [- 2 * LU * q0 * (abs(1- 2*q1+eps) - abs(1- 2*q1))/eps ...
						- 4 * lq1 * U2 * (3 * q1^2 - 4 * q1 + 1);
				2 * LU * (q0 - 1) * (abs(1- 2*q2+eps) - abs(1- 2*q2))/eps ...
						+ 4 * lq2 * U2 * q2 * (2-3*q2)];
end
