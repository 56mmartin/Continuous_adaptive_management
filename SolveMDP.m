function [MDPPolicy, discreteX, discreteQ, discreteA] = SolveMDP(...
	T, nTimeSteps, duration, x0, nXInterval, ...
		nQInterval, nAInterval, minU, maxU, verbose)

%% Solve MDP policy.

	if verbose >= 2, fprintf('\nBuild MDP...\n'); end

	
	%% Build states and actions
	% 3 standard deviations contains 99.7% of the values.
	maxX = x0 + 2 * sqrt(T) + maxU * T; 
	minX = x0 - 2 * sqrt(T) + minU * T;
	nX = length(x0);  % number of different x's to manage (~= nStates !!)
	for iX = 1 : nX
		discreteX(iX, :) = linspace(minX(iX), maxX(iX), nXInterval);		% discretise X
		discreteQ = linspace(0, 1, nQInterval);					% discretise Q
	end
	discreteA = linspace(minU, maxU, nAInterval);		% discretise actions		

	
	nXStates = nXInterval ^ nX;		% # of states x
	nStates = nXStates * nQInterval;		% # of states q and x
	nAction = nAInterval ^ nX;		% # of actions 
	
	
	%% Initialise transitions & reward matrices
	P = cell(1, nAction);
	R = zeros(nStates, nAction);
	nSample = 10^2;
	mu = 0; sigma = sqrt(duration); % Wiener process parameters (Gaussian)
	sample = normrnd(mu,sigma,1,nSample);	% hist(sample, 20);
	
	% First, fill matrices P and R
	for iAction = 1 : nAction				% action number
		row = []; col = []; val = []; 
% 		as = dec2base(iAction-1,nAInterval, nX)-'0' + 1;
		as = dec2bigbase(iAction-1,nAInterval, nX) + 1;
		for iState = 1 : nXStates			% sub-state X number
			xs = dec2bigbase(iState-1,nXInterval, nX) + 1;
			for iQ = 1 : nQInterval				% sub-state Q number
				q = discreteQ(iQ);				% sub-state Q value
				
				binX = 1;		
				x = zeros(1, nX);
				a = zeros(1, nX);
				for iX = 1 : nX
					x(iX) = discreteX(iX, xs(iX));				% sub-state X value		
					a(iX) = discreteA(as(iX));				% action value
				end
				xPosAverage = x + a * duration + mu * duration;  % average future x if realQ = 1
				xNegAverage = x - a * duration + mu * duration;  % average future x if realQ = -1

				for iX = 1 : nX
					% Approximate future x's by their nearest discretised x.
					discreteXMatrix = discreteX(iX, :)' * ones(1, nSample);
					 % Bin in discretised values if realQ = 1, then  if realQ = -1
					xMatrix = ones(nXInterval, 1) * (xPosAverage(iX) + sample);
					[~, xPosProjected] = min((xMatrix - discreteXMatrix).^2);
					xMatrix = ones(nXInterval, 1) * (xNegAverage(iX) + sample);					
					[~, xNegProjected] = min((xMatrix - discreteXMatrix).^2);

					% Count each bin
					binXPos = histcounts(xPosProjected, 0.5:(nXInterval+0.5));		
					binXNeg = histcounts(xNegProjected, 0.5:(nXInterval+0.5));	
					binOneX = (q * binXPos + (1-q) * binXNeg) / nSample;
					assert(abs(sum(binX)-1)< 1);
% 						binX
					binX = binOneX' * binX';
					binX = binX(:);
					if size(binX, 1) == 1, binX = binX'; end
				end

				displ = zeros(length(binX), nX + 3);
				for iFuture = 1 : length(binX)
					if binX(iFuture) ~= 0
						futXStates = dec2bigbase(iFuture-1,nXInterval, nX) + 1;
						for iX = 1 : nX
							futureX(iX) = discreteX(iX, futXStates(iX));				% sub-state X value			
						end
						if q == 1
% 							futureQ = ones(size(discreteX)); % 1 whatever future x is.
							futureQ = 1; % 1 whatever future x is.
						else
							exponent = ((futureX - xPosAverage) * (futureX - xPosAverage)' - ...
												 (futureX - xNegAverage) * (futureX - xNegAverage)') / ...
												 (2 * sigma ^2);
							futureQ = 1 ./ (1 + (1-q) / q * exp(exponent));							
						end
								
						% Approximate future q's by their nearest discretised q.
% 						discreteQMatrix = discreteQ' * ones(1, length(futureQ));
% 						qMatrix = ones(length(discreteQ), 1) * futureQ;
% 						[~, qRounded] = min((qMatrix - discreteQMatrix).^2);
						[~, qRounded] = min((futureQ - discreteQ).^2);
						futureState = (iFuture-1) * nQInterval + qRounded;		% future state number	
						currentState = (iState-1) * nQInterval + iQ;    % current state number
						
						displ(iFuture, :) = [futXStates qRounded futureState currentState];
						row = [row currentState]; col = [col futureState]; val = [val binX(iFuture)]; 
					end
				end
				
	% 		Fill R	
				R(currentState, iAction) = duration * (a * a' + x * x');
			end
		end
% 		Fill P	
		P{iAction} = sparse(row, col, val);
	end
	
	if verbose >= 2, fprintf('\nSolve MDP\n'); end
	% Then, solve MDP with backwards induction
	V = zeros(nStates, 1);
	
					
	MDPPolicy = zeros(nStates, nTimeSteps);
	MDPValue = zeros(nStates, nTimeSteps);
	Q = zeros(nStates, nAction);
	for t = nTimeSteps : - 1 : 1
		Vprev = V;
		
		for iAction = 1 : nAction				% action number			
			Q(:,iAction) = R(:,iAction) + P{iAction}*Vprev;  % Bellman's equation
		end
    [V, policy] = min(Q,[],2);								% Select optimal action
		MDPPolicy(:, t) = policy;			% Store optimal action
		MDPValue(:, t) = V;												% Store optimal value action
		
	end

	if verbose >= 2
		fprintf('\nMDP solved.\n');
	end

% [MDPPolicy(:, 1) MDPPolicy(:, end) ]


end