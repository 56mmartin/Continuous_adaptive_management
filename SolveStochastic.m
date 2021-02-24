function SolveStochastic(nSimu, T, nTimeSteps, maxIter, ...
	Intervals, x0, q0, maxU, solver, verbose)

% For printing figure: something like
% SolveStochastic(50, 10, 50, 20, [30,20,15], [0 5], 0.5, 1, 1, 2)

%% Initialisations
format('shortG')
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

tv = linspace(0,T, nTimeSteps + 1);
duration = tv(2) - tv(1);    % duration of one time step
nSimuDisplay = min(6, nSimu); %  #simulations to display at the end.
minU = -maxU;
UBounds = [minU-(maxU-minU)*0.1,maxU+(maxU-minU)*0.1];   % bounds for displaying control
% 	nStateX = 50; 
% 	nStateQ = 25;
% 	nAction = 50; 
	nXInterval = Intervals(1); 
	nQInterval = Intervals(2);
	nAInterval = Intervals(3); 
nX = length(x0);  % number of different x's to manage

if verbose >= 1
	close all; clc; profile clear; profile on;
	fprintf('\nT = %.1f, x0 = [', T); 
	fprintf(1, '%.1f; ', x0); 
	fprintf('], q0 = %.2f, maxU = %.1f\n', q0, maxU); 
	if solver >= 1
		fprintf('Solved with MDP\n'); 
		fprintf('# X intervals = %i, # Q intervals = %i, # A intervals = %i\n', ...
			nXInterval, nQInterval, nAInterval); 		
	else 
		fprintf('Solved with Dual control\n'); 
	end
end

%% Solve MDP on disretised problem and store policy
if solver >= 1     
	[MDPPolicy, discreteX, discreteQ, discreteA] = SolveMDP(...
		T, nTimeSteps, duration, x0, nXInterval, ...
		nQInterval, nAInterval, minU, maxU, verbose);	

end

%%	Find first control for optimal control approach
if solver == 0 || solver == 2 
	omega = 0.15;					% rate of exponential decay (the higher the faster decay).
	[~, initialControl] = SolveDeterministic(T, x0, q0, ...
		omega, 1e-5, 100, maxU, verbose-2);
end

%% Start simulations
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
	
	%% prepare online display of simulations
	if verbose >= 3   
		figure; 
		set(gcf,'Position',[7 0.5 8 10]);  % standard is [7 0 8 6]
		
		for iX = 1 : nX	   % states x and controls
			subplot(2*nX+1, 1,iX*2-1); grid on; 
			ylabel(sprintf('State (x%i)', iX)); 
			xlim([tv(1),tv(end)]); hold on; pplot(iX*3-2) = plot(0,0);
			
			subplot(2*nX+1, 1,iX*2); grid on; hold on; 
			ylabel(sprintf('Control (u%i)', iX)); 
			xlim([tv(1),tv(end)]); ylim(UBounds); hold on;
			pplot(iX*3-1) = plot(0,0);
			if solver == 2
				pplot(iX*3) = plot(0,0); 
				legend('Optimal Control', 'Dynamic Programming');
			end
		end
		
		% state q.
		subplot(2*nX+1, 1,2*nX+1);grid on; ylabel('Knowledge (y)'); 
			xlim([tv(1),tv(end)]); ylim([0,1]);	hold on; pplot(3*nX+1) = plot(0,0);
		if realQ == 1
			plot(tv, ones(1, size(tv, 2)), '-o');
		else
			plot(tv, zeros(1, size(tv, 2)), '-o');			
		end
		legend('q', 'real q');
	end
	
	
	%% Loop on time steps
	x = x0';
	q = q0;
	for timeStep = 1 : nTimeSteps
		ts = tv(timeStep);
		if solver == 0 || solver == 2 
			%% 1) Based on optimal control
			if ts == 0, control = initialControl; else
				[~, control] = SolveDeterministic(T - ts, x, q, ...
					omega, 1e-5, maxIter, maxU, 0, control);
			end
			u = control(:, 1);										% only use the first control. 
		end
		if solver == 1 || solver == 2
			%% 2) Based on MDPs - apply best action from the nearest discretised state
			xRounded = [];
			for iX = 1 : nX	
				[~, xR] = min((x(iX) - discreteX(iX, :)).^2);
				xRounded = [xRounded xR];				% rounded values of x
			end
			[~, qRounded] = min((q - discreteQ).^2);	% rounded value of q
			xRounded = bigbase2dec(xRounded-1,nXInterval);
			s = xRounded * nQInterval + qRounded;		% state number
			uTemp = (dec2bigbase(MDPPolicy(s, timeStep)-1,nAInterval, nX) + 1);	
			for iX = 1 : nX	
				uMDP(iX, 1) = discreteA(uTemp(iX));  % extract policy
			end
			
			if solver == 1, u = uMDP; end
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
% 	if iSimu <= nSimuDisplay
		uStoreAll = [uStoreAll; uStore]; xStoreAll = [xStoreAll; xStore]; qStoreAll = [qStoreAll; qStore];
% 	end
	
end

% if solver == 0
% 	xlswrite('StoredPolicies.xlsx', storedU);
% end

if verbose >= 1
pm=char(177);
bounds = 1.96 * std(score) / sqrt(nSimu);
fprintf('\n Average cost: %.1f %c %.1f \n', mean(score), pm, bounds);
end

if verbose >= 2
	
	% Plot average plots with errors	
	plotNumber = 0;	
	for iPlot = 1 : (2 * nX	+ 1)
		for q = [1 -1]
			subplot(2 * nX	+ 1, 1,iPlot); grid on; hold on; 
			style = '-';
			if q == 1, plotColor = [0 0 1]; 
			else, plotColor = [1 0 0]; end
					
			if iPlot == 2 * nX	+ 1   % state q
				data = qStoreAll(logical(realQStore == q), :);	
				bounds = [0 1]; ylim(bounds);
				if q == 1, style = '--'; end
				title(''); xlabel('Time (t)', 'FontWeight', 'bold');
				ylabel('(c) Knowledge y(t)', 'FontWeight', 'bold');
			else   % x or u, multiple states
				iX = floor((iPlot-1)/2) + 1;  % state number
				indices = logical(mod((0:nSimu*nX-1), nX) + 1 == iX);
				if mod(iPlot, 2) == 1   % state x
					if q == -1, continue; end
					data = xStoreAll(indices, :);
					plotColor = [0 0 0];	
					ylabel(sprintf('State (x%i)', iX), 'FontWeight', 'bold');  
				else  % control u
					ylim(UBounds);	
					if q == 1, style = '--'; end   % dashed blue line for paper
					ylabel(sprintf('Control (u%i)', iX), 'FontWeight', 'bold');
% 					indices
% 					logical(realQStore == 1)
					data = uStoreAll(indices, :);
% 					data
% 					size()
					data = data(logical(realQStore == q), :);
				end
			end
			if isempty(data), continue; end

			medianData = prctile(data, 50, 1);
			yLow = prctile(data, 5, 1);
			yHigh = prctile(data, 95, 1);		
			f=patch([tv fliplr(tv)],[yHigh fliplr(yLow)],plotColor);% shaded area
			set(f,'EdgeColor','none');  % no edge
			alpha(0.12);		% transparency
			plotNumber = plotNumber + 1;
			h(plotNumber) = plot(tv, medianData, style, 'color', plotColor, 'LineWidth', 2);  % mean

			for iSimu = 1 : nSimuDisplay    % display simulations	
% 				if mod(iPlot, 2) == 1 && iSimu > nSimuDisplay / 2  % only half display for u and q.
% 					break;
% 				end
				try %#ok<ALIGN>
						plot(tv, data(iSimu, :), style, 'color', (1+plotColor)/2);
				catch, continue; end
	% 			end
			end
		end
	end

	legend(h([end-1:end]),{'When real q = 1','When real q = 0'},'Location','Best');
	

	for iSimu = 1 : 0 * nSimu    % display simulations	
		figure; 
		subplot(3, 1,1);
		plot(tv, xStoreAll(iSimu, :));
		grid on;
		title('State (x)');
		subplot(3, 1,2);
		plot(tv, uStoreAll(iSimu, :));
		grid on;
		title('Control (u)');
		ylim(UBounds);
		subplot(3, 1,3);
		plot(tv, qStoreAll(iSimu, :)); hold on;
		grid on;
		if realQStore(iSimu) == 1
			plot(tv, ones(1, size(tv, 2)), '-o');
		else
			plot(tv, zeros(1, size(tv, 2)), '-o');			
		end
		legend('q', 'real q');
		ylim([0,1]);
		title('Knowledge (q)');
		if waitforbuttonpress == 1, break, end
	end
end

% % Save figure as pdf without margins
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MySavedFile','-dpdf');

% Change labels
fig = gcf;
subplot(3, 1,1);
xlabel(' ');
ylabel('State {\itx(t)}');
subplot(3, 1,2);
ylabel('Control {\itu(t)}');
xlabel(' ');
subplot(3, 1,3);
ylabel('Knowledge {\ity(t)}');
xlabel('Time{\it t}');

% Change labels
fig = gcf;
subplot(5, 1,1);
set(gca,'FontSize',10);
xlabel(' ');
ylabel('State {\itx_1(t)}');
subplot(5, 1,2);
set(gca,'FontSize',10);
ylabel('Control {\itu_1(t)}');
subplot(5, 1,3);
set(gca,'FontSize',10);
xlabel(' ');
ylabel('State {\itx_2(t)}');
subplot(5, 1,4);
set(gca,'FontSize',10);
ylabel('Control {\itu_2(t)}');
xlabel(' ');
subplot(5, 1,5);
set(gca,'FontSize',10);
ylabel('Knowledge {\ity(t)}');
xlabel('Time{\it t}');


if verbose >= 1
	profile viewer;
end

end
