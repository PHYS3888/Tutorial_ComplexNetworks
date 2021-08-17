function [expFittedData,A,n] = fitExponential(binCenters,fData)

% Do the fitting:
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.15,5]);
f = fittype('A*exp(-n*x)','options',s);
[c,stats] = fit(binCenters,fData,f);

% The two fitted parameters are c.A and c.n
% We can extract them as:
A = c.A; % Amplitude, A
n = c.n; % Decay rate, n (/mm)

% So we can define a function for the fitted exponential:
f_handle = @(x) A*exp(-n*x);

% Then apply the function to the input bin positions:
expFittedData = f_handle(binCenters);

end
