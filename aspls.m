function [detrended_signal, baseline] = aspls(X, lambda, order, k, itermax, epsilon)
%     Performs Adaptive smoothness penalized least squares smoothing
%     (asPLS) which is a type of detrending which tries to preserve the
%     features of interest in the data (i.e. peaks).
% 
%     Inputs:
%       - X : array-like, The y-values of the measured data, with t data points. 
%           Must not contain missing data (NaN) or Inf.
%       - lambda : float, The smoothing parameter. Larger values will create smoother baselines.
%           Default is 1e5.
%       - order : int, The order of the differential matrix. Must be greater than 0. Default is 2
%           (second order differential matrix). Typical values are 2 or 1.
%       - k : int, coefficient of hardness for the decay of the weights as
%           a function of distance to baseline. A higher k means a harder
%           decay which will push the baseline towards lower value of the
%           input data. 0.5 is empirically recommended (Default).
%       - itermax : int, The max number of fit iterations. Default is 20.
%       - epsilon : float, The exit criteria. Default is 1e-4.
% 
%     Returns:
%       - detrended_signal : (array) The input signal - the calculated baseline.
%       - baseline : (array) The calculated baseline.
% 
%     References:
%       Zhang, F., et al. Baseline correction for infrared spectra using
%       adaptive smoothness parameter penalized least squares method.
%       Spectroscopy Letters, 2020, 53(3), 222-233.
%       
%     Author: Corentin Nelias, Agarwal Lab 2023

if nargin < 6
    epsilon = 10e-4;
    if nargin < 5
        itermax=20;
        if nargin < 4
            k = 0.5;
            if nargin < 3
                order=2;
                if nargin < 2
                    lambda=10e5;
                    if nargin < 1
                        error('Missing input data.');
                    end
                end
            end
        end
    end
end

[features, t] =size(X);
D = diff(speye(t), order);
DD = lambda*D'*D;
for i=1:features %iterating over features
    w = ones(t, 1);
    a = ones(t, 1);
    x=X(i,:);
    for j=1:itermax
        W=spdiags(w, 0, t, t);
        A=spdiags(a, 0, t, t);
%         C = chol(W + A*DD);
%         z = (C\(C'\(w .* x')))';
        z = (W + lambda*A*D'*D)\(W*X');
        d = x-z;
        a = abs(d')/max(abs(d));
        sigma_ = std(d(d < 0)); %standard deviation of events below fit.
        w_next = ones(1, t)./( 1 + exp(k*(d - sigma_)/sigma_) )  ;
        if sum(abs(w - w_next))/sum(abs(w)) < epsilon %testing if weights have converged
            break;
        end
        w = w_next'; %if weights did not converge, update.
    end
    Z(i,:)=z;
end
detrended_signal=X-Z;
baseline = Z;









