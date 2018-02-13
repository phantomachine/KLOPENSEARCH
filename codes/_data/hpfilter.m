function [cycle,trend]=hpfilter(rawdata,LAMBDA)

% HPFILTER.M
% ------------------------------------------------------------
% Implements fast solution to the minimization problem in 
% Hodrick and Prescott's filter.
%
% ------------------------------------------------------------
% (c)-left T.Kam, 2004-
% ------------------------------------------------------------

if size(rawdata,1) < size(rawdata,2)
   rawdata = rawdata';
end

t=size(rawdata,1);
a = 6*LAMBDA + 1;
b = -4*LAMBDA;
c = LAMBDA;
d = [c,b,a];
d = ones(t,1)*d;
m = diag(d(:,3)) + diag(d(1:t-1,2),1) + diag(d(1:t-1,2),-1);
m = m + diag(d(1:t-2,1),2) + diag(d(1:t-2,1),-2);

m(1,1) = 1+LAMBDA;       m(1,2) = -2*LAMBDA;
m(2,1) = -2*LAMBDA;      m(2,2) = 5*LAMBDA+1;
m(t-1,t-1) = 5*LAMBDA+1; m(t-1,t) = -2*LAMBDA;
m(t,t-1) = -2*LAMBDA;    m(t,t) = 1+LAMBDA;

trend = m\rawdata;

cycle = rawdata - trend;


