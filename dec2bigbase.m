
function s = dec2bigbase(d,base,n)
%DEC2BIGBASE Convert decimal integer to base B vector.
% DEC2BIGBASE(D,B) returns the representation of D as a vector of
% digits in base B. D must be a non-negative integer smaller than 2^52
% and B must be an integer greater than 1.
%
% DEC2BIGBASE(D,B,N) produces a representation with at least N digits.
%
% Examples
% dec2bigbase(23,3) returns [2 1 2]
% dec2bigbase(23,3,5) returns [0 0 2 1 2]
%
% See also DEC2BASE, BASE2DEC, DEC2HEX, DEC2BIN.

% written by Douglas M. Schwarz
% Eastman Kodak Company (on leave until 4 Jan 1999)
% schwarz@kodak.com, schwarz@servtech.com
% 1 October 1998

if size(d,2) ~= 1, d = d(:); end

base = floor(base);
if base < 2, error('B must be greater than 1.'); end
if base == 2
  [~,nreq] = log2(max(d));
else
  nreq = ceil(log2(max(d) + 1)/log2(base));
end

if nargin == 3
    nreq = max(nreq,1);
    n = max(n,nreq);
    last = n - nreq + 1;
else
    n = max(nreq,1);
    last = 1;
end

s(:,n) = rem(d,base);
while n ~= last
    n = n - 1;
    d = floor(d/base);
    s(:,n) = rem(d,base);
end