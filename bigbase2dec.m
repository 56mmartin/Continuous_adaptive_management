
function d = bigbase2dec(v,b)
%BIGBASE2DEC Convert base B vector to decimal integer.
% BIGBASE2DEC(V,B) converts the vector V of base B digits into its
% decimal (base 10) equivalent. B must be an integer greater than 1.
%
% If V is a matrix, each row is interpreted as a base B vector.
%
% Example
% base2dec([2 1 2],3) returns 23
%
% See also DEC2BIGBASE, BASE2DEC, DEC2BASE, HEX2DEC, BIN2DEC.

% written by Douglas M. Schwarz
% Eastman Kodak Company (on leave until 4 Jan 1999)
% schwarz@kodak.com, schwarz@servtech.com
% 2 October 1998

if isempty(v), d = []; return, end

[m,n] = size(v);
b = [ones(m,1),cumprod(b(ones(m,n-1)),2)];
d = sum((b .* fliplr(v)),2);