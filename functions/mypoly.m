function S=POLY(n, NumofPixels)
%function POLY(n, NumofPixels)
%
%This function creates polynomial spectra of orders 0 through n.  There
% are NumofPixels pixels in each spectrum, and the values range from 1 to
% 2^n.
%
%The output is a spectral matrix with each spectrum occupying one ROW.

%get the step size correct
Step = 1/(NumofPixels-1);

S=[];

%loop to create the spectra
for i=1:Step:2,
	NewColumn =[];
	%loop to create new column
	for j=0:1:n,
		NewColumn=[NewColumn; i^j];
		end; %of j loop
	S=[S, NewColumn];
   	end; %of i loop

		

