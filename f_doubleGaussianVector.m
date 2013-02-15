% Returns the sum of squared error between Zdata and a double Gaussian
% function parameterized by par. bkgnd gives the constant offset of the
% double Gaussian, while ii and jj are matrices over the domain

% Version history
% 2010-05-17    mlew        Revised fitting parameters
% 2009-07-06    mlew        Changed 8 parameter Gaussian fit to full 11
%                           parameter fit (with different sigma for x and y
%                           and 1st and 2nd lobes
% 2009-06-30    Majid       Replaced for loops with matrix commands
%                           to speed up
% 2009-06-29    Majid       Got it from Mike



function [err, jacobian] = f_doubleGaussianVector(par,Zdata,bkgnd,ii,jj)
  Zmodel = par(1).*exp( -((ii-par(3)).^2+(jj-par(4)).^2) / (2*par(7).^2)) ...
          +par(2).*exp( -((ii-par(5)).^2+(jj-par(6)).^2) / (2*par(8).^2)) ...
          +bkgnd;
 
  err = reshape(Zmodel-Zdata,1,[]);
  
  if nargout > 1
      jacobian = f_doubleGaussianVectorJacobian(par,ii,jj);
  end
