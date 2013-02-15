function [filt_bg,flag] = f_waveletBackground(I)

% decompose image into wavelet basis
% N = 4 represents amount of detail in background
% db3 is type of basis, Daubechies "3"
[C,S] = wavedec2(I,4,'db3');
C_filt = zeros(1,length(C));
k = S(1,1)*S(1,2);      % pull out approx coeff
%n = k+3*S(2,1)*S(2,2);
C_filt(1:k) = C(1:k);
%x = 1:S(1);

% search for sharp peaks in approx coefficents, average out those with
% neighboring coefficients
% corresponds to bright fiducials
med = median(C_filt(1:k));
flag_vec = zeros(k,1);

for ii=2:k-2;
    if abs(C_filt(ii)-med)<1000;
        C_filt(ii) = C_filt(ii);
    elseif abs(C_filt(ii+1)-med)<1000;
        C_filt(ii) = mean([C_filt(ii-1),C_filt(ii+1)]);
        flag_vec(ii) = -1;
    else
        C_filt(ii) = mean([C_filt(ii-1),C_filt(ii+2)]);
        flag_vec(ii) = -1;
    end
end

%report how many times the above loop went into alternate if
%statements
flag = sum(flag_vec);


I_filt = waverec2(C_filt,S,'db3');
% smooth out wavelet reconstruction with Gaussian blur
h = fspecial('gaussian',100,15);
filt_bg = imfilter(I_filt,h,'replicate');


%             imagesc(filt_bg,[0 1000])
%             axis image
%             figure
%             imagesc(I,[0 500])
%             axis image
%             figure
%             imagesc((I-filt_bg),[0 300])
%             axis image

end