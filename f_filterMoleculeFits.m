% Copyright (c)2013, The Board of Trustees of The Leland Stanford Junior
% University. All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
% Neither the name of the Leland Stanford Junior University nor the names 
% of its contributors may be used to endorse or promote products derived 
% from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [lobeDistBounds] = ...
    f_filterMoleculeFits(fitFilePrefix,lobeDistBounds,conversionGain,EMGain)

scrsz = get(0,'ScreenSize');

conversionFactor = conversionGain/EMGain;
ampRatioLimit = 0.5;
% sigmaRatioLimit = 0.4;

%% load localization data
numFitsInFiles = zeros(1,length(fitFilePrefix));
for fileNum=1:length(fitFilePrefix)
    load([fitFilePrefix{fileNum} 'molecule fits.mat'],'totalPSFfits');
    numFitsInFiles(fileNum) = size(totalPSFfits,1);

    if fileNum == 1
        tempPSFfits = totalPSFfits;
    else
        tempPSFfits = [tempPSFfits; totalPSFfits];
    end
end
totalPSFfits = tempPSFfits;
clear tempPSFfits;
    
%% plot various fit distributions to determine filters
figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w');
subplot(2,2,1)
detectedPhotons = totalPSFfits(:,21);
detectedPhotons = detectedPhotons(detectedPhotons >= 0);
hist(detectedPhotons,50);
xlabel('Detected photons');
title('Distributions of all molecule fits');
subplot(2,2,2)
hist(totalPSFfits(:,23),50);
xlabel('Ratio of DH lobe amplitudes');
subplot(2,2,3)
hist(totalPSFfits(:,22),50);
xlabel('Distance between DH lobes (pixels)');
subplot(2,2,4)
normResiduals = totalPSFfits(:,16)*conversionFactor./totalPSFfits(:,21);
normResiduals = normResiduals(normResiduals >= 0 & normResiduals <= 10);
hist(normResiduals,50);
xlabel('Normalized fit residual');

%% ask user for new filter parameters
temp = inputdlg({'Range of detected photons (ex. [min max])', ...
                 'Max DH lobe amplitude', ...
                 'Range of distances between DH lobes', ...
                 'Range of normalized fit residuals'},...
                 'Input new filter parameters (pass-band)',...
                1,...
                {['[500 ' num2str(max(detectedPhotons)) ']'], ...
                 num2str(ampRatioLimit), ...
                 mat2str(lobeDistBounds), ...
                 '[0 3]'});        
if isempty(temp)
    return;
end

photonRange = str2num(temp{1});
ampRatioLimit = str2double(temp{2});
lobeDistBounds = str2num(temp{3});
residualBounds = str2num(temp{4});
clear temp;

%% apply filters to SM localization data
% run filter loop multiple times just in case error code status changes
% filters will only operate on good fits or their own error codes b/c
% we don't want these filters altering other error statuses than their own

for a=1:2
    eligibleFits = totalPSFfits(:,17)>0 | totalPSFfits(:,17)==-1005;
    % interlobe distance within bounds?
    good = totalPSFfits(:,22) >= lobeDistBounds(1) & totalPSFfits(:,22) <= lobeDistBounds(2);
    totalPSFfits(eligibleFits & ~good,17) = -1005;
    totalPSFfits(eligibleFits & good,17) = 10;      % a value that's unique and positive

    eligibleFits = totalPSFfits(:,17)>0 | totalPSFfits(:,17)==-1006;
    % amplitude ratio of lobes less than limit?
    good = totalPSFfits(:,23) <= ampRatioLimit;
    totalPSFfits(eligibleFits & ~good,17) = -1006;
    totalPSFfits(eligibleFits & good,17) = 10;
    
    eligibleFits = totalPSFfits(:,17)>0 | totalPSFfits(:,17)==-1007;
    % normalized error within limit?
    good = totalPSFfits(:,16)*conversionFactor./totalPSFfits(:,21) >= residualBounds(1) ...
        & totalPSFfits(:,16)*conversionFactor./totalPSFfits(:,21) <= residualBounds(2);
    totalPSFfits(eligibleFits & ~good,17) = -1007;
    totalPSFfits(eligibleFits & good,17) = 10;

    eligibleFits = totalPSFfits(:,17)>0 | totalPSFfits(:,17)==-1010;
    % detected photons within bounds?
    good = totalPSFfits(:,21) >= photonRange(1) & totalPSFfits(:,21) <= photonRange(2);
    totalPSFfits(eligibleFits & ~good,17) = -1010;
    totalPSFfits(eligibleFits & good,17) = 10;
end

%% replot distributions with filters applied
goodFits = totalPSFfits(:,17)>0;
figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w');
subplot(2,2,1)
detectedPhotons = totalPSFfits(goodFits,21);
hist(detectedPhotons,50);
xlabel('Detected photons');
title('Distributions of filtered molecule fits');
subplot(2,2,2)
hist(totalPSFfits(goodFits,23),50);
xlabel('Ratio of DH lobe amplitudes');
subplot(2,2,3)
hist(totalPSFfits(goodFits,22),50);
xlabel('Distance between DH lobes (pixels)');
subplot(2,2,4)
normResiduals = totalPSFfits(goodFits,16)*conversionFactor./totalPSFfits(goodFits,21);
hist(normResiduals,50);
xlabel('Normalized fit residual');

%% save filtered data
tempPSFfits = totalPSFfits;
for fileNum=1:length(fitFilePrefix)
    if fileNum == 1
        totalPSFfits = tempPSFfits(1:numFitsInFiles(fileNum),:);
    else
        totalPSFfits = tempPSFfits((sum(numFitsInFiles(1:fileNum-1))+1):sum(numFitsInFiles(1:fileNum)),:);
    end
    % save SM fits .mat file with filter limits
    save([fitFilePrefix{fileNum} 'molecule fits.mat'],'totalPSFfits',...
        'lobeDistBounds','ampRatioLimit','residualBounds','photonRange',...
        '-append');
end

end

