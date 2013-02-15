function f_debugMoleculeFits(totalPSFfits)
        [debugFile, debugPath] = uiputfile({'*.csv';'*.txt'},'Specify a path and name for debug histogram and full data output');
        if isequal(debugFile,0)
            return;
        end
        %makes sure output has an extension
        if debugFile(end-3) ~= '.'
            debugFile = [debugFile '.csv'];
        end
        edges = [-1007:-1000,-3,1,inf];
        vector = histc(totalPSFfits(:,17),edges);
        scrsz = get(0,'ScreenSize');
        hErrors=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');

        bar(vector(1:end-1))
        set(gca,'XTickLabel',{'LS error', 'amp ratio', 'lobe dist',...
                              'sig ratio', 'sig size', 'pks out',...
                              'amp< 0', 'guess out','fit err',...
                              'good fit'})
        set(gca,'FontSize',8);
                          
        print(hErrors,'-dpng',[debugPath debugFile(1:end-4) ' outcomes 1.png']);
        saveas(hErrors,[debugPath debugFile(1:end-4) ' outcomes.png']);
        % open a file for writing
        [fid,message] = fopen([debugPath debugFile], 'w');
        if ~isempty(message)
            error([debugPath debugFile ': ' message]);
            %return;
        end
        % print a title, followed by a blank line
        fprintf(fid, ['frame num,fit flag,SM idx in frame,template x (pix),template y (pix), template idx,' ...
            'match strength,amp1,amp2,peak1 x (pix),peak1 y (pix),' ...
            'peak2 x (pix),peak2 y (pix), sigma1 (pix),sigma2 (pix),mean bkgnd photons,'...
            'fit error,molecule x (nm),molecule y (nm),DHPSF angle,' ...
            'num photons,interlobe distance,amplitude ratio,sigma ratio,x fid-corrected (nm),y fid-corrected (nm), z fid-corrected (nm),'...
            'photons detected,mean background photons\n']);
        fclose(fid);

        dlmwrite([debugPath debugFile],totalPSFfits(:,[1 17 2:16 18:end]),'-append');
        disp('Debug files written successfully.');
end