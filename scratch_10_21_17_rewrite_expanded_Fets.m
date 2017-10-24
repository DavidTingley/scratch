
d  = dir('*201*');
for i=1:length(d)
    cd(d(i).name)
    sessionInfo = bz_getSessionInfo;
    if ~isempty(dir('*clu*'))
    spikes = bz_GetSpikes('noprompt',true);
    if ~isempty(spikes)
    for shank = 1:length(unique(spikes.shankID))
        cd(num2str(shank))
        tkwx = dir('*kwx');
        fets = h5read(tkwx.name,['/channel_groups/' num2str(shank) '/features_masks']);
        fets = double(squeeze(fets(1,:,:)));
        fet = fets';
        fetMultiplier = double(intmax('int32'))/max(abs(fets(:))); % masked klustakwik has small floats that need to be expanded before rounding below
        fet = fet .* fetMultiplier;
        cd ..


        BufSize = inf;
        nFeatures = size(fet, 2);
        formatstring = '%d';
        for ii=2:nFeatures
          formatstring = [formatstring,'\t%d'];
        end
        formatstring = [formatstring,'\n'];

        outputfile = fopen([sessionInfo.FileName '.fet.' num2str(shank)],'w');
        fprintf(outputfile, '%d\n', nFeatures);

        if isinf(BufSize)
            fprintf(outputfile,formatstring,round(fet'));
        else
            nBuf = floor(size(fet,1)/BufSize)+1;

            for i=1:nBuf 
                BufInd = [(i-1)*nBuf+1:min(i*nBuf,size(fet,1))];
                fprintf(outputfile,formatstring,round(fet(BufInd,:)'));
            end
        end
    end
    end
    end
    cd ..
end
