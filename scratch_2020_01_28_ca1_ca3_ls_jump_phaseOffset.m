folders = dir('*201*');
count = 1;

for f = 1:length(folders)
    cd(folders(f).name)
    if exist([folders(f).name '.jumpPiezo.behavior.mat'])
        sessionInfo = bz_getSessionInfo;
        
%         if ~isempty(sessionInfo.ca1)
%             ca1 = bz_GetLFP(sessionInfo.ca1);
%             sessionInfo.recordingDuration = ca1.timestamps(end);
%         else
%             ca1 = [];
%         end
%         if ~isempty(sessionInfo.ca3)
%             ca3 = bz_GetLFP(sessionInfo.ca3);
%             sessionInfo.recordingDuration = ca3.timestamps(end);
%         else
%             ca3 = [];
%         end
%         if ~isempty(sessionInfo.ls)
%            ls = bz_GetLFP(sessionInfo.ls);
%             sessionInfo.recordingDuration = ls.timestamps(end);
%         else
%             ls = [];
%         end
                
        load([folders(f).name '.jumpPiezo.behavior.mat'])
        for cond = 1:length(unique(behavior.events.trialConditions))
            
            trials = find(behavior.events.trialConditions==cond);
            jumpTs = (behavior.events.jumpTime(trials,2));
            
            if ~isempty(sessionInfo.ca1)
            ca1 = bz_GetLFP(sessionInfo.ca1,'intervals',[jumpTs-1 jumpTs+1]);
            end
            if ~isempty(sessionInfo.ca3)
            ca3 =  bz_GetLFP(sessionInfo.ca3,'intervals',[jumpTs-1 jumpTs+1]);
            end
            if ~isempty(sessionInfo.ls)
            ls =  bz_GetLFP(sessionInfo.ls,'intervals',[jumpTs-1 jumpTs+1]);
            end

            for t = 1:length(trials)
                if ~isempty(sessionInfo.ca1) 
                    ca1_traces(count,:) = makelength(single(ca1(t).data),2501);
                else
                    ca1_traces(count,:) = nan(2501,1);
                end
                if ~isempty(sessionInfo.ca3)
                    ca3_traces(count,:) = makelength(single(ca3(t).data),2501);
                else
                    ca3_traces(count,:) = nan(2501,1);
                end
                if ~isempty(sessionInfo.ls)
                    ls_traces(count,:) = makelength(single(ls(t).data),2501);
                else
                    ls_traces(count,:) = nan(2501,1);
                end
                if isnan(ls_traces(count,1)) & isnan(ca1_traces(count,1)) & isnan(ca3_traces(count,1))
                    error
                end
                count = count + 1;
            end
        end
    end
   cd .. 
end


[b a] = butter(3,[5/625 14/625],'bandpass');

for i=1:count-1
   ca1_angles(i,:) = angle(hilbert(filtfilt(b,a,double(ca1_traces(i,:)))));
   ca3_angles(i,:) = angle(hilbert(filtfilt(b,a,double(ca3_traces(i,:)))));
   ls_angles(i,:) = angle(hilbert(filtfilt(b,a,double(ls_traces(i,:))))); 
end

for i=1:count-1
    ca1_ca3_diffs(i,:) = circ_dist(ca1_angles(i,:),ca3_angles(i,:));
    ca1_ls_diffs(i,:) = circ_dist(ca1_angles(i,:),ls_angles(i,:));
    ca3_ls_diffs(i,:) = circ_dist(ca3_angles(i,:),ls_angles(i,:));
end

subplot(3,2,1)
imagesc(ca1_ca3_diffs)
subplot(3,2,2)
imagesc(ca3_ls_diffs)
subplot(3,2,3)
imagesc(ca3_ls_diffs)


