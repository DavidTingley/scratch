files = dir('*fig');
for f = 1:length(files)
%     if ~strcmp(files(f).name(1),'_')
    open(files(f).name)
    s = input('','s');
%     if isempty(s)
%         delete(files(f).name)
%     else
%         movefile(files(f).name,['_' files(f).name])
%     end
    close all
%     end
end