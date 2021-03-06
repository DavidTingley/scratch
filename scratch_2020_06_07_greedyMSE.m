


clear all
warning off
cgm_files = dir('_*mat');
hourHz = 1/(60*60);
dayHz = 1/(24*60*60);
nyquist = 5./(24*60)/2;

% % figure(1)
for ff=1:length(cgm_files)
    data{ff} = load(cgm_files(ff).name,'idx','zt','theta*','count','isa','isig_levels','spSlope','states','emgSig','mov');
    temp = load(cgm_files(ff).name(2:end),'idx');
    data{ff}.idx = temp.idx;
    count = data{ff}.count; 
    isig_levels = data{ff}.isig_levels;
    theta_resid = data{ff}.theta_resid; 
    theta_z = data{ff}.theta_z;
    theta_deltaR = data{ff}.theta_deltaR;
    spSlope = data{ff}.spSlope;
    emgSig = data{ff}.emgSig;
    states = data{ff}.states;
    zeitTime = discretize(data{ff}.zt,24);
        
    if isfield(data{ff},'mov')
        movement = data{ff}.mov;
        bad_movement = 0;
    else
        movement =  rand(1,length(count));%ones(1,length(count));
        bad_movement = 1;
    end
    pred = [nanZscore(theta_deltaR'),nanZscore(theta_z'),nanZscore(theta_resid'),nanZscore(count'),nanZscore(spSlope'),nanZscore(states'),nanZscore(movement'),nanZscore(emgSig')]; %,nanZscore(theta_z'),amps',dur',freq',;
    names = [{'theta_deltaR'},{'theta_z'},{'theta_resid'},{'count'},{'spSlope'},{'states'},{'movement'},{'emgSig'}]; % ,{'theta_z'},'amps','dur','freq',

    for i=1:size(pred(:,1))
        if  ~ismember(i,data{ff}.idx) | isnan(sum(pred(i,:)))
            pred(i,:) = nan;
        end
    end

    isa = [0,diff(data{ff}.isig_levels)];
    glucFlux = nanZscore(isa);
    cc = ccgBinned(count,[(glucFlux)],40);
    [a b] = min(cc);
%     glucFlux = circshift(glucFlux,2);

seq=[];
offs = [];
 for i=1:length(names)
    for j=1:length(names)
        for offset=-3:3
            mdl = fitlm(pred(:,[seq; j]),circshift(glucFlux,offset));
            err(offset+4) = mdl.RMSE;
        end
        [ms(j) offset(j)] = min(err);
    end
    
    ms(seq)=nan;
    [mins(i) b] = min(ms);
    seq = [seq;b];
    offs = [offs; offset(b)];
 end

 subplot(4,2,ff)
 plot(mins)
 title(num2str(seq'))
 
end
 