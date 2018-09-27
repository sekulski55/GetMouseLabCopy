% script to make .tgt files
clear all
close all

namefile = 'extTGT_';
rot = 1.75;
tgt_diam = 6;  % 6 for small; 16 for big

for subnum = 31:32
    
    % determine if CW or CCW
    if mod(subnum,2)==0
        rot_size = rot*-1;
    else
        rot_size = rot;
    end
    
    % different blocks of experiment
    nofb_base = 40;
    fb_base = 40;
    practice = 4;   % inputting as 4 because of # of trials per cycle
    error_clamp = 880; 
    nofb_post = 40;
    fb_post = 40;
     
    trials1 = nofb_base + fb_base + practice + error_clamp + nofb_post + fb_post;
%     trials1 = nofb_base + fb_base + practice + error_clamp + nofb_post;
    baseline_total = nofb_base + fb_base;
    total_trials = trials1;
    
    % add 3 to sum of nofb_base + fb_base
    lastpracticetrial = 83;
    
    startofclamp = 85;
    
    % Constant values
    T.trialnum = (1:total_trials)';
    T.tgt_dist = ones(total_trials, 1)*80;
    T.tgt_angle = ones(total_trials, 1);
    T.rotation = zeros(total_trials, 1);
    T.aiming_targets = ones(total_trials, 1);
    T.aiming_landmarks = zeros(total_trials, 1);
    T.online_fb = ones(total_trials, 1);
    T.endpoint_fb = ones(total_trials, 1);
    T.clamped_fb = zeros(total_trials, 1);
    T.between_blocks = zeros(total_trials, 1);
    T.tgtsize = ones(total_trials,1)*tgt_diam;
     
    probe_tgt_angles = [45 135 225 315];  % target set

    b1 = 10; % no fb baseline
    b2 = b1 + 10; % fb baseline
    pr1 = b2 + 1; % practice
    r1 = pr1 + 220; % perturbation
    p1 = r1 + 10; % 0 clamp post
    p2 = p1 + 10; % fb post
    
    tpe = 4; % trials per epoch - allow for some flexibility for # of tgts
    
    for i = 1:p2;
        
        if i > 0 & i < b1+1 % no fb baseline
            T.tgt_angle((i*tpe)-(tpe-1):(i*tpe),1) = probe_tgt_angles(randperm(tpe));
            T.online_fb((i*tpe)-(tpe-1):(i*tpe),1) = 0;
            T.endpoint_fb((i*tpe)-(tpe-1):(i*tpe),1) = 0;
            
        elseif i > b1 & i < b2+1 % fb baseline
            T.tgt_angle((i*tpe)-(tpe-1):(i*tpe)) = probe_tgt_angles(randperm(tpe));
            T.online_fb((i*tpe)-(tpe-1):(i*tpe)) = 1;
            T.endpoint_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
            
        elseif i > b2 & i < pr1 + 1 % practice
            T.online_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
            T.endpoint_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
            T.tgt_angle((i*tpe)-(tpe-1):(i*tpe),1) = 45; % choose angle for practice
            T.rotation((i*tpe)-(tpe-1):(i*tpe),1) = rot_size;
            T.clamped_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
            
        elseif i > pr1 & i < r1+1 % clamp
            T.online_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
            T.endpoint_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
            T.tgt_angle((i*tpe)-(tpe-1):(i*tpe),1) = probe_tgt_angles(randperm(tpe));
            T.rotation((i*tpe)-(tpe-1):(i*tpe),1) = rot_size;
            T.clamped_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
     
        elseif i > r1 & i < p1+1 % 0 deg clamp 
            T.tgt_angle((i*tpe)-(tpe-1):(i*tpe),1) = probe_tgt_angles(randperm(tpe));
            T.online_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
            T.endpoint_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
            T.rotation((i*tpe)-(tpe-1):(i*tpe),1) = 0;
            T.clamped_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
 
        elseif i > p1 & i < p2+1 % fb post
            T.tgt_angle((i*tpe)-(tpe-1):(i*tpe),1) = probe_tgt_angles(randperm(tpe));
            T.online_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
            T.endpoint_fb((i*tpe)-(tpe-1):(i*tpe),1) = 1;
            T.clamped_fb((i*tpe)-(tpe-1):(i*tpe),1) = 0;
        end 
        
    end
    
    % decide here where you want to place between block messages
    T.between_blocks([b1*tpe], 1) = 1;
    T.between_blocks([b2*tpe], 1) = 5;
    T.between_blocks([b2*tpe+1], 1) = 6;
    T.between_blocks([b2*tpe+2], 1) = 7;
    T.between_blocks(lastpracticetrial, 1) = 2;
    T.between_blocks([r1*tpe-600], 1) = 2;
    T.between_blocks([r1*tpe-320], 1) = 2;
    T.between_blocks([r1*tpe-40], 1) = 2;
    T.between_blocks([p1*tpe], 1) = 3;
    T.between_blocks([p2*tpe], 1) = 4;
    
    % this line is necessary if you are adding practice trials. you must 
    % delete superfluous trials that are artifact of looping by cycles.
    T = structfun(@(x) (x([(1:lastpracticetrial),(startofclamp:end)])), T, 'uniformoutput', 0);
    
    dummy = struct2dataset(T); %%% These two lines are to convert structure into double
    set = double(dummy)
    
    set(:,1) = 1:size(set,1);
    %Add in last Betweenblocks
    %set(end,15) = 1;
    %Variables header file
    header = {'trialnum','tgt_distance','tgt_angle','rotation',...
        'aiming_targets','aiming_landmarks','online_fb', 'endpoint_feedback',...
        'clamped_fb','between_blocks','target_size'};
    
    filename = strcat(namefile,num2str(subnum),'.tgt');
    %If you ever need strings in here use this way
    fid = fopen(filename,'wt');
    [rows,cols] = size(set);
    fprintf(fid,'%s\t',header{:});
    for i = 1:rows
        fprintf(fid,'\n');
        for j = 1:cols
            fprintf(fid,'%3.2f\t',set(i,j));
        end
    end
    fclose(fid)
    
end

    
    
   