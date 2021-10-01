
function GetMouseVMR(name_prefix, tgt_file_name_prefix, tgt_set)
% Easiest way to test the game is to copy and paste: 
% GetMouseVMR('dummy', 'demo', '1') into command window

% Include this line so program can run
Screen('Preference', 'SkipSyncTests', 1);

Priority(1)


% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Initializing for PsychPortAudio
InitializePsychSound(1)

% See BasicAMAndMixScheduleDemo for loading multiple sounds.
% Open the default audio device [], with default mode [] (==Only playback),
% and a required latencyclass of zero 0 == no low-latency mode, as well as
% a frequency of freq and nrchannels sound channels.
% This returns a handle to the audio device:
% Try with the 'freq'uency we wanted:
freq = 22500;
nrchannels = 2; % for stereo, I think
pamaster = PsychPortAudio('Open', [], 1+8, 0, freq, nrchannels, [], 0.015); % last arg is for delay in sound onset
PsychPortAudio('Start', pamaster, 0, 0, 1);
pasound1 = PsychPortAudio('OpenSlave', pamaster, 1);
pasound2 = PsychPortAudio('OpenSlave', pamaster, 1);
pasound3 = PsychPortAudio('OpenSlave', pamaster, 1);

% Load sounds
load ding
load tooslow
load knock_short_quiet

% Fill the audio playback buffer with the audio data 'wavedata':
PsychPortAudio('FillBuffer', pasound1, ding);
PsychPortAudio('FillBuffer', pasound2, tooslow);
PsychPortAudio('FillBuffer', pasound3, knock_short_quiet);

% Clean PTB's pipes by emitting sound
PsychPortAudio('Volume', pasound1, 0.1)
PsychPortAudio('Start', pasound1, 1, 0, 0);
PsychPortAudio('Volume', pasound1, 1.0)

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% make sure monitor refresh is 144 hz
% desired_refresh = 30;
% actual_refresh = Screen('FrameRate',screenNumber);
% actual_refresh
% if actual_refresh < desired_refresh   
%     disp('Fix monitor refresh rate!')
%     return
% end

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
% WinTabMex(0, window); %Initialize tablet driver, connect it to 'win'
ListenChar(2)% -rm- no idea what this is
% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Load cursor image
[cursor_img, ~, cursor_alpha] = imread('cursor.png');
cursor_img(:,:,4) = cursor_alpha(:,:);

% This part depends on whether you want aiming landmarks or not:
% Load aiming landmarks
% I need to figure out best way to preallocate cell arrays. -hk
% % aim_img = cell(20000, 20000);
% % aimtext = cell(10000,1);
% % landmarks_tgt_angle = [60, 65, 70, 75, 80, 85, 90];

% for i = 1:length(landmarks_tgt_angle)
%     aim_img{i} = imread(['landmarks',num2str(landmarks_tgt_angle(i)),'.png']);
%     aimtext{i} = Screen('MakeTexture', window, aim_img{i}); 
% end  

mm2pixel = 3.6137;
cursortext=Screen('MakeTexture', window, cursor_img);
cursor_r = 1.75*mm2pixel;
resx = windowRect(3);
resy = windowRect(4);

% Load target file 
%cd('C:\Users\hyosub\Dropbox\HYO_extSEC\TargetFiles')
cd('TargetFiles')
% name_prefix = 'S1'; % dummy for now
% % pscyhtest = load([tgt_file_name_prefix,tgt_set,'.tgt']);
% tgt_file_name_prefix = 'tgt';
% tgt_set = '2';
tgt_file = dlmread([tgt_file_name_prefix,tgt_set,'.tgt'], '\t', 1, 0); % start reading in from 2nd row, 1st column
trial_num = tgt_file(:,1);
tgt_dist = tgt_file(:,2).*mm2pixel;
tgt_ang = tgt_file(:,3);
rotation = tgt_file(:,4);
aiming_targets = tgt_file(:,5);
aiming_landmarks = tgt_file(:,6);
fixed = [];
ccw = 0;
online_fb = tgt_file(:,7);
endpoint_fb = tgt_file(:,8);
clamped_feedback = tgt_file(:,9);
between_blocks = tgt_file(:,10);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Setup the text type for the window
Screen('TextFont', window, 'Ariel');
Screen('TextSize', window, 28);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Enable alpha blending for anti-aliasing
% For help see: Screen BlendFunction?
% Also see: Chapter 6 of the OpenGL programming guide 
% Shows up as rect instead of smooth dot without this line. -HK
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% some variables
%mm2pixel = 3.6137;  % may need to change depending on monitor 
start_tolerance = 10*mm2pixel;  % eventually change to 10
%startcirclewidth = 6*mm2pixel;
startcirclewidth = 5*mm2pixel;
%rt_dist_thresh = 10*mm2pixel;
rt_dist_thresh = 5*mm2pixel;
%targetsize = 6*mm2pixel;
targetsize = 3*mm2pixel;
%ringsize = 16*mm2pixel;
ringsize = 6*mm2pixel;
red = [255 0 0];
green = [0 255 0]; 
blue = [0 0 255];
insidetime = 0;
curtime = 0;
endptfbtime = 1.0;
pausetime = 0.0;
wait_time = 0.5;
rt = 0;
RTs = [];
mt = 0;
MTs = [];
searchtime = 0;
SearchTimes = [];
fb_time = 0;
hits = 0;
data = [];
hand_vel=0;
prev_vel =0;
%Vels = [];
c=0;


gamephase = 0;
trial = 1;
cursor = [];

% Target location:
tgtx = xCenter + tgt_dist.*cosd(tgt_ang);
tgty = yCenter - tgt_dist.*sind(tgt_ang);
tgtloc = [tgtx tgty];

% aim landmark variables
% fixed = 0;  % default
% increment = 5.625;  % default

% hit will be any part of cursor touching target
hit_tolerance = targetsize./2 + cursor_r;

% We will store cursor position data. -HK
% 1/15/16 update: not storing cursor position this way anymore.
% Artifact of using GetCursor
% cursorx = [];
% cursory = [];

% I think there are 2540 lines per inch (lpi) on tablet
% tablet active area is 19.2 in x 12.0 in
tablet_x_scale = 1/27.625;
tablet_x_offset = -1.1969*2540;
tablet_y_scale = -1/27.625;
tablet_y_offset = 11.724*2540;
% WinTabMex(2); %Empties the packet queue in preparation for collecting actual data

% Set the mouse to the center of the screen to start with
HideCursor;
% 
%     center_tab_x = (xCenter/tablet_x_scale) - tablet_x_offset;
%     center_tab_y = (yCenter/tablet_y_scale) - tablet_y_offset;
%     SetMouse(center_tab_x, center_tab_y, window);
    

% See Matlab for Behavioral Scientists by David A. Rosenbaum. 
desiredSampleRate = 500;
k = 0;
tab_k = 15; 
numtrials = size(tgt_file,1);

% Define the ESC key
KbName('UnifyKeynames');
esc = KbName('ESCAPE');
space = KbName('SPACE');

maxtrialnum = max(numtrials);

hand_angle = nan(maxtrialnum,1);

% Variables that store data -all copied from Ryan's code
MAX_SAMPLES=6e6; %about 1 hour @ 1.6kHz = 60*60*1600
Trials_Samples = 100000;
% timevec=nan(MAX_SAMPLES,1);
% delay_calc_time=nan(MAX_SAMPLES,1);
gamephase_move=nan(MAX_SAMPLES,1);
tablet_queue_length=nan(MAX_SAMPLES,1);
thePoints=nan(MAX_SAMPLES,2);
theHandDist=nan(MAX_SAMPLES,1); %make array for hand dist hopefully
Velocity = zeros(Trials_Samples,1);
Velocity1=[];
Hand_diff = zeros(Trials_Samples,1);
Hand_diff_real = [];
Vels = zeros(Trials_Samples,1);
u=0;
v=1;
w=0;

cursorPoints=nan(MAX_SAMPLES,2);
tabletPoints=uint16(nan(MAX_SAMPLES/8,2)); %reduce # of samples since the tablet is sampled @ 200Hz
tabletTime=nan(MAX_SAMPLES/8,1);
% total_vel=nan(MAX_SAMPLES,1);
% total_displacement=nan(MAX_SAMPLES,1);
% index_of_point_shown = nan(MAX_SAMPLES,1);
% [deltax,deltay]=deal(nan(MAX_SAMPLES,1));
dt_all = nan(MAX_SAMPLES,1);
t = nan(MAX_SAMPLES,1);
trial_time = nan(MAX_SAMPLES,1);
phase2_time = nan(MAX_SAMPLES,1);
trial_move = nan(MAX_SAMPLES,1);
start_x_move = nan(MAX_SAMPLES,1);
start_y_move = nan(MAX_SAMPLES,1);
rotation_move = nan(MAX_SAMPLES,1);

tic;
begintime = GetSecs;
nextsampletime = begintime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Create sine wave rotation
    %rotation = ones(50,1) %new rotation, delete this
    x = 0:0.2513:4*pi;    %makes a plot line for the wave with 50 (maxtrial) intervals
    Sin_rotation = sin(x);   %rotation varible x 10 so I could see
    i=11;
    while i<=50;
        rotation(i)= rotation(i)*Sin_rotation(i-10);
        i=i+1;
    end
    Sine = sin(2*pi*trial/maxtrialnum);
    %A * Sin( 2 * pi * f * n /N)
    
    %creates the negative sine wave gain
    x2 = 0:0.2513:4*pi;    %makes a plot line for the wave with 50 (maxtrial) intervals
    Sin_rotation2 = -sin(x2)+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Loop game until over or ESC key press
while trial <= maxtrialnum;   %     
    
    % Exits experiment when ESC key is pressed. 
    [keyIsDown, secs, keyCode] = KbCheck;
    if keyIsDown
        if keyCode(esc)
%             Screen('CloseAll')
            break
        end
    end
         
    k = k+1;    % just to stay consistent with Ryan's code
    t(k) = GetSecs - begintime;
    dt = toc-curtime;
    dt_all(k) = dt;
    
    if k == 1;
        trial_time(k) = dt;
    else
        trial_time(k) = trial_time(k-1) + dt;
    end
    curtime = toc;
    
    % Flip to the screen
    % last argument - 1: synchronous screen flipping, 2:asynchronous screen flipping
    Screen('Flip', window, 0, 0, 2);
     
    % Record trial number
    trial_move(k) = trial;
    rotation_move(k) = rotation(trial,1);
        
    % Read information from the tablet
%     pkt = WinTabMex(5); % reads the latest data point out of a tablet's event queue
    tablet_queue_length(k) = 0;
    
%     while ~isempty(pkt) % makes sure data are in packet; once pkt is 'grabbed,' then rest of code executes     
%         tabletPoints(tab_k,1:2) = pkt(1:2)';    % placing x and y (pkt rows 1,2) into tabletPoints variable
%         tabletTime(tab_k) = (pkt(6)-tabletTime(16))/1000;   % tab_k initialized to 15; giving a little buffer at start of game?
%         tab_k = tab_k+1;    % now tab_k is just another iterating variable
%         tablet_queue_length(k) = tablet_queue_length(k)+1;  % adding each loop through
%         pkt = WinTabMex(5); % reads the latest data point out of a tablet's event queue
%     end
    
    % HAND COORDINATES
    % x,y coordinates from WinTabMex pkt
    [hX, hY] = GetMouse(window);
    
    %hX=Sin_rotation2(trial)*hX; %multplying the mouse position by the negative sine wave to make the gain
    %hY=Sin_rotation2(trial)*hY;
    %why do values set under .5 trap the cursor in a box?
    
    thePoints(k,:) = [hX hY]; % record full precision points
    
    hand_dist = sqrt((hX-xCenter)^2 + (hY-yCenter)^2);
    %theHandDist(k,:) = hand_dist; 
    
    %makes sure that the distances are evenly spaced apart
    %persistant varibles: remmebered in function workspaces
    %set up logical statment comparing samples that wont store comsecutive
    %distances 
    %would have to change dt as well
    
    
    %Create sine wave rotation
    %rotation = ones(50,1) %new rotation, delete this
    %x = 0:0.2513:4*pi;    %makes a plot line for the wave with 50 (maxtrial) intervals
    %Sin_rotation = 10*sin(x);   %rotation varible x 10 so I could see
    %rotation(trial)= rotation(trial)*Sin_rotation(trial); %makes new sin rotated array
    
    %array = 1:1:50
    
    
    
    % ROTATED CURSOR (including clamp)
    if clamped_feedback(trial,1) == 1
        % Clamped fb location:
        rcX = xCenter + hand_dist.*cosd(tgt_ang(trial,1) + rotation(trial,1)); % may need to subtract rotation in order to make (+) clamp CCW
        %rcX = 2*xCenter + hand_dist.*cosd(tgt_ang(trial,1) + rotation(trial,1));
        rcY = yCenter - hand_dist.*sind(tgt_ang(trial,1) + rotation(trial,1));
        %rcY = 2*yCenter - hand_dist.*sind(tgt_ang(trial,1) + rotation(trial,1));
    else
        [rcX_rotated, rcY_rotated] = rotatexy(round(hX)-xCenter,(round(hY)-yCenter),rotation(trial,1),1);
        %[rcX_rotated, rcY_rotated] = rotatexy(round(hX)-xCenter,(round(hY)-yCenter),rotation(trial,1),Sin_rotation2(trial));
        rcX = rcX_rotated + xCenter;
        rcY = rcY_rotated + yCenter;
    end
   
    % Draw home position.
    if gamephase == 0   % Searching for start location
        searchtime = searchtime + dt;
        SearchTimes(trial) = searchtime;
        Screen('FrameOval', window, white, [xCenter - startcirclewidth/2, yCenter - startcirclewidth/2, xCenter + startcirclewidth/2, yCenter + startcirclewidth/2], 2);

        if hand_dist < start_tolerance*20 
            visible = 1;
%             Screen('FrameOval', window, white, [xCenter - hand_dist, yCenter - hand_dist, xCenter + hand_dist, yCenter + hand_dist], 1);
%         elseif hand_dist < (startcirclewidth/2 + cursor_r)
%             visible = 1;
        else
            visible = 0;
        end
            
        % calculate distance of cursor from start position 
        if hand_dist < startcirclewidth/2
            Screen('DrawDots', window, [xCenter yCenter], startcirclewidth, white, [], 2);
            inside = 1;
            insidetime = insidetime + dt;
        else  
            inside = 0;
            insidetime = 0;
        end
    
        % Starting to use idea of game phases to signify different phases of
        % the game. Similar to use of flags in Pygame scripts.
        if inside ==  1 && insidetime > wait_time
            gamephase = 1;
            insidetime = 0;
        end
    
    elseif gamephase == 1;  % Show target
        visible = online_fb(trial,1);
        %online_fb = ones(50,1);
        %visible = 1;
        %visible=online_fb(trial,1);
        
        if aiming_landmarks(trial);
            Screen('DrawTexture', window, aimtext{aiming_landmarks(trial)},[],[],[],[],[],white);
            %aim_landmarks(window,tgt_dist(trial),xCenter,yCenter,tgt_ang(trial),fixed,ccw,increment,white)
        end
        
%         Screen('DrawDots', window, tgtloc(trial,:), targetsize, green, [], 2);
        Screen('DrawDots', window, tgtloc(trial,:), targetsize, blue, [], 2);
        
        rt = rt + dt;               
        if hand_dist >  rt_dist_thresh 
            RTs(trial) = rt;
            gamephase = 2;
            phase2_time(k) = 0;
        end 
          
    elseif gamephase == 2;  % Moving towards target
        visible = online_fb(trial,1);
        %visible=0;
        mt = mt + dt;
        MTs(trial) = mt;
        
        if aiming_landmarks(trial)
            Screen('DrawTexture', window, aimtext{aiming_landmarks(trial)},[],[],[],[],[],white);
            %aim_landmarks(window,tgt_dist(trial),xCenter,yCenter,tgt_ang(trial),fixed,ccw,increment,white)
        end
        
        Screen('DrawDots', window, tgtloc(trial,:), targetsize, blue, [], 2);
        

        
        theHandDist(k,:) = hand_dist/mm2pixel; %hand dist at current sample
        Hand_diff(k,:) = theHandDist(k) - theHandDist(k-1); 
         
        phase2_time(k) = phase2_time(k-1) + (trial_time(k)-trial_time(k-1)) + dt; %calculates time during phase 2
        Velocity(k) = (theHandDist(k) - theHandDist(k-1)) / (dt); %calculates velocity
        
        %these if statements make sure that if the game is sampling too
        %fast and gives a "false 0" for the velocity it will find the
        %average of two other close velocity points and make it equal to
        %that. This should make sure velocity is only 0 when user not
        %moving
        
         if Velocity(k-3)==0
             Velocity(k-3)= (Velocity(k-2)+Velocity(k-1))/2;
         end 
         
         if Velocity(k-2)==0
             Velocity(k-2)= (Velocity(k-1)+Velocity(k))/2;
         end 
         
         if Velocity(k-1)==0
             Velocity(k-1)= (Velocity(k-2)+Velocity(k))/2;
         end
 
         %if loop if abs hand(k)- hand(k-1) >0
         %handvel is abs hand(k) - hand(k-1) / time
         %prevvel = currentvel
         %else handvel = prev_vel to make sure it doesnt sample too fast
         if abs(theHandDist(k) - theHandDist(k-1)) > 0
             %cur_samptime = Getsecs;
             hand_vel = abs((theHandDist(k) - theHandDist(k-1)/(dt_all(k)-dt_all(k-c))));
             hand_vel = hand_vel / 1000;
             prev_vel = hand_vel;
             c=0;
         else 
             hand_vel = prev_vel;
             c=c+1;
         end
         

         
         Vels(k) = hand_vel;
 
%         
         %5 cm/sec threshhold 
         %understand pixels to time, and convert
         %apporximate pixel to cm
         %make the threshhold start the movement, and end it
         %save as tgt, make the first 5 trials not visible
        
        
        
         if hand_dist >= tgt_dist(trial,1) %had greater than
         
         %New if statement makes sure that the user moves at least slightly
         %during phase 2, and checks to see if both current sample
         %velocity and k-3 velocity are close enough to 0. This should
         %happen only when the user stops moving
         
       %if phase2_time(k) >= 0.5 && Velocity(k)<=1 && Velocity(k-3)<=1
       %if phase2_time(k) >= 0.5 && Vels(k)==Vels(k-3)
            fb_angle = atan2d(rcY-yCenter, rcX-xCenter);
            fb_x = tgt_dist(trial,1)*cosd(fb_angle) + xCenter;
            fb_y = tgt_dist(trial,1)*sind(fb_angle) + yCenter;
            
            hand_angle(trial,1) = atan2d((hY-yCenter)*(-1), hX-xCenter);
            
            if MTs(trial) <= 0.3
                PsychPortAudio('Start', pasound3, 1, 0, 0); % the last input argument '0' instructs Matlab to execute code immediately
            elseif MTs(trial) > 0.3
                PsychPortAudio('Start', pasound2, 1, 0, 0);
            end
            
            if sqrt((fb_x - tgtx).^2 + (fb_y - tgty).^2) <= hit_tolerance
                hits = hits + 1;    % Make sure this is correct
            else
                hits = hits;
            end
            visible = 0;
            gamephase = 3;
            w=0;
        end
             
    elseif gamephase == 3;  % Endpoint feedback 
        %endpoint_fb = ones(50,1); %makes all trials have endpoint feedback
        visible = endpoint_fb(trial,1);
        
        if aiming_landmarks(trial)
            Screen('DrawTexture', window, aimtext{aiming_landmarks(trial)},[],[],[],[],[],blue);
            %aim_landmarks(window,tgt_dist(trial),xCenter,yCenter,tgt_ang(trial),fixed,ccw,increment,white)
        end
            
        if fb_time <= endptfbtime

            if sqrt((fb_x - tgtx).^2 + (fb_y - tgty).^2) <= hit_tolerance
                % Change color or play sound or both
                Screen('DrawDots', window, tgtloc(trial,:), targetsize, blue, [], 2);  % just testing hit thresh for now
            else
                Screen('DrawDots', window, tgtloc(trial,:), targetsize, blue, [], 2);
            end
            
            fb_time = fb_time + dt;
        
        else
            gamephase = 4;
            visible = 0;
        end
        
    elseif gamephase == 4;  % Between Blocks Message
        
        trial_time(k) = 0;
        if between_blocks(trial) == 1
            Screen('DrawText', window, 'Great job!' , xCenter-100, yCenter, white);
            
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyIsDown
                if keyCode(space)
                    gamephase = 0;
                    fb_time = 0;
                    searchtime = 0;
                    rt = 0;
                    mt = 0;
                    trial_time(k) = 0;
                    phase2_time(k)=0;
                    trial = trial + 1;
                end
            end
            
        elseif between_blocks(trial) == 2
            Screen('DrawText', window, 'Ignore the cursor. Move your hand directly to the target.' , xCenter-450, yCenter, white);
            
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyIsDown
                if keyCode(space)
                    gamephase = 0;
                    fb_time = 0;
                    searchtime = 0;
                    rt = 0;
                    mt = 0;
                    trial_time(k) = 0;
                    phase2_time(k)=0;
                    trial = trial + 1;
                end
            end
            
        elseif between_blocks(trial) == 3
            Screen('DrawText', window, 'Move your hand directly to the target.' , xCenter-200, yCenter, white);
            
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyIsDown
                if keyCode(space)
                    gamephase = 0;
                    fb_time = 0;
                    searchtime = 0;
                    rt = 0;
                    mt = 0;
                    trial_time(k) = 0;
                    phase2_time(k)=0;
                    trial = trial + 1;
                end
            end
            
        elseif between_blocks(trial) == 4
            Screen('DrawText', window, 'The experiment is over. Thank you for participating!' , xCenter-430, yCenter, white);
            
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyIsDown
                if keyCode(space)
                    trial = trial + 1;
                    gamephase = 0;
                    fb_time = 0;
                    searchtime = 0;
                    rt = 0;
                    mt = 0;
                    trial_time(k) = 0;
                    phase2_time(k)=0;
                end
            end
            
        else
            gamephase = 0;
            fb_time = 0;
            searchtime = 0;
            rt = 0;
            mt = 0;
            trial_time(k) = 0;
            phase2_time(k)=0;
            trial = trial + 1;
        end        
    end
    
    % Draw Cursor
    if visible
        %hX=2*hX;
        %hY=2*hY;
        if (gamephase == 0)
            cursor = [(hX - cursor_r) (hY - cursor_r) (hX + cursor_r) (hY + cursor_r)];
            %cursor = 2*cursor;
            cursorPoints(k,:) = [hX hY]; % record full precision points
            Screen('DrawTexture', window, cursortext, [], cursor, [], [], [],[255 255 255]);
            
        elseif (gamephase == 1 || gamephase == 2)
            cursor = [(rcX - cursor_r) (rcY - cursor_r) (rcX + cursor_r) (rcY + cursor_r)];
            cursorPoints(k,:) = [rcX rcY]; % record full precision points
            Screen('DrawTexture', window, cursortext, [], cursor, [], [], [],[255 255 255]);
        elseif  (gamephase == 3)
            %cursor = [fb_x - cursor_r, fb_y - cursor_r, fb_x + cursor_r, fb_y + cursor_r];
            cursor = [(hX - cursor_r) (hY - cursor_r) (hX + cursor_r) (hY + cursor_r)];
            cursorPoints(k,:) = [fb_x fb_y]; % record full precision points
            Screen('DrawTexture', window, cursortext, [], cursor, [], [], [],[255 255 255]);
        end
    else
        cursorPoints(k,:) = [NaN NaN];
    end
    
    gamephase_move(k) = gamephase;
    
    sampletime(k) = GetSecs;
    nextsampletime = nextsampletime + 1/desiredSampleRate;
    
    while GetSecs < nextsampletime
    end
    
end
%end

endtime = GetSecs
elapsedTime = endtime - begintime
numberOfSamples = k
actualSampleRate = 1/(elapsedTime / numberOfSamples)
thePoints(1:k,:);

ShowCursor;
% Clear the screen
sca;
% WinTabMex(3); % Stop/Pause data acquisition.
% WinTabMex(1); % Shutdown driver.
ListenChar(0);

% Stop playback:
PsychPortAudio('Stop', pasound1);
PsychPortAudio('Stop', pasound2);

% Close the audio device:
PsychPortAudio('Close')

% timevec = single(timevec(1:k));
% delay_calc_time = single(delay_calc_time(1:k));
% tablet_queue_length = tablet_queue_length(1:k);
% thePoints = single(thePoints(1:k,:));
% index_of_point_shown = index_of_point_shown(1:k);
% tabletPoints = tabletPoints(1:tab_k,:);
% tabletTime = tabletTime(1:tab_k,:);
%timeincirc=single(timeincirc(1:hitcirc_count-1,:));


% Game file
hand_angle = hand_angle;

% Movement file 
trial_move = trial_move;
gamephase_move = gamephase_move;
t = t;
dt_all = dt_all;
rotation_move = rotation_move;
start_x_move = ones(k,1).*xCenter;
start_y_move = ones(k,1).*yCenter;
hand_x(:,1) = thePoints(:,1)-xCenter;
hand_y(:,1) = (thePoints(:,2)-yCenter)*(-1); % adjust for other monitor points
cursor_x(:,1) = cursorPoints(:,1) - xCenter;
cursor_y(:,1) = (cursorPoints(:,2) - yCenter)*(-1); % adjust for other monitor points

% clear sounds
clear ding tooslow aim_img

%
% cd('C:\Dropbox\HYO_extSEC\Data')
cd('Data')

% Save data
name_prefix_all = [tgt_file_name_prefix,'_',name_prefix];
disp('Saving...')
if ~exist([name_prefix_all,'.mat'],'file')
    datafile_name = [name_prefix_all,'.mat'];
  
elseif ~exist([name_prefix_all,'_a.mat'],'file'), datafile_name = [name_prefix_all,'_a.mat'];
elseif ~exist([name_prefix_all,'_b.mat'],'file'), datafile_name = [name_prefix_all,'_b.mat'];
else
    char1='c';
    while exist([name_prefix_all,'_',char1,'.mat'],'file'), char1=char(char1+1); end
    datafile_name = [name_prefix_all,'_',char1,'.mat'];
end
save(datafile_name); 
disp(['Saved ', datafile_name]);
%end


function [rx, ry] = rotatexy(x,y,phi,gain)
% phi is in degrees
phi=phi*pi/180;
[theta r]=cart2pol(x,y);
[rx ry]=pol2cart(theta-phi,gain*r);
return

% function aim_landmarks(window,tgt_dist,xCenter,yCenter,tgt_ang,fixed,ccw,increment,color)
% % This function draws up the numbered aiming landmarks: options for fixed
% % vs rotated and CW vs CCW.
% for i = 1:(360/increment)
%     if ccw
%         aimx(i) = xCenter + tgt_dist*cosd(tgt_ang + i.*increment);
%         aimy(i) = yCenter - tgt_dist*sind(tgt_ang + i.*increment);
%     else
%         aimx(i) = xCenter + tgt_dist*cosd(tgt_ang - i.*increment);
%         aimy(i) = yCenter - tgt_dist*sind(tgt_ang - i.*increment);
%     end
%     
%     if fixed == 1
%         if i < (360./increment)
%             num2display(i) = i;
%             Screen('DrawText', window, num2str(num2display(i)), aimx(i), aimy(i), color);
%         end
%         
%     elseif fixed ~= 1
%         if i < ((360./increment)/2)
%             num2display(i) = i;
%             Screen('DrawText', window, num2str(num2display(i)), aimx(i), aimy(i), color);
%         elseif i == (360./increment./2)
%             Screen('DrawDots', window, [aimx(i)+9, aimy(i)+9], 18, color, [], 2);   % I think position offsets should equal half circle size (eg, aimx + 9 ..., 18)
%         elseif (i > ((360./increment)./2)) && (i < (360./increment))
%             num2display(i) = i - 360./increment;
%             Screen('DrawText', window, num2str(num2display(i)), aimx(i), aimy(i), color);
%         end
%         
%     end
% end
% 
% return
% 
% 
%  