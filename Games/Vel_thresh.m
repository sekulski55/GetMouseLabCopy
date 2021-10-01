%velocity threshhold
%want game to change from phase 2 to 3 based on velocity of movement
%instant velocity = hand_dist / dt
%make initial vel thrshhold
%when it goes above 5m/s start, when it goes below stop


%Varibles from the code--
%[hX, hY] = GetMouse(window);
%[xCenter, yCenter] = RectCenter(windowRect);
%dt = toc-curtime;
curtime = 0;

%hand_dist = sqrt((hX-xCenter)^2 + (hY-yCenter)^2);

%we are in gamephase 2 --
Vels = [];
velocity_thresh = 3; %mm/sec, arbitrary 
velocity = hand_dist / dt;  %not sure how to check right now
velocity = hand_dist(dt)-hand_dist(dt-1) / dt; %would it be this
%velocity = hand_dist(cur_time) - hand_dist(curtime-1) /dt
%initial vel threshold

        if velocity >  velocity_thresh 
            if velocity < velocity_thresh
                VELs(trial) = Max(velocity);
                gamephase = 3;
            else 
        end 
