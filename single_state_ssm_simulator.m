function r_hat = single_state_ssm_simulator(A,B,r,numtrials,clamp)

%%% This function was written to simulate a single-rate state space model  

%%% Other pointers: numtrials should be input as 1x3 vector (ie,[# of 
%%% baseline trials, # of perturbation trials, # of washout trials]);
%%% for a standard visuomotor rotation, set "clamp", to zero. Copy next 
%%% line into command window to see example of this model learning a 30 deg
%%% visuomotor rotation:
%%% single_state_ssm_simulator(.96,.3,30,[10 50 50],0)

%initialize internal model
r_hat=zeros(1,sum(numtrials));

%perturbation schedule
rot=[zeros(1,numtrials(1)) ones(1,numtrials(2)) zeros(1,numtrials(3))];
clamp = clamp*rot;

%if (r)otation is input as a multi-element vector, then use that schedule;
%otherwise, if (r)otation is a scalar, then the rotation is either ON or
%OFF and does not change size
if length(r) > 1
    r = r;
else    
    r = r*rot;
end

%plot figure
figure()
hold on
xlabel('Movement cycle')
ylabel('Hand angle (deg)')
title('State space model simulation','fontsize',14)
xlim([0 sum(numtrials)+2])
ylim([-5 40])
plot(1:length(r),repmat(0,length(r),1),'linewidth',1,'color','k')
h1 = plot(1:length(r),r,'linewidth',2,'color','k')

for i=1:sum(numtrials)-1    
    %if clamp, then the error is constant; otherwise, error is updated
    %every trial
    if clamp(i)==1
        e(i) = r(i);
    elseif clamp(i) == 0 
        e(i) = r(i) - r_hat(i); %equation 1 in Taylor & Ivry (2014)
    end
    
    %internal model
    r_hat(i+1)=A*r_hat(i)+B*e(i); %equation 2 in Taylor & Ivry (2014)
    
    h2 = plot(r_hat(1:i),'b','LineWidth',3)
 
    drawnow    
end
legend([h1 h2],'Perturbation','Internal model')

