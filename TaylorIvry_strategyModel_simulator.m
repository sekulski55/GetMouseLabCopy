%Taylor & Ivry strategy model (Plos Comp Bio - 2011). Running this script
%will generate something that should like Figure 6A from the paper. Main 
%point of model is that strategy acts on target error and adaptation 
%(internal model/state estimate) updates based on the aiming error, or 
%sensory prediction error (SPE). --Hyosub (9/28/2018)--

clear all; close all

A = .991; %retention factor for adaptation 
B = .012; %learning rate for adaptation
E = .999; %retention factor for strategy
K = .985; %certainty in aiming target location
F = .023; %learning rate for strategy (how quickly aiming gets adjusted)

nBsl = 120;
nRot = 322;
nWO = 80;

s(1:nBsl+2) = [zeros(nBsl+2,1)]; %strategy is zero for baseline and first two rotation
s(nBsl+3) = 45; %participant told to aim to 45 on 3rd trial of rotation                                      
                                         
r = [zeros(nBsl,1); ones(nRot,1)*-45; zeros(nWO,1)]; %45 deg VMR
r_est(1) = 0; %initializing IM (state estimate)

sigma_m = 5; %motor variability of our humanoid subject

for n=1:length(r)-1
    
    %Full model
    e_tgt(n) = s(n) + (r(n)-r_est(n));  %eqn3: the error you see is the sum
                                        %of your strategy + the rotation - 
                                        %state estimate
    
    e_aim(n) = r(n) - r_est(n) + (s(n) - K*s(n));   %eqn6: captures idea of 
                                                    %uncertainty in desired
                                                    %aim location
    
    r_est(n+1) = A*r_est(n) + B*e_aim(n); %eqn2: state update is based on 
                                          %aiming error NOT target error 
                                          %(this is key to model)
    
    if n<=nBsl+2 
        s(n) = s(n);
    elseif n>nBsl+2 & n<nBsl+nRot
        s(n+1) = E*s(n) - F*e_tgt(n); %adjusting strategy begins in 
                                      %perturbation block
    elseif n>=nBsl+nRot
        s(n+1) = 0; %we have to set participant's strategy back to zero 
                    %after perturbation is turned off
    end
    
    humanoid(n) = e_tgt(n) + randn(1)*sigma_m; %I'm just corrupting the model 
                                         %with Gaussian noise so that it
                                         %looks like participant in Fig. 6A
                                
end

figure; hold on
ylim([-50 50])
xlim([0 525])
title('Model Fit for AT','fontsize',14)
ylabel('Target error (deg)','fontsize',12)
xlabel('Movement Number','fontsize',12)
line(xlim,[0 0],'color','k','linewidth',2)
line(xlim,[45 45],'color','k','linewidth',2)
line(xlim,[-45 -45],'color','k','linewidth',2)
line([120 120],ylim,'linestyle','--','color','k')
line([442 442],ylim,'linestyle','--','color','k')
plot(humanoid,'k','linewidth',3)
plot(e_tgt,'b','linewidth',5)


%plot Model Learning for Aiming Target group (Fig. 6D)
figure; hold on
ylim([-50 50])
xlim([0 525])
title('Model Learning','fontsize',14)
ylabel('Angle (deg)','fontsize',12)
xlabel('Movement Number','fontsize',12)
line(xlim,[0 0],'color','k','linewidth',2)
line(xlim,[45 45],'color','k','linewidth',2)
line(xlim,[-45 -45],'color','k','linewidth',2)
line([120 120],ylim,'linestyle','--','color','k')
line([442 442],ylim,'linestyle','--','color','k')
h1=plot(s,'--b','linewidth',2)
h2=plot(r_est,'b','linewidth',3)
legend([h1 h2],{'Strategy','Internal Model'})

