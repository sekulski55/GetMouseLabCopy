%This simulator was inspired by a model of use-dependent learning developed
%by Joern Diedrichsen (see Diedrichsen et al. 2010). I write "inspired by"
%only because I've made some adjustments in order to apply it to a
%visuomotor rotation instead of a force field paradigm and so erroneous
%assumptions and implementations  (all my own) may be present. Most salient
%part of model is that it incorporates PARALLEL operation of two different 
%learning processes, adaptation and UDL. 

%Hyosub Kim (9/23/2018)

clear all; close all

A = .9; B = .2; %parameter values for state-space model of adaptation
E = .95; F = .011;  %parameter values for use-dependent learnig process
D = .5; %in original model, this parameter related to inverse stiffness of
        %arm; here, it represents weight on adaptation process 
w0 = 0;    %this can be thought of as intrinsic bias; allows for "slow
            %drift" back to baseline direction
[w(1) v(1) y(1) e(1)] = deal(0);

nBL = 50; nAcq = 50; nTrans = 200;
nTrials = nBL + nAcq + nTrans;

%let's make error a visuomotor rotation
rot = [zeros(nBL,1); ones(nAcq,1)*30; zeros(nTrans,1)];

for i=1:nTrials-1
    
    %actual movement directions
    y(i) = D*(v(i)) + w(i); %eqn 7
    
    e(i) = rot(i) - v(i); %SPE driving adaptation; in original model, error
                          %is (force field - state estimate), but here the 
                          %error should be related to the actual visual
                          %feedback (so perhaps switch v(i) with y(i)??)
    
    %state estimate
    v(i+1) = A*v(i) + B*(e(i)); %eqn 6

    %use-dependent learning
    w(i+1) = E*w(i) + F*y(i) + (1-E-F)*w0; %eqn 5
    
end
figure; hold on
plot(v,'b','linewidth',5)
plot(w,'g','linewidth',5)
plot(y,'r','linewidth',5)
plot(rot,'--k','linewidth',2)
xlabel('Trial number')
ylabel('Hand angle (deg)')
title('Use-dependent simulation','fontsize',14)
legend('IM','UDL','Actual movement','Perturbation')







