function x = dualRateSSM_simulator(Af,Bf,As,Bs)

%%% This function was written to simulate the dual-rate state-space model 
%%% of Smith et al (2006).  
%%% Copy next line into command window to see example of this model 
%%% learn a two different force fields followed by an error clamp block, as
%%% in Fig. 2 from the paper:
%%% dualRateSSM_simulator(.92,.03,.996,.004)

%initialize internal model
nTrials = 800;
[x,xf,xs] = deal(zeros(nTrials,1));

%perturbation schedule
f=[zeros(20,1);ones(400,1);ones(15,1)*-1;zeros(365,1)];


%plot figure
figure;hold on
xlabel('Trial number','fontsize',14)
ylabel('Adaptation','fontsize',14)
title('Multi-rate model','fontsize',14)
xlim([0 sum(nTrials)+2])
ylim([-1.2 1.2])
plot(1:length(x),repmat(0,length(x),1),'linewidth',1,'color','k')
h1 = plot(1:length(f),f,'linewidth',2,'color','k')

for i=1:sum(nTrials)-1    
    e(i) = f(i)-x(i);
    
    %internal model
    xf(i+1) = Af*xf(i)+Bf*e(i);
    xs(i+1) = As*xs(i)+Bs*e(i);
    x(i+1) = xf(i+1)+xs(i+1);
      
    h2 = plot(xf(1:i),'--g','linewidth',3);
    h3 = plot(xs(1:i),'--b','linewidth',3);
    h4 = plot(x(1:i),'r','linewidth',4);
    
    %drawnow 
    pause(0.001)
end
legend([h1 h2 h3 h4],'Perturbation','Fast state', 'Slow state', 'Net adaptation')
legend('boxoff')
