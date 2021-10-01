x = 0:0.2513:4*pi;
Sin_rotation = 4*sin(x);
plot(x,10*sin(x)), grid on
    %Create sine wave rotation
    x = 0:0.2513:4*pi;
    Sin_rotation = sin(x);
    %rotation(trial,1)= rotation(trial,1)*Sin_rotation(trial);

cd('TargetFiles')
tgt_file = dlmread(['demo','1','.tgt'], '\t', 1, 0); % start reading in from 2nd row, 1st column
%rotation = ones(50,1)
rotation = tgt_file(:,4);
array = ones(50,1);
trial = 1;
i=11;
while i<=50;
    rotation(i)= rotation(i)*Sin_rotation(i-10);
    i=i+1;
end
x = 0:0.2513:4*pi;
Sin_rotation = sin(x);

x2 = 0:0.2513:4*pi;    %makes a plot line for the wave with 50 (maxtrial) intervals
Sin_rotation2 = -sin(x2)+1;

plot(Sin_rotation); hold on; plot(Sin_rotation2);