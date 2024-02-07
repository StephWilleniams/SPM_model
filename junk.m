% Matlab generic junk file

%% Load the files

%filename = 'output_burger1.txt';
%filename = 'output_burger1.txt';
%%filename = ['output_firstOrder.txt'];
close all

for i = 1
filename = ['example1/test_' num2str(i) '.txt'];
A = readmatrix(filename);

plot(A(1,:))
hold on
end

%% Make a video 

filename = ['example1/test_' num2str(1) '.txt'];
A = readmatrix(filename);
filename = ['example2/test_' num2str(1) '.txt'];
B = readmatrix(filename);
%v = VideoWriter("ExampleComparison.avi");
%open(v)

% Dimensionalisation
ds = 0.01;
dt = 0.01;

% Array sizes
[Nt,Ns] = size(A);

ss = linspace(0,ds*Ns,Ns); % Sizes
tt = linspace(0,dt*Nt,Nt); % Times

n = 10;

% Loop through times
for i = 1:10:1000

    plot(ss(1:n:end),A(i,1:n:end),'b-','LineWidth',5)
    hold on;
    plot(ss(1:n:end),B(i,1:n:end),'r-','LineWidth',5)
    %plot(ss,exp(-(ss-5-(i-1)*dt).^2),'r--','LineWidth',2)
    ylim([-0.2,1.2])
    pause(0.3)

    %frame = getframe(gcf);
    %writeVideo(v,frame);

    hold off

end
