%%Assignment 1 
% Jinseng Vanderkloot 101031534

clc
clear all
close all 

%Initial (Question 1) 
m0 = 9.10938215e-31;            % electron mass
mn = 0.26*m0;                   % Effective mass
Temp = 300;                     % Inital Temp (K)
kb = 1.3806504e-23;             % Boltzmann constant
tmn = 0.2e-12;                  % Mean time between collision 

% Region Area 
wArea = 200e-9;
lArea = 100e-9;

%Thermal Velocity (Question 1.A) 
vt=sqrt((2*kb*Temp)/mn);        % Sim in 2D so (2*kb*Temp), 3D is (3*kb*Temp)
fprintf("Thermal Velocity = %d m/s \n", vt);

% Mean free path (Velocity * minimum time between collision) (Question 1.B)
meanFreePath = vt * tmn;
fprintf("Mean free path = %d \n", meanFreePath);

%Electron motion (Question 1.C) 
numElec = 1000; 
numEPlot = 10; 
dt = (lArea*wArea); %Typically 1/100 of region size
stepsTot = 1000; 
tTot= stepsTot*dt;
x = zeros(1,numElec);
y = zeros(1,numElec);
vx = zeros(1,numElec);
vy = zeros(1,numElec);
vtot = zeros(1,numElec);
colors = rand(numElec,3);
eTemp=0;

%Electron Graph 
for cnt = 1:numElec
    x(cnt)=rand()*wArea;
    y(cnt)=rand()*lArea;
    angle = (2*pi*rand());
    vx(cnt)=sqrt(vt^2 /2)*cos(angle); %(randi([-1,1])); % velocity * random direction   
    vy(cnt)=sqrt(vt^2 /2)*sin(angle); %(randi([-1,1])); % velocity * random direction
    vtot(cnt)= sqrt (vx(cnt)^2)+(vy(cnt)^2);
    %elec(cnt) = [x,y,vx,vy,vtot];
end

t=0;
intCNT = 2;
while t < tTot
    t = t + dt; 
    
    oldx = x;
    oldy = y;

    x(1:numElec) = x(1:numElec) + (vx(1:numElec).*dt);
    y(1:numElec) = y(1:numElec) + (vy(1:numElec).*dt);
    
    vtot(1:numElec)= sqrt ((vx(1:numElec).^2)+(vy(1:numElec).^2));

    for check = 1:numElec
        if (y(check)<=0 || y(check)>=lArea)
            vy(check) = -vy(check);
         end
        if(x(check)<=0)
           x(check) = x(check) + wArea;
        end
        if(x(check)>=wArea)
          x(check) = x(check) - wArea;
        end
    end 

    for Eplot = 1:numEPlot
        subplot (2,1,1)
        if abs(oldx(Eplot)-x(Eplot))<(wArea/2)
            p = plot([oldx(Eplot),x(Eplot)],[oldy(Eplot),y(Eplot)]);
        end
        p.Color=colors(Eplot,:);
        axis([0,wArea,0,lArea]);
        title('Electron Model'), xlabel('Position (m)', 'FontSize', 10), ylabel('Position (m)', 'FontSize', 10);
        
        hold on;
    end 
    pause(0.01);
    
    subplot (2,1,2)
    Time(:,intCNT) = t;
    allT = ((vtot(:).^2).*mn)./(kb); %since vector is 1D vtot, do not device by 2*kb
    eTemp(:,intCNT) = mean(allT);
     
    plot(Time,eTemp);
    title('Averge Temp'),xlabel('Time (s)', 'FontSize', 10), ylabel('Temp (K)', 'FontSize', 10), ylim([299,301]); 
    hold on;
    intCNT = intCNT +1; 

end 









