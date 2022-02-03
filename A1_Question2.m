%%Assignment 1 
% Jinseng Vanderkloot 101031534

clc
clear all
close all 

%Initial
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

% Mean free path (Velocity * minimum time between collision) (Question 1.B)
meanFreePath = vt * tmn;

%Electron motion 
numElec = 1000;                 %Number of simulated Electrons 
numEPlot = 10;                  %Number of plotted Electrons 
dt = (lArea*wArea);             %Typically 1/100 of region size
stepsTot = 200;                 %Total amount of steps (1000 was a long simulation) 
tTot= stepsTot*dt;              %Total Simulation time 
x = zeros(1,numElec);           %Inital X matrix          
y = zeros(1,numElec);           %Inital y matrix  
vx = zeros(1,numElec);          %Inital velocity x matrix  
vy = zeros(1,numElec);          %Inital velocity y matrix
vtot = zeros(1,numElec);        %Inital velocity matrix
colors = rand(numElec,3);       %Color assignment for each electro
avgTemp=0;                      %Set average Temp to 0 

%Probability of Scatter 
Pscatter = 1-exp(-dt/tmn);      %Scatter Equation 
tScatter = zeros(1,numElec);    %track scatter for each particle 

%Electron Graph 
for cnt = 1:numElec
    x(cnt)=rand()*wArea;
    y(cnt)=rand()*lArea;
    vx(cnt)=sqrt(vt^2)*randn();  % velocity * Gaussian dist   
    vy(cnt)=sqrt(vt^2)*randn();  % velocity * Gaussian dist 
    %Varience = sqrt(kT/m) - Do we use this? 
    vtot(cnt)= sqrt (vx(cnt)^2)+(vy(cnt)^2);
end

%Main Loop
t=0;
intCNT = 1; %Counter with time
while t < tTot
    t = t + dt; 
    
    %Store old position 
    oldx = x;
    oldy = y;
    
    %Update to new position 
    x(1:numElec) = x(1:numElec) + (vx(1:numElec).*dt);
    y(1:numElec) = y(1:numElec) + (vy(1:numElec).*dt);
    
    
    %Check each electron for scratter 
    for scat = 1:numElec
        if Pscatter > rand()
            vx(scat)=sqrt(vt^2 /2)*randn();  
            vy(scat)=sqrt(vt^2 /2)*randn();
            tScatter(scat)= 0; %If collision, time goes to 0
        else
            tScatter(scat)= tScatter(scat) + dt; %track time increaing while no collision
        end
    end 

    vtot(1:numElec)= sqrt ((vx(1:numElec).^2)+(vy(1:numElec).^2));
    
    %Apply boundary conditions 
    for check = 1:numElec
        %If top or bottom contact, bounce off in opposite direction
        if (y(check)<=0 || y(check)>=lArea)
            vy(check) = -vy(check);
        end
        %if left side of box, come out right side 
        if(x(check)<=0)
           x(check) = x(check) + wArea;
        end
        %if right side of box, come out left side
        if(x(check)>=wArea)
          x(check) = x(check) - wArea;
        end
    end 

    %Plot Boundary and map some electrons
    for Eplot = 1:numEPlot
        figure(1)
        subplot (3,1,1)
        %if the electron went out of sides and back on other side, do not
        %draw line
        if abs(oldx(Eplot)-x(Eplot))<(wArea/2)
            p = plot([oldx(Eplot),x(Eplot)],[oldy(Eplot),y(Eplot)]);
        end
        p.Color=colors(Eplot,:);
        axis([0,wArea,0,lArea]);
        title('Electron Model'), xlabel('Position (m)', 'FontSize', 10), ylabel('Position (m)', 'FontSize', 10);
        
        hold on;
    end 
    pause(0.01);
    
    %Calc Average Temp for all t and Plot 
    subplot (3,1,2)
    Time(:,intCNT) = t;
    allT = ((vtot(:).^2).*mn)./(2*kb); %since vector is 1D vtot, do not device by 2*kb
    avgTemp(:,intCNT) = mean(allT);
     
    plot(Time,avgTemp,"r");
    title('Averge Temp'),xlabel('Time (s)', 'FontSize', 10), ylabel('Temp (K)', 'FontSize', 10), ylim([250,600]); 
    hold on;
    intCNT = intCNT +1; 
    
    %Histogram of velocities over time 
    subplot(3,3,7)
    histogram([vtot(:)],30)
    title('Velocity Histgram'),xlabel('Velocity (m/s)', 'FontSize', 10), ylabel('Number of Particles', 'FontSize', 10);

    %Mean time between collision  
    Time(:,intCNT) = t;
    allScat(:,intCNT) = mean(tScatter(:));
    subplot(3,3,8)
    plot(Time,allScat,'r');
    title('Mean Time between Collision'),xlabel('Time (s)', 'FontSize', 10), ylabel('Time(s)', 'FontSize', 10);
    hold on;
   
    %Mean Free Path over time 
    Time(:,intCNT) = t;
    mfp(:,intCNT) = mean(tScatter(:))*mean(vtot(:));
    subplot(3,3,9)
    plot(Time,mfp,'r');
    title('Mean Free Path '),xlabel('Time (s)', 'FontSize', 10), ylabel('Time(s)', 'FontSize', 10);
    hold on;
    
end 