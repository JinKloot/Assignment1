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
numElec = 400;                 %Number of simulated Electrons 
numEPlot = 40;                  %Number of plotted Electrons 
dt = (lArea*wArea)/2;             %Typically 1/100 of region size
stepsTot = 200;                 %Total amount of steps (1000 was a long simulation) 
tTot= stepsTot*dt;              %Total Simulation time 
x = zeros(1,numElec);           %Inital X matrix          
y = zeros(1,numElec);           %Inital y matrix  
vx = zeros(1,numElec);          %Inital velocity x matrix  
vy = zeros(1,numElec);          %Inital velocity y matrix
vtot = zeros(1,numElec);        %Inital velocity matrix
avgTemp=0;                      %Set average Temp to 0 

%Probability of Scatter 
scatOn = 0;                     %Turn Scatter on (1) or off(0)
Pscatter = 1-exp(-dt/tmn);      %Scatter Equation 
tScatter = zeros(1,numElec);

%Bottle Neck Boundary reduced to prevent bleeding 
% X limit (no particle between X1 and X2 - Cosider Y limits) 
boxX1=80e-9;
boxX2=120e-9;
% Y Limit (no particles between 0 and Y1 and Y2 to Y limit) 
boxY1 = 40e-9;
boxY2 = 60e-9;


%Electron Graph 
for cnt = 1:numElec
    x(cnt)=rand()*wArea;
    y(cnt)=rand()*lArea;
    while (x(cnt)>=boxX1 && x(cnt)<=boxX2 && (y(cnt)<=boxY1 || y(cnt>=boxY2))) %Relocate them if in boundary
        x(cnt)=rand()*wArea;
        y(cnt)=rand()*lArea;
    end
    vx(cnt)=sqrt(vt^2)*randn();  % velocity * Gaussian dist   
    vy(cnt)=sqrt(vt^2)*randn();  % velocity * Gaussian dist 
    %Varience = sqrt(kT/m) - Do we use this? 
    vtot(cnt)= sqrt (vx(cnt)^2)+(vy(cnt)^2);
    colors= rand(numElec,3);
   
end


%Boundary Energy/Velocity loss coefficient
vloss = 0.99; 

% %Change location of particle if places inside of bottle neck 
%  for ctt = 1:numElec
%     if (x(ctt) >= boxX1 && x(ctt) <= boxX2)
%         if (y(ctt)<= boxY1)
%             disDir = randi([1 3]); %Randomize which way to displace particle in boundary
%             if disDir == 1
%                 x(ctt) = randn()*boxX1; % relocate electron on left side of bottleneck
%             elseif disDir == 2
%                 y(ctt) = boxY1 + (randn()*(boxY2-boxY1));
%             else
%                 x(ctt) = boxX2 + (randn()*(wArea - boxX2));
%             end
%         end 
%         if (y(ctt)>= boxY2)
%             disDir = randi([1 3]); %Randomize which way to displace particle in boundary
%             if disDir == 1
%                 x(ctt) = randn()*boxX1; % relocate electron on left side of bottleneck
%             elseif disDir == 2
%                 y(ctt) = boxY1 + (randn()*(boxY2-boxY1));
%             else
%                 x(ctt) = boxX2 + (randn()*(wArea - boxX2));
%             end
%         end
%     end 
% end


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
%     for (scat = 1:numElec)
%         if Pscatter > rand()
%             vx(scat)=sqrt(vt^2 /2)*randn();  
%             vy(scat)=sqrt(vt^2 /2)*randn();
%             tScatter(scat)= 0; %If collision, time goes to 0
%         else
%             tScatter(scat)= tScatter(scat) + dt; %track time increaing while no collision
%         end
%     end 

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

    %Apply bottle neck conditions 
    for check2 = 1:numElec
        %If contact on left walls of boundary (not in Gap)
        if (oldx(check2)<boxX1 && x(check2)>=boxX1 && (y(check2)<= boxY1 || y(check2)>= boxY2))
            vx(check2) = -(vx(check2)*vloss); 
        end
        %If contact on right walls of boundary (not in Gap)
        if (oldx(check2)>boxX2 && x(check2)<=boxX2 && (y(check2)<= boxY1 || y(check2)>= boxY2))
            vx(check2) = -(vx(check2)*vloss); 
        end
        %If contact with bottom boundary in gap
        if (x(check2)>boxX1 && x(check2)< boxX2 && oldy(check2)>boxY1 && y(check2)<= boxY1)
            vy(check2) = -(vy(check2)*vloss); 
        end
        %If contact with top boundary in gap 
        if (x(check2)>boxX1 && x(check2)< boxX2 && oldy(check2)<boxY2 && y(check2)>=boxY2)
            vy(check2) = -(vy(check2)*vloss); 
        end
    end 


    %Plot Boundary and map some electrons
    for Eplot = 1:numEPlot
        figure(1)
        subplot (2,1,1)
        %if the electron went out of sides and back on other side, do not
        %draw line
        if abs(oldx(Eplot)-x(Eplot))<(wArea/2)
            p = plot([oldx(Eplot),x(Eplot)],[oldy(Eplot),y(Eplot)]);
        end
        rectangle('Position',[boxX1 0 (boxX2-boxX1) boxY1])
        rectangle('Position',[boxX1 boxY2 (boxX2-boxX1) (lArea-boxY2)])
        p.Color=colors(Eplot,:);
        axis([0,wArea,0,lArea]);
        title('Electron Model'), xlabel('Position (m)', 'FontSize', 10), ylabel('Position (m)', 'FontSize', 10);
        
        hold on;
    end 
    pause(0.01);

    %Electron Density Map 
    eMapX=linspace(0, wArea, 100);
    eMapY=linspace(0, lArea, 50);
    EDM=histcounts2(y,x,eMapY,eMapX);
    subplot(2,2,3)
    imagesc(eMapY,eMapX,EDM),colorbar,title('Electron Density Map');
    
    %Electron Temprature Map 
    allT = ((vtot(:).^2).*mn)./(2*kb);
    xv = linspace(min(x), max(x),100);
    yv = linspace(min(y), max(y),50);
    [X,Y] = meshgrid(xv,yv);
    ETM=griddata(x,y,allT,X,Y);
    subplot(2,2,4);
    imagesc(xv,yv,ETM),colorbar,title("Temprature Map")
    axis([0, wArea, 0 lArea]);
    
    
end 


