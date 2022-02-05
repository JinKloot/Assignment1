%% Assignment 1 - Simulation 3 - 2
% Enhancements - Injection - 
%  Jinseng Vanderkloot 101031534
%% Initialization of individual electron values

clc
clear all
close all 

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
%% Electrons position and velocity arrays
numElec = 20;                   %Number of simulated Electrons 
numEPlot = 20;                  %Number of plotted Electrons 
dt = (lArea*wArea);             %Typically 1/100 of region size
stepsTot = 200;                 %Total amount of steps (1000 was a long simulation) 
tTot= stepsTot*dt;              %Total Simulation time 
x = zeros(1,numElec);           %Inital X matrix          
y = zeros(1,numElec);           %Inital y matrix  
vx = zeros(1,numElec);          %Inital velocity x matrix  
vy = zeros(1,numElec);          %Inital velocity y matrix
vtot = zeros(1,numElec);        %Inital velocity matrix
avgTemp=0;                      %Set average Temp to 0 

%Electron color assignment
for cnt = 1:numElec
    colors= rand(numElec,3);
    y(cnt) = 50e-9;
end

%Probability of Scatter 
scatOn = 1;                     %Turn Scatter on (1) or off(0)
Pscatter = 1-exp(-dt/tmn);      %Scatter Equation 
tScatter = zeros(1,numElec);
%% Bottle Neck Boundary 
% X limit (no particle between X1 and X2 - Cosider Y limits) 
boxX1=80e-9;
boxX2=120e-9;
% Y Limit (no particles between 0 and Y1 and Y2 to Y limit) 
boxY1 = 40e-9;
boxY2 = 60e-9;

%New Box
boxX3 = 140e-9;
boxX4 = 160e-9;
boxY3 = 40e-9;
boxY4 = 60e-9;

%Boundary Energy/Velocity loss coefficient, when hitting wall, increase
%velocity = increase temp, decrease velocity = decrease temp
vloss = 1.1; 
%% Main loop 
t=0;
intCNT = 1; %Counter with time
while t < tTot
    t = t + dt; 
    
    %Store old position 
    oldx = x;
    oldy = y;
    
    %Inject one elctron for each loop iteration
    if intCNT <= numElec
        %Add a velocity for injected electron 
        vx(intCNT)=sqrt(vt^2)*abs(randn());  % velocity * Gaussian dist   
        vy(intCNT)=sqrt(vt^2)*randn();  % velocity * Gaussian dist 
        vtot(intCNT)= sqrt (vx(intCNT)^2)+(vy(intCNT)^2);
    end 

    %Update to new position 
    x(1:numElec) = x(1:numElec) + (vx(1:numElec).*dt);
    y(1:numElec) = y(1:numElec) + (vy(1:numElec).*dt);  
    vtot(1:numElec)= sqrt ((vx(1:numElec).^2)+(vy(1:numElec).^2));
   
    for check = 1:numElec
        %Scatter 
        if scatOn==1
            if Pscatter > rand()
                if vx ~= 0
                vx(check)=sqrt(vt^2 /2)*randn();
                vy(check)=sqrt(vt^2 /2)*randn();
                end
            end
        end

        %Apply Boundary Conditions
        %If bottom contact, bounce off in opposite direction
        if (y(check)<=0)
            y(check) = 0;
            vy(check) = -vy(check);
        end
        %If top contact, bounce off in opposite directio
        if (y(check)>=lArea)
            y(check) = lArea;
            vy(check) = -vy(check);
        end
        %if left side of box, come out right side 
        if(x(check)<=0)
           x(check) = 0;
           vx(check) = -vx(check);
        end
        %if right side of box, come out left side
        if(x(check)>=wArea)
           x(check) = wArea;
           vx(check) = -vx(check);
        end

        %Apply bottle neck conditions 
        %If contact on left walls of boundary (not in Gap)
        if (oldx(check)<boxX1 && x(check)>=boxX1 && (y(check)<= boxY1 || y(check)>= boxY2))
            x(check)=boxX1;
            vx(check) = -(vx(check)*vloss); 
        end
        %If contact on right walls of boundary (not in Gap)
        if (oldx(check)>boxX2 && x(check)<=boxX2 && (y(check)<= boxY1 || y(check)>= boxY2))
            x(check)=boxX2;
            vx(check) = -(vx(check)*vloss); 
        end
        %If contact with bottom boundary in gap
        if (x(check)>boxX1 && x(check)< boxX2 && oldy(check)>boxY1 && y(check)<= boxY1)
            y(check)= boxY1;
            vy(check) = -(vy(check)*vloss); 
        end
        %If contact with top boundary in gap 
        if (x(check)>boxX1 && x(check)< boxX2 && oldy(check)<boxY2 && y(check)>=boxY2)
            y(check)=boxY2;
            vy(check) = -(vy(check)*vloss); 
        end

        %Apply bottle neck conditions for new box 
        if (oldx(check)<boxX3 && x(check)>=boxX3 && (y(check)<= boxY4 || y(check)>= boxY3))
            x(check)=boxX3;
            vx(check) = -(vx(check)*vloss); 
        end
        %If contact on right walls of boundary (not in Gap)
        if (oldx(check)>boxX4 && x(check)<=boxX4 && (y(check)<= boxY4 || y(check)>= boxY3))
            x(check)=boxX4;
            vx(check) = -(vx(check)*vloss); 
        end
        %If contact with bottom boundary in gap
        if (x(check)>boxX3 && x(check)< boxX4 && oldy(check)>boxY4 && y(check)<= boxY4)
            y(check)= boxY4;
            vy(check) = -(vy(check)*vloss); 
        end
        %If contact with top boundary in gap 
        if (x(check)>boxX3 && x(check)< boxX4 && oldy(check)<boxY3 && y(check)>=boxY3)
            y(check)=boxY3;
            vy(check) = -(vy(check)*vloss); 
        end
    end 

    %Plot Boundary and map some electrons
    for Eplot = 1:numEPlot
        figure(1)
        %if the electron went out of sides and back on other side, do not
        %draw line
        if abs(oldx(Eplot)-x(Eplot))<(wArea/2)
            p = plot([oldx(Eplot),x(Eplot)],[oldy(Eplot),y(Eplot)]);
        end
        rectangle('Position',[boxX1 0 (boxX2-boxX1) boxY1],'FaceColor',[0 0 0])
        rectangle('Position',[boxX1 boxY2 (boxX2-boxX1) (lArea-boxY2)],'FaceColor',[0 0 0])
        rectangle('Position',[boxX3 boxY3 (boxX4-boxX3) (boxY4-boxY3)],'FaceColor',[0 0 0])
        p.Color=colors(Eplot,:);
        axis([0,wArea,0,lArea]);
        ptitle=['Electron Model steps (',intCNT, '/',stepsTot];
        title(ptitle), xlabel('Position (m)', 'FontSize', 10), ylabel('Position (m)', 'FontSize', 10);
        hold on;
    end 
    pause(0.01);    
    intCNT = intCNT + 1; 
 
end 


