
%% MAYSAM SOLTANIAN
%%NO LOAD CONDITION

%% cleaning up

clc;
clear;
close all;

%% experimental data (v-i)
% Vph=109.5;
% I=1.368;
% Wr=2999.9999;
%Vph=[196.5 180.3 171.3 162.3 153.2 144.2 135.2 126.2]/sqrt(3);
Vp = [400 350 300 250 200 150];

Icm = [31.3 29.0 24.9 21.1 16.5 13.6];  %Measured Current

Wm = [400 450 500 550 600 650]*(2*pi/60);
Wp = [2*pi*50 2*pi*50 2*pi*50 2*pi*50 2*pi*50 2*pi*50];
Pp = 2;
Wpr =((Wp)-(Pp*Wm)); %Formular of getting the Speed of rotor

%I=[2.21 1.8 1.63 1.51 1.45 1.4 1.3 1.26]/sqrt(2);
%Wr=[2984 2982 2980 2978 2976 2973 2970 2965];

%% problem definitions
% f=50;
% Pc =4;
Rr=0.1023;
dataNumber=numel(Vp);
% number of variables
variableNum=2;

%% algorithm initializing and definitions

% learning factors 
C1=2; C2=2;
% minimom of inertia weight 
W_MIN=0.4;
% maximum of inertia weight 
W_MAX=0.9;
% initial inertia weight 
W=W_MAX;      
POPULATION_NUMBER=500;
% maximum number of iteration
MAXIMUM_ITERATION=1000;
% initial value of swarm best solution
swarm_BestSolution=inf;
% preallocating variables
particle_Velocity=struct;
particle_Position=struct;
Inn=zeros(1,dataNumber);
particle_Solution=zeros(POPULATION_NUMBER,1);
particle_BestSolution=zeros(POPULATION_NUMBER,1);
particle_BestPosition=struct;
A=zeros(variableNum,1);
Var=zeros(variableNum,1);
% upper and lower boundaries of each parameter
VARIABLE_MAX= [200 1];
VARIABLE_MIN= [0 0];

%% algorithm's initializing loop

for j=1:POPULATION_NUMBER
    % initializing positions
    for k=1:variableNum
    particle_Position(j).variable(k)=unifrnd(VARIABLE_MIN(1,k),VARIABLE_MAX(1,k));
    end
    % initializing velocities
    for k=1:variableNum
    particle_Velocity(j).variable(k)=0;
    end
    
    % evaluating objective function
    for k=1:variableNum
    A(k)=particle_Position(j).variable(k);
    end
      
%     S=((((120*f)/2)-Wr)/((120*f)/2));
    
    for k=1:dataNumber
        
        Icc(1, k) = (((A(1).*Rr.*Vp(1,k))./(Wp(1,k)*Wpr(1,k)))) + ((((A(1).*A(2).*Vp(1,k))./Wp(1,k))).*1i);
        
    end
    
    Ic = abs(Icc);        
    z = sum(abs((Icm.^2)-(Ic.^2))); %objevtive function 1
    %z=sqrt(mean((Icm-In).^2)); 
    
    particle_Solution(j)=z;
    
    % initializing particle's best solution
    particle_BestSolution(j)=particle_Solution(j);
    % initializing particle's best positions
    for k=1:variableNum
    particle_BestPosition(j).variable(k)=particle_Position(j).variable(k);
    end
    % updating swarm's best solution & swarm's best positions
    if particle_BestSolution(j)<swarm_BestSolution
        swarm_BestSolution=particle_BestSolution(j);
        for k=1:variableNum
        swarm_BestPosition.variable(k)=particle_BestPosition(j).variable(k);
        end
    end 
end

%% main loop
for KK=1:MAXIMUM_ITERATION
    
    for j=1:POPULATION_NUMBER
        for k=1:variableNum
        % updating Velocities
        particle_Velocity(j).variable(k)=W*particle_Velocity(j).variable(k)+C1*rand.*(particle_BestPosition(j).variable(k)-particle_Position(j).variable(k))+C2*rand.*(swarm_BestPosition.variable(k)-particle_Position(j).variable(k));
        end
        % updating Positions
        for k=1:variableNum
        Var(k)=particle_Position(j).variable(k)+particle_Velocity(j).variable(k);
        if Var(k)<VARIABLE_MAX(1,k) &&Var(k)>VARIABLE_MIN(1,k)
                particle_Position(j).variable(k)=Var(k);
        end
        end
        
        % evaluating objective function
        for k=1:variableNum
          A(k)=particle_Position(j).variable(k);
        end

        for k=1:dataNumber
            
         Icc(1,k) =  ((((A(1).*Rr.*Vp(1,k))./(Wp(1,k)*Wpr(1,k))))) + ((((A(1).*A(2).*Vp(1,k))./Wp(1,k))).*1i);  
        
        end
        
        Ic=abs(Icc);
        %z=sqrt(mean((Icm-In).^2));
        z = sum(abs((Icm.^2)-(Ic.^2))); %objective function 1
        
        particle_Solution(j)=z;
        
        % updating particle's best positions & particle's best solution 
        if particle_Solution(j)<particle_BestSolution(j)
            particle_BestSolution(j)=particle_Solution(j);
            for k=1:variableNum
                particle_BestPosition(j).variable(k)=particle_Position(j).variable(k);
            end    
            % Updating swarm's best solution and swarm's best positions
            if particle_BestSolution(j)<swarm_BestSolution
               swarm_BestSolution=particle_BestSolution(j);
               for k=1:variableNum
                   swarm_BestPosition.variable(k)=particle_BestPosition(j).variable(k);
               end
            end 
        end   
    end
    % updating inertia weight
    W=W_MAX-(W_MAX-W_MIN)*(KK/MAXIMUM_ITERATION);
end

%% monitoring

% printing answers
disp(['Minimum error is = ' num2str(swarm_BestSolution)])
%for k=1:variableNum
%disp(['The best solution for Variable number ' num2str(k)  ' is = ' num2str(swarm_BestPosition.variable(k))])
%end
disp(['M = ' num2str(swarm_BestPosition.variable(1))])
disp(['Lr = ' num2str(swarm_BestPosition.variable(2))])
% disp(['Rr = ' num2str(swarm_BestPosition.variable(3))])
% disp(['Lr = ' num2str(swarm_BestPosition.variable(4))])
% disp(['Lm = ' num2str(swarm_BestPosition.variable(5))])


%% the end

