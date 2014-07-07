% This is a stochastic simulation algorithm (ssa) script designed to
% simulate SRP binding to ribosomes under single molecule conditions.
% The theory behind the script can be found in this reference:
% Gillespie, D. T. (2007). Stochastic Simulation of Chemical Kinetics. 
% Annual Review of Physical Chemistry, 58(1), 35-55.

%% Description of the Reactions and their associated rates:

% Ribosomes are ribo.
% Assumes two conformations of SRP, SRP(d) and SRP(p). 
% When the fluorescent label photobleaches SRP becomes SRP(d)0 or SRP(p)0. 

% Interconversion between the two conformations in solution:
% 1     SRP(d) > SRP(p)
% 2     SRP(p) > SRP(d)

% SRP(d) binding and unbinding ribosomes:
% 3     SRP(d) + Ribo > SRP(d)*Ribo
% 4     SRP(d)*Ribo > SRP(d) + Ribo 

% SRP (p) binding and unbinding ribosomes:
% 5     SRP(p) + Ribo > SRP(p)*Ribo
% 6     SRP(p)*Ribo > SRP(p) + Ribo 

% Interconversion between the two SRP conformations on a ribosome:
% 7     SRP(d)*Ribo > SRP(p)*Ribo
% 8     SRP(p)*Ribo > SRP(d)*Ribo

% Dye photobleaching. (Can only happen on the ribosome due to TIRF ilumination)
% 9     SRP(d)*Ribo > SRP(d)0*Ribo
% 10	SRP(p)*Ribo > SRP(p)0*Ribo

% Reactions above repeated with population of photobleached SRP
% 11	SRP(d)0 > SRP(p)0
% 12	SRP(p)0 > SRP(d)0
% 13	SRP(d)0 + Ribo > SRP(d)0*Ribo
% 14	SRP(d)0*Ribo > SRP(d)0 + Ribo 
% 15	SRP(p)0 + Ribo > SRP(p)0*Ribo
% 16	SRP(p)0*Ribo > SRP(p)0 + Ribo 
% 17	SRP(d)0*Ribo > SRP(p)0*Ribo
% 18	SRP(p)0*Ribo > SRP(d)0*Ribo

%% Stoichiometry matrix:

% define stoichiometry matrix: nu=stoichiometry matrix, n=number of reactions, m=number of species
% reflects how molecule counts change when each reaction happens (each row
% is a reaction in the same order described above)
nu=[-1	1	0	0	0	0	0	0	0;
1	-1	0	0	0	0	0	0	0;
-1	0	1	0	-1	0	0	0	0;
1	0	-1	0	1	0	0	0	0;
0	-1	0	1	-1	0	0	0	0;
0	1	0	-1	1	0	0	0	0;
0	0	-1	1	0	0	0	0	0;
0	0	1	-1	0	0	0	0	0;
0	0	-1	0	0	0	0	1	0;
0	0	0	-1	0	0	0	0	1;
0	0	0	0	0	-1	1	0	0;
0	0	0	0	0	1	-1	0	0;
0	0	0	0	-1	-1	0	1	0;
0	0	0	0	1	1	0	-1	0;
0	0	0	0	-1	0	-1	0	1;
0	0	0	0	1	0	1	0	-1;
0	0	0	0	0	0	0	-1	1;
0	0	0	0	0	0	0	1	-1]';

[n,m] = size(nu); 

%% Initial conditions, should reflect experimental conditions:
concSRP = 15e-9; % [SRP] in M
volSlide = 10e-6; % Volume of slide in liters
numRibos = 1000; % Number of ribosomes per slide
volRibos = volSlide/numRibos;
avogadro = 6.022e23;
solSRP = concSRP*volRibos*avogadro; %molecules of SRP in solution
ssaRate2onRate = volRibos*avogadro; %This is conversion factor to make conventional on-Rate into /molecules*sec needed for ssa

% T is final time (seconds)
T=240;

% These are for the convertion to SM:

%Frame rate of single molecule experiment (in miliseconds)
SMframerate = 100; 
SMt = [SMframerate/1000: SMframerate/1000: T]; % SMt is the final time vector

%% THIS IS WHERE PARAMETERS ARE MODIFIED TO TEST BEST MATCH WITH OBSERVED EXPERIMENTAL RESULTS. 

% Define initial quantities of reaction components: 

dist2prox = 1; % ratio of SRP(d) to SRP(p) in solution
dist = solSRP*(dist2prox/(dist2prox+1)); % molecules of SRP(d)
moles_dist = dist/avogadro; %moles of SRP(d)
prox = dist/dist2prox; %molecules of SRP(p)
moles_prox = prox/avogadro; %moles of SRP(p)

x0=[dist;   % x(1) SRP(d)
prox;       % x(2) SRP(p)
0;          % x(3) SRP(d)*Ribo
0;          % x(4) SRP(p)*Ribo
1;          % x(5) Ribo
0;          % x(6) SRP(d)0
0;          % x(7) SRP(p)0
0;          % x(8) SRP(d)0*Ribo
0];         % x(9) SRP(p)0*Ribo

% Define rate constants (in  /sec or /molecules*sec):

% SRP(d) > SRP(p)
k1=0;  % In this simulation round SRP interconversion assumed not to happen

% SRP(p) > SRP(d)
k2=0;  % In this simulation round SRP interconversion assumed not to happen

% SRP(d) + Ribo > SRP(d)*Ribo
k3 = 5e-10; %This is the /molec*sec rate needed for the ssa, in this simulation 10X faster than SRP(p) binding   
k3onRate = k3*ssaRate2onRate;
k3time = 1/(k3onRate*(moles_dist/volRibos));

% SRP(d)*Ribo > SRP(d) + Ribo
time = 10; % In this simulation SRP(d) resides on ribos 10X shorter than SPR(p)
k4=1/time;  

% SRP(p) + Ribo > SRP(p)*Ribo
k5 = 5e-9; %This is the /molec*sec rate needed for the ssa, in this simulation 10X slower than SRP(d) binding   
k5onRate = k5*ssaRate2onRate;
k5time = 1/(k5onRate*(moles_prox/volRibos));    

% SRP(p)*Ribo > SRP(p) + Ribo
time =100; % In this simulation SRP(p) resides on ribos 10X longer than SPR(d)
k6 = 1/time;

% SRP(d)*Ribo > SRP(p)*Ribo
k7=0; % In this simulation SRP interconversion on ribosomes assumed not to happen

% SRP(p)*Ribo > SRP(d)*Ribo
k8=0; % In this simulation SRP interconversion on ribosomes assumed not to happen

% SRP(d)*Ribo > SRP(d)0*Ribo
k9=1/83.3;   % In this simulation SRP label photobleaching on ribosomes assumed to happen equaly to both conformations  

% SRP(p)*Ribo > SRP(p)0*Ribo
k10=1/83.3;  % In this simulation SRP label photobleaching on ribosomes assumed to happen equaly to both conformations    

% photobleached SRP assumed to be the same as non-photobleached in terms of
% binding rates.
k11=k1;         % SRP(d)0 > SRP(p)0
k12=k2;         % SRP(p)0 > SRP(d)0
k13=k3;         % SRP(d)0 + Ribo > SRP(d)0*Ribo
k14=k4;         % SRP(d)0*Ribo > SRP(d)0 + Ribo
k15=k5;         % SRP(p)0 + Ribo > SRP(p)0*Ribo
k16=k6;         % SRP(p)0*Ribo > SRP(p)0 + Ribo
k17=k7;         % SRP(d)0*Ribo > SRP(p)0*Ribo
k18=k8;         % SRP(p)0*Ribo > SRP(d)0*Ribo

% single summary of all reaction rates
k = [k1; k2; k3; k4; k5; k6; k7; k8; k9; k10; k11; k12; k13; k14; k15; k16; k17; k18];

%% The simulation engine:

% how many times simulation is run
runs = 100; 

for i=1:1:runs
    
    x = x0;
    finished = 0;
    t=0;
    tarr = zeros([1 10^5]);
    xarr = zeros([n 10^5]);
    jarr = zeros([1 10^5]);
    tarr(1) = t; xarr(:,1) = x;
    tauarr = [];
    step=0;
    
    while(~finished)
        step = step+1;
        
        % define propensity functions , one for each reaction
        a1=k1*x(1);
        a2=k2*x(2);
        a3=k3*x(1)*x(5);
        a4=k4*x(3);
        a5=k5*x(2)*x(5);
        a6=k6*x(4);
        a7=k7*x(3);
        a8=k8*x(4);
        a9=k9*x(3);
        a10=k10*x(4);
        a11=k11*x(6);
        a12=k12*x(7);
        a13=k13*x(6)*x(5);
        a14=k14*x(8);
        a15=k15*x(7)*x(5);
        a16=k16*x(9);
        a17=k17*x(8);
        a18=k18*x(9);
        
        
        a= [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12; a13; a14; a15; a16; a17; a18];
        
        a0 = sum(a);
        if (a0==0)  % No more reactions.
            finished=1;
            tau = T-t;
            delx=0;
        else
            
            % Time for next reaction.
            r1 = rand; r2 = rand;
            tau = (1/a0)*log(1/r1);
            
            % Check terminal condn
            finished = ((t+tau)>=T);
            if ((t+tau)>T)
                % No more reactions.
                tau = T-t;
                delx=0;
            else
                % Find next reaction.
                f = r2*a0;
                j=1; 
                jsum=a(1);
                found = (a(1)>f);
                
                while((~found) & (j<m))
                    j=j+1;
                    jsum = jsum+a(j);
                    found = (jsum>f);
                end
                delx = nu(:,j);   % One reaction 'j' happened.
            end
        end
        
        t = t+tau;
        x = x+delx;
        tarr(step) = t;
        xarr(:,step) = x;
        jarr(step) = j;
        
    end
    
    % samples at end of time simulation
    tarr = tarr(1:step);
    xarr = xarr(:,1:step);
    jarr = jarr(1:step);
    
    
    % This converts the simulation into an SM trace:
    SMtarr = tarr;
    SMtarr = round(SMtarr*(1000/SMframerate));
    SMxarr = xarr(3,:) + xarr(4,:);
    SMlength = length(SMtarr);
    SMdata = zeros(size(SMt)); % SMdata is the final trace assignmnent
    
    for l = 1:1:SMlength
        
        if l < SMlength
            if SMtarr(l) == 0
                SMdata(1:SMtarr(l+1)-1) = SMxarr(l);
            else
                SMdata(SMtarr(l):SMtarr(l+1)-1) = SMxarr(l);
            end
        else
            SMdata(SMtarr(l)) = SMxarr(l);
        end
        
    end
    
    % Eliminates binding events that haven't finished by the end of the
    % experiment:
    
    if SMdata(length(SMdata)) == 1
       transitions = diff(SMdata);
       transitions_index = find(transitions==1);
       SMdata(max(transitions_index)+1:length(SMdata)) = 0;
    end
    
    % summatry array of each simulation
    summary(i).tarr = tarr;
    summary(i).xarr = xarr;
    summary(i).SMdata = SMdata;
    
    % display simulation number to keep track of where the script is
    i
end

save ssa_RUN
