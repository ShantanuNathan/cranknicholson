%% This program plots the graph as soon as the wavefunction is calculated.
% Using Natural Units
hbar = 1;
m = 1;
%====================================================================
% Parameters for solving the problem in the interval 0 < x < L
i = sqrt(-1);
L = 200; % Interval Length
N = 1600; % No of points
x = linspace(-L,L,N)'; % Coordinate vector
dx = x(3) - x(2); % Coordinate step
x0 = 10;
%======================================================================
%propagation
NT = 1000; % No. of time steps
TF = 60; 
T = linspace(0,TF,NT); % Time Uector
dt = T(2)-T(1); % Time step
%======================================================================
%parameters for initlizing wave function
ko  = 2;
sig = 1 ;
%initlizing wave function in position space
psi = (exp(-i*ko*x)).*(exp(-((x-x0).^2)/(4*sig^2)));
psi = psi/(2*pi*(sig^2))^(0.25);
%=====================================================================
%Potential Barrier
avgE = hbar*ko*ko/(2*m);
U0 = avgE;
U = zeros([1 N]);
% Parameters for computing psi(x,T) at different times 0 < T < TF
n1 = 900; %%initially n1:n2 = 900:910
n2 = 910;
U(n1:n2) = U0; %Potential Initialized
%=====================================================================
%defining crank_nicholsan and hamiltonian matrix
z = ((i*hbar*dt)/(2*m*(dx^2)));
A = zeros(N);
B = zeros(N);

% v_group = h*ko/m;

%loop to find 



for m=1:N
    for n=1:N
        if(m==n)
            A(m,n) = 2+2*z+(i*dt*U(n))/hbar ;
            B(m,n) = 2-2*z-(i*dt*U(n))/hbar ;
        else if ( m == (n+1) || n == (m+1) )
                A(m,n) = -z;
                B(m,n) = z;
            end
        end
    end
end

%%
t=1;
C =(inv(A))*B;  %Crank Nicolson Matrix

%Plotting the graph
title({'Quantum Tunneling'},'FontSize',17);
xlabel('Postion (x)','FontSize',17);
ylabel('Probability Density |\phi|^{2}','FontSize',17);
set(gcf, 'Position', get(0, 'Screensize')); %to plot the graph in Full screen
axis([-80,150,0,0.3]);

U_scaled = U/8; % Scaling Potential down to have better visual

wf = animatedline('Color','b'); %Line for the wavefunction
pot = animatedline('Color','r'); %Line for Potential Barrier

addpoints(pot,x,U_scaled);

dim1 = [.8 .5 .3 .3];
dim2 = [.2 .5 .3 .3];

sT = 'T = ';
sR = 'R = ';

probT_max = 0;
probU_max = 0;

max_R = 0;
max_T = 0;

t_1 = 0;
t_2 = 0;

psi_1 = 0;
psi_2 = 0;
psi_3 = 0;

while(t<=NT)   
%     probT = 0;
%     probR = 0;
%     probU = 0;

 
    clearpoints(wf);
   
    addpoints(wf,x,abs(power(psi,2)));
    drawnow;
    psi = C'*psi;
        
%     for k=1:length(psi)
%         if(k>=n2)
%         probT = probT + abs(power(psi(k,1),2))*dx;
%         end
%         if(k<=n1)
%         probR = probR + abs(power(psi(k,1),2))*dx;
%         end
%         if(k>n1 && k<n2)
%             probU = probU + abs(power(psi(k,1),2))*dx;
%         end
%             
%     end
%     
%     delete(findall(gcf,'type','annotation'));
%     
%     annotation('textbox',dim1,'String',strcat(sT,num2str(probT)),'FitBoxToText','on');
%     annotation('textbox',dim2,'String',strcat(sR,num2str(probR)),'FitBoxToText','on');
%    
%     pause(0.5e-02);
    
%   Measuring time of Propogation
    if(abs(power(psi(n2,1),2)) == max(abs(power(psi(n2:length(psi)),2))))
        t_1 = t;
    end
    
    %To take snapshot
%     if(t == 800)
%         psi_3 = psi;
%         break;
%     end
    
    t=t+1;
end

