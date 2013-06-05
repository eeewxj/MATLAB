%NLSEmagic1D
%Program to integrate the one-dimensional Nonlinear Shrödinger Equation:
%i*Ut + a*Uxx - V(x)*U + s*|U|^2*U = 0.
%
%©2011 Ronald M Caplan
%Developed with support from the 
%Computational Science Research Center at 
%San Diego State University.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  VERSION:  012    1/8/2012 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear any previous run variables and close all figures:
close all; clear all;
%------------------Simulation parameters----------------------
endtw        = 50;  %End time of simulation.
numframes    = 50;  %Number of desired frames to plot.
h            = 1/10;%Spatial grid spacing.
k            = 0;   %Time step [Set to 0 to auto-compute smallest stable timestep]
method       = 1;   %Method:  1: CD O(4,2), 2: 2SHOC O(4,4)
cuda         = 0;   %Try to use CUDA code (if installed and compiled)
precision    = 2;   %Single (1) or double (2) precision.
tol          = 10;  %Simulation mod-squared tolerance to detect blowup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%             INITIAL CONDITION PARAMETERS        %%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%
%Bright sech soliton
a    = 1;
s    = 1;   %Focusing/Attractive NL
OM   = 1;   %Frequency
c    = 0.1; %Velocity of soliton
x0   = 0;
%Compute numerical domain to be big enough to avoid BC effects:
xmin = -sqrt(a/OM)*asech(eps*sqrt(s/(2*OM)));
xmax =  sqrt(a/OM)*asech(eps*sqrt(s/(2*OM))) + c*endtw;
BC   = 1;   %Dirchilet BC
%%%%%%%%%                                                 %%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------------------------------------------');
disp('--------------------NLSEmagic1D---------------------');
disp('----------------------------------------------------');

%Adjust ranges so they are in units of h (fixes bugs):
xmin = xmin - rem(xmin,h);
xmax = xmax - rem(xmax,h);

%Set up spatial grid:
xvec = (xmin:h:xmax)';   
xres = length(xvec);

%Set up CUDA info
if(cuda==1)    
    %Check that CUDA is installed on the system:
    cudachk=0;
    str = computer;
    if( strcmp(str, 'PCWIN') || strcmp(str, 'PCWIN64'))    
        if(exist([matlabroot '/bin/nvmex.pl'],'file')~=0)
            cudachk=1;
        end
    else
        if(exist('/dev/nvidia0','file')~=0)
            cudachk=1;
        end
    end
    if(cudachk==1)
        comp_cap = getCudaInfo;
        if(precision==2 && (comp_cap(1)==1 && comp_cap(2)<3))
           precision = 1;
           disp('WARNING:  Precision changed to SINGLE because your CUDA card does not support double precision.');
        end  
        cudablocksize     = 512; 
        sharedmemperblock = (cudablocksize+2)*(3+2*method)*(4*precision)/1000;
        numcudablocks     = (ceil(xres/cudablocksize));
        %For MSD BC, need to check this:
        if(BC==2)
            if(xres - cudablocksize*(numcudablocks-1) == 1)
                disp('MSD CUDA ERROR: N (x) is one cell greater than CUDA block,')
                xmax = xmax-h;   
                xmin = xmin+h; 
                disp(['adjusting xmax to ',num2str(xmax), ' and xmin to ',num2str(xmin),' to compensate.']);  
            end   
            %Now, recompute new grid:
            xvec  = (xmin:h:xmax)';  xres = length(xvec);
            numcudablocks  = (ceil(xres/cudablocksize));
        end        
    else
        disp('Sorry, it seems CUDA is not installed.');
        cuda=0;
    end
end%cuda1

%Initialize solutiona and potential matrices:
U = zeros(xres,1);
V = zeros(size(xvec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%             INITIAL CONDITION FORMULATION       %%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%
V = 0.*ones(size(xvec));   
U = sqrt((2*OM)/s).*sech(sqrt(OM/a).*(xvec - x0)).*exp(1i*((c/(2*a))*xvec));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Linear stability bounds for time step k:
if(k==0)  
   hmin = h;  
   if(method==1)        
      klin = hmin^2/(sqrt(2)*a);
   elseif((method==2))      
      klin = (3/4)*hmin^2/(sqrt(2)*a);
   end
   k = 0.8*klin;   
   disp(['Time step computed!  k: ',num2str(k)]);
end

%Compute total number of time steps for simulation:
endt  = endtw - mod(endtw,k); %<--Make endtime multiple of k
if(mod(endtw,k)~=0), disp(['NOTE! End time adjusted to ',num2str(endt)]); end;
steps = floor(endt/k); %<--Compute number of steps required
if(numframes>=steps), numframes=steps; end;

%Compute number of time steps to compute in mex file:
chunk_size   = ceil(steps/numframes);
extra_steps  = chunk_size*numframes - steps;
if(extra_steps>chunk_size)
    numframes    = ceil(steps/chunk_size);
    disp(['To avoid too much over-time compuations, numframes has been altered to ',num2str(numframes)]);
    extra_steps  = chunk_size*numframes - steps;
end
if(extra_steps>0)
    disp(['NOTE! There will be ',num2str(extra_steps),' extra time steps taken in simulation.']);
    endt = chunk_size*numframes*k;
    steps  = steps + extra_steps;
    disp(['Therefore, true end time will be: ',num2str(endt),' taking ',num2str(steps),' time steps']);
end;
  
%Display simulation info:
disp('----------------------------------------------------');
disp('1D NLSE Parameters:')
disp(['a:          ',num2str(a)]);
disp(['s:          ',num2str(s)]);
disp(['xmin:       ',num2str(xmin)]);
disp(['xmax:       ',num2str(xmax)]);
disp(['Endtime:    ',num2str(endt)]);
disp('Numerical Parameters:');
disp(['h:          ',num2str(h)]);
disp(['k:          ',num2str(k)]);
disp(['GridSize:   ',num2str(xres)]);
disp(['TimeSteps:  ',num2str(steps)]);
disp(['ChunkSize:  ',num2str(chunk_size)]); 
disp(['NumFrames:  ',num2str(numframes)]);
if(precision==1)
    disp('Precision:  Single');
else
    disp('Precision:  Double');
end
if(method==1)
    disp('Method:     RK4+CD     Order(4,2)');
elseif(method==2)
    disp('Method:     RK4+2SHOC  Order(4,4)');
end
if(BC==1)
    disp('Boundary Conditions:  Dirichlet');
elseif(BC==2)
    disp('Boundary Conditions:  MSD');
elseif(BC==3)
    disp('Boundary Conditions:  Uxx = 0');
elseif(BC==4)
    disp('Boundary Conditions:  One-Sided Diff');
end
if(cuda == 1)
    disp( 'CUDA Parameters and Info:');           
    disp(['BlockSize:           ',num2str(cudablocksize)]);
    disp(['CUDAGridSize:        ',num2str(numcudablocks)]);
    disp(['NumBlocks:           ',num2str(numcudablocks)]); 
    disp(['Shared Memory/Block: ',num2str(sharedmemperblock),'KB']);
    disp(['TotalGPUMemReq:      ',num2str(xres*(9 + 2*(method-1))*(4*precision)/1024),'KB']); 
end
%%%%%%%%%                                                 %%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%             PLOTTING SETUP                      %%%%%%%%%
%%%%%%%%%                                                 %%%%%%%%%
fsize       = 20;  %Set font size for figures.
fig_count   = 1;

%Find maximums for plot axis:
maxMod2U  = max(U.*conj(U)); 
maxRealU  = max(real(U));
plotmax   = max([maxMod2U maxRealU]);

%Set up plot figure and plot initial condition:
fig_plot  = figure(fig_count);
fig_count = fig_count+1;
set(fig_plot, 'Name','1D NLSE  t = 0','NumberTitle','off');
plot_real = plot(xvec,real(U),   '--b' ,'LineWidth',4); 
hold on;
plot_v    = plot(xvec,V,   '-g');
plot_imag = plot(xvec,imag(U),   '-.r' ,'LineWidth',4);
plot_mod2 = plot(xvec,U.*conj(U),'-k'  ,'LineWidth',6);
%plot_phase = plot(xvec,angle(U),'--c');
axis([xmin xmax -plotmax (plotmax + 0.1*plotmax)]);
xlabel('x','Fontsize',fsize); 
ylabel('\Psi(x)','Fontsize',fsize);
set(gca,'Fontsize',fsize);
hold off;
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%             BEGIN SIMULATION                    %%%%%%%%%
%%%%%%%%% 
n = 0;
for chunk_count = 1:numframes
 
if(method==1)  
    if(precision==1)
       if(cuda==0)
           U = NLSE1D_TAKE_STEPS_CD_F(U,V,s,a,h^2,BC,chunk_size,k);
       else
           U = NLSE1D_TAKE_STEPS_CD_CUDA_F(U,V,s,a,h^2,BC,chunk_size,k);  
       end
    elseif(precision==2)
       if(cuda==0)
           U = NLSE1D_TAKE_STEPS_CD(U,V,s,a,h^2,BC,chunk_size,k);
       else
           U = NLSE1D_TAKE_STEPS_CD_CUDA_D(U,V,s,a,h^2,BC,chunk_size,k);  
       end
    end           
elseif(method==2)   
    if(precision==1)
       if(cuda==0)
            U = NLSE1D_TAKE_STEPS_2SHOC_F(U,V,s,a,h^2,BC,chunk_size,k);
       else
            U = NLSE1D_TAKE_STEPS_2SHOC_CUDA_F(U,V,s,a,h^2,BC,chunk_size,k);  
       end       
    elseif(precision==2)
       if(cuda==0)
           U = NLSE1D_TAKE_STEPS_2SHOC(U,V,s,a,h^2,BC,chunk_size,k);
       else 
           U = NLSE1D_TAKE_STEPS_2SHOC_CUDA_D(U,V,s,a,h^2,BC,chunk_size,k);  
       end
    end   
end
n = n + chunk_size;
%Set fig titles first, so if sim crashes, we know when it happened:  
set(fig_plot, 'Name',['1D NLSE  t = ',num2str((n*k),'%.2f')]);

%Detect blow-up:
if(max(real(U(:))) > tol || sum(isnan(U(:))) > 0)
  disp('CRAAAAAAAAAAAAAAAASSSSSSSSSSSSSSSSSHHHHHHHHHHHHHHHHHHHHH');
  break;
end

%Plot current time step:
set(plot_mod2,'ydata',U.*conj(U)); 
set(plot_real,'ydata',real(U)); 
set(plot_imag,'ydata',imag(U)); 
%set(plot_phase,'ydata',angle(U));
drawnow;  
end %chunk-count

%%%%%%%%%                                                 %%%%%%%%%
%%%%%%%%%         SIMULATION IS OVER                      %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%