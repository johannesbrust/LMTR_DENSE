% This is a script to minimize functions from the CUTEst test set using 
% the Eig(inf,2) and and L-BFGS method. The intend is to test an
% alternative form of gamma.
% J.B 10/13, This script is modified for the CUTEst test suite.
% J.B 11/28 Multiple CUTEst problems for dense initial matrix.
% J.B 12/01 G4 option.
 
% For details, see:
% O.Burdakov, L. Gong, Y. Yuan, S. Zikrin, On Efficiently Combining Limited
% Memory and Trust-Region Techniques, technical report LiTH-MAT-R-2013/13-SE,
% Department of Mathematics, Link�ping University, 2013.
% http://liu.diva-portal.org/smash/record.jsf?pid=diva2%3A667359

clc;
clear;

parentpath = [cd(cd('..')) '/main'];
addpath(parentpath)

% Set input parameters:
params=struct;
params.m = 5;  % number of L-BFGS updates
params.gtol = 1e-5;  % exit if ||g||_2<gtol*max(1,||x||_2)
params.ranktol = 1e-7;  % tolerance for establishing rank of V
params.dflag = 0;  % display parameter, 1 if to display information per iteration
params.trtol = 0.1;  % exit MS algorithm if abs(||s||-delta)<trtol*delta
params.ftol=1e-11;  % tolerance on relative function reduction





params.maxit = 100000;

CUTEst_init  % initialize CUTEr, see appropriate documentation 
fid =fopen('cutest_list.txt','r');  % read file that contains CUTEr test problems
sfil='TestResults';

%leg={'L-BFGS-MTBT','L-BFGS','EIG(\infty,2)','m=10','m=15','EIG-MS','EIG-MS(2,2)','BWX-MS','D-DOGL'};

leg={'EIG(\infty,2)','EIG(\infty,2)-G1'};

paramsG1            =params;
paramsG1.gamma_opt  = 1;

paramsG2            =params;
paramsG2.gamma_opt  = 2;

paramsG3            =params;
paramsG3.gamma_opt  = 3;

paramsG4            =params;
paramsG4.gamma_opt  = 4;

paramsG5            =params;
paramsG5.gamma_opt  = 5;



% Initialize storage for output information
numRuns=1;
numAlgorithms=length(leg);
numProblems=62;
ex=zeros(numProblems,numAlgorithms);
numf=zeros(numProblems,numAlgorithms);
numg=zeros(numProblems,numAlgorithms);
numit=zeros(numProblems,numAlgorithms);
tcpu=zeros(numProblems,numRuns,numAlgorithms);
t_aver=zeros(numProblems,numAlgorithms);
tract=zeros(numProblems,numAlgorithms);
numrst=zeros(numProblems,numAlgorithms);

p=1;
tline = fgets(fid);
while ischar(tline)     
    tline = fgets(fid);       
    
    if ~strcmp(tline(1),'%')  && ischar(tline)   
        
        eval(['!runcutest -p matlab -D ' tline]);
        %eval(['!runcuter --package mx --decode ' tline]);
        prob = cutest_setup();
        x0 = prob.x; 
                
         s=1;        
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE,x0,paramsG1,numRuns);
        
        s=s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE,x0,paramsG5,numRuns);
        
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
%             runAlgorithm(@LMTR_EIG_inf_2_DENSE,x0,paramsG2,numRuns);  
%         
%         s = s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
%             runAlgorithm(@LMTR_EIG_inf_2_DENSE,x0,paramsG3,numRuns); 
%         
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),numrst(p,s)]=...
%             runAlgorithm(@LMTR_EIG_inf_2,x0,params15,numRuns);
%         
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),numrst(p,s)]=...
%             runAlgorithm(@LMTR_EIG_MS,x0,params,numRuns);        
%                 
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),numrst(p,s)]=...
%             runAlgorithm(@LMTR_EIG_MS_2_2,x0,params,numRuns);        
%         
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
%             runAlgorithm(@LMTR_BWX_MS,x0,params,numRuns);                                        
%                 
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
%             runAlgorithm(@LMTR_DDOGL,x0,params,numRuns);                
        
        % Average CPU time
        if p==1 && numRuns > 2
            for si=1:s
                t_aver(p,si) = sum(tcpu(si,3:numRuns,p))/(numRuns-2);
            end
        else
            for si=1:s
                t_aver(p,si) = sum(tcpu(si,2:numRuns,p))/(numRuns-1);
            end
        end             
        
        cutest_terminate();
        
        p=p+1;            
    end
    save(sfil);
end


%% Plot performance profiles

% Fig.1: EIG(inf,2) for m=5,10,15
% indAlg=[3,4,5];
% leg(3)={'m=5'};
% perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg));
% perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg));

for i = 1:numProblems
    %
    %t_aver(i,:) = tcpu(1:2,:,i)';
    t_aver(i,1:2) = tcpu(1:2,:,i)';

end

% Fig.2: EIG(inf,2) and L-BFGS
indAlg=1:2;
%leg(3)={'EIG(\infty,2)'};
perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg));
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg));

% Fig.3: EIG(inf,2) and L-BFGS in those problems where the step-size one was 
% rejected in, at least, 30% of iterations
% indProb=find(tract(:,2)./numit(:,2)>=0.30);
% perf(ex(indProb,indAlg),numit(indProb,indAlg),leg(indAlg));
% perf(ex(indProb,indAlg),t_aver(indProb,indAlg),leg(indAlg));
% 
% % Fig.4: EIG(inf,2), EIG-MS and EIG-MS(2,2)
% indAlg=[3,6,7];
% perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg));
% perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg));
% 
% % Fig.5: EIG(inf,2) and BWX-MS
% indAlg=[3,8];
% perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg));
% perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg));
% 
% % Fig.6: EIG(inf,2) and D-DOGL
% indAlg=[3,9];
% perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg));
% perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg));
