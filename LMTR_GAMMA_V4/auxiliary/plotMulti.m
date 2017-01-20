% Script to plot results for running different gamma values on a set of
% CUTEst problems.
clear;
clc;
close all;
load('ResultsMulti');

leg         = {'Eig-Inf','Eig-Inf-G1','Eig-Inf-G2','Eig-Inf-G3'};
data_l      = {'numf','numg','numit','tcpu'};

numalg      = length(leg);

inclt       = 1;
if inclt > 0
    for i = 1:numProblems
   
        t_aver(i,:) = tcpu(1:numalg,:,i)';
    
    end
end


ex2         = ex;
ex2(ex<0)   = NaN;


data        = [numf;numg;numit;t_aver];


for i = 1:numalg
    
    idx = ((i-1)*numProblems +1):(i*numProblems);
    h   = plot(data(idx,:).*ex2,'*-' );
    
    legend(leg);
    xlabel('Problem');
    ylabel(data_l(i));

    fnamepdf    = strcat(data_l(i),'.pdf');
    fnamejpg    = strcat(data_l(i),'.jpg');
    
    nm          = data_l(i); 
    
    pathpdf     = strcat(pwd,'/FigMulti/');
    pathjpg     = strcat(pwd,'/FigMulti/');
    
    %pathpdf     = strcat(pwd, strcat('/FigMulti/',fnamepdf));
    %pathjpg     = strcat(pwd, strcat('/FigMulti/',fnamejpg));
    
    saveas(gcf,fullfile(pwd,'FigMulti',data_l{i}),'jpg');
    saveas(gcf,fullfile(pwd,'FigMulti',data_l{i}),'pdf');
    %saveas(gcf,nm,'fig');
    
    %saveas(gcf,pathpdf);
    %saveas(gcf,pathjpg);
    
end



