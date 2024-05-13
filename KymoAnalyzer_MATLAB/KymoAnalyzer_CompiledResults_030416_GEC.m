%% Script to compare different genotypes/conditions from KymoAnalyzer Results
% Written by Sylvia Neumann 10/31/14
% Maintained by George Campbell 03/04/16
% Organization of folders: ParentFolder
%                                Subfolders: Genotype1
%                                   Subfolders:Experiment1
%                                   Subfolders:Experiment2
%                                   Subfolders:ExperimentN
%                                Subfolders: Genotype2
%                                Subfolders: GenotypeN

%                           Example folders:
%                           Figure3
%                               WT
%                                   001
%                                       PooledData
%                                   002
%                                       PooledData
%                               KIF5C_KO
%                                   001
%                                       PooledData
% NAMING RESTRICTIONS: Folder names cannot include special characters, spaces, or dashes. Underscores are OK.
%                      Long names may be truncated in certain variables
%                      that compare two different genotypes (to 63 total
%                      characters)
% Script dependencies:  padcat.m
%                       rotateXLabels.m,
%                       KymoAnalyzer_ResultsAll_111015_GEC.m
%                       multigauss_fitting_110414.m
%                       normality_test_with_R_110414.m
%% Clear Matlab
clc;clear; close hidden all;

%% Choose experimental folder (Folder that contains all movies and pooled data)

directoryMain= uigetdir(pwd,'Select your parent folder (Folder that contains all experiments for all genotypes)');
cd (directoryMain);
temp = dir(directoryMain);
% remove hidden MAC OSX files
for k = length(temp):-1:1
    % remove non-folders
    if ~temp(k).isdir
        temp(k) = [ ];
        continue
    end

    % remove folders starting with .
    fname = temp(k).name;
    if fname(1) == '.'
        temp(k) = [ ];
    end
end

%% global variables for MultiGaussian calculations
global sigma 
global BIC
global classification

%% get genotype folders and subfolders
% list is array with folder names that are within parent folder
list=([]);
j=1;
for i=1:length(temp)
    if strcmp(temp(i,1).name,'Figures')==1 || strcmp(temp(i,1).name,'MatlabResults')==1
   % there is nothing here? 
    else
    list(j,1).name=temp(i,1).name;
    j=j+1;
    end
end

% rewriting order of list with user input
% genotype folders are alphabetically sorted. User can input the order of
% genotypes as theywill be plotted in subsequent graphs. E.g. if the
% genotypes are wt and dyn1, without user input plots will show dyn1 first
% and then wt. User can choose to plot wt first and then dyn1. 
% *** does this work? FIX
index={};
genotype={};
defaultAnswer={};

for i=1:length(list)
    index{i}=num2str(i);
    genotype{i}=strcat('Genotype_',index{i});
    defaultAnswer{i}=char(list(i,:).name);
end

numlines=1;
dlg_title='Choose the order of Genotypes to appear in plots';
answer = inputdlg(genotype,dlg_title,numlines,defaultAnswer);

array=([]);
j=0;
for i=1:length(list)
    fname=answer{i};
    if isempty(fname)
        j=j+1;
        array(j,1)=i;
    end
    if ~isempty(fname)
        list(i,:).name=answer{i};
    end
end
list(array)=[ ];

% read in data for each experiment folder
KymoAnalyzer_ResultsAll_111015_GEC(directoryMain);

mkdir(fullfile(directoryMain,'Figures'));
mkdir(fullfile(directoryMain,'MatlabResults'));

%% loop through all Genotype folders
% list is array with names of genotype folders
% list2 is array with names of experiment folder in genotype folder

for i=1:length(list)
    %if length(list(i,1).name>0
    subfolder=fullfile(directoryMain,list(i,1).name,filesep);
    temp2=dir(subfolder);
    
    for k = length(temp2):-1:1
        % remove hidden Mac files
        % remove non-folders
        if ~temp2(k).isdir
            temp2(k) = [ ];
            continue
        end

        % remove folders starting with .
        fname = temp2(k).name;
        if fname(1) == '.'
            temp2(k) = [ ];
        end
    end
    
    list2=([]);
    k=1;
    for j=1:length(temp2)
        list2(k,1).name=temp2(j,1).name;
        k=k+1;
    end
    
    AllCP.(list(i,1).name).Num=([]);
    AllCP.(list(i,1).name).PCT=([]);
    AllNCP.(list(i,1).name).Num=([]);
    AllNCP.(list(i,1).name).PCT=([]);
    AllFlux.(list(i,1).name).Value=([]);
    AllPD.(list(i,1).name).Value=([]);
    AllPF.(list(i,1).name).Value=([]);
    AllPFperSec.(list(i,1).name).Value=([]);
    AllsplitPD.(list(i,1).name).ante=([]);
    AllsplitPF.(list(i,1).name).ante=([]);
    AllsplitPFperSec.(list(i,1).name).ante=([]);
    AllsplitPD.(list(i,1).name).retro=([]);
    AllsplitPF.(list(i,1).name).retro=([]);
    AllsplitPFperSec.(list(i,1).name).retro=([]);
    AllsplitPD.(list(i,1).name).rev=([]);
    AllsplitPF.(list(i,1).name).rev=([]);
    AllsplitPFperSec.(list(i,1).name).rev=([]);
    AllSF.(list(i,1).name).Value=([]);
    AllrevSFperSec.(list(i,1).name).Value=([]);
    AllrevSF.(list(i,1).name).Value=([]);
    AllSFperSec.(list(i,1).name).Value=([]);
    AllSV.(list(i,1).name).ante=([]);
    AllSV.(list(i,1).name).retro=([]);
    AllNV.(list(i,1).name).ante=([]);
    AllNV.(list(i,1).name).retro=([]);
    AllRL.(list(i,1).name).ante=([]);
    AllRL.(list(i,1).name).retro=([]);
    AllcomRL.(list(i,1).name).ante=([]);
    AllcomRL.(list(i,1).name).retro=([]); 
    AllcomSV.(list(i,1).name).ante=([]); 
    AllcomSV.(list(i,1).name).retro=([]); 
    AllDensity.(list(i,1).name).Value=([]);
    
    for j=1:length(list2)
        fileResults=strcat(list(i,1).name,'_',list2(j,1).name,'.m');
        fileResults=fullfile(directoryMain,list(i,1).name,list2(j,1).name,'MatlabResults',fileResults);
        load(fileResults,'-mat');
        
        % concatenate all arrays to pool all data from each genotype
        AllCP.(list(i,1).name).Num=cat(1,AllCP.(list(i,1).name).Num,CP.Num);
        AllCP.(list(i,1).name).PCT=cat(1,AllCP.(list(i,1).name).PCT,CP.PCT);
        AllNCP.(list(i,1).name).Num=cat(1,AllNCP.(list(i,1).name).Num,NCP.Num);
        AllNCP.(list(i,1).name).PCT=cat(1,AllNCP.(list(i,1).name).PCT,NCP.PCT);
        AllFlux.(list(i,1).name).Value=cat(1, AllFlux.(list(i,1).name).Value,Flux.Value);
        AllPD.(list(i,1).name).Value=cat(1, AllPD.(list(i,1).name).Value,PD.Value);
        AllPF.(list(i,1).name).Value=cat(1, AllPF.(list(i,1).name).Value,PF.Value);
        AllPFperSec.(list(i,1).name).Value=cat(1, AllPFperSec.(list(i,1).name).Value,PFperSec.Value);        
        AllsplitPD.(list(i,1).name).ante=cat(1, AllsplitPD.(list(i,1).name).ante,splitPD.ante);
        AllsplitPF.(list(i,1).name).ante=cat(1, AllsplitPF.(list(i,1).name).ante,splitPF.ante);
        AllsplitPFperSec.(list(i,1).name).ante=cat(1, AllsplitPFperSec.(list(i,1).name).ante,splitPFperSec.ante); 
        AllsplitPD.(list(i,1).name).retro=cat(1, AllsplitPD.(list(i,1).name).retro,splitPD.retro);
        AllsplitPF.(list(i,1).name).retro=cat(1, AllsplitPF.(list(i,1).name).retro,splitPF.retro);
        AllsplitPFperSec.(list(i,1).name).retro=cat(1, AllsplitPFperSec.(list(i,1).name).retro,splitPFperSec.retro);
        AllsplitPD.(list(i,1).name).rev=cat(1, AllsplitPD.(list(i,1).name).rev,splitPD.rev);
        AllsplitPF.(list(i,1).name).rev=cat(1, AllsplitPF.(list(i,1).name).rev,splitPF.rev);
        AllsplitPFperSec.(list(i,1).name).rev=cat(1, AllsplitPFperSec.(list(i,1).name).rev,splitPFperSec.rev);
        AllSF.(list(i,1).name).Value=cat(1, AllSF.(list(i,1).name).Value,SF.Value);
        AllSFperSec.(list(i,1).name).Value=cat(1, AllSFperSec.(list(i,1).name).Value,SFperSec.Value);
        AllrevSF.(list(i,1).name).Value=cat(1, AllrevSF.(list(i,1).name).Value,revSF.Value);
        AllrevSFperSec.(list(i,1).name).Value=cat(1, AllrevSFperSec.(list(i,1).name).Value,revSFperSec.Value);
        AllSV.(list(i,1).name).ante=cat(1, AllSV.(list(i,1).name).ante,SV.ante);
        AllSV.(list(i,1).name).retro=cat(1, AllSV.(list(i,1).name).retro,SV.retro);
        AllNV.(list(i,1).name).ante=cat(1, AllNV.(list(i,1).name).ante,NV.ante);
        AllNV.(list(i,1).name).retro=cat(1, AllNV.(list(i,1).name).retro,NV.retro);
        AllRL.(list(i,1).name).ante=cat(1, AllRL.(list(i,1).name).ante,RL.ante);
        AllRL.(list(i,1).name).retro=cat(1, AllRL.(list(i,1).name).retro,RL.retro);
        AllcomRL.(list(i,1).name).ante=cat(1, AllcomRL.(list(i,1).name).ante,comRL.ante);
        AllcomRL.(list(i,1).name).retro=cat(1, AllcomRL.(list(i,1).name).retro,comRL.retro);
        AllcomSV.(list(i,1).name).ante=cat(1, AllcomSV.(list(i,1).name).ante,comSV.ante);
        AllcomSV.(list(i,1).name).retro=cat(1, AllcomSV.(list(i,1).name).retro,comSV.retro);
        AllDensity.(list(i,1).name).Value=cat(1, AllDensity.(list(i,1).name).Value,Density.Value);
    end
end  
varlist = {'CP','Flux','NCP','NV','PD','PF','PFperSec','splitPD','splitPF','splitPFperSec','RL','comRL','SF','SFperSec','revSF','revSFperSec','SV','directory','figure1','figure2','file','histNV','histRL','histSV','fileResults','i','j',...
    'temp','temp2','subfolder'}; %What is this for? FIX
clear (varlist{:}); clear varlist;

%% Calculating the stats: average, STD and SEM
% list.name is array that contains list of genotypes to compare)

for i=1:length(list)
    
    AllCP.(list(i,1).name).Stats(1,1)=100*mean(AllCP.(list(i,1).name).PCT(:,1));
    AllCP.(list(i,1).name).Stats(1,2)=100*mean(AllCP.(list(i,1).name).PCT(:,2));
    AllCP.(list(i,1).name).Stats(1,3)=100*mean(AllCP.(list(i,1).name).PCT(:,3));
    AllCP.(list(i,1).name).Stats(1,4)=100*mean(AllCP.(list(i,1).name).PCT(:,4));
    
    AllCP.(list(i,1).name).Stats(2,1)=100*std(AllCP.(list(i,1).name).PCT(:,1));
    AllCP.(list(i,1).name).Stats(2,2)=100*std(AllCP.(list(i,1).name).PCT(:,2));
    AllCP.(list(i,1).name).Stats(2,3)=100*std(AllCP.(list(i,1).name).PCT(:,3));
    AllCP.(list(i,1).name).Stats(2,4)=100*std(AllCP.(list(i,1).name).PCT(:,4));
    
    AllCP.(list(i,1).name).Stats(3,1)=100*std(AllCP.(list(i,1).name).PCT(:,1))/sqrt(length(AllCP.(list(i,1).name).PCT(:,1)));
    AllCP.(list(i,1).name).Stats(3,2)=100*std(AllCP.(list(i,1).name).PCT(:,2))/sqrt(length(AllCP.(list(i,1).name).PCT(:,2)));
    AllCP.(list(i,1).name).Stats(3,3)=100*std(AllCP.(list(i,1).name).PCT(:,3))/sqrt(length(AllCP.(list(i,1).name).PCT(:,3)));
    AllCP.(list(i,1).name).Stats(3,4)=100*std(AllCP.(list(i,1).name).PCT(:,4))/sqrt(length(AllCP.(list(i,1).name).PCT(:,4)));
    
    AllNCP.(list(i,1).name).Stats(1,1)=100*mean(AllNCP.(list(i,1).name).PCT(:,1));
    AllNCP.(list(i,1).name).Stats(1,2)=100*mean(AllNCP.(list(i,1).name).PCT(:,2));
    AllNCP.(list(i,1).name).Stats(1,3)=100*mean(AllNCP.(list(i,1).name).PCT(:,3));
    
    AllNCP.(list(i,1).name).Stats(2,1)=100*std(AllNCP.(list(i,1).name).PCT(:,1));
    AllNCP.(list(i,1).name).Stats(2,2)=100*std(AllNCP.(list(i,1).name).PCT(:,2));
    AllNCP.(list(i,1).name).Stats(2,3)=100*std(AllNCP.(list(i,1).name).PCT(:,3));
    
    AllNCP.(list(i,1).name).Stats(3,1)=100*std(AllNCP.(list(i,1).name).PCT(:,1))/sqrt(length(AllNCP.(list(i,1).name).PCT(:,1)));
    AllNCP.(list(i,1).name).Stats(3,2)=100*std(AllNCP.(list(i,1).name).PCT(:,2))/sqrt(length(AllNCP.(list(i,1).name).PCT(:,2)));
    AllNCP.(list(i,1).name).Stats(3,3)=100*std(AllNCP.(list(i,1).name).PCT(:,3))/sqrt(length(AllNCP.(list(i,1).name).PCT(:,3)));
   
    AllDensity.(list(i,1).name).Stats(1,1)=mean(AllDensity.(list(i,1).name).Value(:,1));
    AllDensity.(list(i,1).name).Stats(1,2)=mean(AllDensity.(list(i,1).name).Value(:,2));
    AllDensity.(list(i,1).name).Stats(1,3)=mean(AllDensity.(list(i,1).name).Value(:,3));
    AllDensity.(list(i,1).name).Stats(1,4)=mean(AllDensity.(list(i,1).name).Value(:,4));
    
    AllDensity.(list(i,1).name).Stats(2,1)=std(AllDensity.(list(i,1).name).Value(:,1));
    AllDensity.(list(i,1).name).Stats(2,2)=std(AllDensity.(list(i,1).name).Value(:,2));
    AllDensity.(list(i,1).name).Stats(2,3)=std(AllDensity.(list(i,1).name).Value(:,3));
    AllDensity.(list(i,1).name).Stats(2,4)=std(AllDensity.(list(i,1).name).Value(:,4));
    
    AllDensity.(list(i,1).name).Stats(3,1)=std(AllDensity.(list(i,1).name).Value(:,1))/sqrt(length(AllDensity.(list(i,1).name).Value(:,1)));
    AllDensity.(list(i,1).name).Stats(3,2)=std(AllDensity.(list(i,1).name).Value(:,2))/sqrt(length(AllDensity.(list(i,1).name).Value(:,2)));
    AllDensity.(list(i,1).name).Stats(3,3)=std(AllDensity.(list(i,1).name).Value(:,3))/sqrt(length(AllDensity.(list(i,1).name).Value(:,3)));
    AllDensity.(list(i,1).name).Stats(3,4)=std(AllDensity.(list(i,1).name).Value(:,4))/sqrt(length(AllDensity.(list(i,1).name).Value(:,4)));
    
    AllFlux.(list(i,1).name).Stats(1,1)=mean(AllFlux.(list(i,1).name).Value(:,1));
    AllFlux.(list(i,1).name).Stats(1,2)=mean(AllFlux.(list(i,1).name).Value(:,2));
    AllFlux.(list(i,1).name).Stats(1,3)=mean(AllFlux.(list(i,1).name).Value(:,3));
    
    AllFlux.(list(i,1).name).Stats(2,1)=std(AllFlux.(list(i,1).name).Value(:,1));
    AllFlux.(list(i,1).name).Stats(2,2)=std(AllFlux.(list(i,1).name).Value(:,2));
    AllFlux.(list(i,1).name).Stats(2,3)=std(AllFlux.(list(i,1).name).Value(:,3));
    
    AllFlux.(list(i,1).name).Stats(3,1)=std(AllFlux.(list(i,1).name).Value(:,1))/sqrt(length(AllFlux.(list(i,1).name).Value(:,1)));
    AllFlux.(list(i,1).name).Stats(3,2)=std(AllFlux.(list(i,1).name).Value(:,2))/sqrt(length(AllFlux.(list(i,1).name).Value(:,2)));
    AllFlux.(list(i,1).name).Stats(3,3)=std(AllFlux.(list(i,1).name).Value(:,3))/sqrt(length(AllFlux.(list(i,1).name).Value(:,3)));

    AllPD.(list(i,1).name).Stats(1,1)=mean(AllPD.(list(i,1).name).Value(:,1));
    AllPD.(list(i,1).name).Stats(2,1)=std(AllPD.(list(i,1).name).Value(:,1));
    AllPD.(list(i,1).name).Stats(3,1)=std(AllPD.(list(i,1).name).Value(:,1))/sqrt(length(AllPD.(list(i,1).name).Value(:,1)));
    
    AllPF.(list(i,1).name).Stats(1,1)=mean(AllPF.(list(i,1).name).Value(:,1));
    AllPF.(list(i,1).name).Stats(2,1)=std(AllPF.(list(i,1).name).Value(:,1));
    AllPF.(list(i,1).name).Stats(3,1)=std(AllPF.(list(i,1).name).Value(:,1))/sqrt(length(AllPF.(list(i,1).name).Value(:,1)));
    
    AllPFperSec.(list(i,1).name).Stats(1,1)=mean(AllPFperSec.(list(i,1).name).Value(:,1));
    AllPFperSec.(list(i,1).name).Stats(2,1)=std(AllPFperSec.(list(i,1).name).Value(:,1));
    AllPFperSec.(list(i,1).name).Stats(3,1)=std(AllPFperSec.(list(i,1).name).Value(:,1))/sqrt(length(AllPFperSec.(list(i,1).name).Value(:,1)));
    
    AllsplitPD.(list(i,1).name).Stats(1,1)=mean(AllsplitPD.(list(i,1).name).ante(:,1));
    AllsplitPD.(list(i,1).name).Stats(1,2)=mean(AllsplitPD.(list(i,1).name).retro(:,1));
    AllsplitPD.(list(i,1).name).Stats(1,3)=mean(AllsplitPD.(list(i,1).name).rev(:,1));
    AllsplitPD.(list(i,1).name).Stats(2,1)=std(AllsplitPD.(list(i,1).name).ante(:,1));
    AllsplitPD.(list(i,1).name).Stats(2,2)=std(AllsplitPD.(list(i,1).name).retro(:,1));
    AllsplitPD.(list(i,1).name).Stats(2,3)=std(AllsplitPD.(list(i,1).name).rev(:,1));
    AllsplitPD.(list(i,1).name).Stats(3,1)=std(AllsplitPD.(list(i,1).name).ante(:,1))/sqrt(length(AllsplitPD.(list(i,1).name).ante(:,1)));
    AllsplitPD.(list(i,1).name).Stats(3,2)=std(AllsplitPD.(list(i,1).name).retro(:,1))/sqrt(length(AllsplitPD.(list(i,1).name).retro(:,1)));
    AllsplitPD.(list(i,1).name).Stats(3,3)=std(AllsplitPD.(list(i,1).name).rev(:,1))/sqrt(length(AllsplitPD.(list(i,1).name).rev(:,1)));    
    
    AllsplitPF.(list(i,1).name).Stats(1,1)=mean(AllsplitPF.(list(i,1).name).ante(:,1));
    AllsplitPF.(list(i,1).name).Stats(1,2)=mean(AllsplitPF.(list(i,1).name).retro(:,1));
    AllsplitPF.(list(i,1).name).Stats(1,3)=mean(AllsplitPF.(list(i,1).name).rev(:,1));
    AllsplitPF.(list(i,1).name).Stats(2,1)=std(AllsplitPF.(list(i,1).name).ante(:,1));
    AllsplitPF.(list(i,1).name).Stats(2,2)=std(AllsplitPF.(list(i,1).name).retro(:,1));
    AllsplitPF.(list(i,1).name).Stats(2,3)=std(AllsplitPF.(list(i,1).name).rev(:,1));
    AllsplitPF.(list(i,1).name).Stats(3,1)=std(AllsplitPF.(list(i,1).name).ante(:,1))/sqrt(length(AllsplitPF.(list(i,1).name).ante(:,1)));
    AllsplitPF.(list(i,1).name).Stats(3,2)=std(AllsplitPF.(list(i,1).name).retro(:,1))/sqrt(length(AllsplitPF.(list(i,1).name).retro(:,1)));
    AllsplitPF.(list(i,1).name).Stats(3,3)=std(AllsplitPF.(list(i,1).name).rev(:,1))/sqrt(length(AllsplitPF.(list(i,1).name).rev(:,1)));    
    
    AllsplitPFperSec.(list(i,1).name).Stats(1,1)=mean(AllsplitPFperSec.(list(i,1).name).ante(:,1));
    AllsplitPFperSec.(list(i,1).name).Stats(1,2)=mean(AllsplitPFperSec.(list(i,1).name).retro(:,1));
    AllsplitPFperSec.(list(i,1).name).Stats(1,3)=mean(AllsplitPFperSec.(list(i,1).name).rev(:,1));
    AllsplitPFperSec.(list(i,1).name).Stats(2,1)=std(AllsplitPFperSec.(list(i,1).name).ante(:,1));
    AllsplitPFperSec.(list(i,1).name).Stats(2,2)=std(AllsplitPFperSec.(list(i,1).name).retro(:,1));
    AllsplitPFperSec.(list(i,1).name).Stats(2,3)=std(AllsplitPFperSec.(list(i,1).name).rev(:,1));
    AllsplitPFperSec.(list(i,1).name).Stats(3,1)=std(AllsplitPFperSec.(list(i,1).name).ante(:,1))/sqrt(length(AllsplitPFperSec.(list(i,1).name).ante(:,1)));
    AllsplitPFperSec.(list(i,1).name).Stats(3,2)=std(AllsplitPFperSec.(list(i,1).name).retro(:,1))/sqrt(length(AllsplitPFperSec.(list(i,1).name).retro(:,1)));
    AllsplitPFperSec.(list(i,1).name).Stats(3,3)=std(AllsplitPFperSec.(list(i,1).name).rev(:,1))/sqrt(length(AllsplitPFperSec.(list(i,1).name).rev(:,1)));    
    
    AllSF.(list(i,1).name).Stats(1,1)=mean(AllSF.(list(i,1).name).Value(:,1));
    AllSF.(list(i,1).name).Stats(2,1)=std(AllSF.(list(i,1).name).Value(:,1));
    AllSF.(list(i,1).name).Stats(3,1)=std(AllSF.(list(i,1).name).Value(:,1))/sqrt(length(AllSF.(list(i,1).name).Value(:,1)));
    
    AllSFperSec.(list(i,1).name).Stats(1,1)=mean(AllSFperSec.(list(i,1).name).Value(:,1));
    AllSFperSec.(list(i,1).name).Stats(2,1)=std(AllSFperSec.(list(i,1).name).Value(:,1));
    AllSFperSec.(list(i,1).name).Stats(3,1)=std(AllSFperSec.(list(i,1).name).Value(:,1))/sqrt(length(AllSFperSec.(list(i,1).name).Value(:,1)));
    
    AllrevSF.(list(i,1).name).Stats(1,1)=mean(AllrevSF.(list(i,1).name).Value(:,1));
    AllrevSF.(list(i,1).name).Stats(2,1)=std(AllrevSF.(list(i,1).name).Value(:,1));
    AllrevSF.(list(i,1).name).Stats(3,1)=std(AllrevSF.(list(i,1).name).Value(:,1))/sqrt(length(AllrevSF.(list(i,1).name).Value(:,1)));
    
    AllrevSFperSec.(list(i,1).name).Stats(1,1)=mean(AllrevSFperSec.(list(i,1).name).Value(:,1));
    AllrevSFperSec.(list(i,1).name).Stats(2,1)=std(AllrevSFperSec.(list(i,1).name).Value(:,1));
    AllrevSFperSec.(list(i,1).name).Stats(3,1)=std(AllrevSFperSec.(list(i,1).name).Value(:,1))/sqrt(length(AllrevSFperSec.(list(i,1).name).Value(:,1)));
    
    AllSV.(list(i,1).name).Stats(1,1)=mean(AllSV.(list(i,1).name).ante(:,1));
    AllSV.(list(i,1).name).Stats(1,2)=mean(AllSV.(list(i,1).name).retro(:,1));
    AllSV.(list(i,1).name).Stats(2,1)=std(AllSV.(list(i,1).name).ante(:,1));
    AllSV.(list(i,1).name).Stats(2,2)=std(AllSV.(list(i,1).name).retro(:,1));
    AllSV.(list(i,1).name).Stats(3,1)=std(AllSV.(list(i,1).name).ante(:,1))/sqrt(length(AllSV.(list(i,1).name).ante(:,1)));
    AllSV.(list(i,1).name).Stats(3,2)=std(AllSV.(list(i,1).name).retro(:,1))/sqrt(length(AllSV.(list(i,1).name).retro(:,1)));
    
    AllNV.(list(i,1).name).Stats(1,1)=mean(AllNV.(list(i,1).name).ante(:,1));
    AllNV.(list(i,1).name).Stats(1,2)=mean(AllNV.(list(i,1).name).retro(:,1));
    AllNV.(list(i,1).name).Stats(2,1)=std(AllNV.(list(i,1).name).ante(:,1));
    AllNV.(list(i,1).name).Stats(2,2)=std(AllNV.(list(i,1).name).retro(:,1));
    AllNV.(list(i,1).name).Stats(3,1)=std(AllNV.(list(i,1).name).ante(:,1))/sqrt(length(AllNV.(list(i,1).name).ante(:,1)));
    AllNV.(list(i,1).name).Stats(3,2)=std(AllNV.(list(i,1).name).retro(:,1))/sqrt(length(AllNV.(list(i,1).name).retro(:,1)));
    
    AllRL.(list(i,1).name).Stats(1,1)=mean(AllRL.(list(i,1).name).ante(:,1));
    AllRL.(list(i,1).name).Stats(1,2)=mean(AllRL.(list(i,1).name).retro(:,1));
    AllRL.(list(i,1).name).Stats(2,1)=std(AllRL.(list(i,1).name).ante(:,1));
    AllRL.(list(i,1).name).Stats(2,2)=std(AllRL.(list(i,1).name).retro(:,1));   
    AllRL.(list(i,1).name).Stats(3,1)=std(AllRL.(list(i,1).name).ante(:,1))/sqrt(length(AllRL.(list(i,1).name).ante(:,1)));
    AllRL.(list(i,1).name).Stats(3,2)=std(AllRL.(list(i,1).name).retro(:,1))/sqrt(length(AllRL.(list(i,1).name).retro(:,1)));
    
    AllcomRL.(list(i,1).name).Stats(1,1)=mean(AllcomRL.(list(i,1).name).ante(:,1));
    AllcomRL.(list(i,1).name).Stats(1,2)=mean(AllcomRL.(list(i,1).name).retro(:,1));
    AllcomRL.(list(i,1).name).Stats(2,1)=std(AllcomRL.(list(i,1).name).ante(:,1));
    AllcomRL.(list(i,1).name).Stats(2,2)=std(AllcomRL.(list(i,1).name).retro(:,1));   
    AllcomRL.(list(i,1).name).Stats(3,1)=std(AllcomRL.(list(i,1).name).ante(:,1))/sqrt(length(AllcomRL.(list(i,1).name).ante(:,1)));
    AllcomRL.(list(i,1).name).Stats(3,2)=std(AllcomRL.(list(i,1).name).retro(:,1))/sqrt(length(AllcomRL.(list(i,1).name).retro(:,1)));
    
    AllcomSV.(list(i,1).name).Stats(1,1)=mean(AllcomSV.(list(i,1).name).ante(:,1));
    AllcomSV.(list(i,1).name).Stats(1,2)=mean(AllcomSV.(list(i,1).name).retro(:,1));
    AllcomSV.(list(i,1).name).Stats(2,1)=std(AllcomSV.(list(i,1).name).ante(:,1));
    AllcomSV.(list(i,1).name).Stats(2,2)=std(AllcomSV.(list(i,1).name).retro(:,1));   
    AllcomSV.(list(i,1).name).Stats(3,1)=std(AllcomSV.(list(i,1).name).ante(:,1))/sqrt(length(AllcomSV.(list(i,1).name).ante(:,1)));
    AllcomSV.(list(i,1).name).Stats(3,2)=std(AllcomSV.(list(i,1).name).retro(:,1))/sqrt(length(AllcomSV.(list(i,1).name).retro(:,1)));
end

%% organizing data for bar graphs % 
% list.name is array that contains list of genotypes to compare)
% xx.Mean and xx.STD are arrays that group the data together that will be
% compared. The order is always
   % anterograde
   % retrograde
   % reversal
   % stationary
% barXX variables are arrays that contain the datasets to be plotted

CP=([]);
NCP=([]);
Density=([]);
NV=([]);
SV=([]);
RL=([]);
comRL=([]);
PD=([]);
PF=([]);
PFperSec=([]);
splitPD=([]);
splitPF=([]);
splitPFperSec=([]);
SF=([]);
SFperSec=([]);
revSF=([]);
revSFperSec=([]);
Flux = ([]);
comSV=([]);

%write data for bar graphs
for i=1:length(list)
    CP.Mean(1,i)= AllCP.(list(i,1).name).Stats(1,1);    
    CP.Mean(2,i)= AllCP.(list(i,1).name).Stats(1,2);    
    CP.Mean(3,i)= AllCP.(list(i,1).name).Stats(1,3);    
    CP.Mean(4,i)= AllCP.(list(i,1).name).Stats(1,4);    
    
    CP.STD(1,i)= AllCP.(list(i,1).name).Stats(2,1);
    CP.STD(2,i)= AllCP.(list(i,1).name).Stats(2,2);
    CP.STD(3,i)= AllCP.(list(i,1).name).Stats(2,3);
    CP.STD(4,i)= AllCP.(list(i,1).name).Stats(2,4);
    
    CP.SEM(1,i)= AllCP.(list(i,1).name).Stats(3,1);
    CP.SEM(2,i)= AllCP.(list(i,1).name).Stats(3,2);
    CP.SEM(3,i)= AllCP.(list(i,1).name).Stats(3,3);
    CP.SEM(4,i)= AllCP.(list(i,1).name).Stats(3,4);
    
    NCP.Mean(1,i)= AllNCP.(list(i,1).name).Stats(1,1);
    NCP.Mean(2,i)= AllNCP.(list(i,1).name).Stats(1,2);
    NCP.Mean(3,i)= AllNCP.(list(i,1).name).Stats(1,3);
    NCP.Mean(4,i)= 0;
    
    NCP.STD(1,i)= AllNCP.(list(i,1).name).Stats(2,1);
    NCP.STD(2,i)= AllNCP.(list(i,1).name).Stats(2,2);
    NCP.STD(3,i)= AllNCP.(list(i,1).name).Stats(2,3);
    NCP.STD(4,i)= 0;
    
    NCP.SEM(1,i)= AllNCP.(list(i,1).name).Stats(3,1);
    NCP.SEM(2,i)= AllNCP.(list(i,1).name).Stats(3,2);
    NCP.SEM(3,i)= AllNCP.(list(i,1).name).Stats(3,3);
    NCP.SEM(4,i)= 0;
     
    Density.Mean(1,i)= AllDensity.(list(i,1).name).Stats(1,1);
    Density.Mean(2,i)= AllDensity.(list(i,1).name).Stats(1,2);
    Density.Mean(3,i)= AllDensity.(list(i,1).name).Stats(1,3);
    Density.Mean(4,i)= AllDensity.(list(i,1).name).Stats(1,4);
    
    Density.STD(1,i)= AllDensity.(list(i,1).name).Stats(2,1);
    Density.STD(2,i)= AllDensity.(list(i,1).name).Stats(2,2);
    Density.STD(3,i)= AllDensity.(list(i,1).name).Stats(2,3);
    Density.STD(4,i)= AllDensity.(list(i,1).name).Stats(2,4);
    
    Density.SEM(1,i)= AllDensity.(list(i,1).name).Stats(3,1);
    Density.SEM(2,i)= AllDensity.(list(i,1).name).Stats(3,2);
    Density.SEM(3,i)= AllDensity.(list(i,1).name).Stats(3,3);
    Density.SEM(4,i)= AllDensity.(list(i,1).name).Stats(3,4);
    
    Flux.Mean(1,i)= AllFlux.(list(i,1).name).Stats(1,1);
    Flux.Mean(2,i)= AllFlux.(list(i,1).name).Stats(1,2);
    Flux.Mean(3,i)= AllFlux.(list(i,1).name).Stats(1,3);
    Flux.Mean(4,i)= 0;
    
    Flux.STD(1,i)= AllFlux.(list(i,1).name).Stats(2,1);
    Flux.STD(2,i)= AllFlux.(list(i,1).name).Stats(2,2);
    Flux.STD(3,i)= AllFlux.(list(i,1).name).Stats(2,3);
    Flux.STD(4,i)= 0;
    
    Flux.SEM(1,i)= AllFlux.(list(i,1).name).Stats(3,1);
    Flux.SEM(2,i)= AllFlux.(list(i,1).name).Stats(3,2);
    Flux.SEM(3,i)= AllFlux.(list(i,1).name).Stats(3,3);
    Flux.SEM(4,i)= 0;
    
    SV.Mean(1,i)= AllSV.(list(i,1).name).Stats(1,1); 
    SV.Mean(2,i)= AllSV.(list(i,1).name).Stats(1,2); 
    SV.Mean(3,i)= 0;
    SV.Mean(4,i)= 0;
    
    SV.STD(1,i)= AllSV.(list(i,1).name).Stats(2,1);
    SV.STD(2,i)= AllSV.(list(i,1).name).Stats(2,2);
    SV.STD(3,i)= 0;
    SV.STD(4,i)= 0;
    
    SV.SEM(1,i)= AllSV.(list(i,1).name).Stats(3,1);
    SV.SEM(2,i)= AllSV.(list(i,1).name).Stats(3,2);
    SV.SEM(3,i)= 0;
    SV.SEM(4,i)= 0;
    
    NV.Mean(1,i)= AllNV.(list(i,1).name).Stats(1,1); 
    NV.Mean(2,i)= AllNV.(list(i,1).name).Stats(1,2); 
    NV.Mean(3,i)= 0;
    NV.Mean(4,i)= 0;
    
    NV.STD(1,i)= AllNV.(list(i,1).name).Stats(2,1);
    NV.STD(2,i)= AllNV.(list(i,1).name).Stats(2,2);
    NV.STD(3,i)= 0;
    NV.STD(4,i)= 0;
    
    NV.SEM(1,i)= AllNV.(list(i,1).name).Stats(3,1);
    NV.SEM(2,i)= AllNV.(list(i,1).name).Stats(3,2);
    NV.SEM(3,i)= 0;
    NV.SEM(4,i)= 0;
    
    RL.Mean(1,i)= AllRL.(list(i,1).name).Stats(1,1); 
    RL.Mean(2,i)= AllRL.(list(i,1).name).Stats(1,2); 
    RL.Mean(3,i)= 0;
    RL.Mean(4,i)= 0;
    
    RL.STD(1,i)= AllRL.(list(i,1).name).Stats(2,1);
    RL.STD(2,i)= AllRL.(list(i,1).name).Stats(2,2);
    RL.STD(3,i)= 0;
    RL.STD(4,i)= 0;
        
    RL.SEM(1,i)= AllRL.(list(i,1).name).Stats(3,1);
    RL.SEM(2,i)= AllRL.(list(i,1).name).Stats(3,2);
    RL.SEM(3,i)= 0;
    RL.SEM(4,i)= 0;
    
    comRL.Mean(1,i)= AllcomRL.(list(i,1).name).Stats(1,1); 
    comRL.Mean(2,i)= AllcomRL.(list(i,1).name).Stats(1,2); 
    comRL.Mean(3,i)= 0;
    comRL.Mean(4,i)= 0;
    
    comRL.STD(1,i)= AllcomRL.(list(i,1).name).Stats(2,1);
    comRL.STD(2,i)= AllcomRL.(list(i,1).name).Stats(2,2);
    comRL.STD(3,i)= 0;
    comRL.STD(4,i)= 0;
        
    comRL.SEM(1,i)= AllcomRL.(list(i,1).name).Stats(3,1);
    comRL.SEM(2,i)= AllcomRL.(list(i,1).name).Stats(3,2);
    comRL.SEM(3,i)= 0;
    comRL.SEM(4,i)= 0;
    
    comSV.Mean(1,i)= AllcomSV.(list(i,1).name).Stats(1,1); 
    comSV.Mean(2,i)= AllcomSV.(list(i,1).name).Stats(1,2); 
    comSV.Mean(3,i)= 0;
    comSV.Mean(4,i)= 0;
    
    comSV.STD(1,i)= AllcomSV.(list(i,1).name).Stats(2,1);
    comSV.STD(2,i)= AllcomSV.(list(i,1).name).Stats(2,2);
    comSV.STD(3,i)= 0;
    comSV.STD(4,i)= 0;
        
    comSV.SEM(1,i)= AllcomSV.(list(i,1).name).Stats(3,1);
    comSV.SEM(2,i)= AllcomSV.(list(i,1).name).Stats(3,2);
    comSV.SEM(3,i)= 0;
    comSV.SEM(4,i)= 0;
    
    PD.Mean(1,i)= AllPD.(list(i,1).name).Stats(1,1);  
    PD.Mean(2,i)= 0;
    PD.Mean(3,i)= 0;
    PD.Mean(4,i)= 0;
    
    PD.STD(1,i)= AllPD.(list(i,1).name).Stats(2,1);
    PD.STD(2,i)= 0;
    PD.STD(3,i)= 0;
    PD.STD(4,i)= 0;
    
    PD.SEM(1,i)= AllPD.(list(i,1).name).Stats(3,1);
    PD.SEM(2,i)= 0;
    PD.SEM(3,i)= 0;
    PD.SEM(4,i)= 0;
    
    PF.Mean(1,i)= AllPF.(list(i,1).name).Stats(1,1);  
    PF.Mean(2,i)= 0;
    PF.Mean(3,i)= 0;
    PF.Mean(4,i)= 0;
    
    PF.STD(1,i)= AllPF.(list(i,1).name).Stats(2,1);
    PF.STD(2,i)= 0;
    PF.STD(3,i)= 0;
    PF.STD(4,i)= 0;
    
    PF.SEM(1,i)= AllPF.(list(i,1).name).Stats(3,1);
    PF.SEM(2,i)= 0;
    PF.SEM(3,i)= 0;
    PF.SEM(4,i)= 0;
    
    PFperSec.Mean(1,i)= AllPFperSec.(list(i,1).name).Stats(1,1);  
    PFperSec.Mean(2,i)= 0;
    PFperSec.Mean(3,i)= 0;
    PFperSec.Mean(4,i)= 0;
    
    PFperSec.STD(1,i)= AllPFperSec.(list(i,1).name).Stats(2,1);
    PFperSec.STD(2,i)= 0;
    PFperSec.STD(3,i)= 0;
    PFperSec.STD(4,i)= 0;
    
    PFperSec.SEM(1,i)= AllPFperSec.(list(i,1).name).Stats(3,1);
    PFperSec.SEM(2,i)= 0;
    PFperSec.SEM(3,i)= 0;
    PFperSec.SEM(4,i)= 0;
    
    splitPD.Mean(1,i)= AllsplitPD.(list(i,1).name).Stats(1,1);  
    splitPD.Mean(2,i)= AllsplitPD.(list(i,1).name).Stats(1,2);
    splitPD.Mean(3,i)= AllsplitPD.(list(i,1).name).Stats(1,3);
    splitPD.Mean(4,i)= 0;
    
    splitPD.STD(1,i)= AllsplitPD.(list(i,1).name).Stats(2,1);
    splitPD.STD(2,i)= AllsplitPD.(list(i,1).name).Stats(2,2);
    splitPD.STD(3,i)= AllsplitPD.(list(i,1).name).Stats(2,3);
    splitPD.STD(4,i)= 0;
    
    splitPD.SEM(1,i)= AllsplitPD.(list(i,1).name).Stats(3,1);
    splitPD.SEM(2,i)= AllsplitPD.(list(i,1).name).Stats(3,2);
    splitPD.SEM(3,i)= AllsplitPD.(list(i,1).name).Stats(3,3);
    splitPD.SEM(4,i)= 0;
    
    splitPF.Mean(1,i)= AllsplitPF.(list(i,1).name).Stats(1,1);  
    splitPF.Mean(2,i)= AllsplitPF.(list(i,1).name).Stats(1,2);  
    splitPF.Mean(3,i)= AllsplitPF.(list(i,1).name).Stats(1,3);  
    splitPF.Mean(4,i)= 0;
    
    splitPF.STD(1,i)= AllsplitPF.(list(i,1).name).Stats(2,1);
    splitPF.STD(2,i)= AllsplitPF.(list(i,1).name).Stats(2,2);
    splitPF.STD(3,i)= AllsplitPF.(list(i,1).name).Stats(2,3);
    splitPF.STD(4,i)= 0;
    
    splitPF.SEM(1,i)= AllsplitPF.(list(i,1).name).Stats(3,1);
    splitPF.SEM(2,i)= AllsplitPF.(list(i,1).name).Stats(3,2);
    splitPF.SEM(3,i)= AllsplitPF.(list(i,1).name).Stats(3,3);
    splitPF.SEM(4,i)= 0;
    
    splitPFperSec.Mean(1,i)= AllsplitPFperSec.(list(i,1).name).Stats(1,1);  
    splitPFperSec.Mean(2,i)= AllsplitPFperSec.(list(i,1).name).Stats(1,2);  
    splitPFperSec.Mean(3,i)= AllsplitPFperSec.(list(i,1).name).Stats(1,3);  
    splitPFperSec.Mean(4,i)= 0;
    
    splitPFperSec.STD(1,i)= AllsplitPFperSec.(list(i,1).name).Stats(2,1);
    splitPFperSec.STD(2,i)= AllsplitPFperSec.(list(i,1).name).Stats(2,3);
    splitPFperSec.STD(3,i)= AllsplitPFperSec.(list(i,1).name).Stats(2,3);
    splitPFperSec.STD(4,i)= 0;
    
    splitPFperSec.SEM(1,i)= AllsplitPFperSec.(list(i,1).name).Stats(3,1);
    splitPFperSec.SEM(2,i)= AllsplitPFperSec.(list(i,1).name).Stats(3,2);
    splitPFperSec.SEM(3,i)= AllsplitPFperSec.(list(i,1).name).Stats(3,3);
    splitPFperSec.SEM(4,i)= 0;
    
    SF.Mean(1,i)= AllSF.(list(i,1).name).Stats(1,1);  
    SF.Mean(2,i)= 0;
    SF.Mean(3,i)= 0;
    SF.Mean(4,i)= 0;
    
    SF.STD(1,i)= AllSF.(list(i,1).name).Stats(2,1);
    SF.STD(2,i)= 0;
    SF.STD(3,i)= 0;
    SF.STD(4,i)= 0;
    
    SF.SEM(1,i)= AllSF.(list(i,1).name).Stats(3,1);
    SF.SEM(2,i)= 0;
    SF.SEM(3,i)= 0;
    SF.SEM(4,i)= 0;
    
    SFperSec.Mean(1,i)= AllSFperSec.(list(i,1).name).Stats(1,1);  
    SFperSec.Mean(2,i)= 0;
    SFperSec.Mean(3,i)= 0;
    SFperSec.Mean(4,i)= 0;
    
    SFperSec.STD(1,i)= AllSFperSec.(list(i,1).name).Stats(2,1);
    SFperSec.STD(2,i)= 0;
    SFperSec.STD(3,i)= 0;
    SFperSec.STD(4,i)= 0;
    
    SFperSec.SEM(1,i)= AllSFperSec.(list(i,1).name).Stats(3,1);
    SFperSec.SEM(2,i)= 0;
    SFperSec.SEM(3,i)= 0;
    SFperSec.SEM(4,i)= 0;
    
    revSF.Mean(1,i)= AllrevSF.(list(i,1).name).Stats(1,1);  
    revSF.Mean(2,i)= 0;
    revSF.Mean(3,i)= 0;
    revSF.Mean(4,i)= 0;
    
    revSF.STD(1,i)= AllrevSF.(list(i,1).name).Stats(2,1);
    revSF.STD(2,i)= 0;
    revSF.STD(3,i)= 0;
    revSF.STD(4,i)= 0;
    
    revSF.SEM(1,i)= AllrevSF.(list(i,1).name).Stats(3,1);
    revSF.SEM(2,i)= 0;
    revSF.SEM(3,i)= 0;
    revSF.SEM(4,i)= 0;
    
    revSFperSec.Mean(1,i)= AllrevSFperSec.(list(i,1).name).Stats(1,1);  
    revSFperSec.Mean(2,i)= 0;
    revSFperSec.Mean(3,i)= 0;
    revSFperSec.Mean(4,i)= 0;
    
    revSFperSec.STD(1,i)= AllrevSFperSec.(list(i,1).name).Stats(2,1);
    revSFperSec.STD(2,i)= 0;
    revSFperSec.STD(3,i)= 0;
    revSFperSec.STD(4,i)= 0;
    
    revSFperSec.SEM(1,i)= AllrevSFperSec.(list(i,1).name).Stats(3,1);
    revSFperSec.SEM(2,i)= 0;
    revSFperSec.SEM(3,i)= 0;
    revSFperSec.SEM(4,i)= 0;
end


barCP=[CP.Mean(1,:);CP.Mean(2,:);CP.Mean(3,:);CP.Mean(4,:)];
stdCP=[CP.STD(1,:);CP.STD(2,:);CP.STD(3,:);CP.STD(4,:)];
semCP=[CP.SEM(1,:);CP.SEM(2,:);CP.SEM(3,:);CP.SEM(4,:)];

barNCP=[NCP.Mean(1,:);NCP.Mean(2,:);NCP.Mean(3,:);NCP.Mean(4,:)];
stdNCP=[NCP.STD(1,:);NCP.STD(2,:);NCP.STD(3,:);NCP.STD(4,:)];
semNCP=[NCP.SEM(1,:);NCP.SEM(2,:);NCP.SEM(3,:);NCP.SEM(4,:)];

barDensity=[Density.Mean(1,:);Density.Mean(2,:);Density.Mean(3,:);Density.Mean(4,:)];
stdDensity=[Density.STD(1,:);Density.STD(2,:);Density.STD(3,:);Density.STD(4,:)];
semDensity=[Density.SEM(1,:);Density.SEM(2,:);Density.SEM(3,:);Density.SEM(4,:)];

barFlux=[Flux.Mean(1,:);Flux.Mean(2,:);Flux.Mean(3,:);Flux.Mean(4,:)];
stdFlux=[Flux.STD(1,:);Flux.STD(2,:);Flux.STD(3,:);Flux.STD(4,:)];
semFlux=[Flux.SEM(1,:);Flux.SEM(2,:);Flux.SEM(3,:);Flux.SEM(4,:)];

barSV=[SV.Mean(1,:);SV.Mean(2,:);SV.Mean(3,:);SV.Mean(4,:)];
stdSV=[SV.STD(1,:);SV.STD(2,:);SV.STD(3,:);SV.STD(4,:)];
semSV=[SV.SEM(1,:);SV.SEM(2,:);SV.SEM(3,:);SV.SEM(4,:)];

barNV=[NV.Mean(1,:);NV.Mean(2,:);NV.Mean(3,:);NV.Mean(4,:)];
stdNV=[NV.STD(1,:);NV.STD(2,:);NV.STD(3,:);NV.STD(4,:)];
semNV=[NV.SEM(1,:);NV.SEM(2,:);NV.SEM(3,:);NV.SEM(4,:)];

barRL=[RL.Mean(1,:);RL.Mean(2,:);RL.Mean(3,:);RL.Mean(4,:)];
stdRL=[RL.STD(1,:);RL.STD(2,:);RL.STD(3,:);RL.STD(4,:)];
semRL=[RL.SEM(1,:);RL.SEM(2,:);RL.SEM(3,:);RL.SEM(4,:)];

barcomRL=[comRL.Mean(1,:);comRL.Mean(2,:);comRL.Mean(3,:);comRL.Mean(4,:)];
stdcomRL=[comRL.STD(1,:);comRL.STD(2,:);comRL.STD(3,:);comRL.STD(4,:)];
semcomRL=[comRL.SEM(1,:);comRL.SEM(2,:);comRL.SEM(3,:);comRL.SEM(4,:)];

barcomSV=[comSV.Mean(1,:);comSV.Mean(2,:);comSV.Mean(3,:);comSV.Mean(4,:)];
stdcomSV=[comSV.STD(1,:);comSV.STD(2,:);comSV.STD(3,:);comSV.STD(4,:)];
semcomSV=[comSV.SEM(1,:);comSV.SEM(2,:);comSV.SEM(3,:);comSV.SEM(4,:)];

barPD=[PD.Mean(1,:);PD.Mean(2,:);PD.Mean(3,:);PD.Mean(4,:)];
stdPD=[PD.STD(1,:);PD.STD(2,:);PD.STD(3,:);PD.STD(4,:)];
semPD=[PD.SEM(1,:);PD.SEM(2,:);PD.SEM(3,:);PD.SEM(4,:)];

barPF=[PF.Mean(1,:);PF.Mean(2,:);PF.Mean(3,:);PF.Mean(4,:)];
stdPF=[PF.STD(1,:);PF.STD(2,:);PF.STD(3,:);PF.STD(4,:)];
semPF=[PF.SEM(1,:);PF.SEM(2,:);PF.SEM(3,:);PF.SEM(4,:)];

barPFperSec=[PFperSec.Mean(1,:);PFperSec.Mean(2,:);PFperSec.Mean(3,:);PFperSec.Mean(4,:)];
stdPFperSec=[PFperSec.STD(1,:);PFperSec.STD(2,:);PFperSec.STD(3,:);PFperSec.STD(4,:)];
semPFperSec=[PFperSec.SEM(1,:);PFperSec.SEM(2,:);PFperSec.SEM(3,:);PFperSec.SEM(4,:)];

barsplitPD=[splitPD.Mean(1,:);splitPD.Mean(2,:);splitPD.Mean(3,:);splitPD.Mean(4,:)];
stdsplitPD=[splitPD.STD(1,:);splitPD.STD(2,:);splitPD.STD(3,:);splitPD.STD(4,:)];
semsplitPD=[splitPD.SEM(1,:);splitPD.SEM(2,:);splitPD.SEM(3,:);splitPD.SEM(4,:)];

barsplitPF=[splitPF.Mean(1,:);splitPF.Mean(2,:);splitPF.Mean(3,:);splitPF.Mean(4,:)];
stdsplitPF=[splitPF.STD(1,:);splitPF.STD(2,:);splitPF.STD(3,:);splitPF.STD(4,:)];
semsplitPF=[splitPF.SEM(1,:);splitPF.SEM(2,:);splitPF.SEM(3,:);splitPF.SEM(4,:)];

barsplitPFperSec=[splitPFperSec.Mean(1,:);splitPFperSec.Mean(2,:);splitPFperSec.Mean(3,:);splitPFperSec.Mean(4,:)];
stdsplitPFperSec=[splitPFperSec.STD(1,:);splitPFperSec.STD(2,:);splitPFperSec.STD(3,:);splitPFperSec.STD(4,:)];
semsplitPFperSec=[splitPFperSec.SEM(1,:);splitPFperSec.SEM(2,:);splitPFperSec.SEM(3,:);splitPFperSec.SEM(4,:)];

barSF=[SF.Mean(1,:);SF.Mean(2,:);SF.Mean(3,:);SF.Mean(4,:)];
stdSF=[SF.STD(1,:);SF.STD(2,:);SF.STD(3,:);SF.STD(4,:)];
semSF=[SF.SEM(1,:);SF.SEM(2,:);SF.SEM(3,:);SF.SEM(4,:)];

barSFperSec=[SFperSec.Mean(1,:);SFperSec.Mean(2,:);SFperSec.Mean(3,:);SFperSec.Mean(4,:)];
stdSFperSec=[SFperSec.STD(1,:);SFperSec.STD(2,:);SFperSec.STD(3,:);SFperSec.STD(4,:)];
semSFperSec=[SFperSec.SEM(1,:);SFperSec.SEM(2,:);SFperSec.SEM(3,:);SFperSec.SEM(4,:)];

barrevSF=[revSF.Mean(1,:);revSF.Mean(2,:);revSF.Mean(3,:);revSF.Mean(4,:)];
stdrevSF=[revSF.STD(1,:);revSF.STD(2,:);revSF.STD(3,:);revSF.STD(4,:)];
semrevSF=[revSF.SEM(1,:);revSF.SEM(2,:);revSF.SEM(3,:);revSF.SEM(4,:)];

barrevSFperSec=[revSFperSec.Mean(1,:);revSFperSec.Mean(2,:);revSFperSec.Mean(3,:);revSFperSec.Mean(4,:)];
stdrevSFperSec=[revSFperSec.STD(1,:);revSFperSec.STD(2,:);revSFperSec.STD(3,:);revSFperSec.STD(4,:)];
semrevSFperSec=[revSFperSec.SEM(1,:);revSFperSec.SEM(2,:);revSFperSec.SEM(3,:);revSFperSec.SEM(4,:)];

%% Creating Figures for bar graphs with Standard Deviation
% Figure 1 is figure for all averages
% color is an array with color specifications for TEN datasets. If more
% than TEN datasets ought to be compared, then add more colorcodes to
% this array.

% Figure Positions
% [% of WIDTH,% of HEIGHT,,]
%
% ROW 1 - CP, Dens, velocity, RL
% CP    [0.03,0.8,0.08,0.16]
% NetCP [0.16,0.8,0.08,0.16]
% Dens  [0.29,0.8,0.08,0.16]

% Sv    [0.42,0.8,0.08,0.16]
% NV    [0.55,0.8,0.08,0.16]
% RL    [0.68,0.8,0.08,0.16]
% cRL   [0.81,0.8,0.08,0.16]
% ROW 2 - Pauses
% Pd    [0.03,0.55,0.08,0.16]
% sPd   [0.16,0.55,0.08,0.16]
% PF    [0.29,0.55,0.08,0.16]
% sPF   [0.42,0.55,0.08,0.16]
% PfpS  [0.55,0.55,0.08,0.16]
% sPFpS [0.68,0.55,0.08,0.16]
% ROW 3 - Switches
% SF    [0.03,0.3,0.08,0.16]
% SFpS  [0.16,0.3,0.08,0.16]
% SFr   [0.29,0.3,0.08,0.16]
% SFrpS [0.42,0.3,0.08,0.16]
% Flux  [0.55,0.3,0.125,0.16] 
% 

color=[0 0 0.6; 0 0.4 0.8; 0 0.4 0; 0 0.8 0; 0.6 0 0; 1 0.2 0.2; 1 0.69 0.4;1 0.89 0.8; 1 1 0; 1 1 0.8];
STD=zeros(4,length(list));

% write array for legend, underscores have to be replaced in genotype names, because they
% make next letter a subscript in matlab text
legendNames={};
for i=1:length(list)
    legendNames{i}=list(i,:).name;
    
    %replace underscores in genotype names
    findScore=strfind(legendNames{i},'_');
    if findScore>0
        replaceScore=strrep(legendNames{i},'_','/');
        legendNames{i}=replaceScore;
    end
end

lines=([]);
for i=1:length(list)
    lines=cat(1,lines,i);
end

figure1=figure;
set(gcf,'Position',[0 0 2000 1000]);

%% plots for Cargo Population
subplot('Position',[0.03,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barCP,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
ylim ([0 100]);
title('Cargo Population');
ylabel('Cargo Population (%)');
set(gca,'XTickLabel',{'Ante','Retro','Reversing','Stationary'});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barCP,STD,stdCP,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Net Cargo Population
subplot('Position',[0.16,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barNCP,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
ylim ([0 100]);
title('Net Cargo Population');
ylabel('Net Cargo Population (%)');
set(gca,'XTickLabel',{'Ante','Retro','Stationary',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barNCP,STD,stdNCP,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Density
subplot('Position',[0.29,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barDensity,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Density');
ylabel('Density Analysis (Tracks / um axon)');
set(gca,'XTickLabel',{'Ante','Retro','Reversing','Stationary'});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barDensity,STD,stdDensity,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Segmental Velocities
subplot('Position',[0.42,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barSV,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Segmental Velocities');
ylabel('Segmental Velocities (um/sec)');
set(gca,'XTickLabel',{'Ante','Retro','',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barSV,STD,stdSV,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Net Velocities
subplot('Position',[0.55,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barNV,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Net Velocities');
ylabel('Net Velocities (um/sec)');
set(gca,'XTickLabel',{'Ante','Retro','',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barNV,STD,stdNV,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Run Length
subplot('Position',[0.68,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barRL,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Segmental Run Length');
ylabel('Segmental Run Length (um)');
set(gca,'XTickLabel',{'Ante','Retro','',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barRL,STD,stdRL,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for combined Run Length
subplot('Position',[0.81,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barcomRL,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Segmental Run Length (combined)');
ylabel('Segmental Run Length (um)');
set(gca,'XTickLabel',{'Ante','Retro','',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barcomRL,STD,stdcomRL,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Pause Duration
subplot('Position',[0.03,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barPD,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Pause Duration');
ylabel('Pause Duration (sec)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barPD,STD,stdPD,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Split Pause Duration
subplot('Position',[0.16,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barsplitPD,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Split Pause Duration');
ylabel('Pause Duration (sec)');
set(gca,'XTickLabel',{'Ante','Retro','Reversing',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barsplitPD,STD,stdsplitPD,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Pause Frequency
subplot('Position',[0.29,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barPF,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Pause Frequency');
ylabel('Pause Frequency (Pauses/Track)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barPF,STD,stdPF,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Split Pause Frequency
subplot('Position',[0.42,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barsplitPF,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Split Pause Frequency');
ylabel('Pause Frequency (Pauses/Track)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barsplitPF,STD,stdsplitPF,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Pause Frequency per Second
subplot('Position',[0.55,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barPFperSec,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Pause Frequency per Second');
ylabel('Pause Frequency (Pauses/sec)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barPFperSec,STD,stdPFperSec,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for split Pause Frequency per Second
subplot('Position',[0.68,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barsplitPFperSec,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Split Pause Frequency per Second');
ylabel('Pause Frequency (Pauses/sec)');
set(gca,'XTickLabel',{'Ante','Retro','Reversing',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barsplitPFperSec,STD,stdsplitPFperSec,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for combined Segmental Velocity
subplot('Position',[0.81,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barcomSV,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Segmental Velocity (combined)');
ylabel('Segmental Velocity (um/sec)');
set(gca,'XTickLabel',{'Ante','Retro','',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barcomSV,STD,stdcomSV,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Switch Frequency
subplot('Position',[0.03,0.3,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barSF,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Switch Frequency');
ylabel('Switch Frequency (Switches/Track)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barSF,STD,stdSF,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Switch Frequency per Sec
subplot('Position',[0.16,0.3,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barSFperSec,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Switch Frequency per Second');
ylabel('Switch Frequency (Switches/sec)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barSFperSec,STD,stdSFperSec,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Switch Frequency Reversals
subplot('Position',[0.29,0.3,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barrevSF,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Switch Frequency Reversals');
ylabel('Switch Frequency (Switches/Track)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barrevSF,STD,stdrevSF,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Switch Frequency Reversals per Sec
subplot('Position',[0.42,0.3,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barrevSFperSec,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Switch Frequency per Second Reversals');
ylabel('Switch Frequency (Switches/sec)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barrevSFperSec,STD,stdrevSFperSec,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Flux
subplot('Position',[0.55,0.3,0.125,0.16]); 
set(gca,'fontsize',8);
hb=bar(barFlux,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Flux');
ylabel('Flux Analysis (Tracks / (um axon * sec))');
set(gca,'XTickLabel',{'Ante','Retro','Reversing',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barFlux,STD,stdFlux,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% legend 
legend(legendNames{:},'Location','EastOutside');

% saving
fileFigure=fullfile(directoryMain,'Figures');
fileFigure=strcat(fileFigure,filesep,'AveragesSTD.fig');
saveas(figure1,fileFigure);
fileFigure=fullfile(directoryMain,'Figures');
fileFigure=strcat(fileFigure,filesep,'AveragesSTD.eps');
saveas(figure1,fileFigure);

%% Creating Figures for bar graphs with Standard Error of Mean (SEM)
% Figure 2 is figure for all averages
% color is an array with color specifications for TEN datasets. If more
% than TEN datasets ought to be compared, then add more colorcodes to
% this array.
color=[0 0 0.6; 0 0.4 0.8; 0 0.4 0; 0 0.8 0; 0.6 0 0; 1 0.2 0.2; 1 0.69 0.4;1 0.89 0.8; 1 1 0; 1 1 0.8];
STD=zeros(4,length(list));

% write array for legend, underscores have to be replaced, because they
% make next letter a subscript in matlab

legendNames={};
for i=1:length(list)
    legendNames{i}=list(i,:).name;
    
    %replace underscores in genotype names
    findScore=strfind(legendNames{i},'_');
    if findScore>0
        replaceScore=strrep(legendNames{i},'_','/');
        legendNames{i}=replaceScore;
    end
end

lines=([]);
for i=1:length(list)
    lines=cat(1,lines,i);
end

figure2=figure;
set(gcf,'Position',[0 0 2000 1000]);
%% plots for Cargo Population
subplot('Position',[0.03,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barCP,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
ylim ([0 100]);
title('Cargo Population');
ylabel('Cargo Population (%)');
set(gca,'XTickLabel',{'Ante','Retro','Reversing','Stationary'});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barCP,STD,semCP,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Net Cargo Population
subplot('Position',[0.16,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barNCP,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
ylim ([0 100]);
title('Net Cargo Population');
ylabel('Net Cargo Population (%)');
set(gca,'XTickLabel',{'Ante','Retro','Stationary',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barNCP,STD,semNCP,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Density
subplot('Position',[0.29,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barDensity,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Density');
ylabel('Density Analysis (Tracks / um axon)');
set(gca,'XTickLabel',{'Ante','Retro','Reversing','Stationary'});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barDensity,STD,semDensity,'k.');
set(h(lines),'MarkerEdgeColor','none');



%% plots for Segmental Velocities
subplot('Position',[0.42,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barSV,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Segmental Velocities');
ylabel('Segmental Velocities (um/sec)');
set(gca,'XTickLabel',{'Ante','Retro','',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barSV,STD,semSV,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Net Velocities
subplot('Position',[0.55,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barNV,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Net Velocities');
ylabel('Net Velocities (um/sec)');
set(gca,'XTickLabel',{'Ante','Retro','',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barNV,STD,semNV,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Run Length
subplot('Position',[0.68,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barRL,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Segmental Run Length');
ylabel('Segmental Run Length (um)');
set(gca,'XTickLabel',{'Ante','Retro','',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barRL,STD,semRL,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Run Length combined
subplot('Position',[0.81,0.8,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barcomRL,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Segmental Run Length combined');
ylabel('Segmental Run Length (um)');
set(gca,'XTickLabel',{'Ante','Retro','',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barcomRL,STD,semcomRL,'k.');
set(h(lines),'MarkerEdgeColor','none');


%% plots for Pause Duration
subplot('Position',[0.03,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barPD,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Pause Duration');
ylabel('Pause Duration (sec)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barPD,STD,semPD,'k.');
set(h(lines),'MarkerEdgeColor','none');

% plots for Split Pause Duration
subplot('Position',[0.16,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barsplitPD,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Split Pause Duration');
ylabel('Pause Duration (sec)');
set(gca,'XTickLabel',{'Ante','Retro','Reversing',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barsplitPD,STD,semsplitPD,'k.');
set(h(lines),'MarkerEdgeColor','none');

% plots for Pause Frequency
subplot('Position',[0.29,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barPF,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Pause Frequency');
ylabel('Pause Frequency (Pauses/Track)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barPF,STD,semPF,'k.');
set(h(lines),'MarkerEdgeColor','none');

% plots for Split Pause Frequency
subplot('Position',[0.42,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barsplitPF,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Split Pause Frequency');
ylabel('Pause Frequency (Pauses/Track)');
set(gca,'XTickLabel',{'Ante','Retro','Reversing',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barsplitPF,STD,semsplitPF,'k.');
set(h(lines),'MarkerEdgeColor','none');

% plots for Pause Frequency per Second
subplot('Position',[0.55,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barPFperSec,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Pause Frequency per Second');
ylabel('Pause Frequency (Pauses/sec)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barPFperSec,STD,semPFperSec,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Split Pause Frequency per Second
subplot('Position',[0.68,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barsplitPFperSec,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Split Pause Frequency per Second');
ylabel('Pause Frequency (Pauses/sec)');
set(gca,'XTickLabel',{'Ante','Retro','Reversing',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barsplitPFperSec,STD,semsplitPFperSec,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Segmental Velocity combined
subplot('Position',[0.81,0.55,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barcomSV,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Segmental Velocity combined');
ylabel('Segmental Velocity (um/sec)');
set(gca,'XTickLabel',{'Ante','Retro','',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barcomSV,STD,semcomSV,'k.');
set(h(lines),'MarkerEdgeColor','none');


%% plots for Switch Frequency
subplot('Position',[0.03,0.3,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barSF,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Switch Frequency');
ylabel('Switch Frequency (Switches/Track)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barSF,STD,semSF,'k.');
set(h(lines),'MarkerEdgeColor','none');

% plots for Switch Frequency per Sec
subplot('Position',[0.16,0.3,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barSFperSec,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Switch Frequency per Second');
ylabel('Switch Frequency (Switches/sec)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barSFperSec,STD,semSFperSec,'k.');
set(h(lines),'MarkerEdgeColor','none');

% plots for Switch Frequency Reversals
subplot('Position',[0.29,0.3,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barrevSF,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Switch Frequency Reversals');
ylabel('Switch Frequency (Switches/Track)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barrevSF,STD,semrevSF,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Switch Frequency per Sec Reversals
subplot('Position',[0.42,0.3,0.08,0.16]);
set(gca,'fontsize',8);
hb=bar(barrevSFperSec,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Switch Frequency per Second Reversals');
ylabel('Switch Frequency (Switches/sec)');
set(gca,'XTickLabel','');
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barrevSFperSec,STD,semrevSFperSec,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% plots for Flux
subplot('Position',[0.55,0.3,0.125,0.16]);
set(gca,'fontsize',8);
hb=bar(barFlux,'BarWidth',1.0);
for c=1:length(list)
    set(hb(c),'FaceColor',color(c,:));
end
title('Flux');
ylabel('Flux Analysis (Tracks / (um axon * sec))');
set(gca,'XTickLabel',{'Ante','Retro','Reversing',''});
rotateXLabels(gca,45);
set(gca, 'Ticklength', [0 0])
hold on
barCenters=([]);
for i=1:length(hb)
    xData=get(hb(i),'XData');
    barCenters(i,:)=(xData - 0.375 + i*0.25);  
end
barCentersT=transpose(barCenters);
h=errorbar(barCentersT,barFlux,STD,semFlux,'k.');
set(h(lines),'MarkerEdgeColor','none');

%% legend
legend(legendNames{:},'Location','EastOutside');

% saving
fileFigure=fullfile(directoryMain,'Figures');
fileFigure=strcat(fileFigure,filesep,'AveragesSEM.fig');
saveas(figure2,fileFigure);
fileFigure=fullfile(directoryMain,'Figures');
fileFigure=strcat(fileFigure,filesep,'AveragesSEM.eps');
saveas(figure2,fileFigure);

%% %add com sv? 





%% Finding the max x values for histograms and max y values for histograms
% determine max for x axis
maxSV=([]);
maxNV=([]);
maxRL=([]);
maxcomRL=([]);
for i=1:length(list)
    maxSV.ante(i,1)=max(AllSV.(list(i,1).name).ante(:,1));
    maxSV.retro(i,1)=max(AllSV.(list(i,1).name).retro(:,1));
    maxNV.ante(i,1)=max(AllNV.(list(i,1).name).ante(:,1));
    maxNV.retro(i,1)=max(AllNV.(list(i,1).name).retro(:,1));
    maxRL.ante(i,1)=max(AllRL.(list(i,1).name).ante(:,1));
    maxRL.retro(i,1)=max(AllRL.(list(i,1).name).retro(:,1));
    maxcomRL.ante(i,1)=max(AllcomRL.(list(i,1).name).ante(:,1));
    maxcomRL.retro(i,1)=max(AllcomRL.(list(i,1).name).retro(:,1));
end

maxSVante=max(maxSV.ante);
maxSVretro=max(maxSV.retro);
maxNVante=max(maxNV.ante);
maxNVretro=max(maxNV.retro);
maxRLante=max(maxRL.ante);
maxRLretro=max(maxRL.retro);
maxcomRLante=max(maxcomRL.ante);
maxcomRLretro=max(maxcomRL.retro);

binSize=([0.05, 0.1, 0.2, 0.3, 0.4, 0.5]);
binSizeName={'a','b','c','d','e','f'};

histSV=([]);
histNV=([]);
histRL=([]);
histcomRL=([]);
MultiGaussian=([]);
for m=1:length(binSize)
    %set bins
    %hist(x,xbins) sorts x into bins with intervals or categories determined by the vector xbins (here called xVector).
    %If xbins is a vector of evenly spaced values, then hist uses the values as the bin centers.
    xVectorSVante.(binSizeName{m})= binSize(m)/2:binSize(m):round(maxSVante)+1;
    xVectorSVretro.(binSizeName{m})= binSize(m)/2:binSize(m):round(maxSVretro)+1;
    xVectorNVante.(binSizeName{m})= binSize(m)/2:binSize(m):round(maxNVante)+1;
    xVectorNVretro.(binSizeName{m})= binSize(m)/2:binSize(m):round(maxNVretro)+1;
    xVectorRLante.(binSizeName{m})= binSize(m)/2:binSize(m):round(maxRLante)+1;
    xVectorRLretro.(binSizeName{m})= binSize(m)/2:binSize(m):round(maxRLretro)+1;
    xVectorcomRLante.(binSizeName{m})= binSize(m)/2:binSize(m):round(maxcomRLante)+1;
    xVectorcomRLretro.(binSizeName{m})= binSize(m)/2:binSize(m):round(maxcomRLretro)+1;

    % hist and determine max for y axis
    for i=1:length(list)
        histSV.(list(i,1).name).data.Ante.(binSizeName{m})=hist(AllSV.(list(i,1).name).ante,xVectorSVante.(binSizeName{m}));
        FRCT=rdivide(histSV.(list(i,1).name).data.Ante.(binSizeName{m}),sum(histSV.(list(i,1).name).data.Ante.(binSizeName{m})));
        histSV.(list(i,1).name).pct.Ante.(binSizeName{m})=mtimes(100,FRCT);
        ymaxSV.ante(i,1)=max(histSV.(list(i,1).name).pct.Ante.(binSizeName{m}));

        histSV.(list(i,1).name).data.Retro.(binSizeName{m})=hist(AllSV.(list(i,1).name).retro,xVectorSVretro.(binSizeName{m}));
        FRCT=rdivide(histSV.(list(i,1).name).data.Retro.(binSizeName{m}),sum(histSV.(list(i,1).name).data.Retro.(binSizeName{m})));
        histSV.(list(i,1).name).pct.Retro.(binSizeName{m})=mtimes(100,FRCT);
        ymaxSV.retro(i,1)=max(histSV.(list(i,1).name).pct.Retro.(binSizeName{m}));

        histNV.(list(i,1).name).data.Ante.(binSizeName{m})=hist(AllNV.(list(i,1).name).ante,xVectorNVante.(binSizeName{m}));
        FRCT=rdivide(histNV.(list(i,1).name).data.Ante.(binSizeName{m}),sum(histNV.(list(i,1).name).data.Ante.(binSizeName{m})));
        histNV.(list(i,1).name).pct.Ante.(binSizeName{m})=mtimes(100,FRCT);
        ymaxNV.ante(i,1)=max(histNV.(list(i,1).name).pct.Ante.(binSizeName{m}));

        histNV.(list(i,1).name).data.Retro.(binSizeName{m})=hist(AllNV.(list(i,1).name).retro,xVectorNVretro.(binSizeName{m}));
        FRCT=rdivide(histNV.(list(i,1).name).data.Retro.(binSizeName{m}),sum(histNV.(list(i,1).name).data.Retro.(binSizeName{m})));
        histNV.(list(i,1).name).pct.Retro.(binSizeName{m})=mtimes(100,FRCT);
        ymaxNV.retro(i,1)=max(histNV.(list(i,1).name).pct.Retro.(binSizeName{m}));

        histRL.(list(i,1).name).data.Ante.(binSizeName{m})=hist(AllRL.(list(i,1).name).ante,xVectorRLante.(binSizeName{m}));
        FRCT=rdivide(histRL.(list(i,1).name).data.Ante.(binSizeName{m}),sum(histRL.(list(i,1).name).data.Ante.(binSizeName{m})));
        histRL.(list(i,1).name).pct.Ante.(binSizeName{m})=mtimes(100,FRCT);
        ymaxRL.ante(i,1)=max(histRL.(list(i,1).name).pct.Ante.(binSizeName{m}));

        histRL.(list(i,1).name).data.Retro.(binSizeName{m})=hist(AllRL.(list(i,1).name).retro,xVectorRLretro.(binSizeName{m}));
        FRCT=rdivide(histRL.(list(i,1).name).data.Retro.(binSizeName{m}),sum(histRL.(list(i,1).name).data.Retro.(binSizeName{m})));
        histRL.(list(i,1).name).pct.Retro.(binSizeName{m})=mtimes(100,FRCT);
        ymaxRL.retro(i,1)=max(histRL.(list(i,1).name).pct.Retro.(binSizeName{m}));

        histcomRL.(list(i,1).name).data.Ante.(binSizeName{m})=hist(AllcomRL.(list(i,1).name).ante,xVectorcomRLante.(binSizeName{m}));
        FRCT=rdivide(histcomRL.(list(i,1).name).data.Ante.(binSizeName{m}),sum(histcomRL.(list(i,1).name).data.Ante.(binSizeName{m})));
        histcomRL.(list(i,1).name).pct.Ante.(binSizeName{m})=mtimes(100,FRCT);
        ymaxcomRL.ante(i,1)=max(histcomRL.(list(i,1).name).pct.Ante.(binSizeName{m}));

        histcomRL.(list(i,1).name).data.Retro.(binSizeName{m})=hist(AllcomRL.(list(i,1).name).retro,xVectorcomRLretro.(binSizeName{m}));
        FRCT=rdivide(histcomRL.(list(i,1).name).data.Retro.(binSizeName{m}),sum(histcomRL.(list(i,1).name).data.Retro.(binSizeName{m})));
        histcomRL.(list(i,1).name).pct.Retro.(binSizeName{m})=mtimes(100,FRCT);
        ymaxcomRL.retro(i,1)=max(histcomRL.(list(i,1).name).pct.Retro.(binSizeName{m}));

    end

    ymaxSVante=max(ymaxSV.ante);
    ymaxSVretro=max(ymaxSV.retro);
    ymaxNVante=max(ymaxNV.ante);
    ymaxNVretro=max(ymaxNV.retro);
    ymaxRLante=max(ymaxRL.ante);
    ymaxRLretro=max(ymaxRL.retro);
    ymaxcomRLante=max(ymaxcomRL.ante);
    ymaxcomRLretro=max(ymaxcomRL.retro);

% %% Creating Figures for histograms
% % figure 3 is for all histograms
%     one_over_sqrt_2pi = 1 / sqrt(2 * pi);
% 
%     for i=1:length(list)
%     figure3=figure;
%     set(gcf,'Position',[600 100 1200 600]);
% 
%     %% segmental velocities
%     % anterograde
% 
%     % histSV.(list(i,1).name)=([]);
%     % histSV.(list(i,1).name).data.Ante=hist(AllSV.(list(i,1).name).ante(:,1),xVector);
%     % FRCT=rdivide(histSV.(list(i,1).name).data.Ante,sum(histSV.(list(i,1).name).data.Ante));
%     % histSV.(list(i,1).name).pct.Ante=mtimes(100,FRCT);
% 
%     subplot('Position',[0.05,0.6,0.15,0.3]);
%     set(gca,'fontsize',8);
%     hb=bar(xVectorSVante.(binSizeName{m}),histSV.(list(i,1).name).pct.Ante.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
%     hold on
%     title('Anterograde Segmental Velocities');
%     xlim([-1 round(maxSVante)+1]);
%     if exist('ymaxSVanteNew') == 1
%         ymaxSVante=ymaxSVanteNew;
%     end
%     ylim([0 round(ymaxSVante)+10]);
%     xlabel('Segmental Velocities (um/sec)');
%     ylabel('Frequency (%)');
% 
%     % Gaussian Fitting for anterograde segmental velocities
%     data=AllSV.(list(i,1).name).ante;
%     %run function multigaussian fitting
%     try
%     multigauss_fitting_110414(data);
%     % rewrite variables into arrays for later saving
%     MultiGaussian.(list(i,1).name).anteSV.mu=mu; 
%     MultiGaussian.(list(i,1).name).anteSV.sigma=sigma;
%     MultiGaussian.(list(i,1).name).anteSV.sigmasq=sigmasq;
%     MultiGaussian.(list(i,1).name).anteSV.prob=prob;
%     MultiGaussian.(list(i,1).name).anteSV.clusterNum=clusterNum;
%     MultiGaussian.(list(i,1).name).anteSV.clusterBIC=BIC;
%     MultiGaussian.(list(i,1).name).anteSV.clusterClassification=classification;
% 
%     % calculate distribution
%     xGauss = -1:0.1:round(maxSVante)+1;
%     yGauss=([]);
%     for j = 1 : length(xGauss);
%         yGauss(j) = 0;
%         for k = 1 : clusterNum
%             temp = - 0.5 * (xGauss(j) - mu(k))^2 / sigmasq(k);
%             yGauss(j) = yGauss(j) + prob(k) * one_over_sqrt_2pi / sqrt(sigmasq(k)) * exp(temp); 
%         end
%     end
%     yGaussMax=max(yGauss);
%     MultiGaussian.(list(i,1).name).anteSV.xGauss=xGauss;
%     MultiGaussian.(list(i,1).name).anteSV.yGauss.(binSizeName{m})=yGauss/yGaussMax*max(histSV.(list(i,1).name).pct.Ante.(binSizeName{m}));
%     % plot by maximum matching
%     plot(xGauss,yGauss/yGaussMax*max(histSV.(list(i,1).name).pct.Ante.(binSizeName{m})),'r-');
%     hold on
%     %calculate separate clusters
%     for j=1:clusterNum
%         number=num2str(j);
%         clustNum=strcat('Cluster',number);
%         for k=1:length(xGauss);
%             temp = - 0.5 * (xGauss(k) - mu(j))^2 / sigmasq(j);
%             yGauss(k) = prob(j) * one_over_sqrt_2pi / sqrt(sigmasq(j)) * exp(temp); 
%         end
%         MultiGaussian.(list(i,1).name).anteSV.seperatedCluster.(clustNum).xGauss=xGauss;
%         MultiGaussian.(list(i,1).name).anteSV.seperatedCluster.(clustNum).yGauss.(binSizeName{m})=yGauss/yGaussMax*max(histSV.(list(i,1).name).pct.Ante.(binSizeName{m}));
%     end
% 
%     if yGaussMax>ymaxSVante
%         ymaxSVanteNew=yGaussMax;
%         ylim([0 ymaxSVanteNew]);
%     end
% 
%     %plot single clusters by maximum matching
%     for j=1:clusterNum
%         number=num2str(j);
%         clustNum=strcat('Cluster',number);
%         plot(xGauss,MultiGaussian.(list(i,1).name).anteSV.seperatedCluster.(clustNum).yGauss.(binSizeName{m}),'c-');
%     end
% 
%     catch ME
%         MultiGaussian.(list(i,1).name).anteSV=NaN;
%     end
% 
%     % retrograde
%     % histSV.(list(i,1).name).data.Retro=hist(AllSV.(list(i,1).name).retro(:,1),xVector);
%     % FRCT=rdivide(histSV.(list(i,1).name).data.Retro,sum(histSV.(list(i,1).name).data.Retro));
%     % histSV.(list(i,1).name).pct.Retro=mtimes(100,FRCT);
% 
%     subplot('Position',[0.3,0.6,0.15,0.3]);
%     set(gca,'fontsize',8);
%     hb=bar(xVectorSVretro.(binSizeName{m}),histSV.(list(i,1).name).pct.Retro.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
%     hold on
%     title('Retrograde Segmental Velocities');
%     xlim([-1 round(maxSVretro)+1]);
%     if exist('ymaxSVretroNew')==1
%         ymaxSVretro=ymaxSVretroNew;
%     end
%     ylim([0 round(ymaxSVretro)+10]);
%     xlabel('Segmental Velocities (um/sec)');
%     ylabel('Frequency (%)');
% 
%     % Gaussian Fitting for retrograde segmental velocities
%     data=AllSV.(list(i,1).name).retro;
%     %run function multigaussian fitting
%     try
%     multigauss_fitting_110414(data);
%     % rewrite variables into arrays for later saving
%     MultiGaussian.(list(i,1).name).retroSV.mu=mu; 
%     MultiGaussian.(list(i,1).name).retroSV.sigma=sigma;
%     MultiGaussian.(list(i,1).name).retroSV.sigmasq=sigmasq;
%     MultiGaussian.(list(i,1).name).retroSV.prob=prob;
%     MultiGaussian.(list(i,1).name).retroSV.clusterNum=clusterNum;
%     MultiGaussian.(list(i,1).name).retroSV.clusterBIC=BIC;
%     MultiGaussian.(list(i,1).name).retroSV.clusterClassification=classification;
% 
%     % calculate distribution
%     xGauss= -1:0.1:round(maxSVretro)+1;
%     yGauss=([]);
%     for j = 1 : length(xGauss);
%         yGauss(j) = 0;
%         for k = 1 : clusterNum
%             temp = - 0.5 * (xGauss(j) - mu(k))^2 / sigmasq(k);
%             yGauss(j) = yGauss(j) + prob(k) * one_over_sqrt_2pi / sqrt(sigmasq(k)) * exp(temp); 
%         end
%     end
%     yGaussMax=max(yGauss);
%     MultiGaussian.(list(i,1).name).retroSV.xGauss=xGauss;
%     MultiGaussian.(list(i,1).name).retroSV.yGauss.(binSizeName{m})=yGauss/yGaussMax*max(histSV.(list(i,1).name).pct.Retro.(binSizeName{m}));
%     % plot by maximum matching
%     plot(xGauss,yGauss/yGaussMax*max(histSV.(list(i,1).name).pct.Retro.(binSizeName{m})),'r-');
%     hold on
%     %calculate separate clusters
%     for j=1:clusterNum
%         number=num2str(j);
%         clustNum=strcat('Cluster',number);
%         for k=1:length(xGauss);
%             temp = - 0.5 * (xGauss(k) - mu(j))^2 / sigmasq(j);
%             yGauss(k) = prob(j) * one_over_sqrt_2pi / sqrt(sigmasq(j)) * exp(temp); 
%         end
%         MultiGaussian.(list(i,1).name).retroSV.seperatedCluster.(clustNum).xGauss=xGauss;
%         MultiGaussian.(list(i,1).name).retroSV.seperatedCluster.(clustNum).yGauss.(binSizeName{m})=yGauss/yGaussMax*max(histSV.(list(i,1).name).pct.Retro.(binSizeName{m}));
%     end
% 
%     %plot single clusters by maximum matching
%     for j=1:clusterNum
%         number=num2str(j);
%         clustNum=strcat('Cluster',number);
%         plot(xGauss,MultiGaussian.(list(i,1).name).retroSV.seperatedCluster.(clustNum).yGauss.(binSizeName{m}),'c-');
%     end
% 
%     catch ME
%         MultiGaussian.(list(i,1).name).retroSV=NaN;
%     end
% 
%     if yGaussMax>ymaxSVretro
%         ymaxSVretroNew=yGaussMax;
%         ylim([0 ymaxSVretroNew]);
%     end
% 
%     % anterograde net velocities
%     % histNV.(list(i,1).name)=([]);
%     % histNV.(list(i,1).name).data.Ante=hist(AllNV.(list(i,1).name).ante(:,1),xVector);
%     % FRCT=rdivide(histNV.(list(i,1).name).data.Ante,sum(histNV.(list(i,1).name).data.Ante));
%     % histNV.(list(i,1).name).pct.Ante=mtimes(100,FRCT);
% 
%     subplot('Position',[0.55,0.6,0.15,0.3]);
%     set(gca,'fontsize',8);
%     hb=bar(xVectorNVante.(binSizeName{m}),histNV.(list(i,1).name).pct.Ante.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
%     hold on
%     title('Anterograde Net Velocities');
%     xlim([-1 round(maxNVante)+1]);
% 
%     if exist('ymaxNVanteNew') == 1
%         ymaxNVante=ymaxNVanteNew;
%     end
%     ylim([0 round(ymaxNVante)+10]);
%     xlabel('Net Velocities (um/sec)');
%     ylabel('Frequency (%)');
% 
%     % Gaussian Fitting for anterograde segmental velocities
%     data=AllNV.(list(i,1).name).ante;
%     %run function multigaussian fitting
%     try
%     multigauss_fitting_110414(data);
%     % rewrite variables into arrays for later saving
%     MultiGaussian.(list(i,1).name).anteNV.mu=mu; 
%     MultiGaussian.(list(i,1).name).anteNV.sigma=sigma;
%     MultiGaussian.(list(i,1).name).anteNV.sigmasq=sigmasq;
%     MultiGaussian.(list(i,1).name).anteNV.prob=prob;
%     MultiGaussian.(list(i,1).name).anteNV.clusterNum=clusterNum;
%     MultiGaussian.(list(i,1).name).anteNV.clusterBIC=BIC;
%     MultiGaussian.(list(i,1).name).anteNV.clusterClassification=classification;
% 
%     % calculate distribution
%     xGauss= -1:0.1:round(maxNVante)+1;
%     yGauss=([]);
%     for j = 1 : length(xGauss);
%         yGauss(j) = 0;
%         for k = 1 : clusterNum
%             temp = - 0.5 * (xGauss(j) - mu(k))^2 / sigmasq(k);
%             yGauss(j) = yGauss(j) + prob(k) * one_over_sqrt_2pi / sqrt(sigmasq(k)) * exp(temp); 
%         end
%     end
%     yGaussMax=max(yGauss);
%     MultiGaussian.(list(i,1).name).anteNV.xGauss=xGauss;
%     MultiGaussian.(list(i,1).name).anteNV.yGauss.(binSizeName{m})=yGauss/yGaussMax*max(histNV.(list(i,1).name).pct.Ante.(binSizeName{m}));
%     % plot by maximum matching
%     plot(xGauss,yGauss/yGaussMax*max(histNV.(list(i,1).name).pct.Ante.(binSizeName{m})),'r-');
%     hold on
%     %calculate separate clusters
%     for j=1:clusterNum
%         number=num2str(j);
%         clustNum=strcat('Cluster',number);
%         for k=1:length(xGauss);
%             temp = - 0.5 * (xGauss(k) - mu(j))^2 / sigmasq(j);
%             yGauss(k) = prob(j) * one_over_sqrt_2pi / sqrt(sigmasq(j)) * exp(temp); 
%         end
%         MultiGaussian.(list(i,1).name).anteNV.seperatedCluster.(clustNum).xGauss=xGauss;
%         MultiGaussian.(list(i,1).name).anteNV.seperatedCluster.(clustNum).yGauss.(binSizeName{m})=yGauss/yGaussMax*max(histNV.(list(i,1).name).pct.Ante.(binSizeName{m}));
%     end
% 
%     %plot single clusters by maximum matching
%     for j=1:clusterNum
%         number=num2str(j);
%         clustNum=strcat('Cluster',number);
%         plot(xGauss,MultiGaussian.(list(i,1).name).anteNV.seperatedCluster.(clustNum).yGauss.(binSizeName{m}),'c-');
%     end
% 
%     catch ME
%         MultiGaussian.(list(i,1).name).anteNV=NaN;
%     end
% 
%     if yGaussMax>ymaxNVante
%         ymaxNVanteNew=yGaussMax;
%         ylim([0 ymaxNVanteNew]);
%     end
% 
%     %retrograde net velocities
%     % histNV.(list(i,1).name).data.Retro=hist(AllNV.(list(i,1).name).retro(:,1),xVector);
%     % FRCT=rdivide(histNV.(list(i,1).name).data.Retro,sum(histNV.(list(i,1).name).data.Retro));
%     % histNV.(list(i,1).name).pct.Retro=mtimes(100,FRCT);
% 
%     subplot('Position',[0.8,0.6,0.15,0.3]);
%     set(gca,'fontsize',8);
%     hb=bar(xVectorNVretro.(binSizeName{m}),histNV.(list(i,1).name).pct.Retro.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
%     hold on
%     title('Retrograde Net Velocities');
%     xlim([-1 round(maxNVretro)+1]);
% 
%     if exist('ymaxSVretroNew') == 1
%         ymaxNVretro=ymaxNVretroNew;
%     end
%     ylim([0 round(ymaxNVretro)+10]);
%     xlabel('Net Velocities (um/sec)');
%     ylabel('Frequency (%)');
% 
%     % Gaussian Fitting for retrograde net velocities
%     data=AllNV.(list(i,1).name).retro; % change to AllNV for real data set
% 
%     %run function multigaussian fitting
%     try
%     multigauss_fitting_110414(data);
%     
%     % rewrite variables into arrays for later saving
%     MultiGaussian.(list(i,1).name).retroNV.mu=mu; 
%     MultiGaussian.(list(i,1).name).retroNV.sigma=sigma;
%     MultiGaussian.(list(i,1).name).retroNV.sigmasq=sigmasq;
%     MultiGaussian.(list(i,1).name).retroNV.prob=prob;
%     MultiGaussian.(list(i,1).name).retroNV.clusterNum=clusterNum;
%     MultiGaussian.(list(i,1).name).retroNV.clusterBIC=BIC;
%     MultiGaussian.(list(i,1).name).retroNV.clusterClassification=classification;
% 
%     % calculate distribution
%     xGauss= -1:0.1:round(maxNVretro)+1;
%     yGauss=([]);
%     for j = 1 : length(xGauss);
%         yGauss(j) = 0;
%         for k = 1 : clusterNum
%             temp = - 0.5 * (xGauss(j) - mu(k))^2 / sigmasq(k);
%             yGauss(j) = yGauss(j) + prob(k) * one_over_sqrt_2pi / sqrt(sigmasq(k)) * exp(temp); 
%         end
%     end
%     yGaussMax=max(yGauss);
%     MultiGaussian.(list(i,1).name).retroNV.xGauss=xGauss;
%     MultiGaussian.(list(i,1).name).retroNV.yGauss.(binSizeName{m})=yGauss/yGaussMax*max(histNV.(list(i,1).name).pct.Retro.(binSizeName{m}));
%     % plot by maximum matching
%     plot(xGauss,yGauss/yGaussMax*max(histNV.(list(i,1).name).pct.Retro.(binSizeName{m})),'r-');
%     hold on
%     %calculate separate clusters
%     for j=1:clusterNum
%         number=num2str(j);
%         clustNum=strcat('Cluster',number);
%         for k=1:length(xGauss);
%             temp = - 0.5 * (xGauss(k) - mu(j))^2 / sigmasq(j);
%             yGauss(k) = prob(j) * one_over_sqrt_2pi / sqrt(sigmasq(j)) * exp(temp); 
%         end
%         MultiGaussian.(list(i,1).name).retroNV.seperatedCluster.(clustNum).xGauss=xGauss;
%         MultiGaussian.(list(i,1).name).retroNV.seperatedCluster.(clustNum).yGauss.(binSizeName{m})=yGauss/yGaussMax*max(histNV.(list(i,1).name).pct.Retro.(binSizeName{m}));
%     end
% 
%     %plot single clusters by maximum matching
%     for j=1:clusterNum
%         number=num2str(j);
%         clustNum=strcat('Cluster',number);
%         plot(xGauss,MultiGaussian.(list(i,1).name).retroNV.seperatedCluster.(clustNum).yGauss.(binSizeName{m}),'c-');
%     end
% 
%     catch ME
%         MultiGaussian.(list(i,1).name).retroNV=NaN;
%     end
% 
% 
%     if yGaussMax>ymaxNVretro
%         ymaxNVretroNew=yGaussMax;
%         ylim([0 ymaxNVretroNew]);
%     end
% 
%     %% run length
%     % histRL.(list(i,1).name)=([]);
%     % histRL.(list(i,1).name).data.Ante=hist(AllRL.(list(i,1).name).ante(:,1),xVector);
%     % FRCT=rdivide(histRL.(list(i,1).name).data.Ante,sum(histRL.(list(i,1).name).data.Ante));
%     % histRL.(list(i,1).name).pct.Ante=mtimes(100,FRCT);
%     subplot('Position',[0.05,0.1,0.15,0.3]);
%     set(gca,'fontsize',8);
%     hb =bar(xVectorRLante.(binSizeName{m}),histRL.(list(i,1).name).pct.Ante.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
%     title('Anterograde Segmental Run Length');
%     xlim([-1 round(maxRLante)+1]);
%     ylim([0 round(ymaxRLante)+10]);
%     xlabel('Segmental Run Length (um)');
%     ylabel('Frequency (%)');
% 
%     % histRL.(list(i,1).name).data.Retro=hist(AllRL.(list(i,1).name).retro(:,1),xVector);
%     % FRCT=rdivide(histRL.(list(i,1).name).data.Retro,sum(histRL.(list(i,1).name).data.Retro));
%     % histRL.(list(i,1).name).pct.Retro=mtimes(100,FRCT);
%     subplot('Position',[0.3,0.1,0.15,0.3]);
%     set(gca,'fontsize',8);
%     hb =bar(xVectorRLretro.(binSizeName{m}),histRL.(list(i,1).name).pct.Retro.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
%     title('Retrograde Segmental Run Length');
%     xlim([-1 round(maxRLretro)+1]);
%     ylim([0 round(ymaxRLretro)+10]);
%     xlabel('Segmental Run Length (um)');
%     ylabel('Frequency (%)');
% 
%     %combined run length
%     % histRL.(list(i,1).name)=([]);
%     % histRL.(list(i,1).name).data.Ante=hist(AllRL.(list(i,1).name).ante(:,1),xVector);
%     % FRCT=rdivide(histRL.(list(i,1).name).data.Ante,sum(histRL.(list(i,1).name).data.Ante));
%     % histRL.(list(i,1).name).pct.Ante=mtimes(100,FRCT);
%     subplot('Position',[0.55,0.1,0.15,0.3]);
%     set(gca,'fontsize',8);
%     hb=bar(xVectorcomRLante.(binSizeName{m}),histcomRL.(list(i,1).name).pct.Ante.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
%     title('Anterograde Combined Segmental Run Length');
%     xlim([-1 round(maxcomRLante)+1]);
%     ylim([0 round(ymaxcomRLante)+10]);
%     xlabel('Segmental Run Length (um)');
%     ylabel('Frequency (%)');
% 
%     % histRL.(list(i,1).name).data.Retro=hist(AllRL.(list(i,1).name).retro(:,1),xVector);
%     % FRCT=rdivide(histRL.(list(i,1).name).data.Retro,sum(histRL.(list(i,1).name).data.Retro));
%     % histRL.(list(i,1).name).pct.Retro=mtimes(100,FRCT);
%     subplot('Position',[0.8,0.1,0.15,0.3]);
%     set(gca,'fontsize',8);
%     hb=bar(xVectorcomRLretro.(binSizeName{m}),histcomRL.(list(i,1).name).pct.Retro.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
%     title('Retrograde Combined Segmental Run Length');
%     xlim([-1 round(maxcomRLretro)+1]);
%     ylim([0 round(ymaxcomRLretro)+10]);
%     xlabel('Segmental Run Length (um)');
%     ylabel('Frequency (%)');
% 
%     %% saving
%     name=list(i,1).name;
%     mkdir(fullfile(directoryMain,'Figures','HistogramsMultiGaussian',binSizeName{m}));
%     fileFigure=fullfile(directoryMain,'Figures','HistogramsMultiGaussian',binSizeName{m});
%     fileFigure=strcat(fileFigure,filesep,name,'_HistogramsMultiGaussian.fig');
%     saveas(figure3,fileFigure);
%     fileFigure=fullfile(directoryMain,'Figures','HistogramsMultiGaussian',binSizeName{m});
%     fileFigure=strcat(fileFigure,filesep,name,'_HistogramsMultiGaussian.eps');
%     saveas(figure3,fileFigure);
%     close all;
%     end
    
    %% Histograms without MultiGaussianFitting
    mkdir(fullfile(directoryMain,'Figures','Histograms',binSizeName{m}));
    for i=1:length(list)
        figure4=figure;
        set(gcf,'Position',[600 100 1200 600]);

        %% segmental velocities
        % anterograde

        subplot('Position',[0.05,0.6,0.15,0.3]);
        set(gca,'fontsize',8);
        hb=bar(xVectorSVante.(binSizeName{m}),histSV.(list(i,1).name).pct.Ante.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
        title('Anterograde Segmental Velocities');
        xlim([-1 round(maxSVante)+1]);
        ylim([0 round(ymaxSVante)+10]);
        xlabel('Segmental Velocities (um/sec)');
        ylabel('Frequency (%)');

        % retrograde
        subplot('Position',[0.3,0.6,0.15,0.3]);
        set(gca,'fontsize',8);
        hb=bar(xVectorSVretro.(binSizeName{m}),histSV.(list(i,1).name).pct.Retro.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
        hold on
        title('Retrograde Segmental Velocities');
        xlim([-1 round(maxSVretro)+1]);
        ylim([0 round(ymaxSVretro)+10]);
        xlabel('Segmental Velocities (um/sec)');
        ylabel('Frequency (%)');

        % anterograde net velocities

        subplot('Position',[0.55,0.6,0.15,0.3]);
        set(gca,'fontsize',8);
        hb=bar(xVectorNVante.(binSizeName{m}),histNV.(list(i,1).name).pct.Ante.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
        title('Anterograde Net Velocities');
        xlim([-1 round(maxNVante)+1]);
        ylim([0 round(ymaxNVante)+10]);
        xlabel('Net Velocities (um/sec)');
        ylabel('Frequency (%)');

        %retrograde net velocities

        subplot('Position',[0.8,0.6,0.15,0.3]);
        set(gca,'fontsize',8);
        hb=bar(xVectorNVretro.(binSizeName{m}),histNV.(list(i,1).name).pct.Retro.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
        title('Retrograde Net Velocities');
        xlim([-1 round(maxNVretro)+1]);
        ylim([0 round(ymaxNVretro)+10]);
        xlabel('Net Velocities (um/sec)');
        ylabel('Frequency (%)');

        %% run length
        %anterograde
        subplot('Position',[0.05,0.1,0.15,0.3]);
        set(gca,'fontsize',8);
        hb=bar(xVectorRLante.(binSizeName{m}),histRL.(list(i,1).name).pct.Ante.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
        title('Anterograde Segmental Run Length');
        xlim([-1 round(maxRLante)+1]);
        ylim([0 round(ymaxRLante)+10]);
        xlabel('Segmental Run Length (um)');
        ylabel('Frequency (%)');

        % retrograde
        subplot('Position',[0.3,0.1,0.15,0.3]);
        set(gca,'fontsize',8);
        hb=bar(xVectorRLretro.(binSizeName{m}),histRL.(list(i,1).name).pct.Retro.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
        title('Retrograde Segmental Run Length');
        xlim([-1 round(maxRLretro)+1]);
        ylim([0 round(ymaxRLretro)+10]);
        xlabel('Segmental Run Length (um)');
        ylabel('Frequency (%)');

        %combined run length
        %anterograde
        subplot('Position',[0.55,0.1,0.15,0.3]);
        set(gca,'fontsize',8);
        hb=bar(xVectorcomRLante.(binSizeName{m}),histcomRL.(list(i,1).name).pct.Ante.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
        title('Anterograde Combined Segmental Run Length');
        xlim([-1 round(maxcomRLante)+1]);
        ylim([0 round(ymaxcomRLante)+10]);
        xlabel('Segmental Run Length (um)');
        ylabel('Frequency (%)');

        %retrograde
        subplot('Position',[0.8,0.1,0.15,0.3]);
        set(gca,'fontsize',8);
        hb=bar(xVectorcomRLretro.(binSizeName{m}),histcomRL.(list(i,1).name).pct.Retro.(binSizeName{m}),'Facecolor',color(i,:),'BarWidth',1);
        title('Retrograde Combined Segmental Run Length');
        xlim([-1 round(maxcomRLretro)+1]);
        ylim([0 round(ymaxcomRLretro)+10]);
        xlabel('Segmental Run Length (um)');
        ylabel('Frequency (%)');

        % saving
        name=list(i,1).name;
        fileFigure=fullfile(directoryMain,'Figures','Histograms',binSizeName{m});
        fileFigure=strcat(fileFigure,filesep,name,'_Histograms.fig');
        saveas(figure4,fileFigure);
        fileFigure=fullfile(directoryMain,'Figures','Histograms',binSizeName{m});
        fileFigure=strcat(fileFigure,filesep,name,'_Histograms.eps');
        saveas(figure4,fileFigure);
        close all;
    end
end
%% boxplots
boxAnteSV=([]);
boxRetroSV=([]);
boxAnteNV=([]);
boxRetroNV=([]);
boxAnteRL=([]);
boxRetroRL=([]);
boxAntecomRL=([]);
boxRetrocomRL=([]);

grpAnteSV=([]);
grpRetroSV=([]);
grpAnteNV=([]);
grpRetroNV=([]);
grpAnteRL=([]);
grpRetroRL=([]);
grpAntecomRL=([]);
grpRetrocomRL=([]);

%write data for boxplots
for i=1:length(list)
    boxAnteSV=cat(1,boxAnteSV,AllSV.(list(i,1).name).ante);
    boxRetroSV=cat(1,boxRetroSV,AllSV.(list(i,1).name).retro);
    boxAnteNV=cat(1,boxAnteNV,AllNV.(list(i,1).name).ante);
    boxRetroNV=cat(1,boxRetroNV,AllNV.(list(i,1).name).retro);
    boxAnteRL=cat(1,boxAnteRL,AllRL.(list(i,1).name).ante);
    boxRetroRL=cat(1,boxRetroRL,AllRL.(list(i,1).name).retro);
    boxAntecomRL=cat(1,boxAntecomRL,AllcomRL.(list(i,1).name).ante);
    boxRetrocomRL=cat(1,boxRetrocomRL,AllcomRL.(list(i,1).name).retro);
    
    GRPAnteSV=repmat({list(i,1).name},length(AllSV.(list(i,1).name).ante),1);
    GRPRetroSV=repmat({list(i,1).name},length(AllSV.(list(i,1).name).retro),1);
    GRPAnteNV=repmat({list(i,1).name},length(AllNV.(list(i,1).name).ante),1);
    GRPRetroNV=repmat({list(i,1).name},length(AllNV.(list(i,1).name).retro),1);
    GRPAnteRL=repmat({list(i,1).name},length(AllRL.(list(i,1).name).ante),1);
    GRPRetroRL=repmat({list(i,1).name},length(AllRL.(list(i,1).name).retro),1);
    GRPAntecomRL=repmat({list(i,1).name},length(AllcomRL.(list(i,1).name).ante),1);
    GRPRetrocomRL=repmat({list(i,1).name},length(AllcomRL.(list(i,1).name).retro),1);
    
    grpAnteSV=cat(1,grpAnteSV,GRPAnteSV);
    grpRetroSV=cat(1,grpRetroSV,GRPRetroSV);
    grpAnteNV=cat(1,grpAnteNV,GRPAnteNV);
    grpRetroNV=cat(1,grpRetroNV,GRPRetroNV);
    grpAnteRL=cat(1,grpAnteRL,GRPAnteRL);
    grpRetroRL=cat(1,grpRetroRL,GRPRetroRL);
    grpAntecomRL=cat(1,grpAntecomRL,GRPAntecomRL);
    grpRetrocomRL=cat(1,grpRetrocomRL,GRPRetrocomRL);
    
end

figure3=figure;
set(gcf,'Position',[600 100 1200 600]);

subplot('Position',[0.05,0.6,0.15,0.3]);
h=boxplot(boxAnteSV,grpAnteSV,'labelorientation','inline');
set(gca,'fontsize',8);
title('Segmental Velocities');
ylabel('Anterograde Segmental Velocities (um/sec)');

subplot('Position',[0.3,0.6,0.15,0.3]);
h=boxplot(boxRetroSV,grpRetroSV,'labelorientation','inline');
set(gca,'fontsize',8);
title('Segmental Velocities');
ylabel('Retrograde Segmental Velocities (um/sec)');

subplot('Position',[0.55,0.6,0.15,0.3]);
h=boxplot(boxAnteNV,grpAnteNV,'labelorientation','inline');
set(gca,'fontsize',8);
title('Net Velocities');
ylabel('Anterograde Net Velocities (um/sec)');

subplot('Position',[0.8,0.6,0.15,0.3]);
h=boxplot(boxRetroNV,grpRetroNV,'labelorientation','inline');
set(gca,'fontsize',8);
title('Net Velocities');
ylabel('Retrograde Net Velocities (um/sec)');

subplot('Position',[0.05,0.1,0.15,0.3]);
h=boxplot(boxAnteRL,grpAnteRL,'labelorientation','inline');
set(gca,'fontsize',8);
title('Segmental Run Length');
ylabel('Anterograde Run Length(um)');

subplot('Position',[0.3,0.1,0.15,0.3]);
h=boxplot(boxRetroRL,grpRetroRL,'labelorientation','inline');
set(gca,'fontsize',8);
title('Segmental Run Length');
ylabel('Retrograde Run Length(um)');

subplot('Position',[0.55,0.1,0.15,0.3]);
h=boxplot(boxAntecomRL,grpAntecomRL,'labelorientation','inline');
set(gca,'fontsize',8);
title('Combined Segmental Run Length');
ylabel('Anterograde Run Length(um)');

subplot('Position',[0.8,0.1,0.15,0.3]);
h=boxplot(boxRetrocomRL,grpRetrocomRL,'labelorientation','inline');
set(gca,'fontsize',8);
title('Combined Segmental Run Length');
ylabel('Retrograde Run Length(um)');

% saving
fileFigure=fullfile(directoryMain,'Figures');
fileFigure=strcat(fileFigure,filesep,'Boxplots.fig');
saveas(figure3,fileFigure);
fileFigure=fullfile(directoryMain,'Figures');
fileFigure=strcat(fileFigure,filesep,'Boxplots.eps');
saveas(figure3,fileFigure);
close all;

%% Statistics: calculate normality, ranksum and randomized T-test

%% lillietest
normal=([]);
for i=1:length(list)
    %% CP 
    try
        normal.(list(i,1).name).CPante=lillietest(AllCP.(list(i,1).name).PCT(:,1));
    catch ME
        normal.(list(i,1).name).CPante=NaN;
    end
    
    try
        normal.(list(i,1).name).CPretro=lillietest(AllCP.(list(i,1).name).PCT(:,2));
    catch ME
        normal.(list(i,1).name).CPretro=NaN;
    end
    
    try
       normal.(list(i,1).name).CPrev=lillietest(AllCP.(list(i,1).name).PCT(:,3));
    catch ME
        normal.(list(i,1).name).CPrev=NaN;
    end
    
    try
        normal.(list(i,1).name).CPstat=lillietest(AllCP.(list(i,1).name).PCT(:,4));
    catch ME
        normal.(list(i,1).name).CPstat=NaN;
    end
    %% NCP
    try
        normal.(list(i,1).name).NCPante=lillietest(AllNCP.(list(i,1).name).PCT(:,1));
    catch ME
        normal.(list(i,1).name).NCPante=NaN;
    end
    
    try
        normal.(list(i,1).name).NCPretro=lillietest(AllNCP.(list(i,1).name).PCT(:,2));
    catch ME
        normal.(list(i,1).name).NCPretro=NaN;
    end
    
    try
        normal.(list(i,1).name).NCPrev=lillietest(AllNCP.(list(i,1).name).PCT(:,3));
    catch ME
        normal.(list(i,1).name).NCPrev=NaN;
    end
    %% density
    try
        normal.(list(i,1).name).Densityante=lillietest(AllDensity.(list(i,1).name).Value(:,1));
    catch ME
        normal.(list(i,1).name).Densityante=NaN;
    end
    
    try
        normal.(list(i,1).name).Densityretro=lillietest(AllDensity.(list(i,1).name).Value(:,2));
    catch ME
        normal.(list(i,1).name).Densityretro=NaN;
    end
    
    try
        normal.(list(i,1).name).Densityrev=lillietest(AllDensity.(list(i,1).name).Value(:,3));
    catch ME
        normal.(list(i,1).name).Densityrev=NaN;
    end
    
    try
        normal.(list(i,1).name).Densitystat=lillietest(AllDensity.(list(i,1).name).Value(:,4));
    catch ME
        normal.(list(i,1).name).Densitystat=NaN;
    end
    %% flux
    try
    normal.(list(i,1).name).Fluxante=lillietest(AllFlux.(list(i,1).name).Value(:,1));
    catch ME
        normal.(list(i,1).name).Fluxante=NaN;
    end
    
    try
        normal.(list(i,1).name).Fluxretro=lillietest(AllFlux.(list(i,1).name).Value(:,2));
    catch ME
        normal.(list(i,1).name).Fluxretro=NaN;
    end
    
    try
        normal.(list(i,1).name).Fluxrev=lillietest(AllFlux.(list(i,1).name).Value(:,3));
    catch ME
        normal.(list(i,1).name).Fluxrev=NaN;
    end
  %% SV  
    try
        normal.(list(i,1).name).SVante=lillietest(AllSV.(list(i,1).name).ante);
    catch ME
        normal.(list(i,1).name).SVante=NaN;
    end
    
    try
        normal.(list(i,1).name).SVretro=lillietest(AllSV.(list(i,1).name).retro);
    catch ME
        normal.(list(i,1).name).SVretro=NaN;
    end
    %% NV 
    try
        normal.(list(i,1).name).NVante=lillietest(AllNV.(list(i,1).name).ante);
    catch ME
        normal.(list(i,1).name).NVante=NaN;
    end
    
    try
        normal.(list(i,1).name).NVretro=lillietest(AllNV.(list(i,1).name).retro);
    catch ME
        normal.(list(i,1).name).NVretro=NaN;
    end
    %% RL
    try
        normal.(list(i,1).name).RLante=lillietest(AllRL.(list(i,1).name).ante);
    catch ME
        normal.(list(i,1).name).RLante=NaN;
    end
    
    try
        normal.(list(i,1).name).RLretro=lillietest(AllRL.(list(i,1).name).retro);
    catch ME
        normal.(list(i,1).name).RLretro=NaN;
    end
    %% com RL
    try
        normal.(list(i,1).name).comRLante=lillietest(AllcomRL.(list(i,1).name).ante);
    catch ME
        normal.(list(i,1).name).comRLante=NaN;
    end
    
    try
        normal.(list(i,1).name).comRLretro=lillietest(AllcomRL.(list(i,1).name).retro);
    catch ME
        normal.(list(i,1).name).comRLretro=NaN;
    end
    %% PD
    try
        normal.(list(i,1).name).PD=lillietest(AllPD.(list(i,1).name).Value);
    catch ME
        normal.(list(i,1).name).PD=NaN;
    end
        %% split PD
    try
        normal.(list(i,1).name).splitPDante=lillietest(AllsplitPD.(list(i,1).name).ante);
    catch ME
        normal.(list(i,1).name).splitPDante=NaN;
    end

    try
        normal.(list(i,1).name).splitPDretro=lillietest(AllsplitPD.(list(i,1).name).retro);
    catch ME
        normal.(list(i,1).name).splitPDretro=NaN;
    end
    
    try
        normal.(list(i,1).name).splitPDrev=lillietest(AllsplitPD.(list(i,1).name).rev);
    catch ME
        normal.(list(i,1).name).splitPDrev=NaN;
    end
    %% PF
    try
        normal.(list(i,1).name).PF=lillietest(AllPF.(list(i,1).name).Value);
    catch ME
        normal.(list(i,1).name).PF=NaN;
    end
    %% split PF
    try
        normal.(list(i,1).name).splitPFante=lillietest(AllsplitPF.(list(i,1).name).ante);
    catch ME
        normal.(list(i,1).name).splitPFante=NaN;
    end
    
    try
        normal.(list(i,1).name).splitPFretro=lillietest(AllsplitPF.(list(i,1).name).retro);
    catch ME
        normal.(list(i,1).name).splitPFretro=NaN;
    end
    
    try
        normal.(list(i,1).name).splitPFrev=lillietest(AllsplitPF.(list(i,1).name).rev);
    catch ME
        normal.(list(i,1).name).splitPFrev=NaN;
    end
    %% PF per sec
    try
        normal.(list(i,1).name).PFperSec=lillietest(AllPFperSec.(list(i,1).name).Value);
    catch ME
        normal.(list(i,1).name).PFperSeantec=NaN;
    end
    %% split PF per sec
    try
        normal.(list(i,1).name).splitPFperSecante=lillietest(AllsplitPFperSec.(list(i,1).name).ante);
    catch ME
        normal.(list(i,1).name).splitPFperSecante=NaN;
    end
    
    try
        normal.(list(i,1).name).splitPFperSecretro=lillietest(AllsplitPFperSec.(list(i,1).name).retro);
    catch ME
        normal.(list(i,1).name).splitPFperSecretro=NaN;
    end
    
    try
        normal.(list(i,1).name).splitPFperSecrev=lillietest(AllsplitPFperSec.(list(i,1).name).rev);
    catch ME
        normal.(list(i,1).name).splitPFperSecrev=NaN;
    end
    %% SF
    try
        normal.(list(i,1).name).SF=lillietest(AllSF.(list(i,1).name).Value);
    catch ME
        normal.(list(i,1).name).SF=NaN;
    end
    %% SF per sec
    try
        normal.(list(i,1).name).SFperSec=lillietest(AllSFperSec.(list(i,1).name).Value);
    catch ME
        normal.(list(i,1).name).SFperSec=NaN;
    end
    %% revSF
    try
        normal.(list(i,1).name).revSF=lillietest(AllrevSF.(list(i,1).name).Value);
    catch ME
        normal.(list(i,1).name).revSF=NaN;
    end
    
    try
        normal.(list(i,1).name).revSFperSec=lillietest(AllrevSFperSec.(list(i,1).name).Value);
    catch ME
        normal.(list(i,1).name).revSFperSec=NaN;
    end
end

%% ranksum and randomized ttest
if length(list)>1
    for i=1:length(list)
        listName=list(i,1).name;
        k=1;
        for j=1:length(list)
                if strcmp(listName,list(j,1).name)==0
                    list3(k,1).name=list(j,1).name;
                    k=k+1;
                end
        end

        for j=1:length(list3)
            comp=strcat(list(i,1).name,'_vs_',list3(j,1).name);
            
            try
            [pAnteCP,hAnteCP]=ranksum(AllCP.(list(i,1).name).PCT(:,1),AllCP.(list3(j,1).name).PCT(:,1));
            rnk.(comp).anteCP=[pAnteCP,hAnteCP];
            catch
                rnk.(comp).anteCP=NaN;
            end
            
            try
            [pRetroCP,hRetroCP]=ranksum(AllCP.(list(i,1).name).PCT(:,2),AllCP.(list3(j,1).name).PCT(:,2));
            rnk.(comp).retroCP=[pRetroCP,hRetroCP];
            catch
                rnk.(comp).retroCP=NaN;
            end
            
            try
            [pRevCP,hRevCP]=ranksum(AllCP.(list(i,1).name).PCT(:,3),AllCP.(list3(j,1).name).PCT(:,3));
            rnk.(comp).revCP=[pRevCP,hRevCP];
            catch
                rnk.(comp).revCP=NaN;
            end
            
            try
            [pStatCP,hStatCP]=ranksum(AllCP.(list(i,1).name).PCT(:,4),AllCP.(list3(j,1).name).PCT(:,4));
            rnk.(comp).statCP=[pStatCP,hStatCP];
            catch
                rnk.(comp).statCP=NaN;
            end
            
            try
            [pAnteNCP,hAnteNCP]=ranksum(AllNCP.(list(i,1).name).PCT(:,1),AllNCP.(list3(j,1).name).PCT(:,1));
            rnk.(comp).anteNCP=[pAnteNCP,hAnteNCP];
            catch
                rnk.(comp).anteNCP=NaN;
            end
            
            try
            [pRetroNCP,hRetroNCP]=ranksum(AllNCP.(list(i,1).name).PCT(:,2),AllNCP.(list3(j,1).name).PCT(:,2));
            rnk.(comp).retroNCP=[pRetroNCP,hRetroNCP];
            catch
                rnk.(comp).retroNCP=NaN;
            end
            
            try
            [pStatNCP,hStatNCP]=ranksum(AllNCP.(list(i,1).name).PCT(:,3),AllNCP.(list3(j,1).name).PCT(:,3));
            rnk.(comp).statNCP=[pStatNCP,hStatNCP];
            catch
                rnk.(comp).statNCP=NaN;
            end
            
            %Density 
            try
            [pAnteDensity,hAnteDensity]=ranksum(AllDensity.(list(i,1).name).Value(:,1),AllDensity.(list3(j,1).name).Value(:,1));
            rnk.(comp).anteDensity=[pAnteDensity,hAnteDensity];
            catch
                rnk.(comp).anteDensity=NaN;
            end
            
            try
            [pRetroDensity,hRetroDensity]=ranksum(AllDensity.(list(i,1).name).Value(:,2),AllDensity.(list3(j,1).name).Value(:,2));
            rnk.(comp).retroDensity=[pRetroDensity,hRetroDensity];
            catch
                rnk.(comp).retroDensity=NaN;
            end
            
            try
            [pRevDensity,hRevDensity]=ranksum(AllDensity.(list(i,1).name).Value(:,3),AllDensity.(list3(j,1).name).Value(:,3));
            rnk.(comp).revDensity=[pRevDensity,hRevDensity];
            catch
                rnk.(comp).revDensity=NaN;
            end
            
            try
            [pStatDensity,hStatDensity]=ranksum(AllDensity.(list(i,1).name).Value(:,4),AllDensity.(list3(j,1).name).Value(:,4));
            rnk.(comp).statDensity=[pStatDensity,hStatDensity];
            catch
                rnk.(comp).statDensity=NaN;
            end
            
            %flux        
            try
            [pAnteFlux,hAnteFlux]=ranksum(AllFlux.(list(i,1).name).Value(:,1),AllFlux.(list3(j,1).name).Value(:,1));
            rnk.(comp).anteFlux=[pAnteFlux,hAnteFlux];
            catch
                rnk.(comp).anteFlux=NaN;
            end
            
            try
            [pRetroFlux,hRetroFlux]=ranksum(AllFlux.(list(i,1).name).Value(:,2),AllFlux.(list3(j,1).name).Value(:,2));
            rnk.(comp).retroFlux=[pRetroFlux,hRetroFlux];
            catch
                rnk.(comp).retroFlux=NaN;
            end
            
            try
            [pRevFlux,hRevFlux]=ranksum(AllFlux.(list(i,1).name).Value(:,3),AllFlux.(list3(j,1).name).Value(:,3));
            rnk.(comp).revFlux=[pRevFlux,hRevFlux];
            catch
                rnk.(comp).revFlux=NaN;
            end
            
            %SV
            try
            [pAnteSV,hAnteSV]=ranksum(AllSV.(list(i,1).name).ante,AllSV.(list3(j,1).name).ante);
            rnk.(comp).anteSV=[pAnteSV hAnteSV];
            catch ME
                rnk.(comp).anteSV=NaN;
            end

            try
            [pRetroSV,hRetroSV]=ranksum(AllSV.(list(i,1).name).retro,AllSV.(list3(j,1).name).retro);
            rnk.(comp).retroSV=[pRetroSV hRetroSV];
            catch ME
                rnk.(comp).retroSV=NaN;
            end
            
            % comSV
            try
            [pAntecomSV,hAntecomSV]=ranksum(AllcomSV.(list(i,1).name).ante,AllcomSV.(list3(j,1).name).ante);
            rnk.(comp).antecomSV=[pAntecomSV hAntecomSV];
            catch ME
                rnk.(comp).antecomSV=NaN;
            end

            try       
            [pRetrocomSV,hRetrocomSV]=ranksum(AllcomSV.(list(i,1).name).retro,AllcomSV.(list3(j,1).name).retro);
            rnk.(comp).retrocomSV=[pRetrocomSV hRetrocomSV];
            catch ME
                rnk.(comp).retrocomSV=NaN;
            end
            
            %NV
            try
            [pAnteNV,hAnteNV]=ranksum(AllNV.(list(i,1).name).ante,AllNV.(list3(j,1).name).ante);
            rnk.(comp).anteNV=[pAnteNV hAnteNV];
            catch ME
                rnk.(comp).anteNV=NaN;
            end

            try
            [pRetroNV,hRetroNV]=ranksum(AllNV.(list(i,1).name).retro,AllNV.(list3(j,1).name).retro);
            rnk.(comp).retroNV=[pRetroNV hRetroNV];
            catch ME
                rnk.(comp).retroNV=NaN;
            end

            try
            [pAnteRL,hAnteRL]=ranksum(AllRL.(list(i,1).name).ante,AllRL.(list3(j,1).name).ante);
            rnk.(comp).anteRL=[pAnteRL hAnteRL];
            catch ME
                rnk.(comp).anteRL=NaN;
            end

            try       
            [pRetroRL,hRetroRL]=ranksum(AllRL.(list(i,1).name).retro,AllRL.(list3(j,1).name).retro);
            rnk.(comp).retroRL=[pRetroRL hRetroRL];
            catch ME
                rnk.(comp).retroRL=NaN;
            end
            
            try
            [pAntecomRL,hAntecomRL]=ranksum(AllcomRL.(list(i,1).name).ante,AllcomRL.(list3(j,1).name).ante);
            rnk.(comp).antecomRL=[pAntecomRL hAntecomRL];
            catch ME
                rnk.(comp).antecomRL=NaN;
            end

            try       
            [pRetrocomRL,hRetrocomRL]=ranksum(AllcomRL.(list(i,1).name).retro,AllcomRL.(list3(j,1).name).retro);
            rnk.(comp).retrocomRL=[pRetrocomRL hRetrocomRL];
            catch ME
                rnk.(comp).retrocomRL=NaN;
            end

            try
            [pPD,hPD]=ranksum(AllPD.(list(i,1).name).Value,AllPD.(list3(j,1).name).Value);
            rnk.(comp).PD=[pPD hPD];
            catch ME
                rnk.(comp).PD=NaN;
            end

            try
            [pPD,hPD]=ranksum(AllsplitPD.(list(i,1).name).ante,AllsplitPD.(list3(j,1).name).ante);
            rnk.(comp).splitPDante=[pPD hPD];
            catch ME
                rnk.(comp).splitPDante=NaN;
            end
            
            try
            [pPD,hPD]=ranksum(AllsplitPD.(list(i,1).name).retro,AllsplitPD.(list3(j,1).name).retro);
            rnk.(comp).splitPDretro=[pPD hPD];
            catch ME
                rnk.(comp).splitPDretro=NaN;
            end
            
            try
            [pPD,hPD]=ranksum(AllsplitPD.(list(i,1).name).rev,AllsplitPD.(list3(j,1).name).rev);
            rnk.(comp).splitPDrev=[pPD hPD];
            catch ME
                rnk.(comp).splitPDrev=NaN;
            end
            
            try
            [pPF,hPF]=ranksum(AllPF.(list(i,1).name).Value,AllPF.(list3(j,1).name).Value);
            rnk.(comp).PF=[pPF hPF];
            catch ME
                rnk.(comp).PF=NaN;
            end

            try
            [pPF,hPF]=ranksum(AllsplitPF.(list(i,1).name).ante,AllsplitPF.(list3(j,1).name).ante);
            rnk.(comp).splitPFante=[pPF hPF];
            catch ME
                rnk.(comp).splitPFante=NaN;
            end
            
            try
            [pPF,hPF]=ranksum(AllsplitPF.(list(i,1).name).retro,AllsplitPF.(list3(j,1).name).retro);
            rnk.(comp).splitPFretro=[pPF hPF];
            catch ME
                rnk.(comp).splitPFretro=NaN;
            end
            
            try
            [pPF,hPF]=ranksum(AllsplitPF.(list(i,1).name).rev,AllsplitPF.(list3(j,1).name).rev);
            rnk.(comp).splitPFrev=[pPF hPF];
            catch ME
                rnk.(comp).splitPFrev=NaN;
            end
            
            try
            [pPFperSec,hPFperSec]=ranksum(AllPFperSec.(list(i,1).name).Value,AllPFperSec.(list3(j,1).name).Value);
            rnk.(comp).PFperSec=[pPFperSec hPFperSec];
            catch ME
                rnk.(comp).PFperSec=NaN;
            end

            try
            [pPF,hPF]=ranksum(AllsplitPFperSec.(list(i,1).name).ante,AllsplitPFperSec.(list3(j,1).name).ante);
            rnk.(comp).splitPFperSecante=[pPF hPF];
            catch ME
                rnk.(comp).splitPFperSecante=NaN;
            end
            
            try
            [pPF,hPF]=ranksum(AllsplitPFperSec.(list(i,1).name).retro,AllsplitPFperSec.(list3(j,1).name).retro);
            rnk.(comp).splitPFperSecretro=[pPF hPF];
            catch ME
                rnk.(comp).splitPFperSecretro=NaN;
            end
            
            try
            [pPF,hPF]=ranksum(AllsplitPFperSec.(list(i,1).name).rev,AllsplitPFperSec.(list3(j,1).name).rev);
            rnk.(comp).splitPFperSecrev=[pPF hPF];
            catch ME
                rnk.(comp).splitPFperSecrev=NaN;
            end
            
            try
            [pSF,hSF]=ranksum(AllSF.(list(i,1).name).Value,AllSF.(list3(j,1).name).Value);
            rnk.(comp).SF=[pSF hSF];
            catch ME
                rnk.(comp).SF=NaN;
            end

            try
            [pSFperSec,hSFperSec]=ranksum(AllSFperSec.(list(i,1).name).Value,AllSFperSec.(list3(j,1).name).Value);
            rnk.(comp).SFperSec=[pSFperSec hSFperSec];
            catch ME
                rnk.(comp).SFperSec=NaN;
            end
            
            try
            [prevSF,hrevSF]=ranksum(AllrevSF.(list(i,1).name).Value,AllrevSF.(list3(j,1).name).Value);
            rnk.(comp).revSF=[prevSF hrevSF];
            catch ME
                rnk.(comp).revSF=NaN;
            end

            try
            [prevSFperSec,hrevSFperSec]=ranksum(AllrevSFperSec.(list(i,1).name).Value,AllrevSFperSec.(list3(j,1).name).Value);
            rnk.(comp).revSFperSec=[prevSFperSec hrevSFperSec];
            catch ME
                rnk.(comp).revSFperSec=NaN;
            end

    %% randomized T-Test: for this test all arrays have to be transposed,
    % because of the way they have to be input into Matlab

            AllCP.(list(i,1).name).PCT=transpose(AllCP.(list(i,1).name).PCT);
            AllCP.(list3(j,1).name).PCT=transpose(AllCP.(list3(j,1).name).PCT);
            AllNCP.(list(i,1).name).PCT=transpose(AllNCP.(list(i,1).name).PCT);
            AllNCP.(list3(j,1).name).PCT=transpose(AllNCP.(list3(j,1).name).PCT);
            
            AllDensity.(list(i,1).name).Value=transpose(AllDensity.(list(i,1).name).Value);
            AllDensity.(list3(j,1).name).Value=transpose(AllDensity.(list3(j,1).name).Value);
            
            AllFlux.(list(i,1).name).Value=transpose(AllFlux.(list(i,1).name).Value);
            AllFlux.(list3(j,1).name).Value=transpose(AllFlux.(list3(j,1).name).Value);
            
            AllSV.(list(i,1).name).ante=transpose(AllSV.(list(i,1).name).ante);
            AllSV.(list3(j,1).name).ante=transpose(AllSV.(list3(j,1).name).ante);
            AllSV.(list(i,1).name).retro=transpose(AllSV.(list(i,1).name).retro);
            AllSV.(list3(j,1).name).retro=transpose(AllSV.(list3(j,1).name).retro);

            AllNV.(list(i,1).name).ante=transpose(AllNV.(list(i,1).name).ante);
            AllNV.(list3(j,1).name).ante=transpose(AllNV.(list3(j,1).name).ante);
            AllNV.(list(i,1).name).retro=transpose(AllNV.(list(i,1).name).retro);
            AllNV.(list3(j,1).name).retro=transpose(AllNV.(list3(j,1).name).retro);


            AllRL.(list(i,1).name).ante=transpose(AllRL.(list(i,1).name).ante);
            AllRL.(list3(j,1).name).ante=transpose(AllRL.(list3(j,1).name).ante);
            AllRL.(list(i,1).name).retro=transpose(AllRL.(list(i,1).name).retro);
            AllRL.(list3(j,1).name).retro=transpose(AllRL.(list3(j,1).name).retro);
            
            AllcomRL.(list(i,1).name).ante=transpose(AllcomRL.(list(i,1).name).ante);
            AllcomRL.(list3(j,1).name).ante=transpose(AllcomRL.(list3(j,1).name).ante);
            AllcomRL.(list(i,1).name).retro=transpose(AllcomRL.(list(i,1).name).retro);
            AllcomRL.(list3(j,1).name).retro=transpose(AllcomRL.(list3(j,1).name).retro);

            AllPD.(list(i,1).name).Value=transpose(AllPD.(list(i,1).name).Value);
            AllPD.(list3(j,1).name).Value=transpose(AllPD.(list3(j,1).name).Value);

            AllsplitPD.(list(i,1).name).ante=transpose(AllsplitPD.(list(i,1).name).ante);
            AllsplitPD.(list3(j,1).name).ante=transpose(AllsplitPD.(list3(j,1).name).ante);
            AllsplitPD.(list(i,1).name).retro=transpose(AllsplitPD.(list(i,1).name).retro);
            AllsplitPD.(list3(j,1).name).retro=transpose(AllsplitPD.(list3(j,1).name).retro);
            AllsplitPD.(list(i,1).name).rev=transpose(AllsplitPD.(list(i,1).name).rev);
            AllsplitPD.(list3(j,1).name).rev=transpose(AllsplitPD.(list3(j,1).name).rev);
            
            AllPF.(list(i,1).name).Value=transpose(AllPF.(list(i,1).name).Value);
            AllPF.(list3(j,1).name).Value=transpose(AllPF.(list3(j,1).name).Value);

            AllsplitPF.(list(i,1).name).ante=transpose(AllsplitPF.(list(i,1).name).ante);
            AllsplitPF.(list3(j,1).name).ante=transpose(AllsplitPF.(list3(j,1).name).ante);
            AllsplitPF.(list(i,1).name).retro=transpose(AllsplitPF.(list(i,1).name).retro);
            AllsplitPF.(list3(j,1).name).retro=transpose(AllsplitPF.(list3(j,1).name).retro);
            AllsplitPF.(list(i,1).name).rev=transpose(AllsplitPF.(list(i,1).name).rev);
            AllsplitPF.(list3(j,1).name).rev=transpose(AllsplitPF.(list3(j,1).name).rev);
            
            AllPFperSec.(list(i,1).name).Value=transpose(AllPFperSec.(list(i,1).name).Value);
            AllPFperSec.(list3(j,1).name).Value=transpose(AllPFperSec.(list3(j,1).name).Value);

            AllsplitPFperSec.(list(i,1).name).ante=transpose(AllsplitPFperSec.(list(i,1).name).ante);
            AllsplitPFperSec.(list3(j,1).name).ante=transpose(AllsplitPFperSec.(list3(j,1).name).ante);
            AllsplitPFperSec.(list(i,1).name).retro=transpose(AllsplitPFperSec.(list(i,1).name).retro);
            AllsplitPFperSec.(list3(j,1).name).retro=transpose(AllsplitPFperSec.(list3(j,1).name).retro);
            AllsplitPFperSec.(list(i,1).name).rev=transpose(AllsplitPFperSec.(list(i,1).name).rev);
            AllsplitPFperSec.(list3(j,1).name).rev=transpose(AllsplitPFperSec.(list3(j,1).name).rev);
            
            AllSF.(list(i,1).name).Value=transpose(AllSF.(list(i,1).name).Value);
            AllSF.(list3(j,1).name).Value=transpose(AllSF.(list3(j,1).name).Value);

            AllSFperSec.(list(i,1).name).Value=transpose(AllSFperSec.(list(i,1).name).Value);
            AllSFperSec.(list3(j,1).name).Value=transpose(AllSFperSec.(list3(j,1).name).Value);
            
            AllrevSF.(list(i,1).name).Value=transpose(AllrevSF.(list(i,1).name).Value);
            AllrevSF.(list3(j,1).name).Value=transpose(AllrevSF.(list3(j,1).name).Value);

            AllrevSFperSec.(list(i,1).name).Value=transpose(AllrevSFperSec.(list(i,1).name).Value);
            AllrevSFperSec.(list3(j,1).name).Value=transpose(AllrevSFperSec.(list3(j,1).name).Value);

            % try/catch loops replace "no value" entries with NaN.
            try
            pAnteCP=rndttest(AllCP.(list(i,1).name).PCT(1,:),AllCP.(list3(j,1).name).PCT(1,:));
            rnd.(comp).anteCP=pAnteCP;
            catch
                rnd.(comp).anteCP=NaN;
            end
          
            try
            pRetroCP=rndttest(AllCP.(list(i,1).name).PCT(2,:),AllCP.(list3(j,1).name).PCT(2,:));
            rnd.(comp).retroCP=pRetroCP;
            catch
                rnd.(comp).retroCP=NaN;
            end
            
            try
            pRevCP=rndttest(AllCP.(list(i,1).name).PCT(3,:),AllCP.(list3(j,1).name).PCT(3,:));
            rnd.(comp).revCP=pRevCP;
            catch
                rnd.(comp).revCP=NaN;
            end
            
            try
            pStatCP=rndttest(AllCP.(list(i,1).name).PCT(4,:),AllCP.(list3(j,1).name).PCT(4,:));
            rnd.(comp).statCP=pStatCP;
            catch
                rnd.(comp).statCP=NaN;
            end
            
            try
            pAnteNCP=rndttest(AllNCP.(list(i,1).name).PCT(1,:),AllNCP.(list3(j,1).name).PCT(1,:));
            rnd.(comp).anteNCP=pAnteNCP;
            catch
                rnd.(comp).anteNCP=NaN;
            end
          
            try
            pRetroNCP=rndttest(AllNCP.(list(i,1).name).PCT(2,:),AllNCP.(list3(j,1).name).PCT(2,:));
            rnd.(comp).retroNCP=pRetroNCP;
            catch
                rnd.(comp).retroNCP=NaN;
            end
                        
            try
            pStatNCP=rndttest(AllNCP.(list(i,1).name).PCT(3,:),AllNCP.(list3(j,1).name).PCT(3,:));
            rnd.(comp).statNCP=pStatNCP;
            catch
                rnd.(comp).statNCP=NaN;
            end
            
            try
            pAnteDensity=rndttest(AllDensity.(list(i,1).name).Value(1,:),AllDensity.(list3(j,1).name).Value(1,:));
            rnd.(comp).anteDensity=pAnteDensity;
            catch
                rnd.(comp).anteDensity=NaN;
            end
          
            try
            pRetroDensity=rndttest(AllDensity.(list(i,1).name).Value(2,:),AllDensity.(list3(j,1).name).Value(2,:));
            rnd.(comp).retroDensity=pRetroDensity;
            catch
                rnd.(comp).retroDensity=NaN;
            end
            
            try
            pRevDensity=rndttest(AllDensity.(list(i,1).name).Value(3,:),AllDensity.(list3(j,1).name).Value(3,:));
            rnd.(comp).revDensity=pRevDensity;
            catch
                rnd.(comp).revDensity=NaN;
            end
            
            try
            pStatDensity=rndttest(AllDensity.(list(i,1).name).Value(4,:),AllDensity.(list3(j,1).name).Value(4,:));
            rnd.(comp).statDensity=pStatDensity;
            catch
                rnd.(comp).statDensity=NaN;
            end
            
            %flux
            try
            pAnteFlux=rndttest(AllFlux.(list(i,1).name).Value(1,:),AllFlux.(list3(j,1).name).Value(1,:));
            rnd.(comp).anteFlux=pAnteFlux;
            catch
                rnd.(comp).anteFlux=NaN;
            end
          
            try
            pRetroFlux=rndttest(AllFlux.(list(i,1).name).Value(2,:),AllFlux.(list3(j,1).name).Value(2,:));
            rnd.(comp).retroFlux=pRetroFlux;
            catch
                rnd.(comp).retroFlux=NaN;
            end
            
            try
            pRevFlux=rndttest(AllFlux.(list(i,1).name).Value(3,:),AllFlux.(list3(j,1).name).Value(3,:));
            rnd.(comp).revFlux=pRevFlux;
            catch
                rnd.(comp).revFlux=NaN;
            end
            
            %SV
            try
            pAnteSV=rndttest(AllSV.(list(i,1).name).ante,AllSV.(list3(j,1).name).ante);
            rnd.(comp).anteSV=pAnteSV;
            catch ME
                rnd.(comp).anteSV=NaN;
            end

            try
            pRetroSV=rndttest(AllSV.(list(i,1).name).retro,AllSV.(list3(j,1).name).retro);
            rnd.(comp).retroSV=pRetroSV;
            catch ME
                rnd.(comp).retroSV=NaN;
            end
            
            % comSV
            try
            pAntecomSV=rndttest(AllcomSV.(list(i,1).name).ante,AllcomSV.(list3(j,1).name).ante);
            rnd.(comp).antecomSV=pAntecomSV;
            catch ME
                rnd.(comp).antecomSV=NaN;
            end

            try
            pRetrocomSV=rndttest(AllcomSV.(list(i,1).name).retro,AllcomSV.(list3(j,1).name).retro);
            rnd.(comp).retrocomSV=pRetrocomSV;       
            catch ME
                rnd.(comp).retrocomSV=NaN;
            end
            
            % NV
            try
            pAnteNV=rndttest(AllNV.(list(i,1).name).ante,AllNV.(list3(j,1).name).ante);
            rnd.(comp).anteNV=pAnteNV;
            catch ME
                rnd.(comp).anteNV=NaN;
            end

            try
            pRetroNV=rndttest(AllNV.(list(i,1).name).retro,AllNV.(list3(j,1).name).retro);
            rnd.(comp).retroNV=pRetroNV;
            catch ME
                rnd.(comp).retroNV=NaN;
            end

            try
            pAnteRL=rndttest(AllRL.(list(i,1).name).ante,AllRL.(list3(j,1).name).ante);
            rnd.(comp).anteRL=pAnteRL;
            catch ME
                rnd.(comp).anteRL=NaN;
            end

            try
            pRetroRL=rndttest(AllRL.(list(i,1).name).retro,AllRL.(list3(j,1).name).retro);
            rnd.(comp).retroRL=pRetroRL;       
            catch ME
                rnd.(comp).retroRL=NaN;
            end
            
            try
            pAntecomRL=rndttest(AllcomRL.(list(i,1).name).ante,AllcomRL.(list3(j,1).name).ante);
            rnd.(comp).antecomRL=pAntecomRL;
            catch ME
                rnd.(comp).antecomRL=NaN;
            end

            try
            pRetrocomRL=rndttest(AllcomRL.(list(i,1).name).retro,AllcomRL.(list3(j,1).name).retro);
            rnd.(comp).retrocomRL=pRetrocomRL;       
            catch ME
                rnd.(comp).retrocomRL=NaN;
            end

            try
            pPD=rndttest(AllPD.(list(i,1).name).Value,AllPD.(list3(j,1).name).Value);
            rnd.(comp).PD=pPD;
            catch ME
                rnd.(comp).PD=NaN;
            end

            try
            pPD=rndttest(AllsplitPD.(list(i,1).name).ante,AllsplitPD.(list3(j,1).name).ante);
            rnd.(comp).splitPDante=pPD;
            catch ME
                rnd.(comp).splitPDante=NaN;
            end
            
            try
            pPD=rndttest(AllsplitPD.(list(i,1).name).retro,AllsplitPD.(list3(j,1).name).retro);
            rnd.(comp).splitPDretro=pPD;
            catch ME
                rnd.(comp).splitPDretro=NaN;
            end
            
            try
            pPD=rndttest(AllsplitPD.(list(i,1).name).rev,AllsplitPD.(list3(j,1).name).rev);
            rnd.(comp).splitPDrev=pPD;
            catch ME
                rnd.(comp).splitPDrev=NaN;
            end
            
            try
            pPF=rndttest(AllPF.(list(i,1).name).Value,AllPF.(list3(j,1).name).Value);
            rnd.(comp).PF=pPF;
            catch ME
                rnd.(comp).PF=NaN;
            end

            try
            pPF=rndttest(AllsplitPF.(list(i,1).name).ante,AllsplitPF.(list3(j,1).name).ante);
            rnd.(comp).splitPFante=pPF;
            catch ME
                rnd.(comp).splitPFante=NaN;
            end
            
            try
            pPF=rndttest(AllsplitPF.(list(i,1).name).retro,AllsplitPF.(list3(j,1).name).retro);
            rnd.(comp).splitPFretro=pPF;
            catch ME
                rnd.(comp).splitPFretro=NaN;
            end
            
            try
            pPF=rndttest(AllsplitPF.(list(i,1).name).rev,AllsplitPF.(list3(j,1).name).rev);
            rnd.(comp).splitPFrev=pPF;
            catch ME
                rnd.(comp).splitPFrev=NaN;
            end
            
            try
            pPFperSec=rndttest(AllPFperSec.(list(i,1).name).Value,AllPFperSec.(list3(j,1).name).Value);
            rnd.(comp).PFperSec=pPFperSec;
            catch ME
                rnd.(comp).PFperSec=NaN;
            end

            try
            pPF=rndttest(AllsplitPFperSec.(list(i,1).name).ante,AllsplitPFperSec.(list3(j,1).name).ante);
            rnd.(comp).splitPFperSecante=pPF;
            catch ME
                rnd.(comp).splitPFperSecante=NaN;
            end
            
            try
            pPF=rndttest(AllsplitPFperSec.(list(i,1).name).retro,AllsplitPFperSec.(list3(j,1).name).retro);
            rnd.(comp).splitPFperSecretro=pPF;
            catch ME
                rnd.(comp).splitPFperSecretro=NaN;
            end
            
            try
            pPF=rndttest(AllsplitPFperSec.(list(i,1).name).rev,AllsplitPFperSec.(list3(j,1).name).rev);
            rnd.(comp).splitPFperSecrev=pPF;
            catch ME
                rnd.(comp).splitPFperSecrev=NaN;
            end
            
            try       
            pSF=rndttest(AllSF.(list(i,1).name).Value,AllSF.(list3(j,1).name).Value);
            rnd.(comp).SF=pSF;
            catch ME
                rnd.(comp).SF=NaN;
            end

            try       
            pSFperSec=rndttest(AllSFperSec.(list(i,1).name).Value,AllSFperSec.(list3(j,1).name).Value);
            rnd.(comp).SFperSec=pSFperSec;
            catch ME
                rnd.(comp).SFperSec=NaN;
            end
            
            try       
            prevSF=rndttest(AllrevSF.(list(i,1).name).Value,AllrevSF.(list3(j,1).name).Value);
            rnd.(comp).revSF=prevSF;
            catch ME
                rnd.(comp).revSF=NaN;
            end

            try       
            prevSFperSec=rndttest(AllrevSFperSec.(list(i,1).name).Value,AllrevSFperSec.(list3(j,1).name).Value);
            rnd.(comp).revSFperSec=prevSFperSec;
            catch ME
                rnd.(comp).revSFperSec=NaN;
            end

            %% transpose all Arrays into original form
            AllCP.(list(i,1).name).PCT=transpose(AllCP.(list(i,1).name).PCT);
            AllCP.(list3(j,1).name).PCT=transpose(AllCP.(list3(j,1).name).PCT);
            AllNCP.(list(i,1).name).PCT=transpose(AllNCP.(list(i,1).name).PCT);
            AllNCP.(list3(j,1).name).PCT=transpose(AllNCP.(list3(j,1).name).PCT);
            AllDensity.(list(i,1).name).Value=transpose(AllDensity.(list(i,1).name).Value);
            AllDensity.(list3(j,1).name).Value=transpose(AllDensity.(list3(j,1).name).Value);
            
            AllFlux.(list(i,1).name).Value=transpose(AllFlux.(list(i,1).name).Value);
            AllFlux.(list3(j,1).name).Value=transpose(AllFlux.(list3(j,1).name).Value);
            
            AllSV.(list(i,1).name).ante=transpose(AllSV.(list(i,1).name).ante);
            AllSV.(list3(j,1).name).ante=transpose(AllSV.(list3(j,1).name).ante);
            AllSV.(list(i,1).name).retro=transpose(AllSV.(list(i,1).name).retro);
            AllSV.(list3(j,1).name).retro=transpose(AllSV.(list3(j,1).name).retro);

            AllNV.(list(i,1).name).ante=transpose(AllNV.(list(i,1).name).ante);
            AllNV.(list3(j,1).name).ante=transpose(AllNV.(list3(j,1).name).ante);
            AllNV.(list(i,1).name).retro=transpose(AllNV.(list(i,1).name).retro);
            AllNV.(list3(j,1).name).retro=transpose(AllNV.(list3(j,1).name).retro);


            AllRL.(list(i,1).name).ante=transpose(AllRL.(list(i,1).name).ante);
            AllRL.(list3(j,1).name).ante=transpose(AllRL.(list3(j,1).name).ante);
            AllRL.(list(i,1).name).retro=transpose(AllRL.(list(i,1).name).retro);
            AllRL.(list3(j,1).name).retro=transpose(AllRL.(list3(j,1).name).retro);
            
            AllcomRL.(list(i,1).name).ante=transpose(AllcomRL.(list(i,1).name).ante);
            AllcomRL.(list3(j,1).name).ante=transpose(AllcomRL.(list3(j,1).name).ante);
            AllcomRL.(list(i,1).name).retro=transpose(AllcomRL.(list(i,1).name).retro);
            AllcomRL.(list3(j,1).name).retro=transpose(AllcomRL.(list3(j,1).name).retro);

            AllPD.(list(i,1).name).Value=transpose(AllPD.(list(i,1).name).Value);
            AllPD.(list3(j,1).name).Value=transpose(AllPD.(list3(j,1).name).Value);

            AllsplitPD.(list(i,1).name).ante=transpose(AllsplitPD.(list(i,1).name).ante);
            AllsplitPD.(list3(j,1).name).ante=transpose(AllsplitPD.(list3(j,1).name).ante);
            AllsplitPD.(list(i,1).name).retro=transpose(AllsplitPD.(list(i,1).name).retro);
            AllsplitPD.(list3(j,1).name).retro=transpose(AllsplitPD.(list3(j,1).name).retro);
            AllsplitPD.(list(i,1).name).rev=transpose(AllsplitPD.(list(i,1).name).rev);
            AllsplitPD.(list3(j,1).name).rev=transpose(AllsplitPD.(list3(j,1).name).rev);
            
            AllPF.(list(i,1).name).Value=transpose(AllPF.(list(i,1).name).Value);
            AllPF.(list3(j,1).name).Value=transpose(AllPF.(list3(j,1).name).Value);

            AllsplitPF.(list(i,1).name).ante=transpose(AllsplitPF.(list(i,1).name).ante);
            AllsplitPF.(list3(j,1).name).ante=transpose(AllsplitPF.(list3(j,1).name).ante);
            AllsplitPF.(list(i,1).name).retro=transpose(AllsplitPF.(list(i,1).name).retro);
            AllsplitPF.(list3(j,1).name).retro=transpose(AllsplitPF.(list3(j,1).name).retro);
            AllsplitPF.(list(i,1).name).rev=transpose(AllsplitPF.(list(i,1).name).rev);
            AllsplitPF.(list3(j,1).name).rev=transpose(AllsplitPF.(list3(j,1).name).rev);
            
            AllPFperSec.(list(i,1).name).Value=transpose(AllPFperSec.(list(i,1).name).Value);
            AllPFperSec.(list3(j,1).name).Value=transpose(AllPFperSec.(list3(j,1).name).Value);

            AllsplitPFperSec.(list(i,1).name).ante=transpose(AllsplitPFperSec.(list(i,1).name).ante);
            AllsplitPFperSec.(list3(j,1).name).ante=transpose(AllsplitPFperSec.(list3(j,1).name).ante);
            AllsplitPFperSec.(list(i,1).name).retro=transpose(AllsplitPFperSec.(list(i,1).name).retro);
            AllsplitPFperSec.(list3(j,1).name).retro=transpose(AllsplitPFperSec.(list3(j,1).name).retro);
            AllsplitPFperSec.(list(i,1).name).rev=transpose(AllsplitPFperSec.(list(i,1).name).rev);
            AllsplitPFperSec.(list3(j,1).name).rev=transpose(AllsplitPFperSec.(list3(j,1).name).rev);
            
            AllSF.(list(i,1).name).Value=transpose(AllSF.(list(i,1).name).Value);
            AllSF.(list3(j,1).name).Value=transpose(AllSF.(list3(j,1).name).Value);

            AllSFperSec.(list(i,1).name).Value=transpose(AllSFperSec.(list(i,1).name).Value);
            AllSFperSec.(list3(j,1).name).Value=transpose(AllSFperSec.(list3(j,1).name).Value);
            
            AllrevSF.(list(i,1).name).Value=transpose(AllrevSF.(list(i,1).name).Value);
            AllrevSF.(list3(j,1).name).Value=transpose(AllrevSF.(list3(j,1).name).Value);

            AllrevSFperSec.(list(i,1).name).Value=transpose(AllrevSFperSec.(list(i,1).name).Value);
            AllrevSFperSec.(list3(j,1).name).Value=transpose(AllrevSFperSec.(list3(j,1).name).Value);
            
        end
    end
end
    
%% write txt files for generation of boxplots on: http://boxplot.tyerslab.com/
% write data into arrays for printing into txt files


for i=1:length(list)
    allAnteSV{i}=AllSV.(list(i,:).name).ante;
    allRetroSV{i}=AllSV.(list(i,:).name).retro;
    allAnteNV{i}=AllNV.(list(i,:).name).ante;
    allRetroNV{i}=AllNV.(list(i,:).name).retro;
    allAnteRL{i}=AllRL.(list(i,:).name).ante;
    allRetroRL{i}=AllRL.(list(i,:).name).retro;
    allAntecomRL{i}=AllcomRL.(list(i,:).name).ante;
    allRetrocomRL{i}=AllcomRL.(list(i,:).name).retro;
    allPD{i}=AllPD.(list(i,:).name).Value;
    allsplitPDante{i}=AllsplitPD.(list(i,:).name).ante;
    allsplitPDretro{i}=AllsplitPD.(list(i,:).name).retro;
    allsplitPDrev{i}=AllsplitPD.(list(i,:).name).rev;
    allPF{i}=AllPF.(list(i,:).name).Value;
    allsplitPFante{i}=AllsplitPF.(list(i,:).name).ante;
    allsplitPFretro{i}=AllsplitPF.(list(i,:).name).retro;
    allsplitPFrev{i}=AllsplitPF.(list(i,:).name).rev;
    allPFperSec{i}=AllPFperSec.(list(i,:).name).Value;
    allsplitPFperSecante{i}=AllsplitPFperSec.(list(i,:).name).ante;
    allsplitPFperSecretro{i}=AllsplitPFperSec.(list(i,:).name).retro;
    allsplitPFperSecrev{i}=AllsplitPFperSec.(list(i,:).name).rev;
    allSF{i}=AllSF.(list(i,:).name).Value;
    allSFperSec{i}=AllSFperSec.(list(i,:).name).Value;
    allrevSF{i}=AllrevSF.(list(i,:).name).Value;
    allrevSFperSec{i}=AllrevSFperSec.(list(i,:).name).Value;
end

if length(list)>1
    anteSV=padcat(allAnteSV{:});
    retroSV=padcat(allRetroSV{:});
    anteNV=padcat(allAnteNV{:});
    retroNV=padcat(allRetroNV{:});
    anteRL=padcat(allAnteRL{:});
    retroRL=padcat(allRetroRL{:});
    antecomRL=padcat(allAntecomRL{:});
    retrocomRL=padcat(allRetrocomRL{:});
    pd=padcat(allPD{:});
    pf=padcat(allPF{:});
    pfperSec=padcat(allPFperSec{:});   
    splitpdante=padcat(allsplitPDante{:});
    splitpdretro=padcat(allsplitPDretro{:});
    splitpdrev=padcat(allsplitPDrev{:});
    splitpfante=padcat(allsplitPFante{:});
    splitpfretro=padcat(allsplitPFretro{:});
    splitpfrev=padcat(allsplitPFrev{:});
    splitpfperSecante=padcat(allsplitPFperSecante{:});
    splitpfperSecretro=padcat(allsplitPFperSecretro{:});
    splitpfperSecrev=padcat(allsplitPFperSecrev{:});
    sf=padcat(allSF{:});
    sfperSec=padcat(allSFperSec{:});
    revsf=padcat(allrevSF{:});
    revsfperSec=padcat(allrevSFperSec{:});

else 
    anteSV=allAnteSV{:};
    retroSV=allRetroSV{:};
    anteNV=allAnteNV{:};
    retroNV=allRetroNV{:};
    anteRL=allAnteRL{:};
    retroRL=allRetroRL{:};
    antecomRL=allAntecomRL{:};
    retrocomRL=allRetrocomRL{:};
    pd=allPD{:};
    pf=allPF{:};
    pfperSec=allPFperSec{:};
    splitpdante=allsplitPDante{:};
    splitpdretro=allsplitPDretro{:};
    splitpdrev=allsplitPDrev{:};
    splitpfante=allsplitPFante{:};
    splitpfretro=allsplitPFretro{:};
    splitpfrev=allsplitPFrev{:};
    splitpfperSecante=allsplitPFperSecante{:};
    splitpfperSecretro=allsplitPFperSecretro{:};
    splitpfperSecrev=allsplitPFperSecrev{:};
    sf=allSF{:};
    sfperSec=allSFperSec{:};
    revsf=allrevSF{:};
    revsfperSec=allrevSFperSec{:};
end

categories={};
for i=1:length(list)
    categories(:,i)={(list(i,:).name)};
end

%%  filenames
mkdir(fullfile(directoryMain,'Figures','txtfiles'));
fileName=fullfile(directoryMain,'Figures','txtfiles');
fileName1=strcat(fileName,filesep,'AnteSegVel.txt');
fileName2=strcat(fileName,filesep,'RetroSegVel.txt');
fileName3=strcat(fileName,filesep,'AnteNetVel.txt');
fileName4=strcat(fileName,filesep,'RetroNetVel.txt');
fileName5=strcat(fileName,filesep,'AnteRL.txt');
fileName6=strcat(fileName,filesep,'RetroRL.txt');
fileName7=strcat(fileName,filesep,'AntecomRL.txt');
fileName8=strcat(fileName,filesep,'RetrocomRL.txt');
fileName9=strcat(fileName,filesep,'PD.txt'); 
fileName10=strcat(fileName,filesep,'PF.txt'); 
fileName11=strcat(fileName,filesep,'PFperSec.txt'); 
fileName12=strcat(fileName,filesep,'SF.txt'); 
fileName13=strcat(fileName,filesep,'SFperSec.txt'); 
fileName14=strcat(fileName,filesep,'revSF.txt'); 
fileName15=strcat(fileName,filesep,'revSFperSec.txt'); 
fileName16=strcat(fileName,filesep,'Ns.txt'); 
fileName17=strcat(fileName,filesep,'splitPDante.txt');
fileName18=strcat(fileName,filesep,'splitPDretro.txt');
fileName19=strcat(fileName,filesep,'splitPDrev.txt');
fileName20=strcat(fileName,filesep,'splitPFante.txt');
fileName21=strcat(fileName,filesep,'splitPFretro.txt');
fileName22=strcat(fileName,filesep,'splitPFrev.txt');
fileName23=strcat(fileName,filesep,'splitPFperSecante.txt');
fileName24=strcat(fileName,filesep,'splitPFperSecretro.txt');
fileName25=strcat(fileName,filesep,'splitPFperSecrev.txt');

%% write into files
    % ante SV
    fid1=fopen(fileName1,'w');
    fprintf(fid1,'%s\t',categories{1,:});
    fprintf(fid1,'\n');

    for i=1:size(anteSV,1)
        fprintf(fid1,'%f\t',anteSV(i,:));
        fprintf(fid1,'\n');
    end
    fclose(fid1);

    % retro SV
    fid2=fopen(fileName2,'w');
    fprintf(fid2,'%s\t',categories{1,:});
    fprintf(fid2,'\n');

    for i=1:size(retroSV,1)
        fprintf(fid2,'%f\t',retroSV(i,:));
        fprintf(fid2,'\n');
    end
    fclose(fid2);

    % ante NV
    fid3=fopen(fileName3,'w');
    fprintf(fid3,'%s\t',categories{1,:});
    fprintf(fid3,'\n');

    for i=1:size(anteNV,1)
        fprintf(fid3,'%f\t',anteNV(i,:));
        fprintf(fid3,'\n');
    end
    fclose(fid3);

    % retro NV
    fid4=fopen(fileName4,'w');
    fprintf(fid4,'%s\t',categories{1,:});
    fprintf(fid4,'\n');

    for i=1:size(retroNV,1)
        fprintf(fid4,'%f\t',retroNV(i,:));
        fprintf(fid4,'\n');
    end
    fclose(fid4);

    % ante RL
    fid5=fopen(fileName5,'w');
    fprintf(fid5,'%s\t',categories{1,:});
    fprintf(fid5,'\n');

    for i=1:size(anteRL,1)
        fprintf(fid5,'%f\t',anteRL(i,:));
        fprintf(fid5,'\n');
    end
    fclose(fid5);

    % retro RL
    fid6=fopen(fileName6,'w');
    fprintf(fid6,'%s\t',categories{1,:});
    fprintf(fid6,'\n');

    for i=1:size(retroRL,1)
        fprintf(fid6,'%f\t',retroRL(i,:));
        fprintf(fid6,'\n');
    end
    fclose(fid6);

    % ante combined RL
    fid7=fopen(fileName7,'w');
    fprintf(fid7,'%s\t',categories{1,:});
    fprintf(fid7,'\n');

    for i=1:size(antecomRL,1)
        fprintf(fid7,'%f\t',antecomRL(i,:));
        fprintf(fid7,'\n');
    end
    fclose(fid7);

    % retro combined RL
    fid8=fopen(fileName8,'w');
    fprintf(fid8,'%s\t',categories{1,:});
    fprintf(fid8,'\n');

    for i=1:size(retrocomRL,1)
        fprintf(fid8,'%f\t',retrocomRL(i,:));
        fprintf(fid8,'\n');
    end
    fclose(fid8);

% PD
fid9=fopen(fileName9,'w');
fprintf(fid9,'%s\t',categories{1,:});
fprintf(fid9,'\n');
 
for i=1:size(pd,1)
    fprintf(fid9,'%f\t',pd(i,:));
    fprintf(fid9,'\n');
end
fclose(fid9);

% splitPDante
fid17=fopen(fileName17,'w');
fprintf(fid17,'%s\t',categories{1,:});
fprintf(fid17,'\n');
 
for i=1:size(splitpdante,1)
    fprintf(fid17,'%f\t',splitpdante(i,:));
    fprintf(fid17,'\n');
end
fclose(fid17);

% splitPDretro
fid18=fopen(fileName18,'w');
fprintf(fid18,'%s\t',categories{1,:});
fprintf(fid18,'\n');
 
for i=1:size(splitpdretro,1)
    fprintf(fid18,'%f\t',splitpdretro(i,:));
    fprintf(fid18,'\n');
end
fclose(fid18);

% splitPDrev
fid19=fopen(fileName19,'w');
fprintf(fid19,'%s\t',categories{1,:});
fprintf(fid19,'\n');
 
for i=1:size(splitpdrev,1)
    fprintf(fid19,'%f\t',splitpdrev(i,:));
    fprintf(fid19,'\n');
end
fclose(fid19);

% PF
fid10=fopen(fileName10,'w');
fprintf(fid10,'%s\t',categories{1,:});
fprintf(fid10,'\n');
 
for i=1:size(pf,1)
    fprintf(fid10,'%f\t',pf(i,:));
    fprintf(fid10,'\n');
end
fclose(fid10);

% splitPFante
fid20=fopen(fileName20,'w');
fprintf(fid20,'%s\t',categories{1,:});
fprintf(fid20,'\n');
 
for i=1:size(splitpfante,1)
    fprintf(fid20,'%f\t',splitpfante(i,:));
    fprintf(fid20,'\n');
end
fclose(fid20);

% splitPFretro
fid21=fopen(fileName21,'w');
fprintf(fid21,'%s\t',categories{1,:});
fprintf(fid21,'\n');
 
for i=1:size(splitpfretro,1)
    fprintf(fid21,'%f\t',splitpfretro(i,:));
    fprintf(fid21,'\n');
end
fclose(fid21);

% splitPFrev
fid22=fopen(fileName22,'w');
fprintf(fid22,'%s\t',categories{1,:});
fprintf(fid22,'\n');
 
for i=1:size(splitpfrev,1)
    fprintf(fid22,'%f\t',splitpfrev(i,:));
    fprintf(fid22,'\n');
end
fclose(fid22);

% PF per Sec
fid11=fopen(fileName11,'w');
fprintf(fid11,'%s\t',categories{1,:});
fprintf(fid11,'\n');
 
for i=1:size(pfperSec,1)
    fprintf(fid11,'%f\t',pfperSec(i,:));
    fprintf(fid11,'\n');
end
fclose(fid11);

% splitPFperSecante
fid23=fopen(fileName23,'w');
fprintf(fid23,'%s\t',categories{1,:});
fprintf(fid23,'\n');
 
for i=1:size(splitpfperSecante,1)
    fprintf(fid23,'%f\t',splitpfperSecante(i,:));
    fprintf(fid23,'\n');
end
fclose(fid23);

% splitPFperSecretro
fid24=fopen(fileName24,'w');
fprintf(fid24,'%s\t',categories{1,:});
fprintf(fid24,'\n');
 
for i=1:size(splitpfperSecretro,1)
    fprintf(fid24,'%f\t',splitpfperSecretro(i,:));
    fprintf(fid24,'\n');
end
fclose(fid24);

% splitPFperSecrev
fid25=fopen(fileName25,'w');
fprintf(fid25,'%s\t',categories{1,:});
fprintf(fid25,'\n');
 
for i=1:size(splitpfperSecrev,1)
    fprintf(fid25,'%f\t',splitpfperSecrev(i,:));
    fprintf(fid25,'\n');
end
fclose(fid25);

% SF
fid12=fopen(fileName12,'w');
fprintf(fid12,'%s\t',categories{1,:});
fprintf(fid12,'\n');
 
for i=1:size(sf,1)
    fprintf(fid12,'%f\t',sf(i,:));
    fprintf(fid12,'\n');
end
fclose(fid12);

% SF per Sec
fid13=fopen(fileName13,'w');
fprintf(fid13,'%s\t',categories{1,:});
fprintf(fid13,'\n');
 
for i=1:size(sfperSec,1)
    fprintf(fid13,'%f\t',sfperSec(i,:));
    fprintf(fid13,'\n');
end
fclose(fid13);

% reversal SF
fid14=fopen(fileName14,'w');
fprintf(fid12,'%s\t',categories{1,:});
fprintf(fid12,'\n');
 
for i=1:size(revsf,1)
    fprintf(fid14,'%f\t',revsf(i,:));
    fprintf(fid14,'\n');
end
fclose(fid14);

% reversal SF per Sec
fid15=fopen(fileName15,'w');
fprintf(fid15,'%s\t',categories{1,:});
fprintf(fid15,'\n');
 
for i=1:size(revsfperSec,1)
    fprintf(fid15,'%f\t',revsfperSec(i,:));
    fprintf(fid15,'\n');
end
fclose(fid15);

% print all numbers "n" for each parameter
fid16=fopen(fileName16,'w');
fprintf(fid16,'%s\t',' ');
fprintf(fid16,'%s\t',categories{1,:});
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','movies');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(size(AllCP.(list(i,1).name).Num,1))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','ante SV');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllSV.(list(i,1).name).ante))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','retro SV');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllSV.(list(i,1).name).retro))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','ante NV');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllNV.(list(i,1).name).ante))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','retro NV');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllNV.(list(i,1).name).retro))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','ante RL');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllRL.(list(i,1).name).ante))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','retro RL');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllRL.(list(i,1).name).retro))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','ante comRL');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllcomRL.(list(i,1).name).ante))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','retro comRL');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllcomRL.(list(i,1).name).retro))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','PD');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllPD.(list(i,1).name).Value))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','splitPDante');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllsplitPD.(list(i,1).name).ante))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','splitPDretro');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllsplitPD.(list(i,1).name).retro))); 
end
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','splitPDrev');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllsplitPD.(list(i,1).name).rev))); 
end
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','splitPFante');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllsplitPF.(list(i,1).name).ante))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','splitPFretro');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllsplitPF.(list(i,1).name).retro))); 
end
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','splitPFrev');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllsplitPF.(list(i,1).name).rev))); 
end
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','splitPFperSecante');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllsplitPFperSec.(list(i,1).name).ante))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','splitPFperSecretro');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllsplitPFperSec.(list(i,1).name).retro))); 
end
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','splitPFperSecrev');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllsplitPFperSec.(list(i,1).name).rev))); 
end
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','PF');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllPF.(list(i,1).name).Value))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','PFperSec');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllPFperSec.(list(i,1).name).Value))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','SF');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllSF.(list(i,1).name).Value))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','SFperSec');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllSFperSec.(list(i,1).name).Value))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','revSF');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllrevSF.(list(i,1).name).Value))); 
end 
fprintf(fid16,'\n');
fprintf(fid16,'%s\t','SFperSec');
for i=1:length(list)
   fprintf(fid16,'%s\t',num2str(length(AllrevSFperSec.(list(i,1).name).Value))); 
end 
fclose(fid16);

%% Save all Variables
[pathstr,nameMain] = fileparts(directoryMain);
filename=strcat(nameMain,'.m');
file=fullfile(directoryMain,'MatlabResults',filename);

a = exist('rnk')& exist('rnd');
if a == 1
    save(file,'AllCP','AllDensity','AllNCP','AllNV','AllPD','AllPF','AllPFperSec','AllsplitPD','AllsplitPF','AllsplitPFperSec','AllRL','AllcomRL','AllSF','AllSFperSec','AllrevSF','AllrevSFperSec','AllSV',...
     'CP','Density','NCP','NV','PD','PF','PFperSec','SF','SFperSec','revSF','revSFperSec','SV','RL','comRL','histNV','histRL','histcomRL','histSV',...
     'directoryMain','MultiGaussian','color','normal','rnk','rnd','Flux','AllFlux','AllcomSV');
else
    save(file,'AllCP','AllDensity','AllNCP','AllNV','AllPD','AllPF','AllPFperSec','AllsplitPD','AllsplitPF','AllsplitPFperSec','AllRL','AllcomRL','AllSF','AllSFperSec','AllrevSF','AllrevSFperSec','AllSV',...
        'CP','Density','NCP','NV','PD','PF','PFperSec','splitPD','splitPF','splitPFperSec','SF','SFperSec','revSF','revSFperSec','SV','RL','comRL','histNV','histRL','histcomRL','histSV',...
        'directoryMain','MultiGaussian','color','normal','Flux','AllFlux','AllcomSV');
end

% 

%clear Matlab
%clear; clc;
