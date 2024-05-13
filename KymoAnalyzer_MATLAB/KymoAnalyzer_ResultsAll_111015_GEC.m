function KymoAnalyzer_ResultsAll_111015_GEC(directoryMain)

%% Script to read in data from KymoAnalyzer1.0 and to calculate statistics and make figures
% Written by Sylvia Neumann 10/31/14
% Edited by George Campbell 11/10/15

% Script dependencies:  rotateXLabels.m

%%
% Clear Matlab
%clc;clear;

%% Choose experimental folder (Folder that contains all movies and pooled data
% remove hidden MAC OSX files
%directoryMain= uigetdir(pwd,'Select your parent folder (Folder that contains all experiments for all genotypes)');
    cd (directoryMain);
    temp = dir(directoryMain);

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



    list=([]);
    j=1;
    for i=1:length(temp)
        if strcmp(temp(i,1).name,'Figures')==1 || strcmp(temp(i,1).name,'MatlabResults') == 1
        % ? 
        else
        list(j,1).name=temp(i,1).name;
        j=j+1;
        end
    end

    for i=1:length(list)
        subfolder=fullfile(directoryMain,list(i,1).name,filesep);
        temp2=dir(subfolder);

        for k = length(temp2):-1:1
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

        k=1;
        list2=([]);
        for j=1:length(temp2)
            list2(k,1).name=temp2(j,1).name;
            k=k+1;
        end

        for j=1:length(list2)
            directory=fullfile(directoryMain,list(i,1).name,list2(j,1).name,filesep);
            cd(directory);


    %% Importing DATA        
            %% read in data from CargoPopulation 
            % CP is a structured variable that contain the numeric or the percentage
            % values for cargo population (anterograde, retrograde, reversal and
            % stationary. CP.Stats contains averages (row 1) and standard deviation (row 2).

            CP=([]);
            anteNum=fullfile(directory,'PooledData','CargoPopulation','anteNum.txt');
            CP.Num(:,1)=importdata(anteNum);
            retroNum=fullfile(directory,'PooledData','CargoPopulation','retroNum.txt');
            CP.Num(:,2)=importdata(retroNum);
            reversalNum=fullfile(directory,'PooledData','CargoPopulation','reversalNum.txt');
            CP.Num(:,3)=importdata(reversalNum);
            stationaryNum=fullfile(directory,'PooledData','CargoPopulation','stationaryNum.txt');
            CP.Num(:,4)=importdata(stationaryNum);

            antePCT=fullfile(directory,'PooledData','CargoPopulation','antePCT.txt');
            CP.PCT(:,1)=importdata(antePCT);
            retroPCT=fullfile(directory,'PooledData','CargoPopulation','retroPCT.txt');
            CP.PCT(:,2)=importdata(retroPCT);
            reversalPCT=fullfile(directory,'PooledData','CargoPopulation','reversalPCT.txt');
            CP.PCT(:,3)=importdata(reversalPCT);
            stationaryPCT=fullfile(directory,'PooledData','CargoPopulation','stationaryPCT.txt');
            CP.PCT(:,4)=importdata(stationaryPCT);

            CP.StatsPCT(1,1)=100*mean(CP.PCT(:,1));
            CP.StatsPCT(1,2)=100*mean(CP.PCT(:,2));
            CP.StatsPCT(1,3)=100*mean(CP.PCT(:,3));
            CP.StatsPCT(1,4)=100*mean(CP.PCT(:,4));

            %standard deviation
            CP.StatsPCT(2,1)=100*std(CP.PCT(:,1));
            CP.StatsPCT(2,2)=100*std(CP.PCT(:,2));
            CP.StatsPCT(2,3)=100*std(CP.PCT(:,3));
            CP.StatsPCT(2,4)=100*std(CP.PCT(:,4));

            % standard error of mean
            CP.StatsPCT(3,1)=100*std(CP.PCT(:,1))/sqrt(length(CP.PCT(:,1)));
            CP.StatsPCT(3,2)=100*std(CP.PCT(:,2))/sqrt(length(CP.PCT(:,2)));
            CP.StatsPCT(3,3)=100*std(CP.PCT(:,3))/sqrt(length(CP.PCT(:,3)));
            CP.StatsPCT(3,4)=100*std(CP.PCT(:,4))/sqrt(length(CP.PCT(:,4)));

            %% read in data from Density
            % Density is a structured variable that contain the 
            % values for # of tracks / um axon (anterograde, retrograde, reversal and
            % stationary). Density.Stats contains averages (row 1) and standard deviation (row 2).

            Density=([]);
            anteDensity=fullfile(directory,'PooledData','Density','anteDensity.txt');
            Density.Value(:,1)=importdata(anteDensity);
            retroDensity=fullfile(directory,'PooledData','Density','retroDensity.txt');
            Density.Value(:,2)=importdata(retroDensity);
            reversalDensity=fullfile(directory,'PooledData','Density','reversalDensity.txt');
            Density.Value(:,3)=importdata(reversalDensity);
            stationaryDensity=fullfile(directory,'PooledData','Density','stationaryDensity.txt');
            Density.Value(:,4)=importdata(stationaryDensity);

            Density.Stats(1,1)=mean(Density.Value(:,1));
            Density.Stats(1,2)=mean(Density.Value(:,2));
            Density.Stats(1,3)=mean(Density.Value(:,3));
            Density.Stats(1,4)=mean(Density.Value(:,4));

            Density.Stats(2,1)=std(Density.Value(:,1));
            Density.Stats(2,2)=std(Density.Value(:,2));
            Density.Stats(2,3)=std(Density.Value(:,3));
            Density.Stats(2,4)=std(Density.Value(:,4));

            Density.Stats(3,1)=std(Density.Value(:,1))/sqrt(length(Density.Value(:,1)));
            Density.Stats(3,2)=std(Density.Value(:,2))/sqrt(length(Density.Value(:,2)));
            Density.Stats(3,3)=std(Density.Value(:,3))/sqrt(length(Density.Value(:,3)));
            Density.Stats(3,4)=std(Density.Value(:,4))/sqrt(length(Density.Value(:,4)));

            %% read in data from Flux
            % Flux is a structured variable that contain the 
            % values for # of tracks/ (um axon * sec) (anterograde, retrograde, reversal). 
            % Flux.Stats contains averages (row 1) and standrd deviation (row 2).

            Flux=([]);
            anteFlux=fullfile(directory,'PooledData','Flux','anteFlux.txt');
            Flux.Value(:,1)=importdata(anteFlux);
            retroFlux=fullfile(directory,'PooledData','Flux','retroFlux.txt');
            Flux.Value(:,2)=importdata(retroFlux);
            reversalFlux=fullfile(directory,'PooledData','Flux','reversalFlux.txt');
            Flux.Value(:,3)=importdata(reversalFlux);

            Flux.Stats(1,1)=mean(Flux.Value(:,1));
            Flux.Stats(1,2)=mean(Flux.Value(:,2));
            Flux.Stats(1,3)=mean(Flux.Value(:,3));

            Flux.Stats(2,1)=std(Flux.Value(:,1));
            Flux.Stats(2,2)=std(Flux.Value(:,2));
            Flux.Stats(2,3)=std(Flux.Value(:,3));

            Flux.Stats(3,1)=std(Flux.Value(:,1))/sqrt(length(Flux.Value(:,1)));
            Flux.Stats(3,2)=std(Flux.Value(:,2))/sqrt(length(Flux.Value(:,2)));
            Flux.Stats(3,3)=std(Flux.Value(:,3))/sqrt(length(Flux.Value(:,3)));

            %% read in data from  NetCargoPopulation
            % NCP is a structured variable that contain the numeric or the percentage
            % values for net cargo population (anterograde, retrograde and
            % stationary. NCP.Stats contains averages (row 1)and standard deviation (row 2).

            NCP=([]);
            NetanteNum=fullfile(directory,'PooledData','NetCargoPopulation','NetanteNum.txt');
            NCP.Num(:,1)=importdata(NetanteNum);
            NetretroNum=fullfile(directory,'PooledData','NetCargoPopulation','NetretroNum.txt');
            NCP.Num(:,2)=importdata(NetretroNum);
            NetstationaryNum=fullfile(directory,'PooledData','NetCargoPopulation','NetstationaryNum.txt');
            NCP.Num(:,3)=importdata(NetstationaryNum);

            NetantePCT=fullfile(directory,'PooledData','NetCargoPopulation','NetantePCT.txt');
            NCP.PCT(:,1)=importdata(NetantePCT);
            NetretroPCT=fullfile(directory,'PooledData','NetCargoPopulation','NetretroPCT.txt');
            NCP.PCT(:,2)=importdata(NetretroPCT);
            NetstationaryPCT=fullfile(directory,'PooledData','NetCargoPopulation','NetstationaryPCT.txt');
            NCP.PCT(:,3)=importdata(NetstationaryPCT);

            NCP.StatsPCT(1,1)=100*mean(NCP.PCT(:,1));
            NCP.StatsPCT(1,2)=100*mean(NCP.PCT(:,2));
            NCP.StatsPCT(1,3)=100*mean(NCP.PCT(:,3));
            NCP.StatsPCT(1,4)=0;

            NCP.StatsPCT(2,1)=100*std(NCP.PCT(:,1));
            NCP.StatsPCT(2,2)=100*std(NCP.PCT(:,2));
            NCP.StatsPCT(2,3)=100*std(NCP.PCT(:,3));
            NCP.StatsPCT(2,4)=0;

            NCP.StatsPCT(3,1)=100*std(NCP.PCT(:,1))/sqrt(length(NCP.PCT(:,1)));
            NCP.StatsPCT(3,2)=100*std(NCP.PCT(:,2))/sqrt(length(NCP.PCT(:,2)));
            NCP.StatsPCT(3,3)=100*std(NCP.PCT(:,3))/sqrt(length(NCP.PCT(:,3)));
            NCP.StatsPCT(3,4)=0;

            %% Importing Data from Net Velocities
            % NV is a structures variable that contains all NV (anterograde and
            % retrograde) for the entire experiment. NV.Stats contains averages (row 1)
            % and standard deviation (row 2) for anterograde (column 1) and retrograde (column 2) velocities

            NV=([]);
            anteNV=fullfile(directory,'PooledData','NetVelocities','AllAnteNV.txt');
            NV.ante(:,1)=importdata(anteNV);
            retroNV=fullfile(directory,'PooledData','NetVelocities','AllRetroNV.txt');
            NV.retro(:,1)=importdata(retroNV);

            NV.Stats=zeros(2,4);
            NV.Stats(1,1)=mean(NV.ante);
            NV.Stats(1,2)=mean(NV.retro);

            NV.Stats(2,1)=std(NV.ante);
            NV.Stats(2,2)=std(NV.retro);

            NV.Stats(3,1)=std(NV.ante)/sqrt(length(NV.ante));
            NV.Stats(3,2)=std(NV.retro)/sqrt(length(NV.ante));

            %% Importing Data from Pause Duration
            % PD is a structures variable that contains all Pause Durations for the entire experiment. PD.Stats contains averages (row 1)
            % and standard deviation (row 2).

            PD=([]);
            pd=fullfile(directory,'PooledData','PauseDuration','AllPD.txt');
            PD.Value(:,1)=importdata(pd);

            PD.Stats=zeros(2,4);
            PD.Stats(1,1)=mean(PD.Value);
            PD.Stats(2,1)=std(PD.Value);
            PD.Stats(3,1)=std(PD.Value)/sqrt(length(PD.Value));

            %% Importing Data from Split Pause Duration
            % PD is a structures variable that contains all Pause Durations for the entire experiment. splitPD.Stats contains averages (row 1)
            % and standard deviation (row 2) and SEM (row 3) for anterograde (column 1), retrograde (column 2) and reversal (column 3).

            splitPD=([]);
            pd=fullfile(directory,'PooledData','SplitPauseDuration','AllantePD.txt');
            splitPD.ante(:,1)=importdata(pd);
            pd=fullfile(directory,'PooledData','SplitPauseDuration','AllretroPD.txt');
            splitPD.retro(:,1)=importdata(pd);
            pd=fullfile(directory,'PooledData','SplitPauseDuration','AllrevPD.txt');
            splitPD.rev(:,1)=importdata(pd);

            splitPD.Stats=zeros(2,4);
            splitPD.Stats(1,1)=mean(splitPD.ante);
            splitPD.Stats(1,2)=mean(splitPD.retro);
            splitPD.Stats(1,3)=mean(splitPD.rev);

            splitPD.Stats(2,1)=std(splitPD.ante);
            splitPD.Stats(2,2)=std(splitPD.retro);
            splitPD.Stats(2,3)=std(splitPD.rev);

            splitPD.Stats(3,1)=std(splitPD.ante)/sqrt(length(splitPD.ante));
            splitPD.Stats(3,2)=std(splitPD.retro)/sqrt(length(splitPD.retro));
            splitPD.Stats(3,3)=std(splitPD.rev)/sqrt(length(splitPD.rev));

            %% Importing Data from Pause Frequency
            % PF is a structures variable that contains all Pause Frequencies for the entire experiment. PF.Stats contains averages (row 1)
            % and standard deviation (row 2).

            PF=([]);
            pf=fullfile(directory,'PooledData','PauseFrequency','AllPF.txt');
            PF.Value(:,1)=importdata(pf);

            PF.Stats=zeros(2,4);
            PF.Stats(1,1)=mean(PF.Value);
            PF.Stats(2,1)=std(PF.Value);
            PF.Stats(3,1)=std(PF.Value)/sqrt(length(PF.Value));

            %% Importing Data from Split Pause Frequency
            % PD is a structures variable that contains all Pause Frequencies for the entire experiment. splitPF.Stats contains averages (row 1)
            % and standard deviation (row 2) and SEM (row 3) for anterograde (column 1), retrograde (column 2) and reversal (column 3).

            splitPF=([]);
            pf=fullfile(directory,'PooledData','SplitPauseFrequency','AllantePF.txt');
            splitPF.ante(:,1)=importdata(pf);
            pf=fullfile(directory,'PooledData','SplitPauseFrequency','AllretroPF.txt');
            splitPF.retro(:,1)=importdata(pf);
            pf=fullfile(directory,'PooledData','SplitPauseFrequency','AllrevPF.txt');
            splitPF.rev(:,1)=importdata(pf);

            splitPF.Stats=zeros(2,4);
            splitPF.Stats(1,1)=mean(splitPF.ante);
            splitPF.Stats(1,2)=mean(splitPF.retro);
            splitPF.Stats(1,3)=mean(splitPF.rev);

            splitPF.Stats(2,1)=std(splitPF.ante);
            splitPF.Stats(2,2)=std(splitPF.retro);
            splitPF.Stats(2,3)=std(splitPF.rev);

            splitPF.Stats(3,1)=std(splitPF.ante)/sqrt(length(splitPF.ante));
            splitPF.Stats(3,2)=std(splitPF.retro)/sqrt(length(splitPF.retro));
            splitPF.Stats(3,3)=std(splitPF.rev)/sqrt(length(splitPF.rev));

            %% Importing Data from Pause Frequency per Second
            % PF is a structures variable that contains all Pause Frequencies for the entire experiment. PF.Stats contains averages (row 1)
            % and standard deviation (row 2).

            PFperSec=([]);
            pf=fullfile(directory,'PooledData','PauseFrequencyPerSec','AllPF.txt');
            PFperSec.Value(:,1)=importdata(pf);

            PFperSec.Stats=zeros(2,4);
            PFperSec.Stats(1,1)=mean(PFperSec.Value);
            PFperSec.Stats(2,1)=std(PFperSec.Value);
            PFperSec.Stats(3,1)=std(PFperSec.Value)/sqrt(length(PFperSec.Value));

            %% Importing Data from Split Pause Frequency per Second
            % PD is a structures variable that contains all Pause Frequencies per Second for the entire experiment. splitPFperSec.Stats contains averages (row 1)
            % and standard deviation (row 2) and SEM (row 3) for anterograde (column 1), retrograde (column 2) and reversal (column 3).

            splitPFperSec=([]);
            pf=fullfile(directory,'PooledData','SplitPauseFrequencyPerSec','AllantePFperSec.txt');
            splitPFperSec.ante(:,1)=importdata(pf);
            pf=fullfile(directory,'PooledData','SplitPauseFrequencyPerSec','AllretroPFperSec.txt');
            splitPFperSec.retro(:,1)=importdata(pf);
            pf=fullfile(directory,'PooledData','SplitPauseFrequencyPerSec','AllrevPFperSec.txt');
            splitPFperSec.rev(:,1)=importdata(pf);

            splitPFperSec.Stats=zeros(2,4);
            splitPFperSec.Stats(1,1)=mean(splitPF.ante);
            splitPFperSec.Stats(1,2)=mean(splitPF.retro);
            splitPFperSec.Stats(1,3)=mean(splitPF.rev);

            splitPFperSec.Stats(2,1)=std(splitPF.ante);
            splitPFperSec.Stats(2,2)=std(splitPF.retro);
            splitPFperSec.Stats(2,3)=std(splitPF.rev);

            splitPFperSec.Stats(3,1)=std(splitPFperSec.ante)/sqrt(length(splitPFperSec.ante));
            splitPFperSec.Stats(3,2)=std(splitPFperSec.retro)/sqrt(length(splitPFperSec.retro));
            splitPFperSec.Stats(3,3)=std(splitPFperSec.rev)/sqrt(length(splitPFperSec.rev));

            %% Importing Data from Run Length
            % RL is a structures variable that contains all RL (anterograde and
            % retrograde) for the entire experiment. RL.Stats contains averages (row 1)
            % and standard deviation (row 2) for anterograde (column 1) and retrograde (column 2) run lengths

            RL=([]);
            anteRL=fullfile(directory,'PooledData','RunLength','AllAnteRL.txt');
            RL.ante(:,1)=importdata(anteRL);
            retroRL=fullfile(directory,'PooledData','RunLength','AllRetroRL.txt');
            RL.retro(:,1)=importdata(retroRL);

            RL.Stats=zeros(2,4);
            RL.Stats(1,1)=mean(RL.ante);
            RL.Stats(1,2)=mean(RL.retro);

            RL.Stats(2,1)=std(RL.ante);
            RL.Stats(2,2)=std(RL.retro);

            RL.Stats(3,1)=std(RL.ante)/sqrt(length(RL.ante));
            RL.Stats(3,2)=std(RL.retro)/sqrt(length(RL.retro));

            %% Importing Data from Combined Run Length
            % comRL is a structures variable that contains all RL (anterograde and
            % retrograde) for the entire experiment. RL.Stats contains averages (row 1)
            % and standard deviation (row 2) for anterograde (column 1) and retrograde (column 2) run lengths

            comRL=([]);
            anteComRL=fullfile(directory,'PooledData','CombinedRunLength','AllAnteRLCombined.txt');
            comRL.ante(:,1)=importdata(anteComRL);
            retroComRL=fullfile(directory,'PooledData','CombinedRunLength','AllRetroRLCombined.txt');
            comRL.retro(:,1)=importdata(retroComRL);

            comRL.Stats=zeros(2,4);
            comRL.Stats(1,1)=mean(comRL.ante);
            comRL.Stats(1,2)=mean(comRL.retro);

            comRL.Stats(2,1)=std(comRL.ante);
            comRL.Stats(2,2)=std(comRL.retro);

            comRL.Stats(3,1)=std(comRL.ante)/sqrt(length(comRL.ante));
            comRL.Stats(3,2)=std(comRL.retro)/sqrt(length(comRL.retro));

            %% Importing Data from Segmental Velocities
            % SV is a structures variable that contains all SV (anterograde and
            % retrograde) for the entire experiment. SV.Stats contains averages (row 1)
            % and standard deviation (row 2) for anterograde (column 1) and retrograde (column 2) run lengths

            SV=([]);
            anteSV=fullfile(directory,'PooledData','SegmentalVelocities','AllAnteSV.txt');
            SV.ante(:,1)=importdata(anteSV);
            retroSV=fullfile(directory,'PooledData','SegmentalVelocities','AllRetroSV.txt');
            SV.retro(:,1)=importdata(retroSV);

            SV.Stats = zeros(2,4);
            SV.Stats(1,1) = mean(SV.ante);
            SV.Stats(1,2) = mean(SV.retro);

            SV.Stats(2,1)=std(SV.ante);
            SV.Stats(2,2)=std(SV.retro);

            SV.Stats(3,1)=std(SV.ante)/sqrt(length(SV.ante));
            SV.Stats(3,2)=std(SV.retro)/sqrt(length(SV.retro));

            %% Importing Data from Combined Segmental Velocity
            % comSV is a structures variable that contains all SV (anterograde and
            % retrograde) for the entire experiment. SV.Stats contains averages (row 1)
            % and standard deviation (row 2) for anterograde (column 1) and retrograde (column 2) run lengths

            comSV=([]);
            anteComSV=fullfile(directory,'PooledData','CombinedSegmentalVelocities','AllcomAnteSV.txt');
            comSV.ante(:,1)=importdata(anteComSV);
            retroComSV=fullfile(directory,'PooledData','CombinedSegmentalVelocities','AllcomRetroSV.txt');
            comSV.retro(:,1)=importdata(retroComSV);

            comSV.Stats=zeros(2,4);
            comSV.Stats(1,1)=mean(comSV.ante);
            comSV.Stats(1,2)=mean(comSV.retro);

            comSV.Stats(2,1)=std(comSV.ante);
            comSV.Stats(2,2)=std(comSV.retro);

            comSV.Stats(3,1)=std(comSV.ante)/sqrt(length(comSV.ante));
            comSV.Stats(3,2)=std(comSV.retro)/sqrt(length(comSV.retro));

            %% Importing Data from Switch Frequency
            % SF is a structures variable that contains all Switch Frequencies for the entire experiment. SF.Stats contains averages (row 1)
            % and standard deviation (row 2).

            SF=([]);
            sf=fullfile(directory,'PooledData','SwitchFrequency','AllSF.txt');
            SF.Value(:,1)=importdata(sf);

            SF.Stats=zeros(2,4);
            SF.Stats(1,1)=mean(SF.Value);
            SF.Stats(2,1)=std(SF.Value);
            SF.Stats(3,1)=std(SF.Value)/sqrt(length(SF.Value));

            %% Importing Data from Switch Frequency per Sec
            % SFperSec is a structures variable that contains all Pause Frequencies for the entire experiment. SFperSec.Stats contains averages (row 1)
            % and standard deviation (row 2).

            SFperSec=([]);
            sf=fullfile(directory,'PooledData','SwitchFrequencyPerSec','AllSFperSec.txt');
            SFperSec.Value(:,1)=importdata(sf);

            SFperSec.Stats=zeros(2,4);
            SFperSec.Stats(1,1)=mean(SFperSec.Value);
            SFperSec.Stats(2,1)=std(SFperSec.Value);
            SFperSec.Stats(3,1)=std(SFperSec.Value)/sqrt(length(SFperSec.Value));

            %% Importing Data from Switch Frequency for Reversals
            % revSF is a structures variable that contains all Switch Frequencies for reversal tracks for the entire experiment. revSF.Stats contains averages (row 1)
            % and standard deviation (row 2) and sem (row 3).

            revSF=([]);
            sf=fullfile(directory,'PooledData','SwitchFrequencyReversals','AllrevSF.txt');
            revSF.Value(:,1)=importdata(sf);

            revSF.Stats=zeros(2,4);
            revSF.Stats(1,1)=mean(revSF.Value);
            revSF.Stats(2,1)=std(revSF.Value);
            revSF.Stats(3,1)=std(revSF.Value)/sqrt(length(revSF.Value));

            %% Importing Data from Switch Frequency for Reversals per Second
            % revSFperSec is a structures variable that contains all Switch Frequencies for reversal tracks for the entire experiment. revSF.Stats contains averages (row 1)
            % and standard deviation (row 2).

            revSFperSec=([]);
            sf=fullfile(directory,'PooledData','SwitchFrequencyReversalsPerSec','AllrevSFperSec.txt');
            revSFperSec.Value(:,1)=importdata(sf);

            revSFperSec.Stats=zeros(2,4);
            revSFperSec.Stats(1,1)=mean(revSFperSec.Value);
            revSFperSec.Stats(2,1)=std(revSFperSec.Value);
            revSFperSec.Stats(3,1)=std(revSFperSec.Value)/sqrt(length(revSFperSec.Value));

    %% Creating Figures for bar graphs with STD

            STD=zeros(4,1);
            X=([1 2 3 4]);
            % Figure 1 is figure for all averages
            figure1=figure;
            set(gcf,'Position',[600 50 1200 900]);

            % plots for Cargo Population
            subplot('Position',[0.1,0.75,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(CP.StatsPCT(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            ylim ([0 100]);
            title('Cargo Population');
            ylabel('Cargo Population (%)');
            set(gca,'XTickLabel',{'Ante','Retro','Reversing','Stationary'});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,CP.StatsPCT(1,:),STD,CP.StatsPCT(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Net Cargo Population
            subplot('Position',[0.32,0.75,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(NCP.StatsPCT(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            ylim ([0 100]);
            title('Net Cargo Population');
            ylabel('Net Cargo Population (%)');
            set(gca,'XTickLabel',{'Ante','Retro','Stationary',''});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,NCP.StatsPCT(1,:),STD,NCP.StatsPCT(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Density
            subplot('Position',[0.56,0.75,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(Density.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 0.5]);
            title('Density');
            ylabel('Density Analysis (Tracks / um axon)');
            set(gca,'XTickLabel',{'Ante','Retro','Reversal','Stationary'});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,Density.Stats(1,:),STD,Density.Stats(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

    %         % plots for Flux
    %         subplot('Position',[0.56,0.75,0.16,0.2]); % need to re-adjust the
    %         %position of the graph
    %         set(gca,'fontsize',8);
    %         bar(Flux.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
    %         %ylim ([0 0.5]);
    %         title('Flux');
    %         ylabel('Flux Analysis (Tracks / um axon)');
    %         set(gca,'XTickLabel',{'Ante','Retro','Reversal'});
    %         rotateXLabels(gca,45);
    %         set(gca, 'Ticklength', [0 0])
    %         hold on
    %         h=errorbar(X,Flux.Stats(1,:),STD,Flux.Stats(2,:),'k.');
    %         set(h(1),'MarkerEdgeColor','none');

            % plots for Segmental Velocities
            subplot('Position',[0.8,0.75,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(SV.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Segmental Velocities');
            ylabel('Segmental Velocities (um/sec)');
            set(gca,'XTickLabel',{'Ante','Retro','',''});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,SV.Stats(1,:),STD,SV.Stats(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Net Velocities
            subplot('Position',[0.1,0.4,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(NV.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Net Velocities');
            ylabel('Net Velocities (um/sec)');
            set(gca,'XTickLabel',{'Ante','Retro','',''});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,NV.Stats(1,:),STD,NV.Stats(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Run Length
            subplot('Position',[0.32,0.4,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(RL.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 25]);
            title('Segmental Run Length');
            ylabel('Segmental Run Length (um)');
            set(gca,'XTickLabel',{'Ante','Retro','',''});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,RL.Stats(1,:),STD,RL.Stats(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Pause Duration
            subplot('Position',[0.56,0.4,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(PD.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 15]);
            title('Pause Duration');
            ylabel('Pause Duration (sec)');
            set(gca,'XTickLabel','');
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,PD.Stats(1,:),STD,PD.Stats(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Pause Frequency
            subplot('Position',[0.8,0.4,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(PF.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Pause Frequency');
            ylabel('Pause Frequency (Pauses/Track)');
            set(gca,'XTickLabel','');
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,PF.Stats(1,:),STD,PF.Stats(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

             % plots for Pause Frequency per Second
            subplot('Position',[0.1,0.05,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(PFperSec.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Pause Frequency per Second');
            ylabel('Pause Frequency (Pauses/sec)');
            set(gca,'XTickLabel','');
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,PFperSec.Stats(1,:),STD,PFperSec.Stats(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Switch Frequency
            subplot('Position',[0.32,0.05,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(SF.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Switch Frequency');
            ylabel('Switch Frequency (Switches/Track)');
            set(gca,'XTickLabel','');
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,SF.Stats(1,:),STD,SF.Stats(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Switch Frequency per Second
            subplot('Position',[0.56,0.05,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(SFperSec.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Switch Frequency Per Second');
            ylabel('Switch Frequency (Switches/sec)');
            set(gca,'XTickLabel','');
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,SFperSec.Stats(1,:),STD,SFperSec.Stats(2,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

    %% Creating Figures for bar graphs with SEM

            STD=zeros(4,1);
            X=([1 2 3 4]);
            % Figure 1 is figure for all averages
            figure2=figure;
            set(gcf,'Position',[600 50 1200 900]);

            % plots for Cargo Population
            subplot('Position',[0.1,0.75,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(CP.StatsPCT(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            ylim ([0 100]);
            title('Cargo Population');
            ylabel('Cargo Population (%)');
            set(gca,'XTickLabel',{'Ante','Retro','Reversing','Stationary'});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,CP.StatsPCT(1,:),STD,CP.StatsPCT(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Net Cargo Population
            subplot('Position',[0.32,0.75,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(NCP.StatsPCT(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            ylim ([0 100]);
            title('Net Cargo Population');
            ylabel('Net Cargo Population (%)');
            set(gca,'XTickLabel',{'Ante','Retro','Stationary',''});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,NCP.StatsPCT(1,:),STD,NCP.StatsPCT(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Density
            subplot('Position',[0.56,0.75,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(Density.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 0.5]);
            title('Density');
            ylabel('Density Analysis (Tracks / um axon)');
            set(gca,'XTickLabel',{'Ante','Retro','Reversal','Stationary'});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,Density.Stats(1,:),STD,Density.Stats(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

    %         % plots for Flux
    %         subplot('Position',[0.56,0.75,0.16,0.2]);
    %         set(gca,'fontsize',8);
    %         bar(Flux.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
    %         %ylim ([0 0.5]);
    %         title('Flux');
    %         ylabel('Flux Analysis (Tracks / um axon)');
    %         set(gca,'XTickLabel',{'Ante','Retro','Reversal','Stationary'});
    %         rotateXLabels(gca,45);
    %         set(gca, 'Ticklength', [0 0])
    %         hold on
    %         h=errorbar(X,Flux.Stats(1,:),STD,Flux.Stats(3,:),'k.');
    %         set(h(1),'MarkerEdgeColor','none');

            % plots for Segmental Velocities
            subplot('Position',[0.8,0.75,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(SV.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Segmental Velocities');
            ylabel('Segmental Velocities (um/sec)');
            set(gca,'XTickLabel',{'Ante','Retro','',''});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,SV.Stats(1,:),STD,SV.Stats(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Net Velocities
            subplot('Position',[0.1,0.4,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(NV.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Net Velocities');
            ylabel('Net Velocities (um/sec)');
            set(gca,'XTickLabel',{'Ante','Retro','',''});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,NV.Stats(1,:),STD,NV.Stats(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Run Length
            subplot('Position',[0.32,0.4,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(RL.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 25]);
            title('Segmental Run Length');
            ylabel('Segmental Run Length (um)');
            set(gca,'XTickLabel',{'Ante','Retro','',''});
            rotateXLabels(gca,45);
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,RL.Stats(1,:),STD,RL.Stats(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Pause Duration
            subplot('Position',[0.56,0.4,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(PD.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 15]);
            title('Pause Duration');
            ylabel('Pause Duration (sec)');
            set(gca,'XTickLabel','');
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,PD.Stats(1,:),STD,PD.Stats(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Pause Frequency
            subplot('Position',[0.8,0.4,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(PF.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Pause Frequency');
            ylabel('Pause Frequency (Pauses/Track)');
            set(gca,'XTickLabel','');
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,PF.Stats(1,:),STD,PF.Stats(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

             % plots for Pause Frequency per Second
            subplot('Position',[0.1,0.05,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(PFperSec.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Pause Frequency per Second');
            ylabel('Pause Frequency (Pauses/sec)');
            set(gca,'XTickLabel','');
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,PFperSec.Stats(1,:),STD,PFperSec.Stats(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Switch Frequency
            subplot('Position',[0.32,0.05,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(SF.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Switch Frequency');
            ylabel('Switch Frequency (Switches/Track)');
            set(gca,'XTickLabel','');
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,SF.Stats(1,:),STD,SF.Stats(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

            % plots for Switch Frequency per Second
            subplot('Position',[0.56,0.05,0.16,0.2]);
            set(gca,'fontsize',8);
            bar(SFperSec.Stats(1,:),'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.9);
            %ylim ([0 10]);
            title('Switch Frequency Per Second');
            ylabel('Switch Frequency (Switches/sec)');
            set(gca,'XTickLabel','');
            set(gca, 'Ticklength', [0 0])
            hold on
            h=errorbar(X,SFperSec.Stats(1,:),STD,SFperSec.Stats(3,:),'k.');
            set(h(1),'MarkerEdgeColor','none');

    %% Creating Figures for histograms
            % figure 3 is for all histograms

            figure3=figure;
            set(gcf,'Position',[600 100 500 800]);

            % anterograde segmental velocities
            histSV=([]);
            histSV.dataAnte=histc(SV.ante(:,1),0:1:20);
            FRCT=rdivide(histSV.dataAnte,sum(histSV.dataAnte));
            histSV.pctAnte=mtimes(100,FRCT);
            subplot('Position',[0.1,0.75,0.28,0.175]);
            set(gca,'fontsize',8);
            bar(histSV.pctAnte,'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.6);
            title('Anterograde Segmental Velocities');
            xlim([0 20]);
            xlabel('Segmental Velocities (um/sec)');
            ylabel('Frequency (%)');

            % retrograde segmental velocities
            histSV.dataRetro=histc(SV.retro(:,1),0:1:20);
            FRCT=rdivide(histSV.dataRetro,sum(histSV.dataRetro));
            histSV.pctRetro=mtimes(100,FRCT);
            subplot('Position',[0.6,0.75,0.28,0.175]);
            set(gca,'fontsize',8);
            bar(histSV.pctRetro,'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.6);
            title('Retrograde Segmental Velocities');
            xlim([0 20]);
            xlabel('Segmental Velocities (um/sec)');
            ylabel('Frequency (%)');

            % anterograde net velocities
            histNV=([]);
            histNV.dataAnte=histc(NV.ante(:,1),0:1:20);
            FRCT=rdivide(histNV.dataAnte,sum(histNV.dataAnte));
            histNV.pctAnte=mtimes(100,FRCT);
            subplot('Position',[0.1,0.45,0.28,0.175]);
            set(gca,'fontsize',8);
            bar(histNV.pctAnte,'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.6);
            title('Anterograde Net Velocities');
            xlim([0 20]);
            xlabel('Net Velocities (um/sec)');
            ylabel('Frequency (%)');

            % retrograde net velocities
            histNV.dataRetro=histc(NV.retro(:,1),0:1:20);
            FRCT=rdivide(histNV.dataRetro,sum(histNV.dataRetro));
            histNV.pctRetro=mtimes(100,FRCT);
            subplot('Position',[0.6,0.45,0.28,0.175]);
            set(gca,'fontsize',8);
            bar(histNV.pctRetro,'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.6);
            title('Retrograde Net Velocities');
            xlim([0 20]);
            xlabel('Segmental Velocities (um/sec)');
            ylabel('Frequency (%)');

            % anterograde run length
            histRL=([]);
            histRL.dataAnte=histc(RL.ante(:,1),0:1:30);
            FRCT=rdivide(histRL.dataAnte,sum(histRL.dataAnte));
            histRL.pctAnte=mtimes(100,FRCT);
            subplot('Position',[0.1,0.15,0.28,0.175]);
            set(gca,'fontsize',8);
            bar(histRL.pctAnte,'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.6);
            title('Anterograde Segmental Run Length');
            xlim([0 30]);
            xlabel('Run Length (um)');
            ylabel('Frequency (%)');

            % retrograde run length
            histRL.dataRetro=histc(RL.retro(:,1),0:1:30);
            FRCT=rdivide(histRL.dataRetro,sum(histRL.dataRetro));
            histRL.pctRetro=mtimes(100,FRCT);
            subplot('Position',[0.6,0.15,0.28,0.175]);
            set(gca,'fontsize',8);
            bar(histRL.pctRetro,'Facecolor',[0 0.51 0.25],'EdgeColor',[0 0.51 0.25],'BarWidth',0.6);
            title('Retrograde Segmental Run Length');
            xlim([0 30]);
            xlabel('Run Length (um)');
            ylabel('Frequency (%)');

    %% Saving
            [pathstr,name] = fileparts(list2(j,1).name);

            %saving figures
            mkdir(fullfile(directory,'Figures'));
            fileFigure=fullfile(directory,'Figures');
            fileFigure=strcat(fileFigure,filesep,list(i,1).name,'_',name,'_AveragesSTD.fig');
            saveas(figure1,fileFigure);
            fileFigure=fullfile(directory,'Figures');
            fileFigure=strcat(fileFigure,filesep,list(i,1).name,'_',name,'_AveragesSTD.eps');
            saveas(figure1,fileFigure);

            fileFigure=fullfile(directory,'Figures');
            fileFigure=strcat(fileFigure,filesep,list(i,1).name,'_',name,'_AveragesSEM.fig');
            saveas(figure2,fileFigure);
            fileFigure=fullfile(directory,'Figures');
            fileFigure=strcat(fileFigure,filesep,list(i,1).name,'_',name,'_AveragesSEM.eps');
            saveas(figure2,fileFigure);

            fileFigure=fullfile(directory,'Figures');
            fileFigure=strcat(fileFigure,filesep,list(i,1).name,'_',name,'_Histograms.fig');
            saveas(figure3,fileFigure);
            fileFigure=fullfile(directory,'Figures');
            fileFigure=strcat(fileFigure,filesep,list(i,1).name,'_',name,'_Histograms.eps');
            saveas(figure3,fileFigure);
            close all;

            % Clearing variables
            varlist = {'FRCT','h','NetanteNum','NetretroNum','NetstationaryNum','NetantePCT','NetretroPCT','NetstationaryPCT',...
                'anteNum','retroNum','reversalNum','stationaryNum','antePCT','retroPCT','reversalPCT','stationaryPCT',...
                'anteFlux','retroFlux','reversalFlux',...
                'anteNV','retroNV','anteSV','retroSV','anteRL','retroRL','sf','pd','pf','fileFigure'};
            clear (varlist{:}); clear varlist;


            %saving variables (.mat file)
            filename=strcat(list(i,1).name,'_',name,'.m');
            mkdir(fullfile(directory,'MatlabResults'));
            file=fullfile(directory,'MatlabResults',filename);
            varlist = {'pathstr','name','filename',};
            clear (varlist{:}); clear varlist;
            save(file,'CP','Flux','Density','NCP','NV','PD','PF','PFperSec','splitPD','splitPF','splitPFperSec','SF','SFperSec','revSF','revSFperSec','SV','RL','comRL',...
                'histNV','histRL','histSV','comSV');
        end
    end
end
 
    
     


