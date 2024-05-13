//This macro pools all data of i) each movie and ii) each experiment folder
//1. Pooled Net Velocities
//2A. CargoPopulation numeric
//2B. CargoPopulation percentage
//3A. NetCardoPopulation numeric
//3B. NetCardoPopulation percentage
//4a. Segmental Velocities
//4b. Combined Segmental Velocities
//5a. Pause Duration
//5b. Split Pause Duration
//6a. Pooled Pause Frequency
//6b. Pooled Split Pause Frequency
//6c. Pooled Pause Frequency/sec
//6d. Pooled Split Pause Frequency/sec
//7a. Pooled Flux
//7b. Pooled Densities
//8a. Pooled Segmental Run Length
//8b. Pooled Combined Segmental Run Length
//9a. Pooled Switch Frequency
//9b. Pooled Switch Frequency/sec
//9c. Pooled Switch Frequency for reversals
//9c. Pooled Switch Frequency for reversals/sec
//10. Make Montage of Kymographs and overlay CP,NCP,Segments


run("Close All");
roiManager("Reset");

input=getDirectory("Choose Experiment Folder");
File.makeDirectory(input+"PooledData");
ExperimentList=getFileList(input);
inputPooled=input+"PooledData/";
File.makeDirectory(inputPooled+"NetVelocities");
File.makeDirectory(inputPooled+"CargoPopulation");
File.makeDirectory(inputPooled+"NetCargoPopulation");
File.makeDirectory(inputPooled+"SegmentalVelocities");
File.makeDirectory(inputPooled+"CombinedSegmentalVelocities");
File.makeDirectory(inputPooled+"PauseDuration");
File.makeDirectory(inputPooled+"SplitPauseDuration");
File.makeDirectory(inputPooled+"PauseFrequency");
File.makeDirectory(inputPooled+"PauseFrequencyPerSec");
File.makeDirectory(inputPooled+"SplitPauseFrequency");
File.makeDirectory(inputPooled+"SplitPauseFrequencyPerSec");
File.makeDirectory(inputPooled+"Flux");
File.makeDirectory(inputPooled+"Density");
File.makeDirectory(inputPooled+"RunLength");
File.makeDirectory(inputPooled+"CombinedRunLength");
File.makeDirectory(inputPooled+"SwitchFrequency");
File.makeDirectory(inputPooled+"SwitchFrequencyPerSec");
File.makeDirectory(inputPooled+"SwitchFrequencyReversals");
File.makeDirectory(inputPooled+"SwitchFrequencyReversalsPerSec");
File.makeDirectory(inputPooled+"Kymographs");
File.makeDirectory(inputPooled+"PercentTimeMotion");

outputNV=inputPooled+"NetVelocities/";
outputCP=inputPooled+"CargoPopulation/";
outputNCP=inputPooled+"NetCargoPopulation/";
outputSegVel=inputPooled+"SegmentalVelocities/";
outputcomSV=inputPooled+"CombinedSegmentalVelocities/";
outputPD=inputPooled+"PauseDuration/";
outputsplitPD=inputPooled+"SplitPauseDuration/";
outputPF=inputPooled+"PauseFrequency/";
outputPFperSec=inputPooled+"PauseFrequencyPerSec/";
outputsplitPF=inputPooled+"SplitPauseFrequency/";
outputsplitPFperSec=inputPooled+"SplitPauseFrequencyPerSec/";
outputFlux=inputPooled+"Flux/";
outputDensity=inputPooled+"Density/";
outputRL=inputPooled+"RunLength/";
outputcombinedRL=inputPooled+"CombinedRunLength/";
outputSF=inputPooled+"SwitchFrequency/";
outputSFperSec=inputPooled+"SwitchFrequencyPerSec/";
outputrevSF=inputPooled+"SwitchFrequencyReversals/";
outputrevSFperSec=inputPooled+"SwitchFrequencyReversalsPerSec/";
outputKymo=inputPooled+"Kymographs/";
outputPM=inputPooled+"PercentTimeMotion/";

//--------------------------------------------- 1. Pooled Net Velocities ---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){

	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputNetVel=input+ExperimentList[i]+"Net_Velocities/";
		NetVelList=getFileList(inputNetVel);
		anteNV=newArray();
		retroNV=newArray();
		for(k=0;k<NetVelList.length;k++){
		string = File.openAsString(inputNetVel+NetVelList[k]); 
		lines = split(string, "\n");   
		NV=newArray(); 
		
		for (j=0; j<lines.length; j++){	
			NV = Array.concat(NV,parseFloat(lines[j]));
		}
		if (NV[0]>0)
		anteNV=Array.concat(anteNV,abs(NV[0]));
		if (NV[0]<0)
		retroNV=Array.concat(retroNV,abs(NV[0]));
		}
		AnteNV=File.open(outputNV+"anteNV_"+replace(ExperimentList[i],"/",".txt"));
		for (j=0;j<anteNV.length;j++){
			print(AnteNV,anteNV[j]);
		}
		File.close(AnteNV);
		RetroNV=File.open(outputNV+"retroNV_"+replace(ExperimentList[i],"/",".txt"));
		for (j=0;j<retroNV.length;j++){
			print(RetroNV,retroNV[j]);
		}
		File.close(RetroNV);
	}
}

AllNVList=getFileList(outputNV);
AllAnteNV=newArray();
AllRetroNV=newArray();
for (i=0;i<AllNVList.length;i++){
	if(startsWith(AllNVList[i],"ante")){
		string = File.openAsString(outputNV+AllNVList[i]); 
		lines = split(string, "\n");   
		NV=newArray(); 
		for (j=0; j<lines.length; j++){
		NV = Array.concat(NV,parseFloat(lines[j]));	
		}
		for (j=0;j<NV.length;j++){
		AllAnteNV=Array.concat(AllAnteNV,abs(NV[j]));
		}
		allAnteNV=File.open(outputNV+"AllAnteNV"+".txt");
		for (j=0;j<AllAnteNV.length;j++){
			print(allAnteNV,AllAnteNV[j]);
		}
		File.close(allAnteNV);
	}
	if(startsWith(AllNVList[i],"retro")){
		string = File.openAsString(outputNV+AllNVList[i]); 
		lines = split(string, "\n");   
		NV=newArray(); 
		for (j=0; j<lines.length; j++){
		NV = Array.concat(NV,parseFloat(lines[j]));	
		}
		for (j=0;j<NV.length;j++){
		AllRetroNV=Array.concat(AllRetroNV,abs(NV[j]));
		}
		allRetroNV=File.open(outputNV+"AllRetroNV"+".txt");
		for (j=0;j<AllRetroNV.length;j++){
			print(allRetroNV,AllRetroNV[j]);
		}
		File.close(allRetroNV);
	}
		
}
//--------------------------------------------- 2A. Pooled Cargo Population Numeric---------------------------------------------------
anteNum=newArray();
retroNum=newArray();
reversalNum=newArray();
stationaryNum=newArray();
for (i=0;i<ExperimentList.length;i++){	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputCP=input+ExperimentList[i]+"CargoPopulation/";
		CPList=getFileList(inputCP);
		string = File.openAsString(inputCP+CPList[0]); 
		lines = split(string, "\n");   
		CPNum=newArray(); 
		for (j=0; j<lines.length; j++){	
			CPNum = Array.concat(CPNum,parseFloat(lines[j]));
		}
		
		anteNum=Array.concat(anteNum,abs(CPNum[0]));
		retroNum=Array.concat(retroNum,abs(CPNum[1]));
		reversalNum=Array.concat(reversalNum,abs(CPNum[2]));
		stationaryNum=Array.concat(stationaryNum,abs(CPNum[3]));
	}
		Ante=File.open(outputCP+"anteNum"+".txt");
		for (j=0;j<anteNum.length;j++){
			print(Ante,anteNum[j]);
		}
		File.close(Ante);
		Retro=File.open(outputCP+"retroNum"+".txt");
		for (j=0;j<retroNum.length;j++){
			print(Retro,retroNum[j]);
		}
		File.close(Retro);
		Reversal=File.open(outputCP+"reversalNum"+".txt");
		for (j=0;j<reversalNum.length;j++){
			print(Reversal,reversalNum[j]);
		}
		File.close(Reversal);
		Stationary=File.open(outputCP+"stationaryNum"+".txt");
		for (j=0;j<stationaryNum.length;j++){
			print(Stationary,stationaryNum[j]);
		}
		File.close(Stationary);
}
//--------------------------------------------- 2B. Pooled Cargo Population Percentage---------------------------------------------------
antePCT=newArray();
retroPCT=newArray();
reversalPCT=newArray();
stationaryPCT=newArray();
for (i=0;i<ExperimentList.length;i++){	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputCP=input+ExperimentList[i]+"CargoPopulation/";
		CPList=getFileList(inputCP);
		string = File.openAsString(inputCP+CPList[1]); 
		lines = split(string, "\n");   
		CPPCT=newArray(); 
		for (j=0; j<lines.length; j++){	
			CPPCT = Array.concat(CPPCT,parseFloat(lines[j]));
		}
		
		antePCT=Array.concat(antePCT,abs(CPPCT[0]));
		retroPCT=Array.concat(retroPCT,abs(CPPCT[1]));
		reversalPCT=Array.concat(reversalPCT,abs(CPPCT[2]));
		stationaryPCT=Array.concat(stationaryPCT,abs(CPPCT[3]));
	}
		Ante=File.open(outputCP+"antePCT"+".txt");
		for (j=0;j<antePCT.length;j++){
			print(Ante,antePCT[j]);
		}
		File.close(Ante);
		Retro=File.open(outputCP+"retroPCT"+".txt");
		for (j=0;j<retroPCT.length;j++){
			print(Retro,retroPCT[j]);
		}
		File.close(Retro);
		Reversal=File.open(outputCP+"reversalPCT"+".txt");
		for (j=0;j<reversalPCT.length;j++){
			print(Reversal,reversalPCT[j]);
		}
		File.close(Reversal);
		Stationary=File.open(outputCP+"stationaryPCT"+".txt");
		for (j=0;j<stationaryPCT.length;j++){
			print(Stationary,stationaryPCT[j]);
		}
		File.close(Stationary);
}
//--------------------------------------------- 3A. Pooled Net Cargo Population Numeric---------------------------------------------------
anteNum=newArray();
retroNum=newArray();
stationaryNum=newArray();
for (i=0;i<ExperimentList.length;i++){	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputCP=input+ExperimentList[i]+"CargoPopulation/";
		CPList=getFileList(inputCP);
		string = File.openAsString(inputCP+CPList[2]); 
		lines = split(string, "\n");   
		CPNum=newArray(); 
		for (j=0; j<lines.length; j++){	
			CPNum = Array.concat(CPNum,parseFloat(lines[j]));
		}
		
		anteNum=Array.concat(anteNum,abs(CPNum[0]));
		retroNum=Array.concat(retroNum,abs(CPNum[1]));
		stationaryNum=Array.concat(stationaryNum,abs(CPNum[2]));
	}
		Ante=File.open(outputNCP+"NetanteNum"+".txt");
		for (j=0;j<anteNum.length;j++){
			print(Ante,anteNum[j]);
		}
		File.close(Ante);
		Retro=File.open(outputNCP+"NetretroNum"+".txt");
		for (j=0;j<retroNum.length;j++){
			print(Retro,retroNum[j]);
		}
		File.close(Retro);
		Stationary=File.open(outputNCP+"NetstationaryNum"+".txt");
		for (j=0;j<stationaryNum.length;j++){
			print(Stationary,stationaryNum[j]);
		}
		File.close(Stationary);
}
//--------------------------------------------- 3B. Net Pooled Cargo Population Percentage---------------------------------------------------
antePCT=newArray();
retroPCT=newArray();
stationaryPCT=newArray();
for (i=0;i<ExperimentList.length;i++){	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputCP=input+ExperimentList[i]+"CargoPopulation/";
		CPList=getFileList(inputCP);
		string = File.openAsString(inputCP+CPList[3]); 
		lines = split(string, "\n");   
		CPPCT=newArray(); 
		for (j=0; j<lines.length; j++){	
			CPPCT = Array.concat(CPPCT,parseFloat(lines[j]));
		}
		
		antePCT=Array.concat(antePCT,abs(CPPCT[0]));
		retroPCT=Array.concat(retroPCT,abs(CPPCT[1]));
		stationaryPCT=Array.concat(stationaryPCT,abs(CPPCT[2]));
	}
		Ante=File.open(outputNCP+"NetantePCT"+".txt");
		for (j=0;j<antePCT.length;j++){
			print(Ante,antePCT[j]);
		}
		File.close(Ante);
		Retro=File.open(outputNCP+"NetretroPCT"+".txt");
		for (j=0;j<retroPCT.length;j++){
			print(Retro,retroPCT[j]);
		}
		File.close(Retro);
		File.close(Reversal);
		Stationary=File.open(outputNCP+"NetstationaryPCT"+".txt");
		for (j=0;j<stationaryPCT.length;j++){
			print(Stationary,stationaryPCT[j]);
		}
		File.close(Stationary);
}

//--------------------------------------------- 4a. Pooled Segmental Velocities---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){

	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputSegVel=input+ExperimentList[i]+"Segmental_Velocities/";
		SegVelList=getFileList(inputSegVel);
		anteSV=newArray();
		retroSV=newArray();
			for(k=0;k<SegVelList.length;k++){
			string = File.openAsString(inputSegVel+SegVelList[k]); 
			lines = split(string, "\n");   
			SV=newArray(); 
			
			for (j=0; j<lines.length; j++){	
				SV = Array.concat(SV,parseFloat(lines[j]));
			}
	
			for (j=0;j<SV.length;j++){
			if (SV[j]>0){
			anteSV=Array.concat(anteSV,abs(SV[j]));
			}
			if (SV[j]<0){
			retroSV=Array.concat(retroSV,abs(SV[j]));
			}
			}
			
			AnteSV=File.open(outputSegVel+"anteSV_"+replace(ExperimentList[i],"/",".txt"));
			for (j=0;j<anteSV.length;j++){
				print(AnteSV,anteSV[j]);
			}
			File.close(AnteSV);
			RetroSV=File.open(outputSegVel+"retroSV_"+replace(ExperimentList[i],"/",".txt"));
			for (j=0;j<retroSV.length;j++){
				print(RetroSV,retroSV[j]);
			}
			File.close(RetroSV);
		}
}
}

AllSVList=getFileList(outputSegVel);
AllAnteSV=newArray();
AllRetroSV=newArray();
for (i=0;i<AllSVList.length;i++){
	if(startsWith(AllSVList[i],"ante")){
		string = File.openAsString(outputSegVel+AllSVList[i]); 
		lines = split(string, "\n");   
		SV=newArray(); 
		for (j=0; j<lines.length; j++){
		SV = Array.concat(SV,parseFloat(lines[j]));	
		}
		for (j=0;j<SV.length;j++){
		AllAnteSV=Array.concat(AllAnteSV,abs(SV[j]));
		}
		allAnteSV=File.open(outputSegVel+"AllAnteSV"+".txt");
		for (j=0;j<AllAnteSV.length;j++){
			print(allAnteSV,AllAnteSV[j]);
		}
		File.close(allAnteSV);
	}
	if(startsWith(AllSVList[i],"retro")){
		string = File.openAsString(outputSegVel+AllSVList[i]); 
		lines = split(string, "\n");   
		SV=newArray(); 
		for (j=0; j<lines.length; j++){
		SV = Array.concat(SV,parseFloat(lines[j]));	
		}
		for (j=0;j<SV.length;j++){
		AllRetroSV=Array.concat(AllRetroSV,abs(SV[j]));
		}
		allRetroSV=File.open(outputSegVel+"AllRetroSV"+".txt");
		for (j=0;j<AllRetroSV.length;j++){
			print(allRetroSV,AllRetroSV[j]);
		}
		File.close(allRetroSV);
	}
		
}

//--------------------------------------------- 4b. Pooled Combined Segmental Velocities---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){

	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputSegVel=input+ExperimentList[i]+"combinedSV/";
		SegVelList=getFileList(inputSegVel);
		anteSV=newArray();
		retroSV=newArray();
			for(k=0;k<SegVelList.length;k++){
			string = File.openAsString(inputSegVel+SegVelList[k]); 
			lines = split(string, "\n");   
			SV=newArray(); 
			
			for (j=0; j<lines.length; j++){	
				SV = Array.concat(SV,parseFloat(lines[j]));
			}
	
			for (j=0;j<SV.length;j++){
			if (SV[j]>0){
			anteSV=Array.concat(anteSV,abs(SV[j]));
			}
			if (SV[j]<0){
			retroSV=Array.concat(retroSV,abs(SV[j]));
			}
			}
			
			AnteSV=File.open(outputcomSV+"comAnteSV_"+replace(ExperimentList[i],"/",".txt"));
			for (j=0;j<anteSV.length;j++){
				print(AnteSV,anteSV[j]);
			}
			File.close(AnteSV);
			RetroSV=File.open(outputcomSV+"comRetroSV_"+replace(ExperimentList[i],"/",".txt"));
			for (j=0;j<retroSV.length;j++){
				print(RetroSV,retroSV[j]);
			}
			File.close(RetroSV);
		}
}
}

AllSVList=getFileList(outputcomSV);
AllAnteSV=newArray();
AllRetroSV=newArray();
for (i=0;i<AllSVList.length;i++){
	if(startsWith(AllSVList[i],"comAnte")){
		string = File.openAsString(outputcomSV+AllSVList[i]); 
		lines = split(string, "\n");   
		SV=newArray(); 
		for (j=0; j<lines.length; j++){
		SV = Array.concat(SV,parseFloat(lines[j]));	
		}
		for (j=0;j<SV.length;j++){
		AllAnteSV=Array.concat(AllAnteSV,abs(SV[j]));
		}
		allAnteSV=File.open(outputcomSV+"AllcomAnteSV"+".txt");
		for (j=0;j<AllAnteSV.length;j++){
			print(allAnteSV,AllAnteSV[j]);
		}
		File.close(allAnteSV);
	}
	if(startsWith(AllSVList[i],"comRetro")){
		string = File.openAsString(outputcomSV+AllSVList[i]); 
		lines = split(string, "\n");   
		SV=newArray(); 
		for (j=0; j<lines.length; j++){
		SV = Array.concat(SV,parseFloat(lines[j]));	
		}
		for (j=0;j<SV.length;j++){
		AllRetroSV=Array.concat(AllRetroSV,abs(SV[j]));
		}
		allRetroSV=File.open(outputcomSV+"AllcomRetroSV"+".txt");
		for (j=0;j<AllRetroSV.length;j++){
			print(allRetroSV,AllRetroSV[j]);
		}
		File.close(allRetroSV);
	}
		
}


//--------------------------------------------- 5. Pooled Pause Duration ---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
	movPD=File.open(outputPD+"PD_"+replace(ExperimentList[i],"/",".txt"));
		if (startsWith(ExperimentList[i],"PooledData")){	
		} else {
			inputPD=input+ExperimentList[i]+"PD/";
			listPD=getFileList(inputPD);
			
			for(k=0;k<listPD.length;k++){
				string = File.openAsString(inputPD+listPD[k]); 
				lines = split(string, "\n");   
				PD=newArray(); 
				for (j=0; j<lines.length; j++){	
					PD = Array.concat(PD,parseFloat(lines[j]));
				}
				AllPD=newArray();
				for(j=0;j<PD.length;j++){
					AllPD=Array.concat(AllPD,PD[j]);
				}
			
			for (j=0;j<AllPD.length;j++){
				print(movPD,AllPD[j]);
				}
			
			
			}	
		}
		File.close(movPD);
	}
}

PDlist=getFileList(outputPD);
allPD=File.open(outputPD+"AllPD"+".txt");
for (i=0;i<PDlist.length;i++){
	
	if (startsWith(PDlist[i],"All")){
	} else {
		
		string = File.openAsString(outputPD+PDlist[i]); 
			lines = split(string, "\n");   
			PD=newArray(); 
			for (j=0; j<lines.length; j++){	
				PD = Array.concat(PD,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPD,PD[j]);
			}
	}

}
File.close(allPD);

//--------------------------------------------- 5A. Pooled Split Pause Duration ---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		movPD=File.open(outputsplitPD+"antePD_"+replace(ExperimentList[i],"/",".txt"));
			if (startsWith(ExperimentList[i],"PooledData")){	
			} else {
				inputPD=input+ExperimentList[i]+"splitPD/";
				listPD=getFileList(inputPD);
				
				for(k=0;k<listPD.length;k++){
					if (startsWith(listPD[k],"ante")){
						string = File.openAsString(inputPD+listPD[k]); 
						lines = split(string, "\n");   
						PD=newArray(); 
						for (j=0; j<lines.length; j++){	
							PD = Array.concat(PD,parseFloat(lines[j]));
						}
						for (j=0;j<PD.length;j++){
							print(movPD,PD[j]);
						}
					}
				}	
			}
		File.close(movPD);
		movPD=File.open(outputsplitPD+"retroPD_"+replace(ExperimentList[i],"/",".txt"));
			if (startsWith(ExperimentList[i],"PooledData")){	
			} else {
				inputPD=input+ExperimentList[i]+"splitPD/";
				listPD=getFileList(inputPD);
				
				for(k=0;k<listPD.length;k++){
					if (startsWith(listPD[k],"retro")){
						string = File.openAsString(inputPD+listPD[k]); 
						lines = split(string, "\n");   
						PD=newArray(); 
						for (j=0; j<lines.length; j++){	
							PD = Array.concat(PD,parseFloat(lines[j]));
						}
						for (j=0;j<PD.length;j++){
							print(movPD,PD[j]);
						}
					}
				}	
			}
		File.close(movPD);
		movPD=File.open(outputsplitPD+"revPD_"+replace(ExperimentList[i],"/",".txt"));
			if (startsWith(ExperimentList[i],"PooledData")){	
			} else {
				inputPD=input+ExperimentList[i]+"splitPD/";
				listPD=getFileList(inputPD);
				
				for(k=0;k<listPD.length;k++){
					if (startsWith(listPD[k],"rev")){
						string = File.openAsString(inputPD+listPD[k]); 
						lines = split(string, "\n");   
						PD=newArray(); 
						for (j=0; j<lines.length; j++){	
							PD = Array.concat(PD,parseFloat(lines[j]));
						}
						for (j=0;j<PD.length;j++){
							print(movPD,PD[j]);
						}
					}
				}	
			}
		File.close(movPD);
	}
}

PDlist=getFileList(outputsplitPD);
allPD=File.open(outputsplitPD+"AllantePD"+".txt");
for (i=0;i<PDlist.length;i++){
	if (startsWith(PDlist[i],"All")){
	} else {
		if (startsWith(PDlist[i],"ante")){
			string = File.openAsString(outputsplitPD+PDlist[i]); 
			lines = split(string, "\n");   
			PD=newArray(); 
			for (j=0; j<lines.length; j++){	
				PD = Array.concat(PD,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPD,PD[j]);
			}
		}
	}

}
File.close(allPD);

allPD=File.open(outputsplitPD+"AllretroPD"+".txt");
for (i=0;i<PDlist.length;i++){
	if (startsWith(PDlist[i],"All")){
	} else {
		if (startsWith(PDlist[i],"retro")){
			string = File.openAsString(outputsplitPD+PDlist[i]); 
			lines = split(string, "\n");   
			PD=newArray(); 
			for (j=0; j<lines.length; j++){	
				PD = Array.concat(PD,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPD,PD[j]);
			}
		}
	}

}
File.close(allPD);

allPD=File.open(outputsplitPD+"AllrevPD"+".txt");
for (i=0;i<PDlist.length;i++){
	if (startsWith(PDlist[i],"All")){
	} else {
		if (startsWith(PDlist[i],"rev")){
			string = File.openAsString(outputsplitPD+PDlist[i]); 
			lines = split(string, "\n");   
			PD=newArray(); 
			for (j=0; j<lines.length; j++){	
				PD = Array.concat(PD,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPD,PD[j]);
			}
		}
	}

}
File.close(allPD);

//--------------------------------------------- 6a. Pooled Pause Frequency ---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
	movPF=File.open(outputPF+"PF_"+replace(ExperimentList[i],"/",".txt"));
		if (startsWith(ExperimentList[i],"PooledData")){	
		} else {
			inputPF=input+ExperimentList[i]+"PF/";
			listPF=getFileList(inputPF);
			
			for(k=0;k<listPF.length;k++){
				string = File.openAsString(inputPF+listPF[k]); 
				lines = split(string, "\n");   
				PF=newArray(); 
				for (j=0; j<lines.length; j++){	
					PF = Array.concat(PF,parseFloat(lines[j]));
				}
				AllPF=newArray();
				for(j=0;j<PF.length;j++){
					AllPF=Array.concat(AllPF,PF[j]);
				}
			
			for (j=0;j<AllPF.length;j++){
				print(movPF,AllPF[j]);
				}
			
			
			}	
		}
		File.close(movPF);
	}
}

PFlist=getFileList(outputPF);
allPF=File.open(outputPF+"AllPF"+".txt");
for (i=0;i<PFlist.length;i++){
	
	if (startsWith(PFlist[i],"All")){
	} else {
		
		string = File.openAsString(outputPF+PFlist[i]); 
			lines = split(string, "\n");   
			PF=newArray(); 
			for (j=0; j<lines.length; j++){	
				PF = Array.concat(PF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPF,PF[j]);
			}
	}

}
File.close(allPF);

//--------------------------------------------- 6b. Pooled Split Pause Duration ---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		movPF=File.open(outputsplitPF+"antePF_"+replace(ExperimentList[i],"/",".txt"));
			if (startsWith(ExperimentList[i],"PooledData")){	
			} else {
				inputPF=input+ExperimentList[i]+"splitPF/";
				listPF=getFileList(inputPF);
				
				for(k=0;k<listPF.length;k++){
					if (startsWith(listPF[k],"ante")){
						string = File.openAsString(inputPF+listPF[k]); 
						lines = split(string, "\n");   
						PF=newArray(); 
						for (j=0; j<lines.length; j++){	
							PF = Array.concat(PF,parseFloat(lines[j]));
						}
						for (j=0;j<PF.length;j++){
							print(movPF,PF[j]);
						}
					}
				}	
			}
		File.close(movPF);
		movPF=File.open(outputsplitPF+"retroPF_"+replace(ExperimentList[i],"/",".txt"));
			if (startsWith(ExperimentList[i],"PooledData")){	
			} else {
				inputPF=input+ExperimentList[i]+"splitPF/";
				listPF=getFileList(inputPF);
				
				for(k=0;k<listPF.length;k++){
					if (startsWith(listPF[k],"retro")){
						string = File.openAsString(inputPF+listPF[k]); 
						lines = split(string, "\n");   
						PF=newArray(); 
						for (j=0; j<lines.length; j++){	
							PF = Array.concat(PF,parseFloat(lines[j]));
						}
						for (j=0;j<PF.length;j++){
							print(movPF,PF[j]);
						}
					}
				}	
			}
		File.close(movPF);
		movPF=File.open(outputsplitPF+"revPF_"+replace(ExperimentList[i],"/",".txt"));
			if (startsWith(ExperimentList[i],"PooledData")){	
			} else {
				inputPF=input+ExperimentList[i]+"splitPF/";
				listPF=getFileList(inputPF);
				
				for(k=0;k<listPF.length;k++){
					if (startsWith(listPF[k],"rev")){
						string = File.openAsString(inputPF+listPF[k]); 
						lines = split(string, "\n");   
						PF=newArray(); 
						for (j=0; j<lines.length; j++){	
							PF = Array.concat(PF,parseFloat(lines[j]));
						}
						for (j=0;j<PF.length;j++){
							print(movPF,PF[j]);
						}
					}
				}	
			}
		File.close(movPF);
	}
}

PFlist=getFileList(outputsplitPF);
allPF=File.open(outputsplitPF+"AllantePF"+".txt");
for (i=0;i<PFlist.length;i++){
	if (startsWith(PFlist[i],"All")){
	} else {
		if (startsWith(PFlist[i],"ante")){
			string = File.openAsString(outputsplitPF+PFlist[i]); 
			lines = split(string, "\n");   
			PF=newArray(); 
			for (j=0; j<lines.length; j++){	
				PF = Array.concat(PF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPF,PF[j]);
			}
		}
	}

}
File.close(allPF);

allPF=File.open(outputsplitPF+"AllretroPF"+".txt");
for (i=0;i<PFlist.length;i++){
	if (startsWith(PFlist[i],"All")){
	} else {
		if (startsWith(PFlist[i],"retro")){
			string = File.openAsString(outputsplitPF+PFlist[i]); 
			lines = split(string, "\n");   
			PF=newArray(); 
			for (j=0; j<lines.length; j++){	
				PF = Array.concat(PF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPF,PF[j]);
			}
		}
	}

}
File.close(allPF);

allPF=File.open(outputsplitPF+"AllrevPF"+".txt");
for (i=0;i<PFlist.length;i++){
	if (startsWith(PFlist[i],"All")){
	} else {
		if (startsWith(PFlist[i],"rev")){
			string = File.openAsString(outputsplitPF+PFlist[i]); 
			lines = split(string, "\n");   
			PF=newArray(); 
			for (j=0; j<lines.length; j++){	
				PF = Array.concat(PF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPF,PF[j]);
			}
		}
	}

}
File.close(allPF);

//--------------------------------------------- 6c. Pooled Pause Frequency Per Second---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
	movPF=File.open(outputPFperSec+"PFperSec_"+replace(ExperimentList[i],"/",".txt"));
		if (startsWith(ExperimentList[i],"PooledData")){	
		} else {
			inputPFperSec=input+ExperimentList[i]+"PFperSec/";
			listPF=getFileList(inputPFperSec);
			
			for(k=0;k<listPF.length;k++){
				string = File.openAsString(inputPFperSec+listPF[k]); 
				lines = split(string, "\n");   
				PF=newArray(); 
				for (j=0; j<lines.length; j++){	
					PF = Array.concat(PF,parseFloat(lines[j]));
				}
				AllPF=newArray();
				for(j=0;j<PF.length;j++){
					AllPF=Array.concat(AllPF,PF[j]);
				}
			
			for (j=0;j<AllPF.length;j++){
				print(movPF,AllPF[j]);
				}
			
			
			}	
		}
		File.close(movPF);
	}
}

PFlist=getFileList(outputPFperSec);
allPF=File.open(outputPFperSec+"AllPF"+".txt");
for (i=0;i<PFlist.length;i++){
	
	if (startsWith(PFlist[i],"All")){
	} else {
		
		string = File.openAsString(outputPFperSec+PFlist[i]); 
			lines = split(string, "\n");   
			PF=newArray(); 
			for (j=0; j<lines.length; j++){	
				PF = Array.concat(PF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPF,PF[j]);
			}
	}

}
File.close(allPF);

//--------------------------------------------- 6d. Pooled Split Pause Frequency per Second ---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		movPF=File.open(outputsplitPFperSec+"antePFperSec_"+replace(ExperimentList[i],"/",".txt"));
			if (startsWith(ExperimentList[i],"PooledData")){	
			} else {
				inputPF=input+ExperimentList[i]+"splitPFperSec/";
				listPF=getFileList(inputPF);
				
				for(k=0;k<listPF.length;k++){
					if (startsWith(listPF[k],"ante")){
						string = File.openAsString(inputPF+listPF[k]); 
						lines = split(string, "\n");   
						PF=newArray(); 
						for (j=0; j<lines.length; j++){	
							PF = Array.concat(PF,parseFloat(lines[j]));
						}
						for (j=0;j<PF.length;j++){
							print(movPF,PF[j]);
						}
					}
				}	
			}
		File.close(movPF);
		movPF=File.open(outputsplitPFperSec+"retroPFperSec_"+replace(ExperimentList[i],"/",".txt"));
			if (startsWith(ExperimentList[i],"PooledData")){	
			} else {
				inputPF=input+ExperimentList[i]+"splitPFperSec/";
				listPF=getFileList(inputPF);
				
				for(k=0;k<listPF.length;k++){
					if (startsWith(listPF[k],"retro")){
						string = File.openAsString(inputPF+listPF[k]); 
						lines = split(string, "\n");   
						PF=newArray(); 
						for (j=0; j<lines.length; j++){	
							PF = Array.concat(PF,parseFloat(lines[j]));
						}
						for (j=0;j<PF.length;j++){
							print(movPF,PF[j]);
						}
					}
				}	
			}
		File.close(movPF);
		movPF=File.open(outputsplitPFperSec+"revPFperSec_"+replace(ExperimentList[i],"/",".txt"));
			if (startsWith(ExperimentList[i],"PooledData")){	
			} else {
				inputPF=input+ExperimentList[i]+"splitPFperSec/";
				listPF=getFileList(inputPF);
				
				for(k=0;k<listPF.length;k++){
					if (startsWith(listPF[k],"rev")){
						string = File.openAsString(inputPF+listPF[k]); 
						lines = split(string, "\n");   
						PF=newArray(); 
						for (j=0; j<lines.length; j++){	
							PF = Array.concat(PF,parseFloat(lines[j]));
						}
						for (j=0;j<PF.length;j++){
							print(movPF,PF[j]);
						}
					}
				}	
			}
		File.close(movPF);
	}
}

PFlist=getFileList(outputsplitPFperSec);
allPF=File.open(outputsplitPFperSec+"AllantePFperSec"+".txt");
for (i=0;i<PFlist.length;i++){
	if (startsWith(PFlist[i],"All")){
	} else {
		if (startsWith(PFlist[i],"ante")){
			string = File.openAsString(outputsplitPFperSec+PFlist[i]); 
			lines = split(string, "\n");   
			PF=newArray(); 
			for (j=0; j<lines.length; j++){	
				PF = Array.concat(PF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPF,PF[j]);
			}
		}
	}

}
File.close(allPF);

allPF=File.open(outputsplitPFperSec+"AllretroPFperSec"+".txt");
for (i=0;i<PFlist.length;i++){
	if (startsWith(PFlist[i],"All")){
	} else {
		if (startsWith(PFlist[i],"retro")){
			string = File.openAsString(outputsplitPFperSec+PFlist[i]); 
			lines = split(string, "\n");   
			PF=newArray(); 
			for (j=0; j<lines.length; j++){	
				PF = Array.concat(PF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPF,PF[j]);
			}
		}
	}

}
File.close(allPF);

allPF=File.open(outputsplitPFperSec+"AllrevPFperSec"+".txt");
for (i=0;i<PFlist.length;i++){
	if (startsWith(PFlist[i],"All")){
	} else {
		if (startsWith(PFlist[i],"rev")){
			string = File.openAsString(outputsplitPFperSec+PFlist[i]); 
			lines = split(string, "\n");   
			PF=newArray(); 
			for (j=0; j<lines.length; j++){	
				PF = Array.concat(PF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allPF,PF[j]);
			}
		}
	}

}
File.close(allPF);

//--------------------------------------------- 7a. Pooled Flux ---------------------------------------------------
anteFlux=newArray();
retroFlux=newArray();
reversalFlux=newArray();
for (i=0;i<ExperimentList.length;i++){	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputFlux=input+ExperimentList[i]+"Flux/";
		FluxList=getFileList(inputFlux);
		string = File.openAsString(inputFlux+FluxList[0]); 
		lines = split(string, "\n");   
		Flux=newArray(); 
		for (j=0; j<lines.length; j++){	
			Flux = Array.concat(Flux,parseFloat(lines[j]));
		}
		
		anteFlux=Array.concat(anteFlux,abs(Flux[0]));
		retroFlux=Array.concat(retroFlux,abs(Flux[1]));
		reversalFlux=Array.concat(reversalFlux,abs(Flux[2]));
	}
		Ante=File.open(outputFlux+"anteFlux"+".txt");
		for (j=0;j<anteFlux.length;j++){
			print(Ante,anteFlux[j]);
		}
		File.close(Ante);
		Retro=File.open(outputFlux+"retroFlux"+".txt");
		for (j=0;j<retroFlux.length;j++){
			print(Retro,retroFlux[j]);
		}
		File.close(Retro);
		Reversal=File.open(outputFlux+"reversalFlux"+".txt");
		for (j=0;j<reversalFlux.length;j++){
			print(Reversal,reversalFlux[j]);
		}
		File.close(Reversal);
}

//--------------------------------------------- 7b. Pooled Densities ---------------------------------------------------
anteDensity=newArray();
retroDensity=newArray();
reversalDensity=newArray();
stationaryDensity=newArray();
for (i=0;i<ExperimentList.length;i++){	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputDensity=input+ExperimentList[i]+"Density/";
		DensityList=getFileList(inputDensity);
		string = File.openAsString(inputDensity+DensityList[0]); 
		lines = split(string, "\n");   
		Density=newArray(); 
		for (j=0; j<lines.length; j++){	
			Density = Array.concat(Density,parseFloat(lines[j]));
		}
		
		anteDensity=Array.concat(anteDensity,abs(Density[0]));
		retroDensity=Array.concat(retroDensity,abs(Density[1]));
		reversalDensity=Array.concat(reversalDensity,abs(Density[2]));
		stationaryDensity=Array.concat(stationaryDensity,abs(Density[3]));
	}
		Ante=File.open(outputDensity+"anteDensity"+".txt");
		for (j=0;j<anteDensity.length;j++){
			print(Ante,anteDensity[j]);
		}
		File.close(Ante);
		Retro=File.open(outputDensity+"retroDensity"+".txt");
		for (j=0;j<retroDensity.length;j++){
			print(Retro,retroDensity[j]);
		}
		File.close(Retro);
		Reversal=File.open(outputDensity+"reversalDensity"+".txt");
		for (j=0;j<reversalDensity.length;j++){
			print(Reversal,reversalDensity[j]);
		}
		File.close(Reversal);
		Stationary=File.open(outputDensity+"stationaryDensity"+".txt");
		for (j=0;j<stationaryDensity.length;j++){
			print(Stationary,stationaryDensity[j]);
		}
		File.close(Stationary);
}


//--------------------------------------------- 8a. Pooled Run Length ---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){

	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputRL=input+ExperimentList[i]+"RL/";
		RLList=getFileList(inputRL);
		anteRL=newArray();
		retroRL=newArray();
			for(k=0;k<RLList.length;k++){
			string = File.openAsString(inputRL+RLList[k]); 
			lines = split(string, "\n");   
			RL=newArray(); 
			
			for (j=0; j<lines.length; j++){	
				RL = Array.concat(RL,parseFloat(lines[j]));
			}
	
			for (j=0;j<RL.length;j++){
			if (RL[j]>0){
			anteRL=Array.concat(anteRL,abs(RL[j]));
			}
			if (RL[j]<0){
			retroRL=Array.concat(retroRL,abs(RL[j]));
			}
			}
			
			AnteRL=File.open(outputRL+"anteRL_"+replace(ExperimentList[i],"/",".txt"));
			for (j=0;j<anteRL.length;j++){
				print(AnteRL,anteRL[j]);
			}
			File.close(AnteRL);
			RetroRL=File.open(outputRL+"retroRL_"+replace(ExperimentList[i],"/",".txt"));
			for (j=0;j<retroRL.length;j++){
				print(RetroRL,retroRL[j]);
			}
			File.close(RetroRL);
		}
	}
}

AllRLList=getFileList(outputRL);
AllAnteRL=newArray();
AllRetroRL=newArray();
for (i=0;i<AllRLList.length;i++){
	if(startsWith(AllRLList[i],"ante")){
		string = File.openAsString(outputRL+AllRLList[i]); 
		lines = split(string, "\n");   
		RL=newArray(); 
		for (j=0; j<lines.length; j++){
		RL = Array.concat(RL,parseFloat(lines[j]));	
		}
		for (j=0;j<RL.length;j++){
		AllAnteRL=Array.concat(AllAnteRL,abs(RL[j]));
		}
		allAnteRL=File.open(outputRL+"AllAnteRL"+".txt");
		for (j=0;j<AllAnteRL.length;j++){
			print(allAnteRL,AllAnteRL[j]);
		}
		File.close(allAnteRL);
	}
	if(startsWith(AllRLList[i],"retro")){
		string = File.openAsString(outputRL+AllRLList[i]); 
		lines = split(string, "\n");   
		RL=newArray(); 
		for (j=0; j<lines.length; j++){
		RL = Array.concat(RL,parseFloat(lines[j]));	
		}
		for (j=0;j<RL.length;j++){
		AllRetroRL=Array.concat(AllRetroRL,abs(RL[j]));
		}
		allRetroRL=File.open(outputRL+"AllRetroRL"+".txt");
		for (j=0;j<AllRetroRL.length;j++){
			print(allRetroRL,AllRetroRL[j]);
		}
		File.close(allRetroRL);
	}
		
}

//--------------------------------------------- 8b. Pooled Combined Run Length ---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){

	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputRL=input+ExperimentList[i]+"combinedRL/";
		RLList=getFileList(inputRL);
		anteRL=newArray();
		retroRL=newArray();
			for(k=0;k<RLList.length;k++){
			string = File.openAsString(inputRL+RLList[k]); 
			lines = split(string, "\n");   
			RL=newArray(); 
			
			for (j=0; j<lines.length; j++){	
				RL = Array.concat(RL,parseFloat(lines[j]));
			}
	
			for (j=0;j<RL.length;j++){
			if (RL[j]>0){
			anteRL=Array.concat(anteRL,abs(RL[j]));
			}
			if (RL[j]<0){
			retroRL=Array.concat(retroRL,abs(RL[j]));
			}
			}
			
			AnteRL=File.open(outputcombinedRL+"anteRLCombined_"+replace(ExperimentList[i],"/",".txt"));
			for (j=0;j<anteRL.length;j++){
				print(AnteRL,anteRL[j]);
			}
			File.close(AnteRL);
			RetroRL=File.open(outputcombinedRL+"retroRLCombined_"+replace(ExperimentList[i],"/",".txt"));
			for (j=0;j<retroRL.length;j++){
				print(RetroRL,retroRL[j]);
			}
			File.close(RetroRL);
		}
	}
}

AllRLList=getFileList(outputcombinedRL);
AllAnteRL=newArray();
AllRetroRL=newArray();
for (i=0;i<AllRLList.length;i++){
	if(startsWith(AllRLList[i],"ante")){
		string = File.openAsString(outputcombinedRL+AllRLList[i]); 
		lines = split(string, "\n");   
		RL=newArray(); 
		for (j=0; j<lines.length; j++){
		RL = Array.concat(RL,parseFloat(lines[j]));	
		}
		for (j=0;j<RL.length;j++){
		AllAnteRL=Array.concat(AllAnteRL,abs(RL[j]));
		}
		allAnteRL=File.open(outputcombinedRL+"AllAnteRLCombined"+".txt");
		for (j=0;j<AllAnteRL.length;j++){
			print(allAnteRL,AllAnteRL[j]);
		}
		File.close(allAnteRL);
	}
	if(startsWith(AllRLList[i],"retro")){
		string = File.openAsString(outputcombinedRL+AllRLList[i]); 
		lines = split(string, "\n");   
		RL=newArray(); 
		for (j=0; j<lines.length; j++){
		RL = Array.concat(RL,parseFloat(lines[j]));	
		}
		for (j=0;j<RL.length;j++){
		AllRetroRL=Array.concat(AllRetroRL,abs(RL[j]));
		}
		allRetroRL=File.open(outputcombinedRL+"AllRetroRLCombined"+".txt");
		for (j=0;j<AllRetroRL.length;j++){
			print(allRetroRL,AllRetroRL[j]);
		}
		File.close(allRetroRL);
	}
		
}

//--------------------------------------------- 9a. Pooled Switch Frequency ---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
	movSF=File.open(outputSF+"SF_"+replace(ExperimentList[i],"/",".txt"));
		if (startsWith(ExperimentList[i],"PooledData")){	
		} else {
			inputSF=input+ExperimentList[i]+"SF/";
			listSF=getFileList(inputSF);
			
			for(k=0;k<listSF.length;k++){
				string = File.openAsString(inputSF+listSF[k]); 
				lines = split(string, "\n");   
				SF=newArray(); 
				for (j=0; j<lines.length; j++){	
					SF = Array.concat(SF,parseFloat(lines[j]));
				}
				AllSF=newArray();
				for(j=0;j<SF.length;j++){
					AllSF=Array.concat(AllSF,SF[j]);
				}
			
			for (j=0;j<AllSF.length;j++){
				print(movSF,AllSF[j]);
				}
			
			
			}	
		}
		File.close(movSF);
	}
}

SFlist=getFileList(outputSF);
allSF=File.open(outputSF+"AllSF"+".txt");
for (i=0;i<SFlist.length;i++){
	
	if (startsWith(SFlist[i],"All")){
	} else {
		
		string = File.openAsString(outputSF+SFlist[i]); 
			lines = split(string, "\n");   
			SF=newArray(); 
			for (j=0; j<lines.length; j++){	
				SF = Array.concat(SF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allSF,SF[j]);
			}
	}

}
File.close(allSF);

//--------------------------------------------- 9b. Pooled Switch Frequency per Second---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
	movSF=File.open(outputSFperSec+"SFperSec_"+replace(ExperimentList[i],"/",".txt"));
		if (startsWith(ExperimentList[i],"PooledData")){	
		} else {
			inputSFperSec=input+ExperimentList[i]+"SFperSec/";
			listSF=getFileList(inputSFperSec);
			
			for(k=0;k<listSF.length;k++){
				string = File.openAsString(inputSFperSec+listSF[k]); 
				lines = split(string, "\n");   
				SF=newArray(); 
				for (j=0; j<lines.length; j++){	
					SF = Array.concat(SF,parseFloat(lines[j]));
				}
				AllSF=newArray();
				for(j=0;j<SF.length;j++){
					AllSF=Array.concat(AllSF,SF[j]);
				}
			
			for (j=0;j<AllSF.length;j++){
				print(movSF,AllSF[j]);
				}
			
			
			}	
		}
		File.close(movSF);
	}
}

SFlist=getFileList(outputSFperSec);
allSF=File.open(outputSFperSec+"AllSFperSec"+".txt");
for (i=0;i<SFlist.length;i++){
	
	if (startsWith(SFlist[i],"All")){
	} else {
		
		string = File.openAsString(outputSFperSec+SFlist[i]); 
			lines = split(string, "\n");   
			SF=newArray(); 
			for (j=0; j<lines.length; j++){	
				SF = Array.concat(SF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allSF,SF[j]);
			}
	}

}
File.close(allSF);
//--------------------------------------------- 9c. Pooled Switch Frequency Reverals---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
	movSF=File.open(outputrevSF+"revSF_"+replace(ExperimentList[i],"/",".txt"));
		if (startsWith(ExperimentList[i],"PooledData")){	
		} else {
			inputSF=input+ExperimentList[i]+"SF/";
			listSF=getFileList(inputSF);
			
			for(k=0;k<listSF.length;k++){
				string = File.openAsString(inputSF+listSF[k]); 
				lines = split(string, "\n");   
				SF=newArray(); 
				for (j=0; j<lines.length; j++){	
					if (parseFloat(lines[j])!=0){
						SF = Array.concat(SF,parseFloat(lines[j]));
					}
				}
				AllSF=newArray();
				for(j=0;j<SF.length;j++){
					AllSF=Array.concat(AllSF,SF[j]);
				}
			
			for (j=0;j<AllSF.length;j++){
				print(movSF,AllSF[j]);
				}
			
			
			}	
		}
		File.close(movSF);
	}
}

SFlist=getFileList(outputrevSF);
allSF=File.open(outputrevSF+"AllrevSF"+".txt");
for (i=0;i<SFlist.length;i++){
	
	if (startsWith(SFlist[i],"All")){
	} else {
		
		string = File.openAsString(outputrevSF+SFlist[i]); 
			lines = split(string, "\n");   
			SF=newArray(); 
			for (j=0; j<lines.length; j++){	
				SF = Array.concat(SF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allSF,SF[j]);
			}
	}

}
File.close(allSF);

//--------------------------------------------- 9d. Pooled Switch Frequency per Second Reverals---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
	movSF=File.open(outputrevSFperSec+"revSFperSec_"+replace(ExperimentList[i],"/",".txt"));
		if (startsWith(ExperimentList[i],"PooledData")){	
		} else {
			inputSFperSec=input+ExperimentList[i]+"SFperSec/";
			listSF=getFileList(inputSFperSec);
			
			for(k=0;k<listSF.length;k++){
				string = File.openAsString(inputSFperSec+listSF[k]); 
				lines = split(string, "\n");   
				SF=newArray(); 
				for (j=0; j<lines.length; j++){	
					if (parseFloat(lines[j])!=0){
						SF = Array.concat(SF,parseFloat(lines[j]));
					}
				}
				AllSF=newArray();
				for(j=0;j<SF.length;j++){
					AllSF=Array.concat(AllSF,SF[j]);
				}
			
			for (j=0;j<AllSF.length;j++){
				print(movSF,AllSF[j]);
				}
			
			
			}	
		}
		File.close(movSF);
	}
}

SFlist=getFileList(outputrevSFperSec);
allSF=File.open(outputrevSFperSec+"AllrevSFperSec"+".txt");
for (i=0;i<SFlist.length;i++){
	
	if (startsWith(SFlist[i],"All")){
	} else {
		
		string = File.openAsString(outputrevSFperSec+SFlist[i]); 
			lines = split(string, "\n");   
			SF=newArray(); 
			for (j=0; j<lines.length; j++){	
				SF = Array.concat(SF,parseFloat(lines[j]));
			}
			for (j=0; j<lines.length;j++){
				print(allSF,SF[j]);
			}
	}

}
File.close(allSF);

//--------------------------------------------- 10. Pooled Percent Time in Motion (PM)---------------------------------------------------
allantePM=newArray();
allretroPM=newArray();
allpausesPM=newArray();
for (i=0;i<ExperimentList.length;i++){	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		antePM=newArray();
		retroPM=newArray();
		pausesPM=newArray();
		inputPM=input+ExperimentList[i]+"PM/";
		PMList=getFileList(inputPM);
		for (k=0; k<PMList.length; k++){
			string = File.openAsString(inputPM+PMList[k]); 
			lines = split(string, "\n");   
			PM=newArray(); 
			for (j=0; j<lines.length; j++){	
				PM = Array.concat(PM,parseFloat(lines[j]));
			}
			
			antePM=Array.concat(antePM,abs(PM[0]));
			retroPM=Array.concat(retroPM,abs(PM[1]));
			pausesPM=Array.concat(pausesPM,abs(PM[2]));

			allantePM=Array.concat(allantePM,abs(PM[0]));
			allretroPM=Array.concat(allretroPM,abs(PM[1]));
			allpausesPM=Array.concat(allpausesPM,abs(PM[2]));
		}
	
		Ante=File.open(outputPM+"antePM"+replace(ExperimentList[i],"/",".txt"));
		for (j=0;j<antePM.length;j++){
			print(Ante,antePM[j]);
		}
		File.close(Ante);
		Ante=File.open(outputPM+"AllantePM"+".txt");
		for (j=0;j<allantePM.length;j++){
			print(Ante,allantePM[j]);
		}
		File.close(Ante);
		Retro=File.open(outputPM+"retroPM"+replace(ExperimentList[i],"/",".txt"));
		for (j=0;j<retroPM.length;j++){
			print(Retro,retroPM[j]);
		}
		File.close(Retro);
		Retro=File.open(outputPM+"AllretroPM"+".txt");
		for (j=0;j<allretroPM.length;j++){
			print(Retro,allretroPM[j]);
		}
		File.close(Retro);
		Pauses=File.open(outputPM+"pausesPM"+replace(ExperimentList[i],"/",".txt"));
		for (j=0;j<pausesPM.length;j++){
			print(Pauses,pausesPM[j]);
		}
		File.close(Pauses);
		Pauses=File.open(outputPM+"AllpausesPM"+".txt");
		for (j=0;j<allpausesPM.length;j++){
			print(Pauses,allpausesPM[j]);
		}
		File.close(Pauses);
	}
}

//--------------------------------------------- 10a. Make Montage Kymographs---------------------------------------------------
setBatchMode(true);
for (i=0;i<ExperimentList.length;i++){
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputKymo=input+ExperimentList[i]+"Kymograph/";
		listKymo=getFileList(inputKymo);
		open(inputKymo+listKymo[0]);
	}
}

if (nImages>1){
	run("Images to Stack", "method=[Copy (center)] name=Stack title=[] use");
	string="columns=1 rows="+nSlices+" scale=1 first=1 last="+nSlices+" increment=1 border=0 font=20 label";
	run("Make Montage...", string);
	save(outputKymo+"Montage_Kymographs.tif");
}
if (nImages==1){
	save(outputKymo+"Montage_Kymographs.tif");
}
run("Close");
run("Close");

//--------------------------------------------- 10b. Make Montage Kymographs with CargoPopulation---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputKymo=input+ExperimentList[i]+"Kymograph/";
		listKymo=getFileList(inputKymo);
		inputROI=input+ExperimentList[i]+"CP_ROIs/";
		listROI=getFileList(inputROI);
		open(inputKymo+listKymo[0]);
		ID=getImageID;
		name=File.nameWithoutExtension;
		roiManager("Open",inputROI+listROI[0]);
		roiManager("Show All");
		run("Flatten");
		selectImage(ID);
		run("Close");
		rename(name);	
		roiManager("Reset");
	}
}

if (nImages>1){
	run("Images to Stack", "method=[Copy (center)] name=Stack title=[] use");
	string="columns=1 rows="+nSlices+" scale=1 first=1 last="+nSlices+" increment=1 border=0 font=20 label";
	run("Make Montage...", string);
	save(outputKymo+"Montage_Kymographs_CP.tif");
}
if (nImages==1){
	save(outputKymo+"Montage_Kymographs_CP.tif");
}
run("Close");
run("Close");

//--------------------------------------------- 10b. Make Montage Kymographs with Segments---------------------------------------------------
for (i=0;i<ExperimentList.length;i++){
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		inputKymo=input+ExperimentList[i]+"Kymograph/";
		listKymo=getFileList(inputKymo);
		inputROI=input+ExperimentList[i]+"Segment_ROIs/";
		listROI=getFileList(inputROI);
		open(inputKymo+listKymo[0]);
		ID=getImageID;
		name=File.nameWithoutExtension;
		roiManager("Open",inputROI+listROI[0]);
		array=newArray();
		for (n=0; n<roiManager("Count");n++){
			roiManager("Select",n);
			if (startsWith(Roi.getName,"Segment")){
				
			} else {
				array=Array.concat(array,n);
			}
		}
		roiManager("Select",array);
		roiManager("Delete");
		roiManager("Show All");
		run("Flatten");
		for (n=0;n<roiManager("Count");n++){
			roiManager("Select",n);
			getSelectionCoordinates(x,y);
			for (m=0;m<x.length;m++){
				makeRectangle(x[m],y[m],2,2);
				run("Colors...", "foreground=black background=black selection=black");
				run("Fill", "slice");
			}
		}
		selectImage(ID);
		run("Close");
		rename(name);	
		roiManager("Reset");
	}
}
if (nImages>1){
	run("Images to Stack", "method=[Copy (center)] name=Stack title=[] use");
	string="columns=1 rows="+nSlices+" scale=1 first=1 last="+nSlices+" increment=1 border=0 font=20 label";
	run("Make Montage...", string);
	save(outputKymo+"Montage_Kymographs_Segments.tif");
}
if (nImages==1){
	save(outputKymo+"Montage_Kymographs_Segments.tif");
}
run("Close");
run("Close");