// This macro calculates the cargo population and switch frequency
// Tracks are assigned directions and labeled accordingly:
// green: anterograde
// red: retrograde
// yellow: reversing
// blue: stationary
// The cargo population is calculated either numerically or as percentage of the number of tracks and saved
// in the folder CP as txt files either CPNum (numeric) or CPPCT (percent).
// This enables the clicking and calculation of segmental properties for only partial tracks. The thresholds of 10 frames was chosen to 
// account for unprecise clicking at the beginning and the end of the movie.
//
// The Switch frequency is either calculated as switches/track (SF) or switches per second (SFperSec). 

//----------------------------BEGIN MACRO------------------------------------------------------------------
// resetting imageJ
print("\\Clear");			
roiManager("Reset");

// reading in the imaging parameters
Dialog.create("ImagingParameters");
Dialog.addNumber("Pixel Size in µm",0.16);
Dialog.addNumber("Frame Rate in /sec", 10);
Dialog.show;
PxlSize=Dialog.getNumber();
Frame=Dialog.getNumber();

//choose input folder: Folder that contains all movies for an experiment
input = getDirectory("Choose Experiment folder (the folder that contains all your movies)"); 
//----------------------------process the input folder-----------------------------------------------------

processFolder(input);

function processFolder(input){
	list=getFileList(input);
		
	for (i=0;i<list.length;i++){
		if (startsWith(list[i], "PooledData")){
			
		} else {

			// name of the current movie
			movie=list[i];
			movie=substring(movie,0,lengthOf(movie)-1);
			print(movie);
			
			// make output files
			File.makeDirectory(input+list[i]+"CP_ROIs");
			File.makeDirectory(input+list[i]+"X_coor");
			File.makeDirectory(input+list[i]+"Y_coor");
			File.makeDirectory(input+list[i]+"CargoPopulation");
			File.makeDirectory(input+list[i]+"SF");
			File.makeDirectory(input+list[i]+"SFperSec");
			File.makeDirectory(input+list[i]+"RevPoint");
			File.makeDirectory(input+list[i]+"Info");
			File.makeDirectory(input+list[i]+"Flux");
			File.makeDirectory(input+list[i]+"Density");
			
			//set outputs
			outputCP_ROI=input+list[i]+"CP_ROIs/";
			outputX_coor=input+list[i]+"X_coor/";
			outputY_coor=input+list[i]+"Y_coor/";
			outputCPNum=input+list[i]+"CargoPopulation/";
			outputSF=input+list[i]+"SF/";
			outputSFperSec=input+list[i]+"SFperSec/";
			outputRevPoint=input+list[i]+"RevPoint/";
			outputInfo=input+list[i]+"Info/";
			outputFlux=input+list[i]+"Flux/";
			outputDensity=input+list[i]+"Density/";
			// record imaging parameters
			parameters=File.open(outputInfo+"ImageParameters_"+movie+".txt");
			print(parameters,PxlSize);
			print(parameters,Frame);
			File.close(parameters);
			
			//set inputs for specific files
			KymoInput=input+list[i]+"Kymograph/";
			RoiInput=input+list[i]+"ROIs/";
			//get file lists for above input folders
			Kymolist=getFileList(KymoInput);
			Roilist=getFileList(RoiInput);
			//run function CargoPopulation
			CargoPopulation(input, Kymolist[0], Roilist[0]);
		}
	}
}

//-------------FUNCTION TO CALCULATE THE CARGO POPULATION-----------------
function CargoPopulation(input,file1,file2){
	roiManager("Reset");
	// open kymograph
	open(input+list[i]+"Kymograph/"+file1);
	//open ROIs: clicked tracks
	roiManager("Open",input+list[i]+"ROIs/"+file2);
	setOption("Show All",true);
	roiCount=roiManager("count");

	// get XY coordinates for each track. to get higher resolution, coordinates are 
	// of the polyline are multiplied by 10, recorded and then divided by 10, so as to 
	// obtain one decimal for the XY coordinates
	for (n=0;n<roiCount;n++){
		roiManager("Select",n);	
		name=Roi.getName;
		setBatchMode(true);
		getSelectionCoordinates(xTemp, yTemp);
		j=yTemp.length;
		x=newArray();
		y=newArray();
		if (yTemp[0] > yTemp[yTemp.length-1]){
			for (i=0;i<yTemp.length;i++){
				j=j-1;
				x=Array.concat(x,xTemp[j]);
				y=Array.concat(y,yTemp[j]);	
			}
		} else {
			for (i=0;i<yTemp.length;i++){
				x=Array.concat(x,xTemp[i]);
				y=Array.concat(y,yTemp[i]);
			}
		}
		/*
		//calculate length of track and minimum length
		length=(y[y.length-1]-y[0])/Frame;
		length2=((y[y.length-1]-y[0])+1)/Frame;
		h=getHeight;
		minLength=(h-10)/Frame;
		*/
		
		length2=((y[y.length-1]-y[0])+1)/Frame;
		
		//calculate coordinates with decimal points
		xcoor=newArray();
		ycoor=newArray();
		for (i=0;i<x.length-1;i++){
			m=(y[i+1]-y[i])/(x[i+1]-x[i]);
			b=y[i]- (m*x[i]);
			Y=y[i]-0.1;
			for (j=y[i]*10;j<y[i+1]*10;j++){
				if (x[i+1]==x[i]){
					Y=Y+0.1;
					xcoor=Array.concat(xcoor,x[i]);
					ycoor=Array.concat(ycoor,Y);	
				} else {
					Y=Y+0.1;
					X=(Y-b)/m;
					xcoor=Array.concat(xcoor,X);
					ycoor=Array.concat(ycoor,Y);
				}
			}
		}
		xcoor=Array.concat(xcoor,x[x.length-1]);
		ycoor=Array.concat(ycoor,y[y.length-1]);

		// print XY coordinates into txt file
		fxcoor=File.open(outputX_coor+"X_Coor_"+name+".txt");
		for (i=0;i<xcoor.length;i++){
			print(fxcoor,xcoor[i]);
		}
		File.close(fxcoor);
		
		fycoor=File.open(outputY_coor+"Y_Coor_"+name+".txt");
		for (i=0;i<xcoor.length;i++){
			print(fycoor,ycoor[i]);
		}
		File.close(fycoor);
		
		// calculating the track center of the track
		sumX=(xcoor[0]);
		for (i=1;i<xcoor.length;i++){
			sumX=sumX+xcoor[i];
		}					
		trackCenter=sumX/(xcoor.length);
		Array.getStatistics(xcoor, min, max);
		dXtrackCenterMin=PxlSize*(abs(min-trackCenter));
		dXtrackCenterMax=PxlSize*(abs(max-trackCenter));
		//print(name,"\t", trackCenter,"\t", min,"\t", max,"\t", dXtrackCenterMin,"\t",dXtrackCenterMax);
		
		//test if stationary
		//threshold of cmin=0.35 µm, if a particle translocates more than 0.5 µm from its track center it will be regarded as mobile
		cmin=0.35;				
		dXtrackCenter=newArray();
		SF=0;					//default value for switch frequency, will be changed if track is reversal track
		SFperSec=0;				//default value for switch frequency, will be changed if track is reversal track
		
		for (i=0;i<xcoor.length-1;i++){
			dXtrackCenter=Array.concat(dXtrackCenter,PxlSize*abs(xcoor[i]-trackCenter));

			// testing if mobile tracks are reversals 
			if (dXtrackCenter[i]>cmin){
				roiManager("Rename","mobile_"+name);
		//------------------------------------------ define reversals----------------------------------------------------				
				next=0;
				prev=0;
				limit=0;
				// cminRV=0.35 µm. This is a good threshold for vesicles, for mitochondria may use larger threshold.
				// cminRV is the distance that two points have to have to be counted as reversals if their velocity vectors 
				// have opposite directions
				cminRV=0.35;				
				index=newArray();
				direction=newArray();
				RevPointX=newArray();	
				RevPointY=newArray();	
					
				//check if track is reversing
				for (i=1;i<xcoor.length-1;i++){
					done=false;
					prev=i-1;
					next=i+1;

					// distance between two points at frame i
					dXprev=PxlSize*(parseFloat(xcoor[i])-parseFloat(xcoor[prev]));				
					dXnext=PxlSize*(parseFloat(xcoor[i])-parseFloat(xcoor[next]));
					if(prev==0){
						done=true;
					}
					if(next==xcoor.length-1){
						done=true;
					}

					// limit is the frame at which a previous reversal has occured
					if(prev==limit){
						done=true;
					}
					while(abs(dXprev)<=cminRV && abs(dXnext)<=cminRV && !done){
						prev--;
						next++;
						dXprev=PxlSize*(parseFloat(xcoor[i])-parseFloat(xcoor[prev]));
						dXnext=PxlSize*(parseFloat(xcoor[i])-parseFloat(xcoor[next]));
						if(next==xcoor.length-1){
							done=true;
						}
						if(prev==0){
							done=true;
						}
						if(prev==limit){
							done=true;
						}
					}
						
					if(abs(dXprev)>=cminRV && abs(dXnext)<=cminRV){
						while(abs(dXnext)<=cminRV && !done){
							next++;
							dXnext=PxlSize*(parseFloat(xcoor[i])-parseFloat(xcoor[next]));
							if(next==xcoor.length-1){
								done=true;
							}
						}		
					}
					if(abs(dXprev)<=cminRV && abs(dXnext)>=cminRV){
						while(abs(dXprev)<=cminRV && !done){
							prev--;
							dXprev=PxlSize*(parseFloat(xcoor[i])-parseFloat(xcoor[prev]));
							if(prev==0){
							done=true;
							}
							if(prev==limit){
								done=true;
							}	
						}
					}		
					
					if(abs(dXprev)>=cminRV && abs(dXnext)>=cminRV){
						if(dXprev<0 && dXnext<0){
							//waitForUser;
							done=true;
							begin=i;
							end=next;

							index=Array.concat(index,i);
							index=Array.concat(index,next);
							roiManager("Rename","reversal_"+name);
							roiManager("Set Color", "yellow");
							roiManager("Set Line Width", 2);
	
							// calculate reversal point
							xcoorTemp=newArray();
							ycoorTemp=newArray();
							RevPointTemp=newArray();
							for (k=begin; k<end; k++){
								xcoorTemp=Array.concat(xcoorTemp,xcoor[k]);
								ycoorTemp=Array.concat(ycoorTemp,ycoor[k]);
							}
									
							xcoorTemp=Array.concat(xcoorTemp,abs(xcoor[end]));
							ycoorTemp=Array.concat(ycoorTemp,ycoor[end]);		
							if (dXprev<0 && dXnext<0){
								Array.getStatistics(xcoorTemp,min,max);
								RevPointTemp=Array.concat(RevPointTemp,min);
							}
							if (dXprev>0 && dXnext>0){
								Array.getStatistics(xcoorTemp,min,max);
								RevPointTemp=Array.concat(RevPointTemp,max);
							}
				
							RevPointTemp2=newArray();
							frameTemp=newArray();
							frame=newArray();
							for (l=0;l<xcoorTemp.length;l++){
								if (RevPointTemp[0]==xcoorTemp[l]){
									RevPointTemp2=Array.concat(RevPointTemp2,ycoorTemp[l]);
									frameTemp=Array.concat(frameTemp,begin+l);
								}
							}
				
							// if reversal happens at pause (dX == 0) then record first and last XY coordinates
							if (RevPointTemp2.length>1){
								RevPointY=Array.concat(RevPointY,RevPointTemp2[0]);
								RevPointY=Array.concat(RevPointY,RevPointTemp2[RevPointTemp2.length-1]);
								frame=Array.concat(frame, frameTemp[frameTemp.length-1]);
								RevPointX=Array.concat(RevPointX,RevPointTemp[0]);
								RevPointX=Array.concat(RevPointX,RevPointTemp[0]);
								i=frameTemp[0];
								limit=i;
							}
				
							if (RevPointTemp2.length==1){
								RevPointY=Array.concat(RevPointY,RevPointTemp2[0]);
								RevPointX=Array.concat(RevPointX,RevPointTemp[0]);
								frame=Array.concat(frame, frameTemp[0]);
								i=frameTemp[0];
								limit=i;
								}
							}
								
							if(dXprev>0 && dXnext>0){
								//waitForUser;
								done=true;
								begin=i;
								end=next;
								index=Array.concat(index,i);
								index=Array.concat(index,next);
								roiManager("Rename","reversal_"+name);
								roiManager("Set Color", "yellow");
								roiManager("Set Line Width", 2);
							
								// calculate reversal point
								xcoorTemp=newArray();
								ycoorTemp=newArray();
								RevPointTemp=newArray();
								for (k=begin; k<end; k++){
									xcoorTemp=Array.concat(xcoorTemp,xcoor[k]);
									ycoorTemp=Array.concat(ycoorTemp,ycoor[k]);
								}
								
								xcoorTemp=Array.concat(xcoorTemp,abs(xcoor[end]));
								ycoorTemp=Array.concat(ycoorTemp,ycoor[end]);		
								if (dXprev<0 && dXnext<0){
									Array.getStatistics(xcoorTemp,min,max);
									RevPointTemp=Array.concat(RevPointTemp,min);
								}
								if (dXprev>0 && dXnext>0){
									Array.getStatistics(xcoorTemp,min,max);
									RevPointTemp=Array.concat(RevPointTemp,max);
								}
				
								RevPointTemp2=newArray();
								frameTemp=newArray();
								frame=newArray();
								for (l=0;l<xcoorTemp.length;l++){
									if (RevPointTemp[0]==xcoorTemp[l]){
										RevPointTemp2=Array.concat(RevPointTemp2,ycoorTemp[l]);
										frameTemp=Array.concat(frameTemp,begin+l);
									}
								}
				
								// if reversal happens at pause (dX == 0) then record first and last XY coordinates
								if (RevPointTemp2.length>1){
									RevPointY=Array.concat(RevPointY,RevPointTemp2[0]);
									RevPointY=Array.concat(RevPointY,RevPointTemp2[RevPointTemp2.length-1]);
									frame=Array.concat(frame, frameTemp[frameTemp.length-1]);
									RevPointX=Array.concat(RevPointX,RevPointTemp[0]);
									RevPointX=Array.concat(RevPointX,RevPointTemp[0]);
									i=frameTemp[0];
									limit=i;
								}
				
								if (RevPointTemp2.length==1){
									RevPointY=Array.concat(RevPointY,RevPointTemp2[0]);
									RevPointX=Array.concat(RevPointX,RevPointTemp[0]);
									frame=Array.concat(frame, frameTemp[0]);
									i=frameTemp[0];
									limit=i;
								}
							}
						} //parenthesis for if(abs(dXprev)>=cminRV && abs(dXnext)>=cminRV)
					} //parenthesis for checking for reversals for that track 

					//Array.show(index,RevPointX,RevPointY);
					//waitForUser;
					
					// writing txt files to record the reversal points
					if (RevPointX.length>0){
						rpX=File.open(outputRevPoint+"RevPointX_"+name+".txt");
						for (j=0;j<RevPointX.length;j++){
							print(rpX,RevPointX[j]);
						}
						File.close(rpX);
	
						rpY=File.open(outputRevPoint+"RevPointY_"+name+".txt");
						for (j=0;j<RevPointY.length;j++){
							print(rpX,RevPointY[j]);
						}
						File.close(rpY);		
					}

					// recording switch frequency
					//if (length>minLength){
					SF=index.length/2;
					SFperSec=SF/length2;
					sf=File.open(outputSF+"SF_"+name+".txt");
					print(sf,SF);
					File.close(sf);	
					sf=File.open(outputSFperSec+"SFperSec_"+name+".txt");
					print(sf,SFperSec);
					File.close(sf);	
					//}
				
			}//parenthesis for mobile track: if (dXtrackCenter[i]>cmin)
				
	//------------------------------------------define anterograde and retrograde----------------------------------------------------
					
			if(startsWith(Roi.getName,"mobile")){
				dXmobile=parseFloat(xcoor[xcoor.length-1])-parseFloat(xcoor[0]);
				if (dXmobile>0){
					roiManager("Rename","anterograde_"+name);
					roiManager("Set Color", "green");
					roiManager("Set Line Width", 2);
				} 
				if (dXmobile<0) {
					roiManager("Rename","retrograde_"+name);
					roiManager("Set Color", "red");
					roiManager("Set Line Width", 2);		
				}

				if (dXmobile==0) {
					roiManager("Rename","stationary_"+name);
					roiManager("Set Color", "blue");
					roiManager("Set Line Width", 2);		
				}

				// recording switch frequency for non reversing tracks (SF=0)
				//if (length>minLength){
				sf=File.open(outputSF+"SF_"+name+".txt");
				print(sf,SF);
				File.close(sf);	
				sf=File.open(outputSFperSec+"SFperSec_"+name+".txt");
				print(sf,SFperSec);
				File.close(sf);	
				//}
			}//parenthesis for mobile tracks	
		} //parenthesis for xcoor/ one track for (i=0;i<xcoor.length-1;i++)
	}//parenthesis for finishing all tracks for (n=0;n<roiCount;n++) 
	setBatchMode(false);
//------------------------------------------ save different populations----------------------------------------------------
	
	// counting the number of anterograde, retrograde, reversing or stationary tracks
	// for tracks that have minimum length
	anterograde=newArray();
	retrograde=newArray();
	reversal=newArray();
	stationary=newArray();
	for (n=0;n<roiCount;n++){
		roiManager("Select",n);
		getSelectionCoordinates(x, y);
		//length=abs(y[y.length-1]-y[0]);
		//minLength=h-10;
		//if (length>minLength){
			if(startsWith(Roi.getName,"anterograde")){
				anterograde=Array.concat(anterograde,n);
			}
			if(startsWith(Roi.getName,"retrograde")){
				retrograde=Array.concat(retrograde,n);
			}
			if(startsWith(Roi.getName,"reversal")){
				reversal=Array.concat(reversal,n);
			}
			if(startsWith(Roi.getName,"track")){
				name=Roi.getName;
				stationary=Array.concat(stationary,n);
				roiManager("Rename","stationary_"+name);
				roiManager("Set Color", "blue");
				roiManager("Set Line Width", 2);
				
			}
		/*} else {
			if(startsWith(Roi.getName,"track")){
				name=Roi.getName;
				roiManager("Rename","stationary_"+name);
				roiManager("Set Color", "blue");
				roiManager("Set Line Width", 2);
			}
			// print tracks that are excluded and their length
			print("track excluded from cargo population analysis: ", Roi.getName, "\t (",length/Frame," sec)");
		}
		*/
	}

	// number of tracks in each category: numeric CP
	Ante=anterograde.length;
	Retro=retrograde.length;
	Reverse=reversal.length;
	Stationary=stationary.length;
	tracksNum=Ante+Retro+Reverse+Stationary;

	// CP as percentage
	AntePCT=anterograde.length/tracksNum;
	RetroPCT=retrograde.length/tracksNum;
	ReversePCT=reversal.length/tracksNum;
	StationaryPCT=stationary.length/tracksNum;

	// Calculating Density as tracks per um axon
	w=PxlSize*getWidth();
	AnteDensity=anterograde.length/w;
	RetroDensity=retrograde.length/w;
	ReverseDensity=reversal.length/w;
	StationaryDensity=stationary.length/w;

	// Calculating Flux as tracks per um axon and sec
	w=PxlSize*getWidth();
	h=getHeight()/Frame;
	AnteFlux=anterograde.length/(w*h);
	RetroFlux=retrograde.length/(w*h);
	ReverseFlux=reversal.length/(w*h);

	// printing into txt files
	CargoPop = newArray(Ante,Retro,Reverse,Stationary);
	CPNum=File.open(outputCPNum +"CP_Num_"+replace(file1,".tif",".txt"));
	for(i=0;i<CargoPop.length;i++){
		print(CPNum,CargoPop[i]);
	}
	File.close(CPNum);
	
	CargoPopulationPCT = newArray(AntePCT,RetroPCT,ReversePCT,StationaryPCT);
	CPPCT=File.open(outputCPNum +"CP_PCT_"+replace(file1,".tif",".txt"));
	for(i=0;i<CargoPopulationPCT.length;i++){
		print(CPPCT,CargoPopulationPCT[i]);
	}
	File.close(CPPCT);

	Density = newArray(AnteDensity,RetroDensity,ReverseDensity,StationaryDensity);
	density=File.open(outputDensity +"Density_"+replace(file1,".tif",".txt"));
	for(i=0;i<Density.length;i++){
		print(density,Density[i]);
	}
	File.close(density);

	Flux = newArray(AnteFlux,RetroFlux,ReverseFlux);
	flux=File.open(outputFlux +"Flux_"+replace(file1,".tif",".txt"));
	for(i=0;i<Flux.length;i++){
		print(flux,Flux[i]);
	}
	File.close(flux);
	
	// saving assigned tracks CPROIs
	count=roiManager("count"); 
	array=newArray(count); 
	for(i=0; i<count;i++) { 
		array[i] = i; 
	} 
		
	roiManager("Select", array); 
	roiManager("Save", outputCP_ROI + "CP_"+replace(getInfo("image.filename"), ".tif",".zip"));
	roiManager("Reset");
	
	/*
	// saving info about excluded tracks
	if (getInfo("Log")==""){			
	} else {
		selectWindow("Log");
		save(outputInfo+"ExcludedTracks_CP_"+movie+".txt");
		run("Close");	// Closes Log Window
	}
	*/
	run("Close");
} // parenthesis function CargoPopulation(input,file1,file2)
run ("Close All");


