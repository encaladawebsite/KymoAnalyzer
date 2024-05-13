// this macro calculates the cargo population based on net velocities and the net velocities

//----------------------------Variables-------------------------------------------------------------------
// list: is the array that contains all movie folders of the experiment folder

//----------------------------BEGIN MACRO------------------------------------------------------------------
print("\\Clear");			
roiManager("Reset");
Dialog.create("ImagingParameters");
Dialog.addNumber("Pixel Size in µm",0.16);
Dialog.addNumber("Frame Rate in /sec", 10);
Dialog.show;
PxlSize=Dialog.getNumber();
Frame=Dialog.getNumber();
input = getDirectory("Choose Experiment folder (the folder that contains all your movies)"); 

//----------------------------process the input folder-----------------------------------------------------
processFolder(input);

function processFolder(input){
	list=getFileList(input);
		
	for (i=0;i<list.length;i++){
		if (startsWith(list[i], "PooledData")){
			
		} else {	
			movie=list[i];
			movie=substring(movie,0,lengthOf(movie)-1);
			print(movie);
			// make output files
			File.makeDirectory(input+list[i]+"CP_ROIs");
			File.makeDirectory(input+list[i]+"X_coor");
			File.makeDirectory(input+list[i]+"Y_coor");
			File.makeDirectory(input+list[i]+"CargoPopulation");
			File.makeDirectory(input+list[i]+"Net_Velocities");
			File.makeDirectory(input+list[i]+"Info");
			
			//set outputs
			outputCP_ROI=input+list[i]+"CP_ROIs/";
			outputX_coor=input+list[i]+"X_coor/";
			outputY_coor=input+list[i]+"Y_coor/";
			outputCPNum=input+list[i]+"CargoPopulation/";
			outputNV=input+list[i]+"Net_Velocities/";
			outputInfo=input+list[i]+"Info/";

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
			NetCargoPopulation(input, Kymolist[0], Roilist[0]);
		}
	}
}

//-------------FUNCTION TO CALCULATE THE CARGO POPULATION-----------------
function NetCargoPopulation(input,file1,file2){
	// resets ROI Manager
	roiManager("Reset");
	// open kymograph and saved ROIs (tracks)
	open(input+list[i]+"Kymograph/"+file1);
	roiManager("Open",input+list[i]+"ROIs/"+file2);
	setOption("Show All",true);
	roiCount=roiManager("count");
	// for each track calculate the XY coordinates with one decimal
	for (n=0;n<roiCount;n++){
		roiManager("Select",n);	
		name=Roi.getName;
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
		length=(abs(y[y.length-1]-y[0]))/Frame;
		h=getHeight;
		minLength=(h-10)/Frame;
		*/
		
		//calculate coordinates with decimal points
		xcoor=newArray();
		ycoor=newArray();
		for (i=0;i<x.length-1;i++){
			m=(y[i+1]-y[i])/(x[i+1]-x[i]);
			b=y[i]- (m*x[i]);
			Y=y[i]-0.1;
			for (j=y[i]*10;j<y[i+1]*10;j++){
				
				if (startsWith(m,"Infinity") || startsWith(m,"-Infinity")){
					Y=Y+0.1;
					xcoor=Array.concat(xcoor,x[i]);
					ycoor=Array.concat(ycoor,Y);	
				} else {
					Y=Y+0.1;
					X=(Y-b)/m;
					//print(X,Y);	
					xcoor=Array.concat(xcoor,X);
					ycoor=Array.concat(ycoor,Y);
				}
			}
		}
		xcoor=Array.concat(xcoor,x[x.length-1]);
		Array.getStatistics(xcoor, min, max);
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
		//Array.getStatistics(xcoor, min, max);
		dXtrackCenterMin=PxlSize*(abs(min-trackCenter));
		dXtrackCenterMax=PxlSize*(abs(max-trackCenter));
		dXtrackCenterFirst=PxlSize*abs(xcoor[0]-trackCenter);
		dXtrackCenterLast=PxlSize*abs(xcoor[xcoor.length-1]-trackCenter);
		//print(name,"\t", trackCenter,"\t",xcoor[0],"\t", xcoor[xcoor.length-1],"\t", dXtrackCenterFirst,"\t", dXtrackCenterLast);
		//test if stationary
		// cmin=3 pixels is 480 nm movement away froom track center (on Artemis) Reis et al 2012: 550 nm cutoff 
		cmin=0.35;
		if (dXtrackCenterMin>cmin || dXtrackCenterMax>cmin){

				i=xcoor.length-1;
				//dXmobile: net translocation from fisrt to last frame
				dXmobile=PxlSize*(parseFloat(xcoor[xcoor.length-1])-parseFloat(xcoor[0]));
				dYmobile=((parseFloat(ycoor[ycoor.length-1])-parseFloat(ycoor[0]))+1)/Frame;
				NV=dXmobile/dYmobile;
				
				//reversal tracks could be net stationary, therefore, set cutoff OF 500 nm total translocation from track center

				
				
				if (abs(dXmobile)>0.35){
					if (dXmobile>0){
						roiManager("Rename","anterograde_"+name);
						roiManager("Set Color", "green");
						roiManager("Set Line Width", 2);
					} else {
						roiManager("Rename","retrograde_"+name);
						roiManager("Set Color", "red");
						roiManager("Set Line Width", 2);		
					}
					
					//length=dYmobile*Frame;
					//minLength=h-10;
					//if (length>minLength){
					netVel=File.open(outputNV+"NV_"+name+".txt");
					print(netVel,NV);
					File.close(netVel);
					//}
				}
			}//parenthesis for mobile tracks
				
			//parenthesis for xcoor
		}//parenthesis for finishing one track
		
//------------------------------------------ save different populations----------------------------------------------------
		anterograde=newArray();
		retrograde=newArray();
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
			if(startsWith(Roi.getName,"track")){
				name=Roi.getName;
				stationary=Array.concat(stationary,n);
				roiManager("Rename","stationary_"+name);
				roiManager("Set Color", "blue");
				roiManager("Set Line Width", 2);
				
			}
			/*
			} else {
				if(startsWith(Roi.getName,"track")){
					name=Roi.getName;
					roiManager("Rename","stationary_"+name);
					roiManager("Set Color", "blue");
					roiManager("Set Line Width", 2);
				}
				print("track excluded from net cargo population analysis: ", Roi.getName, "\t (",length/Frame," sec)");
			}
			*/
		}
			
		Ante=anterograde.length;
		Retro=retrograde.length;
		Stationary=stationary.length;
		
		tracksNum=Ante+Retro+Stationary;
		
		AntePCT=anterograde.length/tracksNum;
		RetroPCT=retrograde.length/tracksNum;
		StationaryPCT=stationary.length/tracksNum;

		CargoPop = newArray(Ante,Retro,Stationary);
		//fycoor=File.open(outputY_coor+"Y_Coor_"+name+".txt");
		CPNum=File.open(outputCPNum +"Net_CP_Num_"+replace(file1,".tif",".txt"));
		for(i=0;i<CargoPop.length;i++){
			print(CPNum,CargoPop[i]);
		}
		File.close(CPNum);
		
		CargoPopulationPCT=newArray(AntePCT,RetroPCT,StationaryPCT);
		CPPCT=File.open(outputCPNum +"Net_CP_PCT_"+replace(file1,".tif",".txt"));
		for(i=0;i<CargoPopulationPCT.length;i++){
			print(CPPCT,CargoPopulationPCT[i]);
		}
		File.close(CPPCT);
		
		count=roiManager("count"); 
		array=newArray(count); 
		for(i=0; i<count;i++) { 
        	array[i] = i; 
	} 
roiManager("Select", array); 
roiManager("Save", outputCP_ROI + "NetCP_"+replace(getInfo("image.filename"), ".tif",".zip"));
roiManager("Reset");
/*
if (getInfo("Log")==""){			
} else {
selectWindow("Log");
save(outputInfo+"ExcludedTracks_NCP_"+movie+".txt");
run("Close");	// Closes Log Window
}
*/
run("Close");
}
run ("Close All");


