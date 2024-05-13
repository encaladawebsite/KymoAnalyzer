//This is a macro that determines the segments of all tracks saved for a movie

// Choose the directory to save the track coordinates and ROI managers;
print("\\Clear");
roiManager("Reset");

// define variabe used in both functions
var x=newArray();
var y=newArray();
var Xsegment=newArray();
var Ysegment=newArray();
var index=newArray();
var trackName

//input pixel size, frame rate and experiment folder
Dialog.create("ImagingParameters");
Dialog.addNumber("Pixel Size in µm",0.16);
Dialog.addNumber("Frame Rate in /sec", 10);
Dialog.show;
PxlSize=Dialog.getNumber();
Frame=Dialog.getNumber();
Factor=10/Frame;

input = getDirectory("Choose Experiment folder"); 
ExperimentList=getFileList(input);
File.makeDirectory(input+"PooledData");
File.makeDirectory(input+"PooledData/DataPerKymograph");

//----------------------------------makes a output folder to save coordinates of segments--------------------------------------------------------------------
for (i=0; i<ExperimentList.length;i++){
	
	if (startsWith(ExperimentList[i], "PooledData")){
			
		} else {
		File.makeDirectory(input+ExperimentList[i]+"Segment_Xcoor");
		File.makeDirectory(input+ExperimentList[i]+"Segment_Ycoor");
		File.makeDirectory(input+ExperimentList[i]+"Segment_ROIs");
		File.makeDirectory(input+ExperimentList[i]+"Segmental_Velocities");
		File.makeDirectory(input+ExperimentList[i]+"PD");
		File.makeDirectory(input+ExperimentList[i]+"PF");
		File.makeDirectory(input+ExperimentList[i]+"PFperSec");
		File.makeDirectory(input+ExperimentList[i]+"splitPD");
		File.makeDirectory(input+ExperimentList[i]+"splitPF");
		File.makeDirectory(input+ExperimentList[i]+"splitPFperSec");
		File.makeDirectory(input+ExperimentList[i]+"RL");
		File.makeDirectory(input+ExperimentList[i]+"combinedRL");
		File.makeDirectory(input+ExperimentList[i]+"combinedSV");
		File.makeDirectory(input+ExperimentList[i]+"Info");
		File.makeDirectory(input+ExperimentList[i]+"DataPerKymograph");
		File.makeDirectory(input+ExperimentList[i]+"PM");
	}
}

//----------------------------------processes all tracks for one movie to assign Segments--------------------------------------------------------------------
processFolder1(input);
function processFolder1(input) {
	for (i = 0; i < ExperimentList.length; i++) {
		if (startsWith(ExperimentList[i], "PooledData")){
			noMovies=ExperimentList.length-1;
		} else {
			// output files for assigning segments
			outputSegmentX=input+ExperimentList[i]+"Segment_Xcoor/";
			outputSegmentY=input+ExperimentList[i]+"Segment_Ycoor/";
			outputSegment_ROI=input+ExperimentList[i]+"Segment_ROIs/";
			outputInfo=input+ExperimentList[i]+"Info/";

			//output files for parameters
			noMovies=ExperimentList.length;
			outputSVs=input+ExperimentList[i]+"Segmental_Velocities/";
			outputNV=input+ExperimentList[i]+"Net_Velocities/";
			outputPD=input+ExperimentList[i]+"PD/";
			outputPF=input+ExperimentList[i]+"PF/";
			outputPFperSec=input+ExperimentList[i]+"PFperSec/";
			outputsplitPD=input+ExperimentList[i]+"splitPD/";
			outputsplitPF=input+ExperimentList[i]+"splitPF/";
			outputsplitPFperSec=input+ExperimentList[i]+"splitPFperSec/";
			outputRL=input+ExperimentList[i]+"RL/";
			outputcomRL=input+ExperimentList[i]+"combinedRL/";
			outputcomSV=input+ExperimentList[i]+"combinedSV/";
			outputPM=input+ExperimentList[i]+"PM/";

			// open kymograph and ROI Manager
			inputKymo=input+ExperimentList[i]+"Kymograph/";
			KymoList=getFileList(inputKymo);
			open(input+ExperimentList[i]+"Kymograph/"+KymoList[0]);
			name=File.nameWithoutExtension;
			print(name);
			inputROIs=input+ExperimentList[i]+"CP_ROIs/";
			CP_ROIList=getFileList(input+ExperimentList[i]+"CP_ROIs/");
			roiManager("Open",input+ExperimentList[i]+"CP_ROIs/"+CP_ROIList[0]);
			
			count=roiManager("Count");
			trackCount=count;
			Index=newArray();
			Name=newArray();
			Ante=newArray();
			Retro=newArray();
			Reversal=newArray();
			Stationary=newArray();
			for(k=0;k<count;k++){
				roiManager("Select",k);			
				if(startsWith(Roi.getName,"anterograde")){
				Index=Array.concat(Index,k);	
				Name=Array.concat(Name,Roi.getName);
				}	
				if(startsWith(Roi.getName,"retrograde")){
				Index=Array.concat(Index,k);	
				Name=Array.concat(Name,Roi.getName);
				}
				if(startsWith(Roi.getName,"reversal")){
				Index=Array.concat(Index,k);
				Name=Array.concat(Name,Roi.getName);
				}	
				if(startsWith(Roi.getName,"stationary")){	
				}		
			}
											
			for (j=0; j<Index.length; j++){
				AssignSegments(Index[j]);
				SegmentalVelocities(Name[j]);
			}

			//saving ROI Manager
			count=roiManager("Count");
			array=newArray(count);
			for (j=0;j<count;j++){
				array[j]=j;
			}
			roiManager("Select",array);
			roiManager("Save", outputSegment_ROI + "SegmentROI_"+replace(getInfo("image.filename"), ".tif",".zip"));
			
			if (getInfo("Log")==""){
				
			} else {
			selectWindow("Log");
			save(outputInfo+"CombinedPauses_"+replace(ExperimentList[i],"/",".txt"));
			run("Close");	// Closes Log Window
			}
		}
		roiManager("Reset");
		
	}
}

//----------------------------------processes all tracks for one movie to calculate segmental velocities and net velocities--------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//----------------------------------FUNCTION TO ASSIGN SEGMENTS-----------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//**********************************new function**************************************************************************************************
function AssignSegments(track){
	roiManager("Select",track);
	trackName=substring(Roi.getName,lengthOf(Roi.getName)-8,lengthOf(Roi.getName));
	
	// xTemp and yTemp are arrays with the original xy coordinates of the clicked polyline
	getSelectionCoordinates(xTemp,yTemp);
	
	// rewriting coordinates if track was clicked from the bottom
	// xTemp2 and yTemp2 is new array of xy coordinates that start at the beginning of the movie
	// x and y are the same arrays as xTemp2 and yTemp2 for use in the second function
	j=yTemp.length;
	xTemp2=newArray();
	yTemp2=newArray();
	x=newArray();
	y=newArray();
	if (yTemp[0] > yTemp[yTemp.length-1]){
		for (i=0;i<yTemp.length;i++){
			j=j-1;
			xTemp2=Array.concat(xTemp2,xTemp[j]);
			yTemp2=Array.concat(yTemp2,yTemp[j]);	
			x=Array.concat(x,xTemp[j]);
			y=Array.concat(y,yTemp[j]);	
		}
	} else {
		for (i=0;i<yTemp.length;i++){
			xTemp2=Array.concat(xTemp2,xTemp[i]);
			yTemp2=Array.concat(yTemp2,yTemp[i]);
			x=Array.concat(x,xTemp[i]);
			y=Array.concat(y,yTemp[i]);
		}
	}

	// calculating velocities of segments as defined by polyline to find pauses
	SVTemp=newArray();
	for (i=0;i<xTemp2.length-1;i++){
		distance=PxlSize*(xTemp2[i+1]-xTemp2[i]);
		time=(yTemp2[i+1]-yTemp2[i]+1)/Frame;
		SV=distance/time;
		SVTemp=Array.concat(SVTemp,SV);	
	}
	
	// finding pauses that are neighboring
	index=newArray();
	for (i=0;i<SVTemp.length-1;i++){
		// removing neighboring pauses
		if (abs(SVTemp[i])<0.1/Factor && abs(SVTemp[i+1])<0.1/Factor){
				index=Array.concat(index,i+1);
		}
	}

	for (i=0;i<index.length;i++){
		xTemp2[index[i]]="NaN";
		yTemp2[index[i]]="NaN";
	}

	//combining neighboring pauses
	// xTemp3 and yTemp3 are arrays with the xy coordinates after pauses are combined
	xTemp3=newArray();
	yTemp3=newArray();
	for (i=0;i<xTemp2.length;i++){
		if (startsWith(xTemp2[i],"NaN")){			
		} else {
			xTemp3=Array.concat(xTemp3,xTemp2[i]);	
		}
		if (startsWith(yTemp2[i],"NaN")){		
		} else {
			yTemp3=Array.concat(yTemp3,yTemp2[i]);	
		}
	}

	// printing info for combined pauses
	if (index.length>0){	
	print("Number of neighboring pauses combined: ",trackName,":", index.length+1);
	}

	//calculating segmental velocities to combine segments with similar velocity
	SVTemp=newArray();
	for (i=0;i<xTemp3.length-1;i++){
		distance=PxlSize*(xTemp3[i+1]-xTemp3[i]);
		time=(yTemp3[i+1]-yTemp3[i]+1)/Frame;
		SV=distance/time;
		SVTemp=Array.concat(SVTemp,SV);	
	}
	//Array.show(xTemp3,yTemp3,SVTemp);
	//waitForUser;
	
	// finding segments with similar velocity
	// XsegmentTemp and YsegmentTemp are arrays with xy coordinates of defined segments 
	index=newArray();
	
	ThreshFast=0.05/Factor;
	ThreshSlow=0.2/Factor;
	InfoSeg=File.open(outputInfo+"SegmentThresholds.txt");
	print(InfoSeg, "Threshold for SV>= 3 µm/sec: " + ThreshFast);
	print(InfoSeg, "Threshold for SV< 3 µm/sec: " + ThreshSlow);
	File.close(InfoSeg);
	for (i=0;i<SVTemp.length-1;i++){
		
		if (SVTemp[i]>=3/Factor){
			a=SVTemp[i]-ThreshFast*(abs(SVTemp[i]));			//this is a threshold in which the slope can change, should be refined based on data
			b=SVTemp[i]+ThreshFast*(abs(SVTemp[i]));			//this is a threshold in which the slope can change, should be refined based on data
			
			if (SVTemp[i]>=0 && a<0){
				a=0;
			}
			if (SVTemp[i]<=0 && b>0){
				b=0;
			}
		} else {
			a=SVTemp[i]-ThreshSlow*(abs(SVTemp[i]));			//this is a threshold in which the slope can change, should be refined based on data
			b=SVTemp[i]+ThreshSlow*(abs(SVTemp[i]));			//this is a threshold in which the slope can change, should be refined based on data
				
			if (SVTemp[i]>=0 && a<0){
				a=0;
			}
			if (SVTemp[i]<=0 && b>0){
				b=0;
			}
		}

		if (SVTemp[i+1]>a && SVTemp[i+1]<b){
		} else {
			index=Array.concat(index,i+1);
		}
	}
	
	XsegmentTemp=newArray();
	YsegmentTemp=newArray();
	XsegmentTemp=Array.concat(XsegmentTemp,xTemp3[0]);
	YsegmentTemp=Array.concat(YsegmentTemp,yTemp3[0]);
	for (i=0;i<index.length;i++){
		XsegmentTemp=Array.concat(XsegmentTemp,xTemp3[index[i]]);
		YsegmentTemp=Array.concat(YsegmentTemp,yTemp3[index[i]]);
	}
	XsegmentTemp=Array.concat(XsegmentTemp,xTemp3[xTemp3.length-1]);
	YsegmentTemp=Array.concat(YsegmentTemp,yTemp3[yTemp3.length-1]);

	// repeat segment definition to combine segments with similar "average" velocity
	SVTemp=newArray();
	for (i=0;i<XsegmentTemp.length-1;i++){
		distance=PxlSize*(XsegmentTemp[i+1]-XsegmentTemp[i]);
		time=(YsegmentTemp[i+1]-YsegmentTemp[i]+1)/Frame;
		SV=distance/time;
		SVTemp=Array.concat(SVTemp,SV);	
	}
		
	// finding segments with similar velocity
	// Xsegment and Ysegment are arrays with xy coordinates of defined segments 
	index=newArray();
	for (i=0;i<SVTemp.length-1;i++){
		
		if (SVTemp[i]>=3/Factor){
			a=SVTemp[i]-ThreshFast*(abs(SVTemp[i]));			//this is a threshold in which the slope can change, should be refined based on data
			b=SVTemp[i]+ThreshFast*(abs(SVTemp[i]));			//this is a threshold in which the slope can change, should be refined based on data
			
			if (SVTemp[i]>=0 && a<0){
				a=0;
			}
			if (SVTemp[i]<=0 && b>0){
				b=0;
			}
		} else {
			a=SVTemp[i]-ThreshSlow*(abs(SVTemp[i]));			//this is a threshold in which the slope can change, should be refined based on data
			b=SVTemp[i]+ThreshSlow*(abs(SVTemp[i]));			//this is a threshold in which the slope can change, should be refined based on data
				
			if (SVTemp[i]>=0 && a<0){
				a=0;
			}
			if (SVTemp[i]<=0 && b>0){
				b=0;
			}
		}

		if (SVTemp[i+1]>a && SVTemp[i+1]<b){
		} else {
			index=Array.concat(index,i+1);
		}
	}

	Xsegment=newArray();
	Ysegment=newArray();
	Xsegment=Array.concat(Xsegment,XsegmentTemp[0]);
	Ysegment=Array.concat(Ysegment,YsegmentTemp[0]);
	for (i=0;i<index.length;i++){
		Xsegment=Array.concat(Xsegment,XsegmentTemp[index[i]]);
		Ysegment=Array.concat(Ysegment,YsegmentTemp[index[i]]);
	}
	Xsegment=Array.concat(Xsegment,XsegmentTemp[XsegmentTemp.length-1]);
	Ysegment=Array.concat(Ysegment,YsegmentTemp[YsegmentTemp.length-1]);

//----------------Writes Result Files--------------------------------------------------------
	Xseg=File.open(outputSegmentX+"SegmentsX_"+trackName+".txt");
	for(i=0;i<Xsegment.length;i++){
		print(Xseg,Xsegment[i]);
	}
	File.close(Xseg);
	
	Yseg=File.open(outputSegmentY+"SegmentsY_"+trackName+".txt");
	for(i=0;i<Ysegment.length;i++){
		print(Yseg,Ysegment[i]);
	}
	File.close(Yseg);
	
//-----------------Shows Assigned Segments----------------------------------------------------
	makeSelection("polyline", Xsegment,Ysegment);
	roiManager("Add");
	trackCount=roiManager("Count");
	roiManager("Select",trackCount-1);
	roiManager("Rename","Segment_"+trackName);
	roiManager("Set Color","magenta");
	roiManager("Set Line Width", 1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//----------------------------------FUNCTION TO CALCULATE SEGMENTAL VELOCITIES--------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function SegmentalVelocities(name){
	//waitForUser;
	//add distance for each polyline segment for one segment
	distance=newArray();
	for (i=0;i<Xsegment.length-1;i++){
		dX=PxlSize*(Xsegment[i+1]-Xsegment[i]);
		distance=Array.concat(distance,dX);	
	}

	// calculate the duration of one segment
	// has to be dY+1 (frame). dY is differences between time stamps of when the capture started, however, 
	// the duration of the movie is the difference of when the first picture was taken until the last picture 
	// was taken plus the duration of the last frame. 

	time=newArray();
	for (i=0;i<Ysegment.length-1;i++){
		dY=(Ysegment[i+1]-Ysegment[i]+1)/Frame;
		time=Array.concat(time,dY);
	}

	SegVel=newArray();
	for (i=0;i<time.length;i++){
		SegVel=Array.concat(SegVel,distance[i]/time[i]);
	}
	
	//calculating segmental run length
	RL=newArray();
	for (i=0;i<SegVel.length;i++){
		if (abs(SegVel[i])>0.1/Factor){
			RL=Array.concat(RL,distance[i]);
		}
	}

	// calculating %time in motion
	// distance... array with RL
	// time... array with number of frames+1/framerate
	// for %time in motion, find all ante/retro segments, add number of frames
	sumAnte=newArray();
	sumRetro=newArray();
	sumPause=newArray();
	sumTotal=newArray();
	for (i=0; i<distance.length;i++){
		sumTotal=Array.concat(sumTotal, time[i]);
		if (abs(distance[i]/time[i])>0.1/Factor){
			if (distance[i]>0){
				sumAnte=Array.concat(sumAnte,time[i]);
			}
			if (distance[i]<0){
				sumRetro=Array.concat(sumRetro,time[i]);
			}
		}
		if (abs(distance[i]/time[i])<=0.1/Factor){
			sumPause=Array.concat(sumPause, time[i]);
		}
	}
	
	anteTime=0;
	retroTime=0;
	pauseTime=0;
	totalTime=0;
	for (i=0; i<sumAnte.length; i++){
		anteTime=anteTime+sumAnte[i];
	}
	for (i=0; i<sumRetro.length; i++){
		retroTime=retroTime+sumRetro[i];
	}
	for (i=0; i<sumPause.length; i++){
		pauseTime=pauseTime+sumPause[i];
	}
	for (i=0; i<sumTotal.length; i++){
		totalTime=totalTime+sumTotal[i];
	}

	anteMotion=100*anteTime/totalTime;
	retroMotion=100*retroTime/totalTime;
	pauseMotion=100*pauseTime/totalTime;
		
	// calculating combined run length
	//tempRL is an array that also contains run length of pauses in order to calculate combined RL
	tempRL=newArray();
	for (i=0;i<SegVel.length;i++){	
		tempRL=Array.concat(tempRL,distance[i]);
	}

	tempFrames=newArray();
	for (i=0;i<Ysegment.length-1;i++){
		dY=Ysegment[i+1]-Ysegment[i];
		tempFrames=Array.concat(tempFrames,dY);
	}
	
	//calculating combined segmental run length
	comRL=newArray();
	comTime=newArray;
	combinedRL=0;
	combinedFrames=0;
	for (i=0;i<SegVel.length;i++){
		done=false;
		while (SegVel[i]>0.1/Factor && done!=true){
			combinedRL=combinedRL+tempRL[i];
			combinedFrames=combinedFrames+tempFrames[i];
			
			if (i==SegVel.length-1){
				done=true;
				i--;
			}
			i++;
			
			if(SegVel[i]<=0.1/Factor){
				i--;
				done=true;
			}
		}
		if (combinedRL!=0){
			comRL=Array.concat(comRL,combinedRL);
			combinedTime=(combinedFrames+1)/Frame;
			comTime=Array.concat(comTime,combinedTime);
		}
		combinedRL=0;
		combinedFrames=0;
		
		while (SegVel[i]<-0.1/Factor && done!=true){
			combinedRL=combinedRL+tempRL[i];
			combinedFrames=combinedFrames+tempFrames[i];
			if (i==SegVel.length-1){
				done=true;
				i--;
			}
			i++;
			if(SegVel[i]>=-0.1/Factor){
				i--;
				done=true;
			}
		}
		if (combinedRL!=0){
			comRL=Array.concat(comRL,combinedRL);
		}
		if (combinedFrames!=0){
			combinedTime=(combinedFrames+1)/Frame;
			comTime=Array.concat(comTime,combinedTime);
		}
		combinedRL=0;
		combinedFrames=0;
	}

	comSV=newArray();
	for (i=0;i<comRL.length;i++){
		comSV=Array.concat(comSV,comRL[i]/comTime[i]);
	}

	segVel=File.open(outputSVs+"SV_"+trackName+".txt");
	for(i=0;i<SegVel.length;i++){
		if (abs(SegVel[i])>0.1/Factor){
			print(segVel,SegVel[i]);
		}
	}
	File.close(segVel);
	
	segRL=File.open(outputRL+"RL_"+trackName+".txt");
	for(i=0;i<RL.length;i++){
		print(segRL,RL[i]);
	}
	File.close(segRL);

	percentMotion=File.open(outputPM+"PM_"+trackName+".txt");
	print(percentMotion,anteMotion);
	print(percentMotion,retroMotion);
	print(percentMotion,pauseMotion);
	File.close(percentMotion);

	combinedRL=File.open(outputcomRL+"comRL_"+trackName+".txt");
	for(i=0;i<comRL.length;i++){
		print(combinedRL,comRL[i]);
	}
	File.close(combinedRL);

	combinedSV=File.open(outputcomSV+"comSV_"+trackName+".txt");
	for(i=0;i<comSV.length;i++){
		print(combinedSV,comSV[i]);
	}
	File.close(combinedSV);

	PD=newArray();							//a pause is defined as a segment with a velocity < 0.1 µm/sec
	for (i=0; i<SegVel.length;i++){
		if (abs(SegVel[i])<=0.1/Factor){
			PD=Array.concat(PD,time[i]);			//pause duration in sec
		}
	}
	pauses=File.open(outputPD+"PD_"+trackName+".txt");
	for(i=0;i<PD.length;i++){
		print(pauses,PD[i]);
	}
	File.close(pauses);

	//h=getHeight();
	//length=y[y.length-1]-y[0];
	//minLength=h-10;
	//if (length>minLength){
	duration=(y[y.length-1]-y[0]+1)/Frame;
	PFperSec=PD.length/duration;	
	pauses=File.open(outputPF+"PF_"+trackName+".txt");
	print(pauses,PD.length);
	File.close(pauses);

	pauses=File.open(outputPFperSec+"PFperSec_"+trackName+".txt");
	print(pauses,PFperSec);
	File.close(pauses);
	//}

	//split pauses into anterograde and retrograde pauses
	if (startsWith(name,"anterograde")){
		PD=newArray();							//a pause is defined as a segment with a velocity < 0.1 µm/sec
		for (i=0; i<SegVel.length;i++){
			if (abs(SegVel[i])<=0.1/Factor){
				PD=Array.concat(PD,time[i]);			//pause duration in sec
			}
		}
		pauses=File.open(outputsplitPD+"antePD_"+trackName+".txt");
		for(i=0;i<PD.length;i++){
			print(pauses,PD[i]);
		}
		File.close(pauses);
	
		//h=getHeight();
		//length=y[y.length-1]-y[0];
		//minLength=h-10;
		//if (length>minLength){
		duration=(y[y.length-1]-y[0]+1)/Frame;
		PFperSec=PD.length/duration;	
		pauses=File.open(outputsplitPF+"antePF_"+trackName+".txt");
		print(pauses,PD.length);
		File.close(pauses);
	
		pauses=File.open(outputsplitPFperSec+"antePFperSec_"+trackName+".txt");
		print(pauses,PFperSec);
		File.close(pauses);
		//}
	}

	if (startsWith(name,"retrograde")){
		PD=newArray();							//a pause is defined as a segment with a velocity < 0.1 µm/sec
		for (i=0; i<SegVel.length;i++){
			if (abs(SegVel[i])<=0.1/Factor){
				PD=Array.concat(PD,time[i]);			//pause duration in sec
			}
		}
		pauses=File.open(outputsplitPD+"retroPD_"+trackName+".txt");
		for(i=0;i<PD.length;i++){
			print(pauses,PD[i]);
		}
		File.close(pauses);
	
		//h=getHeight();
		//length=y[y.length-1]-y[0];
		//minLength=h-10;
		//if (length>minLength){
		duration=(y[y.length-1]-y[0]+1)/Frame;
		PFperSec=PD.length/duration;	
		pauses=File.open(outputsplitPF+"retroPF_"+trackName+".txt");
		print(pauses,PD.length);
		File.close(pauses);
	
		pauses=File.open(outputsplitPFperSec+"retroPFperSec_"+trackName+".txt");
		print(pauses,PFperSec);
		File.close(pauses);
		//}
	}

	if (startsWith(name,"reversal")){
		//waitForUser;
		antePD=newArray();
		retroPD=newArray();								//a pause is defined as a segment with a velocity < 0.1 µm/sec
		revPD=newArray();
		if (abs(SegVel[0])<=0.1/Factor && SegVel[1]>0.1/Factor){
			antePD=Array.concat(antePD,time[0]);					//pause duration in sec
		}

		//a pause is defined as a segment with a velocity < 0.1 µm/sec
		if (abs(SegVel[0])<=0.1/Factor && SegVel[1]<-0.1/Factor){
			retroPD=Array.concat(retroPD,time[0]);					//pause duration in sec
		}
		
		for (i=1; i<SegVel.length-1;i++){
			if (abs(SegVel[i])<=0.1/Factor){
				if (SegVel[i-1]>0.1/Factor && SegVel[i+1]>0.1/Factor){
					antePD=Array.concat(antePD,time[i]);			//pause duration in sec
				}
						
				if (SegVel[i-1]<-0.1/Factor && SegVel[i+1]<-0.1/Factor){
					retroPD=Array.concat(retroPD,time[i]);			//pause duration in sec
				}
				if (SegVel[i-1]<-0.1/Factor && SegVel[i+1]>0.1/Factor){
					revPD=Array.concat(revPD,time[i]);			//pause duration in sec
				}
				if (SegVel[i-1]>0.1/Factor && SegVel[i+1]<-0.1/Factor){
					revPD=Array.concat(revPD,time[i]);			//pause duration in sec
				}
			}
		}

		if (abs(SegVel[SegVel.length-1])<=0.1/Factor && SegVel[SegVel.length-2]>0.1/Factor){
			antePD=Array.concat(antePD,time[time.length-1]);					//pause duration in sec
		}

		//a pause is defined as a segment with a velocity < 0.1 µm/sec
		if (abs(SegVel[SegVel.length-1])<=0.1/Factor && SegVel[SegVel.length-2]<-0.1/Factor){
			retroPD=Array.concat(retroPD,time[time.length-1]);					//pause duration in sec
		}

		pauses=File.open(outputsplitPD+"antePD_"+trackName+".txt");
		for(i=0;i<antePD.length;i++){
			print(pauses,antePD[i]);
		}
		File.close(pauses);
		pauses=File.open(outputsplitPD+"retroPD_"+trackName+".txt");
		for(i=0;i<retroPD.length;i++){
			print(pauses,retroPD[i]);
		}
		File.close(pauses);
		pauses=File.open(outputsplitPD+"revPD_"+trackName+".txt");
		for(i=0;i<revPD.length;i++){
			print(pauses,revPD[i]);
		}
		File.close(pauses);
		//h=getHeight();
		//length=y[y.length-1]-y[0];
		//minLength=h-10;
		//if (length>minLength){
		duration=(y[y.length-1]-y[0]+1)/Frame;
		antePFperSec=antePD.length/duration;	
		retroPFperSec=retroPD.length/duration;	
		revPFperSec=revPD.length/duration;	
		pauses=File.open(outputsplitPF+"antePF_"+trackName+".txt");
		print(pauses,antePD.length);
		File.close(pauses);
	
		pauses=File.open(outputsplitPFperSec+"antePFperSec_"+trackName+".txt");
		print(pauses,antePFperSec);
		File.close(pauses);
	
		pauses=File.open(outputsplitPF+"retroPF_"+trackName+".txt");
		print(pauses,retroPD.length);
		File.close(pauses);
	
		pauses=File.open(outputsplitPFperSec+"retroPFperSec_"+trackName+".txt");
		print(pauses,retroPFperSec);
		File.close(pauses);

		pauses=File.open(outputsplitPF+"revPF_"+trackName+".txt");
		print(pauses,revPD.length);
		File.close(pauses);
	
		pauses=File.open(outputsplitPFperSec+"revPFperSec_"+trackName+".txt");
		print(pauses,revPFperSec);
		File.close(pauses);
		//}
	}
	
}// parenthesis function SegmentalVelocities

run("Close All");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Pool Data per Kymograph															//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
for (i=0;i<ExperimentList.length;i++){	
	if (startsWith(ExperimentList[i],"PooledData")){	
	} else {
		name=replace(ExperimentList[i],"Analysis_","");
		name=replace(name,"/","");
					
		ante=newArray();
		retro=newArray();
		reversal=newArray();
		stationary=newArray();

		//read in files CP_Num
		inputCP=input+ExperimentList[i]+"CargoPopulation/";
		CPList=getFileList(inputCP);
		string = File.openAsString(inputCP+CPList[0]); 
		lines = split(string, "\n");   
		CPNum=newArray(); 
		for (j=0; j<lines.length; j++){	
			CPNum = Array.concat(CPNum,parseFloat(lines[j]));
		}
		
		ante=Array.concat(ante,abs(CPNum[0]));
		retro=Array.concat(retro,abs(CPNum[1]));
		reversal=Array.concat(reversal,abs(CPNum[2]));
		stationary=Array.concat(stationary,abs(CPNum[3]));

		//read in files CP_PCT
		inputCP=input+ExperimentList[i]+"CargoPopulation/";
		CPList=getFileList(inputCP);
		string = File.openAsString(inputCP+CPList[1]); 
		lines = split(string, "\n");   
		CPPCT=newArray(); 
		for (j=0; j<lines.length; j++){	
			CPPCT = Array.concat(CPPCT,parseFloat(lines[j]));
		}
		
		ante=Array.concat(ante,abs(CPPCT[0]));
		retro=Array.concat(retro,abs(CPPCT[1]));
		reversal=Array.concat(reversal,abs(CPPCT[2]));
		stationary=Array.concat(stationary,abs(CPPCT[3]));

		//read in files NCP_Num
		inputCP=input+ExperimentList[i]+"CargoPopulation/";
		CPList=getFileList(inputCP);
		string = File.openAsString(inputCP+CPList[2]); 
		lines = split(string, "\n");   
		CPNum=newArray(); 
		for (j=0; j<lines.length; j++){	
			CPNum = Array.concat(CPNum,parseFloat(lines[j]));
		}
		
		ante=Array.concat(ante,abs(CPNum[0]));
		retro=Array.concat(retro,abs(CPNum[1]));
		reversal=Array.concat(reversal,NaN);
		stationary=Array.concat(stationary,abs(CPNum[2]));

		//read in files NCP_PCT
		inputCP=input+ExperimentList[i]+"CargoPopulation/";
		CPList=getFileList(inputCP);
		string = File.openAsString(inputCP+CPList[3]); 
		lines = split(string, "\n");   
		CPPCT=newArray(); 
		for (j=0; j<lines.length; j++){	
			CPPCT = Array.concat(CPPCT,parseFloat(lines[j]));
		}
		
		ante=Array.concat(ante,abs(CPPCT[0]));
		retro=Array.concat(retro,abs(CPPCT[1]));
		reversal=Array.concat(reversal,NaN);
		stationary=Array.concat(stationary,abs(CPPCT[2]));

		//read in files Flux
		inputFlux=input+ExperimentList[i]+"Flux/";
		FluxList=getFileList(inputFlux);
		string = File.openAsString(inputFlux+FluxList[0]); 
		lines = split(string, "\n");   
		Flux=newArray(); 
		for (j=0; j<lines.length; j++){	
			Flux = Array.concat(Flux,parseFloat(lines[j]));
		}
		
		ante=Array.concat(ante,abs(Flux[0]));
		retro=Array.concat(retro,abs(Flux[1]));
		reversal=Array.concat(reversal,abs(Flux[2]));
		stationary=Array.concat(stationary,NaN);

		//read in files Density
		inputDensity=input+ExperimentList[i]+"Density/";
		DensityList=getFileList(inputDensity);
		string = File.openAsString(inputDensity+DensityList[0]); 
		lines = split(string, "\n");   
		Density=newArray(); 
		for (j=0; j<lines.length; j++){	
			Density = Array.concat(Density,parseFloat(lines[j]));
		}
		
		ante=Array.concat(ante,abs(Density[0]));
		retro=Array.concat(retro,abs(Density[1]));
		reversal=Array.concat(reversal,abs(Density[2]));
		stationary=Array.concat(stationary,abs(Density[3]));

		//read in files PM
		inputPM=input+ExperimentList[i]+"PM/";
		PMList=getFileList(inputPM);
		PM_ante=newArray(); 
		PM_retro=newArray(); 
		PM_pauses=newArray(); 
		for (k=0; k<PMList.length; k++){
			string = File.openAsString(inputPM+PMList[k]); 
			lines = split(string, "\n");   
			PM=newArray();
			for (j=0; j<lines.length; j++){	
				PM = Array.concat(PM,parseFloat(lines[j]));
			}
			PM_ante=Array.concat(PM_ante,abs(PM[0])); 
			PM_retro=Array.concat(PM_retro,abs(PM[1])); 
			PM_pauses=Array.concat(PM_pauses,abs(PM[2])); 
		}
		
		//read in NV
		inputNetVel=input+ExperimentList[i]+"Net_Velocities/";
		NetVelList=getFileList(inputNetVel);
		anteNV=newArray();
		retroNV=newArray();
		NV_ante=newArray();
		NV_retro=newArray();
		NV=newArray();
		for(k=0;k<NetVelList.length;k++){
		string = File.openAsString(inputNetVel+NetVelList[k]); 
		lines = split(string, "\n");   
		NVTemp=newArray(); 
		for (j=0; j<lines.length; j++){	
			NVTemp = Array.concat(NVTemp,parseFloat(lines[j]));
		}
		NV=Array.concat(NV,NVTemp);
		if (NVTemp[0]>0)
		anteNV=Array.concat(anteNV,abs(NVTemp[0]));
		if (NVTemp[0]<0)
		retroNV=Array.concat(retroNV,abs(NVTemp[0]));
		}
		NV_ante=Array.concat(NV_ante,anteNV);
		NV_retro=Array.concat(NV_retro,retroNV);

		//read in files SV
		inputSegVel=input+ExperimentList[i]+"Segmental_Velocities/";
		SegVelList=getFileList(inputSegVel);
		anteSV=newArray();
		retroSV=newArray();
		SV_ante=newArray();
		SV_retro=newArray();
		SV=newArray();
		
		for(k=0;k<SegVelList.length;k++){
			string = File.openAsString(inputSegVel+SegVelList[k]); 
			lines = split(string, "\n");   
			SVTemp=newArray(); 
			
			for (j=0; j<lines.length; j++){	
				SVTemp = Array.concat(SVTemp,parseFloat(lines[j]));
			}
			
			for (j=0;j<SVTemp.length;j++){
				if (abs(SVTemp[j])>0.1){
				SV=Array.concat(SV,SVTemp[j]);
				}
				if (SVTemp[j]>0.1){
				anteSV=Array.concat(anteSV,abs(SVTemp[j]));
				}
				if (SVTemp[j]<-0.1){
				retroSV=Array.concat(retroSV,abs(SVTemp[j]));
				}
			}
		}
		SV_ante=Array.concat(SV_ante,anteSV);
		SV_retro=Array.concat(SV_retro,retroSV);
		

		//read in files RL
		inputRL=input+ExperimentList[i]+"RL/";
		RLList=getFileList(inputRL);
		anteRL=newArray();
		retroRL=newArray();
		RL_ante=newArray();
		RL_retro=newArray();
		RL=newArray();
		
		for(k=0;k<RLList.length;k++){
			string = File.openAsString(inputRL+RLList[k]); 
			lines = split(string, "\n");   
			RLTemp=newArray(); 
			
			for (j=0; j<lines.length; j++){	
				RLTemp = Array.concat(RLTemp,parseFloat(lines[j]));
			}
			RL=Array.concat(RL,RLTemp);
			
			for (j=0;j<RLTemp.length;j++){
				if (RLTemp[j]>0){
				anteRL=Array.concat(anteRL,abs(RLTemp[j]));
				}
				if (RLTemp[j]<0){
				retroRL=Array.concat(retroRL,abs(RLTemp[j]));
				}
			}
		}
		RL_ante=Array.concat(RL_ante,anteRL);
		RL_retro=Array.concat(RL_retro,retroRL);

		//read in files cSV
		inputSegVel=input+ExperimentList[i]+"combinedSV/";
		SegVelList=getFileList(inputSegVel);
		anteSV=newArray();
		retroSV=newArray();
		cSV_ante=newArray();
		cSV_retro=newArray();
		cSV=newArray();
		for(k=0;k<SegVelList.length;k++){
			string = File.openAsString(inputSegVel+SegVelList[k]); 
			lines = split(string, "\n");   
			SVTemp=newArray(); 
			
			for (j=0; j<lines.length; j++){	
				SVTemp = Array.concat(SVTemp,parseFloat(lines[j]));
			}
	
			for (j=0;j<SVTemp.length;j++){
				if (abs(SVTemp[j])>0.1){
					cSV=Array.concat(cSV,SVTemp[j]);
				}
				if (SVTemp[j]>0.1){
					anteSV=Array.concat(anteSV,abs(SVTemp[j]));
				}
				if (SVTemp[j]<-0.1){
					retroSV=Array.concat(retroSV,abs(SVTemp[j]));
				}
			}
		}
		cSV_ante=Array.concat(cSV_ante,anteSV);
		cSV_retro=Array.concat(cSV_retro,retroSV);

		// read in files cRL
		inputRL=input+ExperimentList[i]+"combinedRL/";
		RLList=getFileList(inputRL);
		anteRL=newArray();
		retroRL=newArray();
		cRL_ante=newArray();
		cRL_retro=newArray();
		cRL=newArray();
		for(k=0;k<RLList.length;k++){
			string = File.openAsString(inputRL+RLList[k]); 
			lines = split(string, "\n");   
			RLTemp=newArray(); 
						
			for (j=0; j<lines.length; j++){	
				RLTemp = Array.concat(RLTemp,parseFloat(lines[j]));
			}
			cRL=Array.concat(cRL,RLTemp);
			for (j=0;j<RLTemp.length;j++){
				if (RLTemp[j]>0){
					anteRL=Array.concat(anteRL,abs(RLTemp[j]));
				}
				if (RLTemp[j]<0){
					retroRL=Array.concat(retroRL,abs(RLTemp[j]));
				}
			} 
		}
		cRL_ante=Array.concat(cRL_ante,anteRL);
		cRL_retro=Array.concat(cRL_retro,retroRL);

		//read in files PD
		inputPD=input+ExperimentList[i]+"PD/";
		listPD=getFileList(inputPD);
		PD=newArray;
		for(k=0;k<listPD.length;k++){
			string = File.openAsString(inputPD+listPD[k]); 
			lines = split(string, "\n");   
			PDTemp=newArray(); 
			for (j=0; j<lines.length; j++){	
				PDTemp = Array.concat(PDTemp,parseFloat(lines[j]));
			}
			PD=Array.concat(PD,PDTemp);
		}

		// read in files splitPD
		inputPD=input+ExperimentList[i]+"splitPD/";
		listPD=getFileList(inputPD);
		sPD_ante=newArray();
		sPD_retro=newArray();
		sPD_rev=newArray();
		for(k=0;k<listPD.length;k++){
			if (startsWith(listPD[k],"ante")){
				string = File.openAsString(inputPD+listPD[k]); 
				lines = split(string, "\n");   
				PDTemp=newArray(); 
				for (j=0; j<lines.length; j++){	
					PDTemp = Array.concat(PDTemp,parseFloat(lines[j]));
				}	
			}
			sPD_ante=Array.concat(sPD_ante,PDTemp);
			PDTemp=newArray();
			
			if (startsWith(listPD[k],"retro")){
				string = File.openAsString(inputPD+listPD[k]); 
				lines = split(string, "\n");   
				PDTemp=newArray(); 
				for (j=0; j<lines.length; j++){	
					PDTemp = Array.concat(PDTemp,parseFloat(lines[j]));
				}	
			}
			sPD_retro=Array.concat(sPD_retro,PDTemp);
			PDTemp=newArray();
			
			if (startsWith(listPD[k],"rev")){
				string = File.openAsString(inputPD+listPD[k]); 
				lines = split(string, "\n");   
				PDTemp=newArray(); 
				for (j=0; j<lines.length; j++){	
					PDTemp = Array.concat(PDTemp,parseFloat(lines[j]));
				}	
			}
			sPD_rev=Array.concat(sPD_rev,PDTemp);
			PDTemp=newArray();
		}	

		//read in files PF
		inputPF=input+ExperimentList[i]+"PF/";
		listPF=getFileList(inputPF);
		PF=newArray;
		for(k=0;k<listPF.length;k++){
			string = File.openAsString(inputPF+listPF[k]); 
			lines = split(string, "\n");   
			PFTemp=newArray(); 
			for (j=0; j<lines.length; j++){	
				PFTemp = Array.concat(PFTemp,parseFloat(lines[j]));
			}
			PF=Array.concat(PF,PFTemp);
		}
	
		// read in files splitPF
		inputPF=input+ExperimentList[i]+"splitPF/";
		listPF=getFileList(inputPF);
		sPF_ante=newArray();
		sPF_retro=newArray();
		sPF_rev=newArray();
		for(k=0;k<listPF.length;k++){
			if (startsWith(listPF[k],"ante")){
				string = File.openAsString(inputPF+listPF[k]); 
				lines = split(string, "\n");   
				PFTemp=newArray(); 
				for (j=0; j<lines.length; j++){	
					PFTemp = Array.concat(PFTemp,parseFloat(lines[j]));
				}	
			}
			sPF_ante=Array.concat(sPF_ante,PFTemp);
			PFTemp=newArray();
			
			if (startsWith(listPF[k],"retro")){
				string = File.openAsString(inputPF+listPF[k]); 
				lines = split(string, "\n");   
				PFTemp=newArray(); 
				for (j=0; j<lines.length; j++){	
					PFTemp = Array.concat(PFTemp,parseFloat(lines[j]));
				}	
			}
			sPF_retro=Array.concat(sPF_retro,PFTemp);
			PFTemp=newArray();
			
			if (startsWith(listPF[k],"rev")){
				string = File.openAsString(inputPF+listPF[k]); 
				lines = split(string, "\n");   
				PFTemp=newArray(); 
				for (j=0; j<lines.length; j++){	
					PFTemp = Array.concat(PFTemp,parseFloat(lines[j]));
				}	
			}
			sPF_rev=Array.concat(sPF_rev,PFTemp);
			PFTemp=newArray();
		}	

		//read in files PF per Sec
		inputPF=input+ExperimentList[i]+"PFperSec/";
		listPF=getFileList(inputPF);
		PFperSec=newArray;
		for(k=0;k<listPF.length;k++){
			string = File.openAsString(inputPF+listPF[k]); 
			lines = split(string, "\n");   
			PFTemp=newArray(); 
			for (j=0; j<lines.length; j++){	
				PFTemp = Array.concat(PFTemp,parseFloat(lines[j]));
			}
			PFperSec=Array.concat(PFperSec,PFTemp);
		}

		// read in files splitPFperSec
		inputPF=input+ExperimentList[i]+"splitPFperSec/";
		listPF=getFileList(inputPF);
		sPFperSec_ante=newArray();
		sPFperSec_retro=newArray();
		sPFperSec_rev=newArray();
		for(k=0;k<listPF.length;k++){
			if (startsWith(listPF[k],"ante")){
				string = File.openAsString(inputPF+listPF[k]); 
				lines = split(string, "\n");   
				PFTemp=newArray(); 
				for (j=0; j<lines.length; j++){	
					PFTemp = Array.concat(PFTemp,parseFloat(lines[j]));
				}	
			}
			sPFperSec_ante=Array.concat(sPFperSec_ante,PFTemp);
			PFTemp=newArray();
			
			if (startsWith(listPF[k],"retro")){
				string = File.openAsString(inputPF+listPF[k]); 
				lines = split(string, "\n");   
				PFTemp=newArray(); 
				for (j=0; j<lines.length; j++){	
					PFTemp = Array.concat(PFTemp,parseFloat(lines[j]));
				}	
			}
			sPFperSec_retro=Array.concat(sPFperSec_retro,PFTemp);
			PFTemp=newArray();
			
			if (startsWith(listPF[k],"rev")){
				string = File.openAsString(inputPF+listPF[k]); 
				lines = split(string, "\n");   
				PFTemp=newArray(); 
				for (j=0; j<lines.length; j++){	
					PFTemp = Array.concat(PFTemp,parseFloat(lines[j]));
				}	
			}
			sPFperSec_rev=Array.concat(sPFperSec_rev,PFTemp);
			PFTemp=newArray();
		}

		//read in files SF
		inputSF=input+ExperimentList[i]+"SF/";
		listSF=getFileList(inputSF);
		SF=newArray;
		SF_rev=newArray();
		for(k=0;k<listSF.length;k++){
			string = File.openAsString(inputSF+listSF[k]); 
			lines = split(string, "\n");   
			SFTemp=newArray(); 
			revSF=newArray();
			for (j=0; j<lines.length; j++){	
				SFTemp = Array.concat(SFTemp,parseFloat(lines[j]));
				if (parseFloat(lines[j])!=0){
					revSF = Array.concat(revSF,parseFloat(lines[j]));
				}
			}
			SF=Array.concat(SF,SFTemp);
			SF_rev=Array.concat(SF_rev,revSF);
		}
		
		//read in files SF per Sec
		inputSF=input+ExperimentList[i]+"SFperSec/";
		listSF=getFileList(inputSF);
		SFperSec=newArray;
		SFperSec_rev=newArray();
		for(k=0;k<listSF.length;k++){
			string = File.openAsString(inputSF+listSF[k]); 
			lines = split(string, "\n");   
			SFTemp=newArray(); 
			revSF=newArray();
			for (j=0; j<lines.length; j++){	
				SFTemp = Array.concat(SFTemp,parseFloat(lines[j]));
				if (parseFloat(lines[j])!=0){
					revSF = Array.concat(revSF,parseFloat(lines[j]));
				}				
			}
			SFperSec=Array.concat(SFperSec,SFTemp);
			SFperSec_rev=Array.concat(SFperSec_rev,revSF);			
		}

		blank=newArray();			
		CargoParameters=newArray("CP_Num","CP_PCT","NCP_Num","NCP_PCT","Flux","Density");
		Array.show("ResultsperKymograph", CargoParameters, ante, retro, reversal, stationary, NV, NV_ante, NV_retro,SV, SV_ante,SV_retro, RL, RL_ante, RL_retro, cSV, cSV_ante, cSV_retro,cRL, cRL_ante, cRL_retro, PD, sPD_ante, sPD_retro, sPD_rev, PF, sPF_ante, sPF_retro, sPF_rev, PFperSec, sPFperSec_ante, sPFperSec_retro, sPFperSec_rev, SF, SFperSec, SF_rev, SFperSec_rev, PM_ante, PM_retro, PM_pauses);
		selectWindow("ResultsperKymograph");
		saveAs(input+ExperimentList[i]+"DataPerKymograph/"+name+"_ResultsperKymograph.txt");
		saveAs(input+"PooledData/DataPerKymograph/"+name+"_ResultsperKymograph.txt");
		run("Close"); 
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Pool Track Coordinates and Segment Coordinates														//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
setBatchMode(true);
newImage("Untitled", "8-bit black", 1024, 1024, 1);
for (i=0;i<ExperimentList.length;i++){	
	if (startsWith(ExperimentList[i], "PooledData")){
	} else {
		print("\\Clear");
		name=replace(ExperimentList[i],"Analysis_","");
		name=replace(name,"/","");
		roiManager("reset");
		ROIList=getFileList(input+ExperimentList[i]+"ROIs/");
		roiManager("Open",input+ExperimentList[i]+"ROIs/"+ROIList[0]);
		Segment_ROIList=getFileList(input+ExperimentList[i]+"Segment_ROIs/");
		roiManager("Open",input+ExperimentList[i]+"Segment_ROIs/"+Segment_ROIList[0]);
		
		//print track coordinates
		string="";
		for (j=0; j<roiManager("count"); j++){
			roiManager("Select",j);
			if (startsWith(Roi.getName,"track")){
				string=string+Roi.getName+"_XCoor"+"\t";
				string=string+Roi.getName+"_YCoor"+"\t";
			}
		}
		headers=string;
		print(headers, "\n");
		// getting maximum length of Roi: maximum of n
		for (j=0; j<roiManager("count"); j++){
			roiManager("Select",j);
			Roi.getCoordinates(xpoints,ypoints);
			n=Array.concat(n,lengthOf(xpoints));
		}
		Array.sort(n);
		maximum=n[n.length-1]-1;
		// printing row by row of all columns
		for (j=0; j<maximum; j++){
			string="";
			for (k=0; k<roiManager("count"); k++){
				roiManager("Select",k);
				if (startsWith(Roi.getName,"track")){
					Roi.getCoordinates(xpoints,ypoints);
					done=false;
					while (xpoints.length-1>=j && !done){
						string=string+xpoints[j]+"\t"+ypoints[j]+"\t";
						done=true;
					}
					if (j>xpoints.length-1) {
						string=string+"\t"+"\t";
					}
				}
			}
			string=string+"\n";
			print(string);
		}
		selectWindow("Log");
		save(input+ExperimentList[i]+"DataPerKymograph/"+name+"_TrackCoordinatesPerKymograph.txt");
		save(input+"PooledData/DataPerKymograph/"+name+"_TrackCoordinatesPerKymograph.txt");
	
		// print Segments coordinates
		print("\\Clear");
		string="";
		for (j=0; j<roiManager("count"); j++){
			roiManager("Select",j);
			if (startsWith(Roi.getName,"Segment")){
				string=string+Roi.getName+"_XCoor"+"\t";
				string=string+Roi.getName+"_YCoor"+"\t";
			}
		}
		headers=string;
		print(headers, "\n");
		// getting maximum length of Roi: maximum of n
		for (j=0; j<roiManager("count"); j++){
			roiManager("Select",j);
			Roi.getCoordinates(xpoints,ypoints);
			n=Array.concat(n,lengthOf(xpoints));
		}
		Array.sort(n);
		maximum=n[n.length-1]-1;
		// printing row by row of all columns
		for (j=0; j<maximum; j++){
			string="";
			for (k=0; k<roiManager("count"); k++){
				roiManager("Select",k);
				if (startsWith(Roi.getName,"Segment")){
					Roi.getCoordinates(xpoints,ypoints);
					done=false;
					while (xpoints.length-1>=j && !done){
						string=string+xpoints[j]+"\t"+ypoints[j]+"\t";
						done=true;
					}
					if (j>xpoints.length-1) {
						string=string+"\t"+"\t";
					}
				}
			}
			string=string+"\n";
			print(string);
		}
		selectWindow("Log");
		save(input+ExperimentList[i]+"DataPerKymograph/"+name+"_SegmentCoordinatesPerKymograph.txt");
		save(input+"PooledData/DataPerKymograph/"+name+"_SegmentCoordinatesPerKymograph.txt");
	}
}
selectWindow("Log"); 
run("Close"); 
close();




	

