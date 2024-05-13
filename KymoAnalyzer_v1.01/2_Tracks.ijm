//This is a macro to select multiple tracks

// Choose the directory to save the track coordinates and ROI managers;

input = getDirectory("Choose movie folder (not Experiment Folder)"); 

//----------------------------------makes a output folder to save ROIs and CP_ROIs--------------------------------------------------------------------
File.makeDirectory(input+"ROIs"); 
File.makeDirectory(input+"Info"); 
outputInfo=input+"Info/";
inputKymograph=input+"Kymograph/";
processFolder(inputKymograph);
function processFolder(inputKymograph) {
	
	KymoList=getFileList(inputKymograph);
	for (i = 0; i < KymoList.length; i++) {
		if(endsWith(KymoList[i], ".tif"))
			outputROI=input+"ROIs/";
			processFile(inputKymograph, KymoList[i]);
	}
}

function processFile(inputKymograph,file){
//open image and ROIs that already have been saved;
roiManager("Reset");
open(inputKymograph+file);
setLocation(200,250);
run("In [+]");  // zooms in once
run("In [+]");  // zooms in a second time
run("Scale to Fit");
run ("Brightness/Contrast...");
//Dialog.create("Image Settings");
//Dialog.addNumber("pixel size in µm", 0.16);
//pixelSize=Dialog.getNumber();
//Dialog.addNumber("frame rate in /sec", 0.5);
//frameRate=Dialog.getNumber();
//Dialog.addNumber("line width",5);
//lineWidth=Dialog.getNumber();
//Dialog.show();
//print(pixelSize,frameRate,lineWidth);
prevROI=File.exists(outputROI + replace(getInfo("image.filename"), ".tif",".zip"));
if (prevROI==true){
	roiManager("Open",outputROI + replace(getInfo("image.filename"), ".tif",".zip"));
}

//click polylines along tracks;
run("Enhance Contrast", "Auto");
counterTrack = roiManager("count");
trackNo=0;
if (roiManager("Count")>0){
	roiManager("Select",roiManager("Count")-1)
	trackNo=substring(Roi.getName,lengthOf(Roi.getName)-3,lengthOf(Roi.getName));
	trackNo=parseFloat(trackNo);
}
state = getBoolean("Do you want to add a track?");
while (state==1) {
	trackNo = trackNo + 1;
	setTool("polyline");
	getLine(x1, y1, x2, y2, lineWidth);
	waitForUser("Draw a track and press OK");
	roiManager("add");
	roiManager("Select",roiManager("count")-1);

	if (trackNo<10){
	roiManager("Rename","track" + "00"+toString(trackNo));
	}
	if (trackNo>=10 && counterTrack<100){
	roiManager("Rename","track" + "0"+toString(trackNo));	
	}
	if (trackNo>=100){
	roiManager("Rename","track"+toString(trackNo));	
	}
	state = getBoolean("Do you want to add a track?");
}
if (state==0) {
	setOption("Show All",true);
	inspect=getBoolean("Do you want to remove tracks?");
	if (inspect==1){
		waitForUser("Select a track in ROI manager and press -- Delete -- press OK when done");
	}	
} 

count=roiManager("count"); 
removedTracks=newArray();
for (i=0;i<count-1;i++){
	roiManager("Select",i);
	getSelectionCoordinates(x, y);
	prevTrackX=newArray(x.length);
	prevTrackY=newArray(x.length);
	for (j=0;j<x.length;j++){
			prevTrackX[j]=x[j];
			prevTrackY[j]=y[j];
		}
	n=i+1;
	for (n=i+1;n<count;n++){
		roiManager("Select",n);
		getSelectionCoordinates(x, y);
		nextTrackX=newArray(x.length);
		nextTrackY=newArray(x.length);
		for (j=0;j<x.length;j++){
			nextTrackX[j]=x[j];
			nextTrackY[j]=y[j];
		}
		if (prevTrackX.length==nextTrackX.length){
			j=0;
			while (j<x.length-1 && prevTrackX[j]==nextTrackX[j] && prevTrackY[j]==nextTrackY[j]){
				j++;
			}
			if (j==x.length-1){
				removedTracks=Array.concat(removedTracks,n);
			}
		}
	}
}

if (removedTracks.length>0){
roiManager("Select",removedTracks);
roiManager("Delete");
}

print("Number of duplicate Tracks removed: ",removedTracks.length);
print("\n");
name=File.nameWithoutExtension;
fid=File.open(outputInfo+"CheckTracks_"+name+".txt");
index=newArray();
for (i=0;i<roiManager("Count");i++){
	roiManager("Select",i);
	getSelectionCoordinates(x,y);

	
	for (j=0;j<y.length-1;j++){
		if (y[j+1]==y[j]){
			index=Array.concat(index,Roi.getName);
			j=y.length-1;
		}
	}
	
}

if (index.length>0){
	print(fid,name);
	print(fid, "these tracks contain horizontal lines, please correct these tracks");
	for (j=0;j<index.length;j++){
		print(fid,index[j]);	
	}

	print(name);
	print("these tracks contain horizontal lines, please correct these tracks");
	for (j=0;j<index.length;j++){
		print(index[j]);	
	}
}

if (index.length==0){
	print(fid,"None of the tracks contain horizontal lines");
}
File.close(fid);

count=roiManager("count"); 
array=newArray(count); 
for(i=0; i<count;i++) { 
  array[i] = i; 
}
roiManager("Select", array); 
roiManager("Save", outputROI + replace(getInfo("image.filename"), ".tif",".zip"));
roiManager("Deselect");
setOption("Show All",true);
}
	


	
	
	
	
	