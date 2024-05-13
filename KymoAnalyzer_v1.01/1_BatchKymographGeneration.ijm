// Sylvia Neumann, 04/14/14
//
// Macro template to make kymographs from multiple images in a folder 
//
// *The macro uses the make kymograph macro from: download newest version at www.embl.de/eamnet/html/kymograph.html. In order to use the macro
// the files: MultipleKymograph_.java, MultipleOverlay_.java, StackDifference_.java, WalkingAverage_.java have to be saved into the "Plugins" folder
//
// 1. Choose parent folder in which movies are stored.
// 2. Choose if movies are stored as .nd2 or .tif files
// 3. The macro will generate a subdirectory in parent folder for each movie named "Analysis_NameOfTheMovie"
// 4. Each movie file in the parent folder will be opened (one by one) and a polyline along the axon is clicked by the user, then press "OK". 
//    Be aware of the orientation of the axon!
// 5. A kymograph will be generated with line width "1" as default. If required chooosing the line width by the user can be implemented.
// 6. The kymograph is inverted and saved as a 8-bit TIFF into a subdirectory "Kymograph" for each movie with the name of the original movie. 
//    (in order to facilitate copy pasting into Illustrator). 


//----------------------------------Chooses parent folder ----------------------------------
input = getDirectory("Choose parent folder"); 
list = getFileList(input); 
Dialog.create("File type");
Dialog.addString("File suffix: ", ".nd2", 5);
Dialog.show();
suffix = Dialog.getString();
Dialog.create("");
state=getBoolean("Do you want to generate more than one kymographs for a single movie?");
if (state==1){
	MultipleKymographs();
} else {
	SingleKymograph();
}
function MultipleKymographs(){
	for (i=0; i<list.length; i++){
		if (endsWith (list[i], suffix)){
			open(input+list[i]);
			ID=getImageID();
			roiManager("Reset");
			run("Enhance Contrast", "Auto");
			setTool("polyline");
			getLine(x1, y1, x2, y2, lineWidth);
			NameKymo=newArray("_A","_B","_C","_D","_E","_F","_G","_H","_I","_J","_K","_L","_M","_N","_O","_P","_Q","_R","_S","_T","_U","_V","_W","_X","_Y","_Z");
			countName=0;
			state1 = 1;
			while (state1==1) {
				msg = "Select axon with polyLine tool then click \"OK\".";
				waitForUser(msg);
				roiManager("Reset");
				roiManager("Add");
				roiManager("Select", 0); 
				roiManager("Rename",replace(list[i],suffix,NameKymo[countName]));
				run("MultipleKymograph ", "linewidth=1");
				run("Invert"); 
				File.makeDirectory(input+"Analysis_"+replace(list[i],suffix,NameKymo[countName]));  
				File.makeDirectory(input+"Analysis_"+replace(list[i],suffix,NameKymo[countName])+"/Kymograph"); 
				File.makeDirectory(input+"Analysis_"+replace(list[i],suffix,NameKymo[countName])+"/KymographPolyline");  
				output=input+"Analysis_"+replace(list[i],suffix,NameKymo[countName])+"/Kymograph/";
				output2=input+"Analysis_"+replace(list[i],suffix,NameKymo[countName])+"/KymographPolyline/";
				selectWindow("Kymograph");
				run("8-bit");
				save(output + replace(list[i],suffix,NameKymo[countName])+".tif");
				roiManager("Select", 0); 
				roiManager("Save", output2 + replace(list[i],suffix,NameKymo[countName]+".zip"));
				selectWindow("Kymograph");
				close();
				state1 = getBoolean("Do you want to generate another kymograph from this movie?");
				if (state1==1){
					if (countName<NameKymo.length-1){
						countName=countName+1;
					} else {
						Dialog.create("");
						Dialog.addMessage("The maximum number of kymographs for this movie is exceeded");
						Dialog.show();
						state1=0;
					}
				}
			}
		selectImage(ID);
		close();	
		}
	}
}

function SingleKymograph(){
	for (i=0; i<list.length; i++){
		if (endsWith (list[i], suffix)){
			open(input+list[i]);
			ID=getImageID();
			roiManager("Reset");
			run("Enhance Contrast", "Auto");
			setTool("polyline");
			getLine(x1, y1, x2, y2, lineWidth);
			msg = "Select axon with polyLine tool then click \"OK\".";
			waitForUser(msg);
			roiManager("Reset");
			roiManager("Add");
			roiManager("Select", 0); 
			roiManager("Rename",replace(list[i],suffix,".zip"));
			run("MultipleKymograph ", "linewidth=1");
			run("Invert"); 
			File.makeDirectory(input+"Analysis_"+replace(list[i],suffix,""));  
			File.makeDirectory(input+"Analysis_"+replace(list[i],suffix,"")+"/Kymograph"); 
			File.makeDirectory(input+"Analysis_"+replace(list[i],suffix,"")+"/KymographPolyline");  
			output=input+"Analysis_"+replace(list[i],suffix,"")+"/Kymograph/";
			output2=input+"Analysis_"+replace(list[i],suffix,"")+"/KymographPolyline/";
			selectWindow("Kymograph");
			run("8-bit");
			save(output + replace(list[i],suffix,".tif"));
			roiManager("Select", 0); 
			roiManager("Save", output2 + replace(list[i],suffix,".zip"));
			selectWindow("Kymograph");
			close();
			selectImage(ID);
			close();
		}
	}
}

