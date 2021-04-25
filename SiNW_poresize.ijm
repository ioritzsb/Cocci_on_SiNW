//Ioritz Sorzabal 2020.10.10
//Obtain max. insc. circles in porous surfaces:
/*Macro sequence:
 * 	Open a user defined folder with SEM images to be analysed,
 * 	Create a binary mask by Default threshopld binarisation
 * 	Apply Median Filter and invert image
 * 	Obtain Max inscribed circles within pores https://github.com/BIOP/ijp-max-inscribed-circles
 * 	Measure Circles
 * 	Flatten and save images and calculated data
 * */
 
//SAVE DIRECTORIES
//--File folder
inputDir = getDirectory("Choose Directory");
fileList = getFileList(inputDir);


//--Set Batch mode ON
setBatchMode(true);

minArray = newArray(10,15,20,30);




//--Iterate through images
for (i = 0; i < fileList.length; i++) { 
	
	if (endsWith(fileList[i], ".jpg")){
		

				
		//Create directories one-by-one
		File.makeDirectory(inputDir +"results");
		
		for (k = 0; k < minArray.length; k++) {				

			file = inputDir+fileList[i]; 
	 		open(file);
	 		
		 	//--Remove .jpg from name 
		  	name_full = getTitle();
			dotIndex = indexOf(name_full, ".");
			name = substring(name_full, 0, dotIndex); 						
			rename(name);

			
			outputDir = inputDir +"results" + File.separator + d2s(minArray[k], 0);

			//Create directories one-by-one

			File.makeDirectory(outputDir);					
			outputDir = outputDir + File.separator;
			

			
			///// IMAGE PROCESSING//////

			//----Generate binary image
			setAutoThreshold("Default dark");
			setOption("BlackBackground", true);
			run("Convert to Mask", "method=Triangle background=Dark calculate black");
			//----Median filter
			run("Median...", "radius=2");
			//----Invert image
			run("Invert");

			//----Obtain Maximum Inscribed Circles
			run("Max Inscribed Circles", "minimum="+d2s(minArray[k], 0)+" minimum_0=0.50 closeness=5");

			//----Set measurements
			run("Set Scale...", "distance=250.3333 known=5 unit=µm");
			run("Set Measurements...", "area min centroid perimeter feret's redirect=None decimal=1");
			
			//----ROI Manager Measurements
			roiManager("Measure");

			//----Create Overlay
			run("Invert");
			roiManager("show all");
			run("Flatten");
			rename(name+"_"+k);

			//---Delete ROIs
			roiManager("delete");			

			//----Save image
			selectImage(name+"_"+k);				
			saveAs("tiff", outputDir + name+"_"+k);

			//----Save results
			selectWindow("Results");
			saveAs("Results",  outputDir + name+"_"+ k + ".tsv");

			//--Close All Images
			run("Fresh Start");

			}
		}
	}

