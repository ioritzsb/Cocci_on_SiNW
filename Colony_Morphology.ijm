//Ioritz Sorzabal 2021.10.10
//Get number of cells per colony, Feret's diameter and Circularity. 
//Obtain  cell connectivity matrix
//GPU accelerated image processing using CLIJ2
/*Macro sequence:
 * 	Open a user defined folder with labelled masks (obtained with Stardist 2D)
 * 	Get touch matrix
 * 	Get binary mask from labelled images
 * 	Get colonies based on cell to cell connectivity, apply morphological filters and gwet morphological data.
 * */

//Close all 
run("Fresh Start");

//SAVE DIRECTORIES
//--File folder
inputDir = getDirectory("Choose Directory");
fileList = getFileList(inputDir);

//Create output directory
outputDir = inputDir +"analysis";
File.makeDirectory(outputDir);

//--Set Batch mode ON
setBatchMode(true);

//Define arrays
setOption("ExpandableArrays", true);
title = newArray;
col_label = newArray;
n_cells = newArray;
circ = newArray;
max_feret_d = newArray;
insc_circ_r = newArray;

n=0;
//--Initialise GPU
run("CLIJ2 Macro Extensions", "cl_device=[]");
Ext.CLIJ2_clear();



//--Iterate through images
for (i = 0; i < fileList.length; i++) { 
	
	if (endsWith(fileList[i], ".tif")){			

		file = inputDir+fileList[i]; 
 		open(file);
 		
	 	//--Remove .tif from name 
	  	name_full = getTitle();
		dotIndex = indexOf(name_full, ".");
		name = substring(name_full, 0, dotIndex); 						
		rename(name);
		
		///// Image Analysis//////

		//Touch Matrix
		//--Generate touch matrix
		tm = name + "_touch";
		am = name + "_adjMatrix";
		
		Ext.CLIJ2_push(name);		
		Ext.CLIJ2_generateTouchMatrix(name, tm);
		Ext.CLIJ2_pull(tm);
		//--Save touch matrix
		selectImage(tm);
		saveAs("Text Image", outputDir + File.separator + tm);
		close("*");

		//Get colonies
		//--Clear data		
		Ext.CLIJ2_pull(name);
		Ext.CLIJ2_clear();
		//--Get colonies
		selectImage(name);
		run("Duplicate...", " ");
		colony = name+"_colony";
		rename(colony);		
		setThreshold(1.0000, 1000000000000000000000000000000.0000);
		setOption("BlackBackground", true);
		run("Convert to Mask");
		Ext.CLIJ2_push(colony);
		col_labels = colony + "_label";
		col_labels_copy = col_labels + "_copy";
		Ext.CLIJ2_connectedComponentsLabelingDiamond(colony, col_labels);

		//Copy and dilate colony labels
		col_labels_dil = col_labels + "_dil";
		col_labels_solid = col_labels + "_solid";
		Ext.CLIJ2_copy(col_labels, col_labels_copy);
		Ext.CLIJx_extendLabelsWithMaximumRadius(col_labels_copy, col_labels_dil, 2)
		Ext.CLIJx_morphoLibJFillHoles(col_labels_dil, col_labels_solid);
		
		
		Ext.CLIJ2_pull(col_labels);
		Ext.CLIJ2_pull(col_labels_solid);
		
		Ext.CLIJ2_clear();

		//Get colony label maximum
		selectWindow(col_labels);
		getStatistics(area, mean, min, max);
		Ext.CLIJ2_push(col_labels);

		//Get colony aprox. morphology
		selectWindow(col_labels_solid);
		rename(col_labels_dil);
		run("Set Scale...", "distance=199 known=10 unit=um");
		run("Analyze Regions", "circularity max._feret max._inscribed_disc");
		IJ.renameResults(col_labels_dil+"-Morphometry","Results_dil");
		

		for (j = 1; j < max+1; j++) {

			//Get colony mask
			col_mask = colony +"_"+j;			
			Ext.CLIJ2_labelToMask(col_labels, col_mask, j);
			Ext.CLIJ2_pull(col_mask);
			Ext.CLIJ2_release(col_mask);

			imageCalculator("Multiply create", name,col_mask);
			col_bac = "bac_" + col_mask;
			rename(col_bac);
			run("Analyze Regions", " ");
			IJ.renameResults(col_bac+"-Morphometry","Results");

			//Save data in arrays
			title[n] = name;
			col_label[n] = j;
			n_cells[n] = nResults;
			
			//Get colony morphology info
			IJ.renameResults("Results",col_bac+"-Morphometry");
			wait(1);
			IJ.renameResults("Results_dil","Results");
			wait(1);
			
			selectWindow("Results");
			insc_circ_r[n] = getResult("InscrDisc.Radius", j-1);
			circ[n] = getResult("Circularity", j-1);
			max_feret_d[n] = getResult("MaxFeretDiam", j-1);			

			IJ.renameResults("Results","Results_dil");
			IJ.renameResults(col_bac+"-Morphometry","Results");
			n++;

			//Close all but label image	
			selectWindow(name);
			close("\\Others");
			close("Results");		
			

		}

		close("Results_dil");		
	}

}

//Close all 
run("Fresh Start");

//Create and save new results table
for (i = 0; i < n_cells.length; i++) {
	
	setResult("Image", i, title[i]);
	setResult("Colony", i, col_label[i]);
	setResult("N_cells", i, n_cells[i]);
	setResult("Circularity", i, circ[i]);
	setResult("MaxFeretDiam", i, max_feret_d[i]);
	setResult("InscrCircR", i, insc_circ_r[i]);	
	
	}
updateResults();
saveAs("Results",  outputDir + File.separator + "Colony_cells.tsv");
//Close all 
run("Fresh Start");

