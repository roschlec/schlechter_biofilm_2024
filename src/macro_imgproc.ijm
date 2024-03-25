//	Process biofilm images
	
// 	Set global measurements
run("Fresh Start");
run("Set Measurements...", "area min redirect=None decimal=3");
setForegroundColor(0, 0, 0);
setBackgroundColor(0, 0, 0);

//	Set the directory of the images and an output directory
#@ File (label="Select input directory", style="directory") dir
#@ File (label= "Select ouput directory", style = "directory") out


list = getFileList(dir); // Get a list of all the files that will be analysed

setBatchMode(true);
for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], '.czi')){
		if (matches(list[i], '.*z-stack.*')) {
			print("Analysing image " + list[i]);
			s = "open=["+dir+"/"+list[i]+"]";
			run("Bio-Formats Windowless Importer", s);
			if(nSlices>3){
				run("Z Project...", "projection=[Max Intensity]");
			}
		}else{
			print("Analysing image " + list[i]);
			s = "open=["+dir+"/"+list[i]+"]";
			run("Bio-Formats Windowless Importer", s); 
			}
			
		titleID = getTitle();
		atitle = replace(titleID, ".czi", "");
		run("Split Channels");

		// Channel C1: Sytox Green
		selectWindow("C1-" + titleID);
		c1 = getTitle();
		c1id = getImageID();
		run("Subtract Background...", "rolling=50");
		run("Enhance Contrast...", "saturated=0.1");
			
		// Channel C2: Red
		selectWindow("C2-" + titleID);
		c2 = getTitle();
		c2id = getImageID();
		run("Subtract Background...", "rolling=50");
		run("Enhance Contrast...", "saturated=0.1");

		//	Channel C3: Brightfield 
		if ("C3-" + titleID == 1) {
			selectWindow("C3-" + titleID);
			c3 = getTitle();
			c3id = getImageID();
			selectImage(c3id);
			saveAs("Tiff", out+"/"+atitle+"_C3.tif");
			run("RGB Color");
			saveAs("PNG", out+"/"+atitle+"_C3.png");
		}

		//	Merge channels C1 + C2
		run("Merge Channels...", "c6=" + c2 + " c7=" + c1 + " create keep");
		selectWindow("Composite");
		compid = getImageID();

		//	Saving files
		selectImage(c1id);
		saveAs("Tiff", out+"/"+atitle+"_C1.tif");
		run("RGB Color");
		saveAs("PNG", out+"/"+atitle+"_C1.png");

		selectImage(c2id);
		saveAs("Tiff", out+"/"+atitle+"_C2.tif");
		run("RGB Color");
		saveAs("PNG", out+"/"+atitle+"_C2.png");

		selectImage(compid);
		saveAs("Tiff", out+"/"+atitle+"_C1+C2.tif");
		saveAs("PNG", out+"/"+atitle+"_C1+C2.png");
		}
		while (nImages()>0) {
          selectImage(nImages());  
          run("Close");
    }
}