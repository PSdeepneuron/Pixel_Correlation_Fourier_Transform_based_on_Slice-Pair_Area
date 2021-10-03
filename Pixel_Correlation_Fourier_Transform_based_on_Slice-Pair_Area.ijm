//Set directory to save results in table to use for analysis
#@ File (style="directory") imageFolder;
dir = File.getDefaultDir;
dir = replace(dir,"\\","/");

//Get pixel dimensions
setTool("hand");
waitForUser("Get pixel dimensions","Click on the non-annotated input image to get the pixel dimensions\nthen press OK");
getPixelSize(unit, pixelWidth, pixelHeight);
getVoxelSize(width, height, depth, unit);
pixel_area = pixelWidth*pixelHeight;

waitForUser("Starting frame","Put the annotated image at the slice from where the process should start\nthen press OK");

setTool("multipoint");
waitForUser("Select object and background","click on what you want to measure the area and border length of in the annotated image, then on the background\nand then press OK");
getSelectionCoordinates(xpoints, ypoints);
foreground = getValue(xpoints[0], ypoints[0]);
background = getValue(xpoints[1], ypoints[1]);
setTool("hand");

//Get image dimensions
w = getWidth();
h = getHeight();

sn = getSliceNumber();

end_slice = nSlices;

title = getTitle();

run("Duplicate...", "duplicate");
rename("lag=0");
selectImage("lag=0"); 
run("Slice Remover", "first=end_slice last=end_slice increment=1");

selectImage(title);
run("Duplicate...", "duplicate");
rename("lag=1");
selectImage("lag=1"); 
run("Slice Remover", "first=1 last=1 increment=1");

imageCalculator("Average create stack", "lag=0","lag=1");
rename("lag=overlay");
close("lag=0");
close("lag=1"); 
selectImage("lag=overlay");

slice_array = newArray();
depth_array = newArray();
corr_array = newArray();

sl = sn;
lags = 0;
for (v=sn;v<nSlices+1;v++){
	//Dendrite pixel count
	overlap = 0;
	n = 0;
	//Get coordinates of dendrite pixels
	//Loop for every pixel in current slice
	for (x=0;x<w;x++){
		for (y=0;y<h;y++){
			if (getPixel(x,y) != background){
				n += 1;
			}
			if (getPixel(x,y) == foreground){
				overlap += 1;
			}
		}
	}
//Get mean correlation
corr = overlap/n;
slice_depth = (sl-1)*depth; 
print(corr);
slice_array = Array.concat(slice_array,sl);
depth_array = Array.concat(depth_array,slice_depth);
corr_array = Array.concat(corr_array,corr);
sl += 1;
lags += 1;
run("Next Slice [>]");
}

Array.print(corr_array);
Array.getStatistics(corr_array, min, max, mean, stdDev);
neg_corr_normalised_array = newArray();
for (n=0;n<nSlices-sn;n++) {
	neg_corr_normalised = -1 * (corr_array[n]-mean);
	neg_corr_normalised_array = Array.concat(neg_corr_normalised_array,neg_corr_normalised); 
}

//Loop for every periodicity
periodicity = 1;
periodicity_array = newArray();
periodicity_value = newArray();
//Iterates from a periodicity of >1 frames since 1 frame is supposed to represent an infintismal point (a periodicity of 1 frame is no periodicity and also gives infinite when log() function is applied to enhance the value of the periodicity)
for (f=2;f<nSlices-sn+1;f++) {
	periodicity += 1;
	periodicity_array = Array.concat(periodicity_array,periodicity);
	periodicity_sum = 0;
	occurance = 0;
	for (s=0;s<nSlices-sn;s++) {
	    periodicity_sum += neg_corr_normalised_array[s]*cos((s*2*PI/f)+(2*PI/1));
      //periodicity_sum += (neg_corr_normalised_array[s])*sqrt(pow(cos((s*2*PI/f)+(2*PI/1)),2) + pow(sin((s*2*PI/f)+(2*PI/1)),2));
		occurance += 1;
	}
	periodicity_value = Array.concat(periodicity_value,periodicity_sum);
}

waitForUser("Progress","Done");

Plot.create("Correlation", "slice", "correlation", periodicity_array,periodicity_value);

save_option = getBoolean("Want to save results?");
if (save_option == 1){
//Make a table containing the arrays
Table.create("Correlation_Based_On_First_Slice");
Table.setColumn("period in slices",periodicity_array);
Table.setColumn("deviation from mean correlation",periodicity_value);
Table.save(dir+"Dendrite_Cross-sectional_Area_Correlation_Fourier"+".csv");
}