/*
* For internal use in Aoki lab.
* 2018.11.14 Yohei Kondo
* This macro produces IMD-like visualization from
* ratio images and background-subtracted donor images.
* Note that the output is 8-bit RGB image, 
* and thus it does not have FRET ratio information.
*/
Dialog.create("IMD parameters");
Dialog.addNumber("Ratio Max", 1.6);
Dialog.addNumber("Ratio Min", 0.5);
Dialog.addNumber("Donor Max", 4000);
Dialog.addNumber("Donor Min",700);
Dialog.addString("Ratio stack name", "ratio");
Dialog.addString("Donor stack name", "CFP");
Dialog.show();
rmax = Dialog.getNumber();
rmin = Dialog.getNumber();
dmax = Dialog.getNumber();
dmin = Dialog.getNumber();
rname = Dialog.getString();
dname = Dialog.getString();
rrange = rmax - rmin
drange = dmax - dmin
run("Conversions...", " ");// disable intensity scaling which is undireable in this case
selectWindow(dname);
run("Duplicate...", "duplicate");
rename("mask");
run("32-bit");
run("Subtract...", "value=&dmin stack");
run("Divide...", "value=&drange stack");
run("Min...", "value=0 stack");
run("Max...", "value=1 stack");
selectWindow(rname);
run("Duplicate...", "duplicate");
rename("tempratio");
run("32-bit");
changeValues(NaN, NaN, 0);
run("Subtract...", "value=&rmin stack");
run("Divide...", "value=&rrange stack");
run("Min...", "value=0 stack");
run("Max...", "value=1 stack");
run("Multiply...", "value=255 stack");
run("8-bit");
//run("Enhance Contrast", "saturated=0.35");
//run("Duplicate...", "duplicate");
//selectWindow("tempratio");
run("physics");
run("RGB Color");
run("RGB Stack");
run("Split Channels");
selectWindow("C1-tempratio");
run("32-bit");
selectWindow("C2-tempratio");
run("32-bit");
selectWindow("C3-tempratio");
run("32-bit");
imageCalculator("Multiply stack", "C1-tempratio", "mask");
selectWindow("C1-tempratio");
rename("red");
run("8-bit");
imageCalculator("Multiply stack", "C2-tempratio", "mask");
selectWindow("C2-tempratio");
rename("green");
run("8-bit");
imageCalculator("Multiply stack", "C3-tempratio", "mask");
selectWindow("C3-tempratio");
rename("blue");
run("8-bit");
//run("Merge Channels...", "c1=red c2=green c3=blue keep");
run("Merge Channels...", "c1=red c2=green c3=blue");
run("Conversions...", "scale");
fn="IMD_"+rmax+"_"+rmin+"_"+dmax+"_"+dmin;
rename(fn);

