/*
This file is part of Entropy-based Super-resolution Imaging (ESI).

ESI is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

ESI is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ESI.  If not, see <http://www.gnu.org/licenses/>
*/
package de.bio_photonics.esi;

import ij.plugin.PlugIn;	    // ImageJ plugin interface
import ij.ImageStack;		    // ImageStack
import ij.ImagePlus;		    // ImagePlus to show images
import ij.IJ;			    // static ImageJ control/log functions

import ij.process.FloatProcessor;   // Image with float-valued pixels
import ij.process.ImageProcessor;   // General (abstract) image
import ij.io.Opener;		    // Image opener for command line mode
import ij.gui.GenericDialog;		// to get the parameters
import ij.plugin.filter.GaussianBlur;	// blur filter


import java.io.InputStream;	    // for command-line reading of images
import java.io.FileInputStream;	    // for command-line reading of images

/**
 * Implementation of the ESI image analysis algorithm.
 * <p> 
 * Please cite our publication: ACS Photonics, 
 * "Entropy-based Super-resolution Imaging (ESI): From Disorder to Fine Detail"
 * <p> 
 * ESI algorithm by:	    Idir Yahiatene <br>
 * JAVA implementation by:  Marcel Mueller <br>
 * <p>
 * ESI is licensed under GPLv2 or later (see source code)
 * */
public class ESI_Analysis implements PlugIn {

    /** Default parameters for the analysis, can be changed in GUI */
    int   nrResImages=100;	// Images in result stack
    int	  nrBins=100;		// #bins/states in entropy calculation
    int   esi_order=4;		// ESI order
    boolean doMultiCore=true;	// run multi-threaded 
    boolean normOutput=true;	// normalize the output images

    /** This runs all steps of the analysis, given an ImageStack as input. */
    void runAnalysis( ImageStack stck, float pxMin, float pxMax  ) {

	// helpers
	Timing t1 = new Timing();
	GaussianBlur gb = new GaussianBlur();
	final int imgPerResult  = stck.getSize() / nrResImages;

	IJ.log("ESI: input images per resulting images: "+imgPerResult);

	// stacks to store the results
	ImageStack rec      = 
	    new ImageStack(2*stck.getWidth(), 2*stck.getHeight());

	ImagePlus interm = null;

	// image to show the summed-up result
	FloatProcessor summedImg = 
	    new FloatProcessor( 2*stck.getWidth(), 2*stck.getHeight());
	ImagePlus summedImgP = 
	    new ImagePlus("summed ESI result", summedImg) ;
	summedImgP.show();

	// Array to hold one TraceStack per Image Substack
	TraceStack [] trStArray = new TraceStack[ nrResImages ];

	// loop over subimages
	for( int k=0; k < nrResImages ; k++) {
	    t1.start();

	    // generate and normalize the TraceStack
	    TraceStack trSt = new TraceStack( stck , k*imgPerResult , (k+1)*imgPerResult );
	    trSt.normalizeFrom( pxMin, pxMax );
	    trSt.createNormBinning( nrBins ); 

	    // run the analysis
	    FloatProcessor reconstruction;
	    if (doMultiCore)  {
		MTP mp = new MTP(this);
		reconstruction = mp.doESI(trSt, esi_order);
	    } else {
		reconstruction = doESI(trSt, esi_order);
	    }

	    gb.blurFloat( reconstruction , 0.8, 0.8, 0.01);
	    if (normOutput)
		normalize( reconstruction );

	    // add these to the result stacks
	    rec.addSlice(reconstruction);

	    // add the slice to the current sum
	    ESI_Analysis.addFP( reconstruction, summedImg);
	    summedImg.resetMinAndMax();
	    summedImgP.updateAndDraw();
	    
	    // present some intermediate result
	    if ( (rec.getSize()>1)||(nrResImages==1))
	    if ( interm == null) {
		    interm= new ImagePlus("ESI result", rec);
		    interm.show();
	    } else {
		    interm.setPosition( rec.getSize() );
	    }

	    t1.stop();
	    IJ.log("ESI: SubImg "+k+"/"+nrResImages+" took "+t1 );
	}

    }

    /** The entropy-based reconstruction. It needs a fully set-up
     *  TraceStack as input. 'lStart' and 'lEnd' are line ranges for
     *  a poor mans multi-threading, see the 'doESI' functions. */
    static void ESI_internal(
	FloatProcessor res, TraceStack pxl, int order, 
	int lStart, int lEnd ) {

	// loop to create the new image
	if (lStart == 0) lStart = 1;
	
	
	for(int y=lStart; y<lEnd; y++)		// loop lines (start to end for this thread)
	for(int x=1; x<pxl.getWidth()  -1; x++)	// loop pixel position in line     
	    for (int i=0; i<2;i++)		// loop res improvement offset in x,y
	    for (int j=0; j<2;j++) {
		
		// on existing pixel: replace by cross-correlation of the 4
		// next-neighbor pixels
		if ( (i==0) && (j==0)) {
		    float tmp =
			PixelTrace.H_cross2( 
			    pxl.get(x-1,y-1) , pxl.get(x+1,y+1),
			    order) +
			PixelTrace.H_cross2( 
			    pxl.get(x-1,y-1) , pxl.get(x+1,y+1),
			    order);
		    res.setf(x*2,y*2, tmp/2.f);
		}
		
		// new pixel, but only offset in x or y, not both:
		// cross correlation between neighbors
		else if (i+j<2) {
		    res.setf(x*2+i,y*2+j , 
			PixelTrace.H_cross2( 
			    pxl.get(x,y) , pxl.get(x+i,y+j),
			    order));
		}
		
		// new pixel, on diagonal:
		// averaged cross-correlation
		else {
		    float tmp =
			PixelTrace.H_cross2( 
			    pxl.get(x,y) , pxl.get(x+i,y+j),order)+
			PixelTrace.H_cross2( 
			    pxl.get(x+i,y) , pxl.get(x,y+j),order);
		    res.setf(x*2+1,y*2+1, tmp/2.f);
		}
	    }
    }



    // --------------------------------------------------
    // plugin 'run' and command line 'main'
    // --------------------------------------------------
    
    
    /** Main methods for test runs from comamnd line. */
    public static void main( String [] args ) throws Exception {

	if (args.length != 1) {
	    System.out.println("Provide a TIFF file");
	    System.exit(-1);
	}

	// open a tiff stack
	InputStream img = new FileInputStream(args[0]);
	Opener op = new Opener();
	ImagePlus ip = op.openTiff( img, args[0]  );
	if (ip == null) {
	    System.out.println("Failed to open: "+args[0]);
	    System.exit(-1);
	}
	ip.show();
	
	// run an ESI analysis on it
	ESI_Analysis es = new ESI_Analysis();
	//TraceStack ptr = new TraceStack( ip.getStack());
	es.run("");
    }

    /** The run method called by ImageJ / Fiji when
     *	starting the plugin. Displays the dialog, then starts
     *	the analysis. */
    @Override 
    public void run(String arg) {
	// get the active image plus instance
	ImagePlus aip = ij.WindowManager.getCurrentImage();
	if (aip == null) {
		IJ.showMessage("No active image stack selected");
		return;
	}

	// currently: check if these are 3 phases
	int numImages = aip.getStack().getSize();
	if (numImages <3) {
		IJ.showMessage("Stack should contain some images");
		return;
	}

	// check the image format
	if (( aip.getType() != ImagePlus.GRAY8 )  &&
		( aip.getType() != ImagePlus.GRAY16 ) &&
		( aip.getType() != ImagePlus.GRAY32 )    ) {
			IJ.showMessage("Stack should be a grayscale image");
			return ;
		}

	// get min and max value
	float [] tmp = ESI_Analysis.getMinMax( aip.getStack() );
	float pxMin = tmp[0], pxMax = tmp[1];

	// ask for the parameters
	GenericDialog gd = new GenericDialog("ESI Parameters");
	gd.addNumericField("#images in output", nrResImages, 0); 
	gd.addNumericField("#bins for entropy", nrBins, 0); 
	gd.addNumericField("Order", esi_order, 0); 
	gd.addNumericField("min px value", pxMin, 2); 
	gd.addNumericField("max px value", pxMax, 2);
	gd.addCheckbox("Multicore",doMultiCore);
	gd.addCheckbox("normalize output",normOutput);
	gd.showDialog();
	if (gd.wasCanceled()) return;

	// store parameters
	nrResImages   =   (int)gd.getNextNumber();
	nrBins	      =   (int)gd.getNextNumber();
	esi_order     =   (int)gd.getNextNumber();
	pxMin	      =   (float)gd.getNextNumber();
	pxMax	      =   (float)gd.getNextNumber();
	doMultiCore   =	  gd.getNextBoolean();
	normOutput    =	  gd.getNextBoolean();


	// start the analysis
	Timing tfull = new Timing();
	tfull.start();

	IJ.log("ESI: Starting ESI analysis (normalizing from "+pxMin+" to "+pxMax+")");
	runAnalysis(aip.getStack(), pxMin, pxMax);
	
	tfull.stop();
	IJ.log("ESI: FINISHED in "+tfull);
    }



    /** Same as MTP.doEBI, but single-threaded, one loop
     *  over all image lines. */
    FloatProcessor doESI( TraceStack pxl, int order ){
	FloatProcessor ret =
	    new FloatProcessor(2*pxl.getWidth(),2*pxl.getHeight());
	ESI_internal( ret, pxl, order, 0, pxl.getHeight()-1);
	if (normOutput)
	    normalize( ret );
	return ret;
    }


    /** Helper class for multi-threading the analysis. 
     *  Quite basic / low-level multithreading,
     *  each thread just calculates a subset of output lines in the image. */
    class MTP implements Runnable {
	
	// for parameter access
	final private ESI_Analysis param;
	MTP( ESI_Analysis p ) { param = p; };
	
	// parameters for each thread to know which lines to calculate
	private int job, start, end, order;
	private TraceStack source;
	private FloatProcessor target;

	/** Starts as many threads as there are CPU cores, each thread calculates
	 *  a subset of output image lines.
	 *  This starts multiple threads, where each runs 'ESI_internal' with
	 *  matching line ranges */
	FloatProcessor doESI( TraceStack trSt , int ord ) {
	    int cores = Runtime.getRuntime().availableProcessors();
	    return doESI( trSt, ord, cores );
	}

	/** Start nrThreads, each calculating a subset of output image lines.
	 *  This starts multiple threads, where each runs 'ESI_internal' with
	 *  matching line ranges */
	FloatProcessor doESI( TraceStack trSt , int ord , int nrThreads ) {

	    // result image
	    FloatProcessor ret = 
		new FloatProcessor(trSt.getWidth()*2, trSt.getHeight()*2);

	    // size of input stack
	    final int depth = trSt.getDepth();
	    final int height = trSt.getHeight();

	    // thread setup and start
	    MTP [] trds = new MTP[nrThreads];
	    Thread [] ts = new Thread[ nrThreads ];
	    for (int i=0; i<nrThreads; i++) {
		trds[i] = new MTP( param );
		trds[i].job=1;
		trds[i].start = (height/nrThreads)*i;
		int tmp = (height/nrThreads)*(i+1);
		trds[i].end   = (tmp>height-1)?(height-1):tmp;
		trds[i].order = ord;
		trds[i].source = trSt;
		trds[i].target = ret;
		ts[i] = new Thread(trds[i]);
		ts[i].start();
	    }

	    // wait for completion
	    try {
	    for ( Thread  i : ts )
		i.join();
	    } catch (java.lang.InterruptedException e) {
		IJ.log(e.toString());
	    }
	
	    // normalize result image
	    if (normOutput)
		ESI_Analysis.normalize(ret);

	    return ret;
	}

	/** Run function, calls into ESI_internal with matching parameters */
	public void run() {
	    if (job==1)
		ESI_internal( target, source, order, start, end);
	};

    }

    
    /** normalize the pixels of a FloatProcessor to a given range. */
    static void normalize( FloatProcessor fp, 
	double low, double high ) {

	float [] data = (float [])fp.getPixels();

	// find min and max
	float min = Float.MAX_VALUE;
	float max = Float.MIN_VALUE;
	for ( float i : data ) {
		if (i<min) min = i;
		if (i>max) max = i;
	    }

	// scale
	for (int i=0; i<data.length; i++) 
	    data[i] = (float)
		((( data[i]-min )/( max-min )) * (high-low) + low);

	// tell it to recalculate its display scaling
	fp.resetMinAndMax();
    }

    /** normalize the pixels of a FloatProcessor to 0..1 */
    static void normalize( FloatProcessor fp ) {
	normalize( fp, 0.0, 1.0);
    }

    /** returns the minimum and maximum value of an ImageStack */
    static float [] getMinMax( ImageStack is ) {

	final int w = is.getWidth();
	final int h = is.getHeight();
	float [] ret = new float[2];
	ret[0]=Float.MAX_VALUE;
	ret[1]=Float.MIN_VALUE;

	for (int z=1;z<=is.getSize();z++) {
	    ImageProcessor ip = is.getProcessor(z);
	    for (int y=0;y<h;y++)
	    for (int x=0;x<w;x++) {
		float val = ip.getf(x,y);
		if (val<ret[0]) ret[0] = val;
		if (val>ret[1]) ret[1] = val;
	    }
	}
    
	return ret;
    }

    /** add one FloatProcessor to a sum. */
    static void addFP( FloatProcessor in, FloatProcessor sum) {

	float [] inD  = (float [])in.getPixels();
	float [] sumD = (float [])sum.getPixels();
	
	for (int i=0;i<sumD.length;i++)
	    sumD[i] += inD[i];

    }


    /** Provides a simple timer for runtime measurement */
    class Timing {
	long start, stop, runtime, outtime;
	Timing() { start =  System.currentTimeMillis(); }
	public void start() { start = System.currentTimeMillis(); };
	public void stop() { 
	    stop = System.currentTimeMillis(); 
	    runtime += stop-start;
	    outtime=runtime;
	    runtime=0;
	    }
	public void hold(){
	    stop = System.currentTimeMillis();
	    runtime += stop-start;
	    outtime  = runtime;
	    start =stop;
	}
	@Override public String toString(){ return("ms: "+(outtime));}

    }

}
