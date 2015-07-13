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

import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.IJ;

/** A class to compute and store 
 *  {@link PixelTrace} from a 3d stack. */
public class TraceStack {

    /** pixel data */
    protected PixelTrace [] data;

    final int width;
    final int height;
    final int nrTraces; // width * height
    final int depth;	// length of each trace, z-dim of stack

    int rStart, rEnd;	// range to analyse

    /** Construct from a full ImageStack. Same as TraceStack( is, 0, is.getSize()). */
    public TraceStack( ImageStack is ) {
	this(is,0,is.getSize());
    }

    /** Construct from an image stack (sub-stack range 'start' to 'end'). 
     * Fiji's ImageProcessor.getf(x,y) is used to obtain float values. */
    public TraceStack( ImageStack is, int start, int end) {
	
	// check sizes
	if (start<0) start = 0;
	if (end>is.getSize()) end = is.getSize();

	// set our size
	width    = is.getWidth();
	height   = is.getHeight();
	depth    = (end-start);
	nrTraces = width*height;

	// allocate the arrays
	data = new PixelTrace[ nrTraces ];
	for ( int y = 0; y < height; y++)
	for ( int x = 0; x < width ; x++) 
		data[y*width+x] = new PixelTrace(depth);

	// copy the data
	initFromStack( is, start );

    }

    /** Inits from stack, internal function called by constructor */
    void initFromStack( ImageStack is, int start) {
	
	// loop the stack in z
	for ( int z=0; z < depth ; z++) {
		// get the current image at z
		IJ.showProgress( z, depth);
		ImageProcessor curImg = is.getProcessor(z+start+1);

		// assign its pixel data into the traces
		for ( int y = 0; y < height; y++)
		for ( int x = 0; x < width ; x++) {
			data[y*width+x].setf(z, curImg.getf(x,y) );
		}
	}
    }


    /** Return the pixel trace at x,y. */
    public PixelTrace get(int x, int y) {
	return get( y * width + x);
    }
    
    /** Return the pixel trace at n=x*width + y. */
    public PixelTrace get(int n) {
	return data[n];
    }
    
    /** Return the width */
    public int getWidth() { return width; }
    /** Return the height */
    public int getHeight() { return height; }
    /** Return the depth */
    public int getDepth() { return depth; }

    /** Normalize the complete TraceStack to 0..1, given a min and max value */
    void normalizeFrom(float curMin, float curMax) {
	for ( PixelTrace i : data ) {
	    i.scale( curMin, curMax ); 
	}
    }

    /** Calculate the probability binning for all traces.
     * Data outside [min,max] map to the lowest/highest bin */
    public void createBinning( int nrBins , float min, float max ) {
	for ( PixelTrace i : data )
	    i.createBinning( nrBins, min, max);
    }

    /** Calculate the probability binning for all traces.
     * Data outside [0..1] map to the lowest/highest bin */
    public void createNormBinning( int nrBins ) {
	createBinning( nrBins, 0, 1);
    }

}
