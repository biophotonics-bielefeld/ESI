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

/** This class holds one trace along a pixel in z-direction.
 *  It is created by {@link TraceStack}. It also provides
 *  various mathematical functions for the ESI analysis. All
 *  these functions are applied to a range [start:end], allowing
 *  to easily analyse substacks.
 *  */
public class PixelTrace {

    /** pixel data */
    final float [] data; 
    /** depth = data.length */
    final int depth;
    /** probability data */
    double [] probs = null;

    /** Create a new pixel trace length d */
    PixelTrace( int d ) {
	depth = d;
	data = new float [depth];
    }

    /** Returns the value of the n'th pixel */
    public float getf( int n ) {
	return data[n];
    }

    /** Sets the value of the n'th pixel */
    public void setf( int n, float val) {
	data[n] = val;
    }

    // ---------------------------------------------------------
    
    /** Creates the probability binning for H_cross2. */
    void createBinning( int nrBins , float min, float max) {
	probs = new double [ nrBins ]; 
	//final double binSize = (max-min)/(nrBins);

	// map our data to the bins
	for ( float d : data ) {
	    // calc the index
	    int idx = (int)( ((d-min)/(max-min)) * nrBins );
    	    // check the index (TODO: should this be an exception?)
	    if (idx < 0) idx = 0;
	    if (idx >= nrBins) idx = nrBins-1;
	    // counter up for the bin
	    probs[idx]++;
	}

	for (int i=0;i<probs.length;i++)
	    probs[i]/=data.length;

    }

    // ---------------------------------------------------------
    
    
    /** Calculate the joint n'th moment of two pixel traces. */
    public static double joint_moment( PixelTrace X, PixelTrace Y, 
	double order) {

	// TODO: check for matching ranges
	final double meanX=X.weightedMean(1);
	final double meanY=Y.weightedMean(1);

	// create storage
	PixelTrace dummy = new PixelTrace(X.depth);

	// loop the range
	for(int i=0;i < X.depth ; i++)
	    dummy.data[i] = (float)(
	        Math.pow(X.data[i]-meanX,order) * 
		Math.pow(Y.data[i]-meanY,order) );

	// TODO: should this also be a weighted mean?
	Double jmom = dummy.weightedMean(1);

	// TODO: should / will this be Nan??
	if(jmom.isNaN()) //Test for NAN
	    throw new RuntimeException("That was NaN");

	return jmom;

    }

    /** Calculating the cross entropy between two pixel traces. 
     *  Eq. 2 and 8 in the paper... */
    public static float H_cross2( PixelTrace X, PixelTrace Y, double order) {

	double h_sum_temp=0;
	
	// check if the probs have been calculated
	if ((X.probs==null)||(Y.probs==null))
	    throw new RuntimeException("Bining not initialized!");

	// Loop over binned probability space, see eq. 2
	for(int i=0; i<X.probs.length; i++)
	{
	    if(Y.probs[i]>0)
		h_sum_temp+= 
		    X.probs[i] *(Math.log(Y.probs[i])/Math.log(2));
	    if(X.probs[i]>0)
		h_sum_temp+= 
		    Y.probs[i] *(Math.log(X.probs[i])/Math.log(2));
	}

	// eq. 8, multiply with joined moment
	return (float)(-h_sum_temp*joint_moment(X,Y,order)/2.);
    }
        
    
    // ---------------------------------------------------------

    /** Returns the mean (for the set range). */
    public float mean() {
	float ret=0;
        for( float i : data )
	    ret+=i;
	return ret/depth;
    }

    /** Summed up data as sum data^order, */
    public float mean( float order) {
	float ret=0;
        for( float i : data )
	    ret+=(float)Math.pow(i,order);
	return ret/depth;
    }
   
    /** Returns the weighted mean. TODO: document the weighting */
    public double weightedMean( float weight)
    {
	double mean=0;
	final double factor=weight/depth;

	for( int i=0; i<depth; i++ ) {
	    // TODO: check and understand what this does
	    weight+=factor;
            mean+=weight*data[i];
        }

	return mean/depth;
    }


    /** Returns the maximum number in the trace */ 
    public float max(){
	float ret = Float.MIN_VALUE;
	for( float j : data ) 
	    if ( j > ret ) ret = j;
	return ret;
    }

    /** Returns the minimum number in the trace */ 
    public float min(){
	float ret = Float.MAX_VALUE;
	for( float j : data ) 
	    if ( j < ret ) ret = j;
	return ret;
    }
   
    /** for scaling, provide min and max, scales to 0..1 */
    void scale( float min, float max ) {
	for ( int i=0; i<depth; i++ )
	    data[i] = (data[i]-min)/(max-min) ;
	    
    }

    
}

