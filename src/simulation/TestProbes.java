package simulation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;

public class TestProbes
{
	public static enum MatchResult{ CORRECT, AMBIVALENT, WRONG, NO_MATCH, NOT_ENOUGH_PROBES };

	public static int[] randomlySample( final ArrayList< CombingProbe > probes, final long length, final int iterations, final long min, final long max )
	{
		final Random rnd = new Random( 35 );
		return randomlySample( probes, length, iterations, min, max, rnd );
	}

	public static int[] randomlySample( final ArrayList< CombingProbe > probes, final long length, final int iterations, long min, long max, final Random rnd )
	{
		Collections.sort( probes );

		// set a linearly increasing ID as the matching itselfs relies on it to determine the correct match (HORRIBLE!!)
		for ( int i = 0; i < probes.size(); ++i )
			probes.get( i ).setTmpId( i + 1 );

		final double error = 10;
		/*
		long min = probes.get( 0 ).start();
		long max = probes.get( 0 ).end();

		for ( final CombingProbe p : probes )
		{
			min = Math.min( Math.min( min, p.start() ), p.end() );
			max = Math.max( Math.max( max, p.start() ), p.end() );
		}
		*/

		// NOT ANYMORE: only consider random segments that contain at least some part of the GMC area
		min -= length/3;
		final long size = max - min + 1 - length + length/3; // only consider probes that lie within the area

		// histogram of how many probes are hit
		//final int[] containsHist = new int[ 100 ];

		// the distances
		final double[] distGMCs = new double[ probes.size() - 1 ];

		for ( int i = 0; i < probes.size() - 1; ++i )
		{
			final CombingProbe p = probes.get( i );
			distGMCs[ i ] = p.distanceTo( probes.get( i + 1 ) ) / CombingProbe.nucleotidesPerPixel();
			//System.out.println( i + ": " + distGMCs[ i ] + " >>> " + (p.id() - 1) + "->" + ( probes.get( i + 1 ).id() - 1 ) );
		}

		// histogram of how many probes are hit
		final int[] resultHist = new int[ MatchResult.values().length ];

		for ( int i = 0; i < iterations; ++i )
		{
			final long from = Math.round( rnd.nextDouble() * size ) + min;
			final long to = from + length;

			final ArrayList< CombingProbe > contained = new ArrayList< CombingProbe >();
	
			for ( final CombingProbe p : probes )
				if ( p.isContained( from, to ) > 0.25 )
					contained.add( p );

			//if ( contained.size() == 0 )
			//	System.out.println( from + " >>> " + to );

			//++containsHist[ contained.size() ];

			final MatchResult m;

			if ( contained.size() > 1 )
			{
				final double[] distDetect = new double[ contained.size() - 1 ];

				for ( int j = 0; j < contained.size() - 1; ++j )
					distDetect[ j ] = contained.get( j ).distanceTo( contained.get( j + 1 ) ) / CombingProbe.nucleotidesPerPixel();

				final double maxError = Math.max( 1, rnd.nextGaussian() * 2 + error );
				m = match( distGMCs, distDetect, maxError, contained.get( 0 ).id() - 1, false ); // probe ID starts with 1, not 0

				/*
				if ( m == MatchResult.AMBIVALENT )
				{
					System.out.println( m );
					System.out.println( from + " >>> " + to + " contains:" );
					for ( final CombingProbe p : contained )
						System.out.println( p );
					match( distGMCs, distDetect, 10, contained.get( 0 ).id() - 1, mirror, true ); // probe ID starts with 1, not 0
					System.exit( 0 );
				}*/
			}
			else
			{
				m = MatchResult.NOT_ENOUGH_PROBES;
			}

			++resultHist[ m.ordinal() ];
		}

		//for ( int i = 0; i < containsHist.length; ++i )
		//	System.out.println( containsHist[ i ] );

		
		//for ( int i = 0; i < resultHist.length; ++i )
		//	System.out.println( /*MatchResult.values()[ i ] + ": " + */ resultHist[ i ] );
		
		//System.out.println();

		for ( int i = 0; i < probes.size(); ++i )
			probes.get( i ).resetTmpId();

		return resultHist;
	}

	public static MatchResult match( final double[] distGMCs, final double distDetect[], final double maxErrorPx, final int correctMatch, final boolean debug )
	{
		int countMatches = 0;
		int minErrorIndex = -1;
		double minError = Double.MAX_VALUE;

		if ( debug )
		{
			System.out.println( "GMCdist:" );
			for ( int i = 0; i < distGMCs.length; ++i )
				System.out.println( i + ":" + distGMCs[ i ] );

			// we should detect nothing in this direction
			{
				System.out.println( "GMCdist (mirror):" );
				final double[] distGMCsMirrored = new double[ distGMCs.length ];
				for ( int i = 0; i < distGMCs.length; ++i )
					distGMCsMirrored[ distGMCs.length - 1 - i ] = distGMCs[ i ];
				for ( int i = 0; i < distGMCsMirrored.length; ++i )
					System.out.println( i + ":" + distGMCsMirrored[ i ] );
			}

			System.out.println( "detect:" );
			for ( int i = 0; i < distDetect.length; ++i )
				System.out.println( i + ":" + distDetect[ i ] );

			System.out.println( "errors: (correct match=" + correctMatch + ")" );
		}

		for ( int i = 0; i < distGMCs.length - distDetect.length + 1; ++i )
		{
			double error = 0;

			for ( int j = 0; j < distDetect.length; ++j )
				error += Math.abs( distGMCs[ i + j ] - distDetect[ j ] );

			if ( debug )
				System.out.println( i + ": " + error );

			if ( error <= maxErrorPx )
			{
				++countMatches;
				if ( error < minError )
				{
					minError = error;
					minErrorIndex = i;
				}
			}
		}

		// we should detect nothing in this direction
		{
			final double[] distGMCsMirrored = new double[ distGMCs.length ];
			for ( int i = 0; i < distGMCs.length; ++i )
				distGMCsMirrored[ distGMCs.length - 1 - i ] = distGMCs[ i ];
	
			for ( int i = 0; i < distGMCsMirrored.length - distDetect.length + 1; ++i )
			{
				double error = 0;
	
				for ( int j = 0; j < distDetect.length; ++j )
					error += Math.abs( distGMCsMirrored[ i + j ] - distDetect[ j ] );
	
				if ( debug )
					System.out.println( i + " (mirror): " + error );
	
				if ( error <= maxErrorPx )
				{
					++countMatches;
					if ( error < minError )
					{
						minError = error;
						minErrorIndex = -1;
					}
				}
			}
		}

		if ( countMatches == 0 )
			return MatchResult.NO_MATCH;
		else if ( countMatches > 1 )
			return MatchResult.AMBIVALENT;
		else if ( minErrorIndex == correctMatch )
			return MatchResult.CORRECT;
		else
			return MatchResult.WRONG;
	}

	public static void drawProbes( final FloatProcessor fp, ArrayList< CombingProbe > probes, final long min, final int y )
	{
		for ( final CombingProbe p : probes )
		{
			int x0 = (int) Math.round( (p.start()-min) / CombingProbe.nucleotidesPerPixel() );
			int x1 = (int) Math.round( (p.end()-min) / CombingProbe.nucleotidesPerPixel() );

			ColorProcessor cp;

			for ( int x = x0; x <= x1; ++x )
			{
				fp.putPixelValue( x, y - 1, 1.0 );
				fp.putPixelValue( x, y, 2.0 );
				fp.putPixelValue( x, y + 1, 1.0 );
			}
		}
	}

	public static void main( String[] args ) throws IOException
	{
		final ArrayList< CombingProbe > allProbesDouble = DesignMorseCode.loadAllProbesets();

		final int combingLength = 350000;

		final long min = DesignMorseCode.min( allProbesDouble );
		final long max = DesignMorseCode.max( allProbesDouble );

		System.out.println( "Probes ranging from " + min + " to " + max + " (size=" + (max-min) + ", equals " + ( Math.round( (max-min)/CombingProbe.nucleotidesPerPixel() ) + 1 ) + " px)" );

		long minDraw = min - combingLength/3;
		long maxDraw = max + combingLength/3;

		System.out.println( "Area drawn ranges from " + minDraw + " to " + maxDraw + " (size=" + (maxDraw-minDraw) + ", equals " + ( Math.round( (maxDraw-minDraw)/CombingProbe.nucleotidesPerPixel() ) + 1 ) + " px)" );

		final ImageStack stack = new ImageStack( (int)( Math.round( (maxDraw-minDraw)/CombingProbe.nucleotidesPerPixel() ) + 1 ), 500 );

		for ( int it = 1; it < 5; ++it )
		{
			
			final FloatProcessor fp = new FloatProcessor( (int)( Math.round( (maxDraw-minDraw)/CombingProbe.nucleotidesPerPixel() ) + 1 ), 500 );
	
			for ( int i = 7; i <= 25; ++i )
			{
				final File f;
	
				//f = new File( "GMC_" + i + ".csv" );
				f = new File( "GMC_" + i + "_design_" + it + ".csv.tmp" );
	
				ArrayList< CombingProbe > probes = CombingProbe.loadFile( f, i );
		
				Collections.sort( probes );
		
				System.out.print( i + "\t" + probes.size() + "\t" );
				
				//for ( final CombingProbe p : probes )
				//	System.out.println( p );
		
				drawProbes( fp, probes, minDraw, 50 + (i-7)*15 );
	
				//System.out.println( randomlySample( probes, combingLength, 100000, min, max )[ 0 ] );
			}

			stack.addSlice( fp );
		}
		
		new ImageJ();
		new ImagePlus( "probes", stack ).show();
	}

}
