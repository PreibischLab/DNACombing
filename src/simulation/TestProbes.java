package simulation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class TestProbes
{
	public static enum MatchResult{ CORRECT, AMBIVALENT, WRONG, NO_MATCH, NOT_ENOUGH_PROBES };

	public static int[] randomlySample( final ArrayList< CombingProbe > probes, final long length, final int iterations, final boolean mirror )
	{
		final Random rnd = new Random( 35 );
		return randomlySample( probes, length, iterations, mirror, rnd );
	}

	public static int[] randomlySample( final ArrayList< CombingProbe > probes, final long length, final int iterations, final boolean mirror, final Random rnd )
	{
		long min = probes.get( 0 ).start();
		long max = probes.get( 0 ).end();

		for ( final CombingProbe p : probes )
		{
			min = Math.min( Math.min( min, p.start() ), p.end() );
			max = Math.max( Math.max( max, p.start() ), p.end() );
		}

		// only consider random segments that contain at least some part of the GMC area
		min -= length;
		final long size = max - min + 1;

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

				final double maxError = Math.max( 1, rnd.nextGaussian() * 2 + 10 );
				m = match( distGMCs, distDetect, maxError, contained.get( 0 ).id() - 1, mirror, false ); // probe ID starts with 1, not 0

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

		return resultHist;
	}

	public static MatchResult match( final double[] distGMCs, final double distDetect[], final double maxErrorPx, final int correctMatch, final boolean mirror, final boolean debug )
	{
		int countMatches = 0;
		int minErrorIndex = -1;
		double minError = Double.MAX_VALUE;

		if ( debug )
		{
			System.out.println( "GMCdist:" );
			for ( int i = 0; i < distGMCs.length; ++i )
				System.out.println( i + ":" + distGMCs[ i ] );

			if ( mirror ) // we should detect nothing in this direction
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

		if ( mirror ) // we should detect nothing in this direction
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

	public static void main( String[] args ) throws IOException
	{
		//System.out.println( new String( "Probe Id;Chromosome;Begin;End;Gap length" ).matches( "Hello" ) );
		for ( int i = 2; i <= 2; ++i )
		{
			ArrayList< CombingProbe > probes = CombingProbe.loadFile( new File( "GMC_" + i + ".csv" ), i );
	
			Collections.sort( probes );
	
			System.out.println( probes.size() + " probes total." );
			
			//for ( final CombingProbe p : probes )
			//	System.out.println( p );
	
			System.out.println( randomlySample( probes, 400000, 10000, true )[ 0 ] );
		}
	}

}
