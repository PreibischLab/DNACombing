package simulation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.Random;

import net.imglib2.KDTree;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPoint;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;

public class DesignMorseCode
{
	final private static boolean containsSameSet( final ArrayList< CombingProbe > selectedProbes, final ArrayList< ArrayList< CombingProbe > > bestProbes )
	{
		for ( final ArrayList< CombingProbe > set : bestProbes )
			if ( equals( selectedProbes, set ) )
				return true;

		return false;
	}

	final private static boolean equals( final ArrayList< CombingProbe > selectedProbes, final ArrayList< CombingProbe > bestProbes )
	{
		if ( selectedProbes.size() != bestProbes.size() )
			return false;

		for ( int i = 0; i < selectedProbes.size(); ++i )
		{
			final CombingProbe p = selectedProbes.get( i );
			final CombingProbe q = bestProbes.get( i );

			if ( p.start() != q.start() || p.end() != q.end() )
				return false;
		}

		return true;
	}

	final private static void sortIntoList( final int r, final ArrayList< CombingProbe > selectedProbes, final ArrayList< Integer > bestInt, final ArrayList< ArrayList< CombingProbe > > bestProbes )
	{
		if ( r > bestInt.get( bestInt.size() - 1 ) )
		{
			for ( int j = bestInt.size() - 1; j > 0; --j )
			{
				if ( r < bestInt.get( j - 1 ) )
				{
					bestInt.add( j, r );
					bestProbes.add( j, selectedProbes );

					bestInt.remove( bestInt.size() -1 );
					bestProbes.remove( bestProbes.size() - 1 );

					return;
				}
			}

			bestInt.add( 0, r );
			bestProbes.add( 0, selectedProbes );

			bestInt.remove( bestInt.size() -1 );
			bestProbes.remove( bestProbes.size() - 1 );

			//System.out.println( i + ", best: " + ( (double)best / (double)testIterations ) * 100.0 + " %" );
		}
	}

	public static Pair< ArrayList< Integer >, ArrayList< ArrayList< CombingProbe > > > designBestEqual(
			final ArrayList< CombingProbe > allProbes,
			final int numProbes,
			final int combingLength,
			final double minDistanceProbes,
			final int numIterations,
			final int testIterations,
			final int nBest,
			final Random rnd )
	{
		final ArrayList< Integer > bestInt = new ArrayList< Integer >();
		final ArrayList< ArrayList< CombingProbe > > bestProbes = new ArrayList< ArrayList<CombingProbe> >();

		final ArrayList< RealLocalizable > positions = new ArrayList< RealLocalizable >();
		final ArrayList< CombingProbe > probes = new ArrayList< CombingProbe >();

		final long start = allProbes.get( 0 ).start();

		for ( final CombingProbe p : allProbes )
		{
			positions.add( new RealPoint( ( p.start() - start ) / CombingProbe.nucleotidesPerPixel() ) );
			positions.add( new RealPoint( ( p.end() - start ) / CombingProbe.nucleotidesPerPixel() ) );
			probes.add( p );
			probes.add( p );

			//System.out.println( (( p.start() - start ) / CombingProbe.nucleotidesPerPixel() ) + " >>> " + ( ( p.end() - start ) / CombingProbe.nucleotidesPerPixel() ) + ": " + p );
		}

		final double avgStep = (allProbes.get( allProbes.size() - 1 ).end() - allProbes.get( 0 ).start()) / numProbes / CombingProbe.nucleotidesPerPixel();

		final KDTree< CombingProbe > tree = new KDTree< CombingProbe >( probes, positions );
		final NearestNeighborSearchOnKDTree< CombingProbe > search = new NearestNeighborSearchOnKDTree< CombingProbe >( tree );

		for ( int i = 0; i < nBest; ++i )
		{
			bestInt.add( -1 );
			bestProbes.add( new ArrayList< CombingProbe >() );
		}

		int i = 0;

		do
		{
			final ArrayList< CombingProbe > selectedProbes = new ArrayList< CombingProbe >();
			final double sigma = 3 + rnd.nextInt( (int)Math.round( avgStep/2 ) - 3 );

			//System.out.println( sigma + " > " + (avgStep/2) );
			for ( int k = 0; k < numProbes; ++k )
			{
				final double position = (rnd.nextGaussian() * sigma + avgStep*k );
				search.search( new RealPoint( position ) );
				//System.out.println( Math.round( k*avgStep ) + ": " + position + " " + search.getSampler().get() );
				selectedProbes.add( search.getSampler().get() );
			}

			//System.out.println( isValid( selectedProbes, minDistanceProbes ) );
			
			if ( !containsSameSet( selectedProbes, bestProbes ) && isValid( selectedProbes, minDistanceProbes ) )
			{
				final int r = TestProbes.randomlySample( selectedProbes, combingLength, testIterations, true, rnd )[ 0 ];
				sortIntoList( r, selectedProbes, bestInt, bestProbes );

			}
		}
		while ( ++i < numIterations );

		return new Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> >( bestInt, bestProbes );
	}

	private static final boolean isValid( final ArrayList< CombingProbe > probes, final double minDistanceProbes )
	{
		Collections.sort( probes );

		for ( int i = 0; i < probes.size() - 1; ++i )
		{
			final CombingProbe p = probes.get( i );
			final CombingProbe q = probes.get( i + 1 );
	
			if ( p.distanceTo( q ) / CombingProbe.nucleotidesPerPixel() < minDistanceProbes )
			{
				//System.out.println( "Fail on:" + p + " <>" + q );
				return false;
			}
		}

		return true;
	}

	public static Pair< ArrayList< Integer >, ArrayList< ArrayList< CombingProbe > > > designBest(
			final ArrayList< CombingProbe > allProbes,
			final int numProbes,
			final int combingLength,
			final double minDistanceProbes,
			final int numIterations,
			final int testIterations,
			final int nBest,
			final Random rnd )
	{
		final ArrayList< Integer > bestInt = new ArrayList< Integer >();
		final ArrayList< ArrayList< CombingProbe > > bestProbes = new ArrayList< ArrayList<CombingProbe> >();

		for ( int i = 0; i < nBest; ++i )
		{
			bestInt.add( -1 );
			bestProbes.add( new ArrayList< CombingProbe >() );
		}

		final ArrayList< CombingProbe > allProbesCopy = new ArrayList< CombingProbe >();

B:		for ( int i = 0; i < numIterations; ++i )
		{
			allProbesCopy.clear();
			allProbesCopy.addAll( allProbes );

			final ArrayList< CombingProbe > selectedProbes = new ArrayList< CombingProbe >();

A:			do
			{
				if ( allProbesCopy.size() == 0 )
					continue B;

				final int random = rnd.nextInt( allProbesCopy.size() );
				final CombingProbe p = allProbesCopy.get( random );
				allProbesCopy.remove( random );

				// too close to one of the probes?
				for ( final CombingProbe q : selectedProbes )
					if ( p.distanceTo( q ) / CombingProbe.nucleotidesPerPixel() < minDistanceProbes )
						continue A;

				selectedProbes.add( p );
			}
			while ( selectedProbes.size() < numProbes );

			final int r = TestProbes.randomlySample( selectedProbes, combingLength, testIterations, true, rnd )[ 0 ];

			sortIntoList( r, selectedProbes, bestInt, bestProbes );
		}

		return new Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> >( bestInt, bestProbes );
	}

	public static Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> > iterateAll(
			final ArrayList< CombingProbe > allProbes,
			final ArrayList<ArrayList<CombingProbe>> probes,
			final int combingLength,
			final double minDistanceProbes,
			final int testIterations,
			final int nBest,
			final Random rnd )

	{
		final ArrayList< Integer > bestInt = new ArrayList< Integer >();
		final ArrayList< ArrayList< CombingProbe > > bestProbes = new ArrayList< ArrayList<CombingProbe> >();

		for ( final ArrayList< CombingProbe > probeset : probes )
		{
			final Random random;
			
			if ( rnd == null )
				random = new Random( 35 );
			else
				random = rnd;

			final Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> > local = iterate( allProbes, probeset, combingLength, minDistanceProbes, testIterations, nBest, random );

			for ( int i = 0; i < local.getA().size(); ++i )
			{
				if ( !contains( bestProbes, local.getB().get( i ) ) )
				{
					bestInt.add( local.getA().get( i ) );
					bestProbes.add( local.getB().get( i ) );
				}
			}
		}

		return new Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> >( bestInt, bestProbes );
	}

	public static boolean contains( final ArrayList< ArrayList< CombingProbe > > all, final ArrayList< CombingProbe > add )
	{
		if ( all.size() == 0 )
			return false;

		for ( int i = 0; i < all.size(); ++i )
		{
			final ArrayList< CombingProbe > oldP = new ArrayList< CombingProbe >();
			final ArrayList< CombingProbe > newP = new ArrayList< CombingProbe >();
	
			oldP.addAll( all.get( i ) );
			newP.addAll( add );

			Collections.sort( oldP );
			Collections.sort( newP );

			for ( int x = 0; x < oldP.size(); ++x )
			{
				final CombingProbe o1 = oldP.get( x );
				final CombingProbe o2 = newP.get( x );

				if ( o1.gmcId() != o2.gmcId() || o1.id() != o2.id() )
					return false;
			}
		}

		return true;
	}

	public static Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> > iterate(
			final ArrayList< CombingProbe > allProbes,
			final ArrayList< CombingProbe > probes,
			final int combingLength,
			final double minDistanceProbes,
			final int testIterations,
			final int nBest,
			final Random rnd )
	{
		final ArrayList< Integer > bestInt = new ArrayList< Integer >();
		final ArrayList< ArrayList< CombingProbe > > bestProbes = new ArrayList< ArrayList<CombingProbe> >();

		for ( int i = 0; i < nBest; ++i )
		{
			bestInt.add( -1 );
			bestProbes.add( new ArrayList< CombingProbe >() );
		}

		for ( int i = 0; i < probes.size(); ++i )
		{
A:			for ( int j = 0; j < allProbes.size(); ++j )
			{
				final ArrayList< CombingProbe > selectedProbes = new ArrayList< CombingProbe >();
				selectedProbes.addAll( probes );
				selectedProbes.remove( i );

				final CombingProbe p = allProbes.get( j );

				// too close to one of the probes?
				for ( final CombingProbe q : selectedProbes )
					if ( p.distanceTo( q ) / CombingProbe.nucleotidesPerPixel() < minDistanceProbes )
						continue A;

				selectedProbes.add( p );

				final int r = TestProbes.randomlySample( selectedProbes, combingLength, testIterations, true, rnd )[ 0 ];

				sortIntoList( r, selectedProbes, bestInt, bestProbes );
			}
		}

		return new Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> >( bestInt, bestProbes );
	}

	public static Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> > exchangeAll(
			final ArrayList< CombingProbe > allProbes,
			final ArrayList<ArrayList<CombingProbe>> probes,
			final int numProbes,
			final int combingLength,
			final double minDistanceProbes,
			final int numIterations,
			final int testIterations,
			final int nBest,
			final Random rnd )

	{
		final ArrayList< Integer > bestInt = new ArrayList< Integer >();
		final ArrayList< ArrayList< CombingProbe > > bestProbes = new ArrayList< ArrayList<CombingProbe> >();

		for ( final ArrayList< CombingProbe > probeset : probes )
		{
			final Random random;
			
			if ( rnd == null )
				random = new Random( 35 );
			else
				random = rnd;

			final Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> > local = exchange( allProbes, probeset, numProbes, combingLength, minDistanceProbes, numIterations, testIterations, nBest, random );

			for ( int i = 0; i < local.getA().size(); ++i )
			{
				if ( !contains( bestProbes, local.getB().get( i ) ) )
				{
					bestInt.add( local.getA().get( i ) );
					bestProbes.add( local.getB().get( i ) );
				}
			}
		}

		return new Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> >( bestInt, bestProbes );
	}

	public static Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> > exchange(
			final ArrayList< CombingProbe > allProbes,
			final ArrayList< CombingProbe > probes,
			final int numProbes,
			final int combingLength,
			final double minDistanceProbes,
			final int numIterations,
			final int testIterations,
			final int nBest,
			final Random rnd )
	{
		final ArrayList< Integer > bestInt = new ArrayList< Integer >();
		final ArrayList< ArrayList< CombingProbe > > bestProbes = new ArrayList< ArrayList<CombingProbe> >();

		for ( int i = 0; i < nBest; ++i )
		{
			bestInt.add( -1 );
			bestProbes.add( new ArrayList< CombingProbe >() );
		}

		for ( int i = 0; i < numIterations; ++i )
		{
			final ArrayList< CombingProbe > selectedProbes = new ArrayList< CombingProbe >();
			selectedProbes.addAll( probes );

			for ( int j = 0; j < numProbes; ++j )
				selectedProbes.remove( rnd.nextInt( selectedProbes.size() ) );

A:			do
			{
				final CombingProbe p = allProbes.get( rnd.nextInt( allProbes.size() ) );

				// too close to one of the probes?
				for ( final CombingProbe q : selectedProbes )
					if ( p.distanceTo( q ) / CombingProbe.nucleotidesPerPixel() < minDistanceProbes )
						continue A;

				selectedProbes.add( p );
			}
			while( selectedProbes.size() < probes.size() );

			final int r = TestProbes.randomlySample( selectedProbes, combingLength, testIterations, true, rnd )[ 0 ];

			sortIntoList( r, selectedProbes, bestInt, bestProbes );
		}

		return new Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> >( bestInt, bestProbes );
	}

	public static Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> > filterBest(
			final Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> > best,
			final int combingLength,
			final int numIterations,
			final int nBest,
			final Random rnd )
	{
		final ArrayList< Integer > bestInt = new ArrayList< Integer >();
		final ArrayList< ArrayList< CombingProbe > > bestProbes = new ArrayList< ArrayList<CombingProbe> >();
		
		for ( int i = 0; i < nBest; ++i )
		{
			bestInt.add( -1 );
			bestProbes.add( new ArrayList< CombingProbe >() );
		}

		for ( int i = 0; i < best.getA().size(); ++i )
		{
			if ( best.getB().get( i ).size() > 0 )
			{
				final Random random;
				if ( rnd == null )
					random = new Random( 35 );
				else
					random = rnd;

				final int r = TestProbes.randomlySample( best.getB().get( i ), combingLength, numIterations, true, random )[ 0 ];
				sortIntoList( r, best.getB().get( i ), bestInt, bestProbes );
			}
		}

		return new Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> >( bestInt, bestProbes );
	}

	public static Pair< Integer, ArrayList< CombingProbe > > optimalProbesFor( final ArrayList< CombingProbe > allProbes, final int numProbes )
	{
		final Random rnd = new Random( 353 );

		final int combingLength = 400000;
		final double minDistanceProbes = 10.0;
		final int numIterations = 10000000;
		final int testIterations = 100;
		final int nBest = 100;

		Pair< ArrayList<Integer>, ArrayList<ArrayList<CombingProbe>> > best = 
				designBestEqual( allProbes, numProbes, combingLength, minDistanceProbes, numIterations, testIterations, nBest, rnd );

//		for ( int i = 0; i < 10; ++i )
//			System.out.println( best.getA().get( i ) );

//		System.out.println( "..." );
//		System.out.println( best.getA().get( best.getA().size() - 1 ) );

		
		// TODO: hack - put it their design and improve on it
		{
			int gmcId = -1;
			
			if ( numProbes == 16 )
				gmcId = 1;
			else if ( numProbes == 17 )
				gmcId = 2;
			else if ( numProbes == 18 )
				gmcId = 3;
			else if ( numProbes == 20 )
				gmcId = 4;
			else if ( numProbes == 21 )
				gmcId = 5;
			else if ( numProbes == 22 )
				gmcId = 6;
			else if ( numProbes == 23 )
				gmcId = 7;
			else if ( numProbes == 24 )
				gmcId = 8;
			else if ( numProbes == 29 )
				gmcId = 9;
			else if ( numProbes == 30 )
				gmcId = 10;

			final ArrayList< CombingProbe > designedProbes = new ArrayList< CombingProbe >();
			
			for ( final CombingProbe p : allProbes )
				if ( p.gmcId() == gmcId )
					designedProbes.add( p );
	
			if ( gmcId > 0 )
			{
				System.out.println( "Adding designed probes gmcId=" + gmcId );
				best.getA().add( 1 );
				best.getB().add( designedProbes );
			}
		}
		
		final int numIterationsFilter = 10000;
		final int nBestFilter = 10;

		best = filterBest( best, combingLength, numIterationsFilter, nBestFilter, rnd );

		//for ( int i = 0; i < nBestFilter; ++i )
		//	System.out.println( best.getA().get( i ) );

		//System.out.println();

		int bestAll = -1;
		ArrayList< CombingProbe > bestProbesAll = null;

		int noChange = 0;
		int noBetter = 0;
		int last = -1;

		for ( int x = 0; x < 10000; ++x )
		{
			//System.out.print( x +": " );

			if ( noChange >= 5 || noBetter >= 50 )
			{
				noChange = 0;
				noBetter = 0;
				final int ex = Math.max(  2, rnd.nextInt(numProbes-1)/2 );
				//System.out.print( "ex=" + ex + " " );
				best = exchangeAll( allProbes, best.getB(), ex, combingLength, minDistanceProbes, 1000, testIterations, nBestFilter, rnd );

				// every third time put the best one back in
				if ( rnd.nextInt( 3 ) == 0 )
				{
					final ArrayList< CombingProbe > bestP = new ArrayList< CombingProbe >();
					bestP.addAll( bestProbesAll );
					best.getA().add( bestAll );
					best.getB().add( bestP );
					//System.out.print( "best back (" + bestAll + ") " );
				}
			}
			else
			{
				best = iterateAll( allProbes, best.getB(), combingLength, minDistanceProbes, testIterations, nBestFilter, rnd );
			}

			//System.out.print( best.getA().get( 0 ) + " >>> " );

			//for ( int i = 0; i < best.getA().size(); ++i )
			//	System.out.println( best.getA().get( i ) );

			best = filterBest( best, combingLength, numIterationsFilter, nBestFilter, null );

			//System.out.print( best.getA().get( 0 ) );

			if ( last == best.getA().get( 0 ) )
			{
				++noChange;
				++noBetter;
			}
			else
			{
				noChange = 0;
				last = best.getA().get( 0 );
			}

			if ( best.getA().get( 0 ) < bestAll )
				++noBetter;

			if ( best.getA().get( 0 ) > bestAll )
			{
				bestAll = best.getA().get( 0 );
				bestProbesAll = new ArrayList< CombingProbe >();
				bestProbesAll.addAll( best.getB().get( 0 ) );
				//System.out.println( ", new best." );
				noBetter = 0;
			}
			else
			{
				//System.out.println();
			}
			
			//if ( x > 0 && x % 1000 == 0 )
			//	System.out.println( TestProbes.randomlySample( bestProbesAll, 400000, 100000, true )[ 0 ] );
		}

		System.out.println( new Date( System.currentTimeMillis() ) + ", " + numProbes + ": "  + TestProbes.randomlySample( bestProbesAll, 400000, 100000, true )[ 0 ] );

		return new Pair< Integer, ArrayList<CombingProbe> >( bestAll, bestProbesAll );
	}
	
	public static void main( String[] args ) throws IOException
	{
		final ArrayList< CombingProbe > allProbes = new ArrayList< CombingProbe >();

		for ( int i = 1; i <= 10; ++i )
		{
			final ArrayList< CombingProbe > probes = CombingProbe.loadFile( new File( "GMC_" + i + ".csv" ), i );

			for ( final CombingProbe p : probes )
			{
				boolean contains = false;

				for ( final CombingProbe q : allProbes )
					if ( q.start() == p.start() && q.end() == p.end() )
						contains = true;

				if ( !contains )
					allProbes.add( p );
			}
		}

		Collections.sort( allProbes );

		System.out.println( allProbes.size() + " probes total." );

		for ( int i = 16; i <= 30; ++i )
		{
			final Pair< Integer, ArrayList< CombingProbe > > best = optimalProbesFor( allProbes, i );
		}
	}

}
