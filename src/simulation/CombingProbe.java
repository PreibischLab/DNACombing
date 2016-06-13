package simulation;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class CombingProbe implements Comparable< CombingProbe >
{
	//P2;chr7;116062000;116074000;50

	final int gmcId;
	final int originalId;
	final int chr;
	final long start, end;

	int tmpId;

	public CombingProbe( final int gmcId, final int id, final int chr, final long start, final long end )
	{
		this.gmcId = gmcId;
		this.originalId = id;
		this.chr = chr;
		this.start = start;
		this.end = end;

		this.tmpId = originalId;
	}

	public void setTmpId( final int id ) { this.tmpId = id; }
	public void resetTmpId() { this.tmpId = originalId; }

	public int gmcId() { return gmcId; }
	public int id() { return tmpId; }
	public int chr() { return chr; }
	public long start() { return start; }
	public long end() { return end; }
	public int length() { return (int)( end()-start() + 1 ); }
	public double inPixels() { return (double)length() / nucleotidesPerPixel(); }

	public CombingProbe copy()
	{
		final CombingProbe p =  new CombingProbe(gmcId, originalId, chr, start, end);
		p.setTmpId( p.id() );
		return p;
	}

	public static double nucleotidesPerPixel() { return 1500.0; }

	public long distanceTo( final CombingProbe p )
	{
		if ( p.start() > end() )
			return p.start() - end();
		else if ( start() > p.end() )
			return start() - p.end();
		else
			return 0; // overlaps
	}

	public double isContained( final long from, final long to )
	{
		if ( end() < from || start() > to )
			return 0.0;
		else if ( start() > from && end() < to )
			return 1.0;
		else if ( start() <= from )
			return (double)( end() - from + 1 ) / (double)length();
		else
			return (double)( to - start() + 1 ) / (double)length();
	}

	@Override
	public int compareTo( final CombingProbe arg0 )
	{
		return (int)(start() - arg0.start());
	}

	@Override
	public String toString() { return "GMC" + gmcId() + ";P" + id() + ";chr" + chr() + ";" + start + ";" + end() + ";l=" + length() + " (" + inPixels() + "px)"; }

	public static ArrayList< CombingProbe > loadFile( final File file, final int gmcId ) throws IOException
	{
		final BufferedReader in = TextFileAccess.openFileRead( file );
		final ArrayList< CombingProbe > probes = new ArrayList< CombingProbe >();

		while ( in.ready() )
		{
			final String line = in.readLine().trim();

			if ( line.matches( "P[0-9]+.*" ) )
			{
				final String[] elements = line.split( ";" );
				
				if ( elements.length < 4 )
					continue;

				final int id = Integer.parseInt( elements[ 0 ].substring( 1, elements[ 0 ].length() ) );
				final int chr = Integer.parseInt( elements[ 1 ].substring( 3, elements[ 1 ].length() ) );
				final long start = Long.parseLong( elements[ 2 ] );
				final long end = Long.parseLong( elements[ 3 ] );

				probes.add( new CombingProbe( gmcId, id, chr, start, end ) );
			}
		}

		in.close();

		return probes;
	}
}
