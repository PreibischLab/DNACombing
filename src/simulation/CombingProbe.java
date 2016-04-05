package simulation;

public class CombingProbe implements Comparable< CombingProbe >
{
	//P2;chr7;116062000;116074000;50

	final int id;
	final int chr;
	final long start, end;

	public CombingProbe( final int id, final int chr, final long start, final long end )
	{
		this.id = id;
		this.chr = chr;
		this.start = start;
		this.end = end;
	}

	public int id() { return id; }
	public int chr() { return chr; }
	public long start() { return start; }
	public long end() { return end; }
	public int length() { return (int)( end()-start() + 1 ); }
	public double inPixels() { return (double)length() / nucleotidesPerPixel(); }

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
	public String toString() { return "P" + id() + ";chr" + chr() + ";" + start + ";" + end() + ";l=" + length() + " (" + inPixels() + "px)"; }
}
