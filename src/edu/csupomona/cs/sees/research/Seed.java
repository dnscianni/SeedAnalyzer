package edu.csupomona.cs.sees.research;

public class Seed implements Comparable<Seed>{
	
	private float length;
	private int startLocation;
	
	public Seed(int sL,float l) {
		startLocation = sL;
		length = l;
	}
	
	/**
	 * @return the length
	 */
	public float getLength() {
		return length;
	}

	/**
	 * @return the startLocation
	 */
	public int getStartLocation() {
		return startLocation;
	}

	@Override
	public int compareTo(Seed o) {
		if(startLocation < o.getStartLocation()) {
			return -1;
		} else if (startLocation > o.getStartLocation()) {
			return 1;
		}
		return 0;
	}
}
