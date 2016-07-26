package net.joshuahughes.jnablaslapack.pointer;

import java.util.ArrayList;

import com.sun.jna.Structure;

public class STRUCTURE extends Structure{
	@Override
	protected ArrayList<String> getFieldOrder() {
		ArrayList<String> list = new ArrayList<String>();
	     list.add("re");
	     list.add("im");
	     return list;
	}
	public static class FloatStructure extends STRUCTURE
	{
		public float re;
		public float im;
	}
	public static class DoubleStructure extends STRUCTURE
	{
		public double re;
		public double im;
	}
}
