package net.joshuahughes.jnablaslapack.pointer;

import com.sun.jna.ptr.FloatByReference;

public class REAL extends FloatByReference{
	public REAL(){}
	public REAL(float value)
	{
		setValue(value);
	}
}
