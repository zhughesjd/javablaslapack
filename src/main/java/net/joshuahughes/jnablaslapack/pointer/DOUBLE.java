package net.joshuahughes.jnablaslapack.pointer;

import com.sun.jna.ptr.DoubleByReference;

public class DOUBLE extends DoubleByReference{
	public DOUBLE(){}
	public DOUBLE(double value)
	{
		setValue(value);
	}
}
