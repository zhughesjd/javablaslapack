package net.joshuahughes.jnablaslapack.pointer;

import com.sun.jna.ptr.IntByReference;

public class INTEGER extends IntByReference{
	public INTEGER(){}
	public INTEGER(int value)
	{
		setValue(value);
	}
}
