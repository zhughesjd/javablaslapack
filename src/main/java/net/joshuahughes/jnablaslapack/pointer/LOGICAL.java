package net.joshuahughes.jnablaslapack.pointer;

import com.sun.jna.ptr.ByteByReference;

public class LOGICAL extends ByteByReference{
	public LOGICAL(){}
	public LOGICAL(boolean value)
	{
		setValue((byte) (value?1:0));
	}
}
