package net.joshuahughes.jnablaslapack;

import com.sun.jna.ptr.ByReference;

public class CHARACTER extends ByReference{
	public CHARACTER() {
		super(Character.BYTES);
	}
	public CHARACTER(char character)
	{
		this();
		this.setValue(character);
	}
	public char getValue()
	{
		return getPointer().getChar(0);
	}
	public void setValue(char value)
	{
		getPointer().setChar(0,value);
	}
}
