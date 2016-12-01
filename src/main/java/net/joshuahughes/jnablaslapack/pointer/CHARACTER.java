package net.joshuahughes.jnablaslapack.pointer;

import com.sun.jna.ptr.ByReference;
/**
 * Blas/Lapack code bases are Fortran 90 code bases:
 * 
 * 		http://www.netlib.org/lapack/
 * 		
 * Characters types are represented by a single byte in fortran 90: 
 * 	
 *		http://www.cs.uwm.edu/~cs151/Bacon/Lecture/HTML/ch06s09.html
 *		http://earth.uni-muenster.de/~joergs/doc/f90/unix-um/dfum_034.html
 *		
 * 
 * @author hughes
 *
 */
public class CHARACTER extends ByReference{
	public CHARACTER() {
		super(Byte.BYTES);
	}
	public CHARACTER(char character)
	{
		this();
		this.setValue(character);
	}
	public char getValue()
	{
		return (char) getPointer().getByte(0);
	}
	public void setValue(char value)
	{
		getPointer().setByte(0,(byte) value);
	}
}
