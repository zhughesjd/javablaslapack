package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackXE extends Library
{

	public static LapackXE instance = (LapackXE) Native.loadLibrary("liblapack",LapackXE.class);

	public void xerbla_(CHARACTER SRNAME,INTEGER INFO);
	public void xerbla_array_(char[] SRNAME_ARRAY,INTEGER SRNAME_LEN,INTEGER INFO);

}