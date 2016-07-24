package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSD extends Library
{

	public static LapackSD instance = (LapackSD) Native.loadLibrary("liblapack",LapackSD.class);

	public void sdisna_(CHARACTER JOB,INTEGER M,INTEGER N,float[] D,float[] SEP,INTEGER INFO);

}