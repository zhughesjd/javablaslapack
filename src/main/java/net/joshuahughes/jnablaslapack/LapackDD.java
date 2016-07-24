package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDD extends Library
{

	public static LapackDD instance = (LapackDD) Native.loadLibrary("liblapack",LapackDD.class);

	public void ddisna_(CHARACTER JOB,INTEGER M,INTEGER N,double[] D,double[] SEP,INTEGER INFO);

}