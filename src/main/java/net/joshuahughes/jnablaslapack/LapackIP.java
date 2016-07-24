package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackIP extends Library
{

	public static LapackIP instance = (LapackIP) Native.loadLibrary("liblapack",LapackIP.class);

	public int iparmq_(INTEGER ISPEC,CHARACTER NAME,CHARACTER OPTS,INTEGER N,INTEGER ILO,INTEGER IHI,INTEGER LWORK);

}