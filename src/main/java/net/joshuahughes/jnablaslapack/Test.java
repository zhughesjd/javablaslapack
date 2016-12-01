package net.joshuahughes.jnablaslapack;

import java.util.Random;

import net.joshuahughes.jnablaslapack.pointer.CHARACTER;
import net.joshuahughes.jnablaslapack.pointer.COMPLEX16;
import net.joshuahughes.jnablaslapack.pointer.DOUBLE;
import net.joshuahughes.jnablaslapack.pointer.INTEGER;



public class Test {
	public static void main(String[] args)
	{
		CHARACTER c = new CHARACTER('c');
		System.out.println(c.getValue());

		int numIter = 1;
		int dim = 500;
		DOUBLE ALPHA = new DOUBLE(1.0);
		DOUBLE BETA  = new DOUBLE(0.0);
		INTEGER DIM1 = new INTEGER(dim);

		Random rand = new Random();
		double[] A = new double[dim*dim];
		double[] B = new double[dim*dim];
		for (int i = 0; i < A.length; ++i)
		{
			A[i] = rand.nextDouble();
			B[i] = rand.nextDouble();
		}
		double[] C = new double[A.length];

		Blas blas = Blas.instance;
		CHARACTER ch_n  = new CHARACTER('n');
		INTEGER dimm = new INTEGER(dim);

		long time = System.currentTimeMillis();
		for (int i = 0; i < numIter; ++i)
		{
			blas.dgemm_(ch_n, ch_n, dimm, dimm, dimm, ALPHA, A, DIM1, B, DIM1, BETA, C, DIM1);
		}
		System.err.println("BLAS time: " + ((System.currentTimeMillis()-time)/1000.0) + " sec");

		COMPLEX16 dc = LapackZL.instance.zladiv_(new double[]{14,32},new double[]{7,8});
		System.out.println(dc.re+":"+dc.im);
		System.out.println((byte)LapackCH.instance.chla_transtype_(new INTEGER(0)));
	}
}
