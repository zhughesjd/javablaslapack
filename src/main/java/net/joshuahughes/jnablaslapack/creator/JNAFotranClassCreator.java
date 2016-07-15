package net.joshuahughes.jnablaslapack.creator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;

import net.joshuahughes.jnablaslapack.CHARACTER;
import net.joshuahughes.jnablaslapack.DOUBLE;
import net.joshuahughes.jnablaslapack.INTEGER;
import net.joshuahughes.jnablaslapack.REAL;

import com.sun.jna.Library;
import com.sun.jna.Native;

import fortran.ofp.FrontEnd;

public class JNAFotranClassCreator {
	public static BufferedWriter writer;
	public static void main(String[] args) throws Exception {
		String dirpath  = "C:/Users/hughes/Desktop/lapack-3.6.1/SRC";
		File validFile = File.createTempFile("fortran", ".f");
		File[] files = new File(dirpath).listFiles();
		String[] packageList = JNAFotranClassCreator.class.getPackage().getName().split("\\.");
		packageList = Arrays.copyOf(packageList, packageList.length-1);
		
		String className = (dirpath.toLowerCase().contains("blas")?"Blas":"Lapack");

		String classPath = System.getProperty("user.dir")+"\\src\\main\\java\\"+String.join("\\",packageList)+"\\"+className+".java";
		writer = new BufferedWriter(new FileWriter(classPath));
		writer.write("package "+String.join(".", packageList)+";\n\n");
		writer.write("import "+Library.class.getName()+";\n");
		writer.write("import "+Native.class.getName()+";\n");
		for(Class<?> clazz : new Class[]{CHARACTER.class,REAL.class,INTEGER.class,DOUBLE.class})
			writer.write("import "+clazz.getName()+";\n");

		writer.write("\npublic interface "+className+" extends Library\n{\n\n");
		writer.write("\tpublic static "+className+" instance = ("+className+") Native.loadLibrary(\"lib"+className.toLowerCase()+"\","+className+".class);\n\n");
		for(File file : files)
		{
			if(file.isFile() && file.getName().endsWith(".f"))
			{
				ArrayList<String> initialComments = new ArrayList<>();
				BufferedReader br = new BufferedReader(new FileReader(file));
				BufferedWriter bw = new BufferedWriter(new FileWriter(validFile));
				boolean initialCommentMode = true;
				String line;
				while((line = br.readLine())!=null)
				{
					if(line.startsWith("*") && initialCommentMode)
						initialComments.add(line);
					else
						initialCommentMode = false;
					bw.write(
							line
							.replace("COMPLEX*16", "DOUBLECOMPLEX")
							+"\n"
							);
				}
				br.close();
				bw.close();
				writer.write("/**\n");
				for(String cLine : initialComments)
					writer.write(cLine+"<br>\n");
				writer.write("*/\n");
				FrontEnd.main(new String[]{"--class",JNAAction.class.getName(),validFile.getAbsolutePath()});
			}
		}
		writer.write("\n}");
		writer.close();
	}
}
