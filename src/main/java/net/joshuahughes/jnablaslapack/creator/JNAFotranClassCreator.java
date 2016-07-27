package net.joshuahughes.jnablaslapack.creator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;

import com.sun.jna.Library;
import com.sun.jna.Native;

import fortran.ofp.FrontEnd;
import net.joshuahughes.jnablaslapack.pointer.REAL;

public class JNAFotranClassCreator {
	public static BufferedWriter writer;
	public static void main(String[] args) throws Exception {
		String dirpath  = args[0];
		boolean insertComments = args.length>=2 && Boolean.parseBoolean(args[1]);
		File validFile = File.createTempFile("fortran", ".f");
		File[] files = new File(dirpath).listFiles();
		String[] packageList = JNAFotranClassCreator.class.getPackage().getName().split("\\.");
		boolean isLapack = !dirpath.toLowerCase().contains("blas");
		int appendageLength = isLapack?2:0;
		String className = isLapack?"Lapack":"Blas";
		packageList = Arrays.copyOf(packageList, packageList.length-1);
		String priorChars = null;
		for(File file : files)
		{
			if(file.isFile() && file.getName().endsWith(".f"))
			{
				String initChars = file.getName().substring(0,appendageLength).toUpperCase();
				if(!initChars.equals(priorChars))
				{
					if(writer!=null)
					{
						writer.write("\n}");
						writer.close();
					}
					writer = create(packageList,className,initChars);
					priorChars = initChars;
				}
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
				if(insertComments)
				{
					writer.write("/**\n");
					for(String cLine : initialComments)
						writer.write(cLine+"<br>\n");
					writer.write("*/\n");
				}
				FrontEnd.main(new String[]{"--class",JNAAction.class.getName(),validFile.getAbsolutePath()});
			}
		}
		writer.write("\n}");
		writer.close();
	}
	private static BufferedWriter create(String[] packageList,String libname,String appendage) throws Exception {
		String className = libname+appendage;
		String classPath = System.getProperty("user.dir")+"\\src\\main\\java\\"+String.join("\\",packageList)+"\\"+className+".java";
		BufferedWriter writer = new BufferedWriter(new FileWriter(classPath));
		writer.write("package "+String.join(".", packageList)+";\n\n");
		for(String importString : new String[]{Library.class.getName(),Native.class.getName(),REAL.class.getPackage().getName()+".*"})
			writer.write("import "+importString+";\n");
		writer.write("\npublic interface "+className+" extends Library\n{\n\n");
		writer.write("\tpublic static "+className+" instance = ("+className+") Native.loadLibrary(\"lib"+libname.toLowerCase()+"\","+className+".class);\n\n");
		return writer;
	}
}
