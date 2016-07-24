package net.joshuahughes.jnablaslapack.creator;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import org.antlr.runtime.Token;
import org.antlr.v4.runtime.misc.Pair;

import fortran.ofp.parser.java.FortranParserActionNull;
import fortran.ofp.parser.java.IFortranParser;
import net.joshuahughes.jnablaslapack.pointer.CHARACTER;
import net.joshuahughes.jnablaslapack.pointer.DOUBLE;
import net.joshuahughes.jnablaslapack.pointer.INTEGER;
import net.joshuahughes.jnablaslapack.pointer.LOGICAL;
import net.joshuahughes.jnablaslapack.pointer.REAL;

public class JNAAction extends FortranParserActionNull {
	LinkedHashMap<String,Pair<String,Boolean>> map = new LinkedHashMap<>();
	boolean endOfArgs = false;
	String latestType;
	LinkedHashSet<String> argNameList = new LinkedHashSet<>();
	private String subprogramName;
	private String subprogramReturnType;
	public JNAAction(String[] arg0, IFortranParser arg1, String arg2) {
		super(arg0, arg1, arg2);
	}
	public void intrinsic_type_spec(Token keyword1, Token keyword2,
			int type, boolean hasKindSelector)
	{
		latestType = keyword1.getText();
	}
	public void generic_name_list_part(Token id)
	{
		if(!endOfArgs)
			argNameList.add(id.getText());
	}
	public void generic_name_list(int count) {
		endOfArgs = true;
	}
	public void function_stmt(Token label, Token keyword, Token name, 
			Token eos, boolean hasGenericNameList, 
			boolean hasSuffix) {
		subprogramName = name.getText();
		subprogramReturnType = latestType;
	}
	public void subroutine_stmt(Token label, Token keyword, Token name, 
			Token eos, boolean hasPrefix, 
			boolean hasDummyArgList, boolean hasBindingSpec,
			boolean hasArgSpecifier) {

		subprogramName = name.getText();
	}
	public void dummy_arg(Token token)
	{
		argNameList.add(token.getText());
	}
	public void entity_decl(Token id, boolean hasArraySpec, boolean hasCoarraySpec,
			boolean hasCharLength, boolean hasInitialization)
	{
		String candidate = id.getText();
		if(argNameList.contains(candidate) && !map.containsKey(candidate))
			map.put(candidate, new Pair<>(latestType,hasArraySpec));
	}
	public void end() {
		try {
			JNAFotranClassCreator.writer.write("\tpublic "+returnByValue(subprogramReturnType)+" "+subprogramName.toLowerCase()+"_(");
			boolean initial = true;
			for(String name : argNameList)
				if(map.containsKey(name)){
					Pair<String,Boolean> pair = map.get(name);
					JNAFotranClassCreator.writer.write((initial?"":",")+passByReference(pair.a,pair.b)+" "+name);
					initial = false;

				}
			JNAFotranClassCreator.writer.write(");\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		subprogramReturnType = null;
		map.clear();
		argNameList.clear();
		endOfArgs = false;
	}
	public void end_subroutine_stmt(Token label, Token keyword1, 
			Token keyword2, Token name, Token eos) {
		end();
	}
	public void end_function_stmt(Token label, Token keyword1, Token keyword2, 
			Token name, Token eos) {
		end();
	}
	public static String returnByValue(String returnType)
	{
		if("LOGICAL".equals(returnType)) return boolean.class.getSimpleName();
		if("INTEGER".equals(returnType)) return int.class.getSimpleName();
		if("CHARACTER".equals(returnType)) return char.class.getSimpleName();
		if("REAL".equals(returnType)) return float.class.getSimpleName();
		if("DOUBLE".equals(returnType)) return double.class.getSimpleName();

		if("COMPLEX".equals(returnType)) return "float[]";
		if("DOUBLECOMPLEX".equals(returnType)) return "double[]";
		return void.class.getSimpleName();
	}
	public static String passByReference(String fortranType, boolean isArray)
	{
		String ref = void.class.getSimpleName();
		if("COMPLEX".equals(fortranType) || ("REAL".equals(fortranType) && isArray))
		{
			isArray = true;
			ref = float.class.getSimpleName();
		}
		else if("DOUBLECOMPLEX".equals(fortranType) || ("DOUBLE".equals(fortranType) && isArray))
		{
			isArray = true;
			ref = double.class.getSimpleName();
		}
		else if(isArray)
		{
			if("LOGICAL".equals(fortranType)) return boolean.class.getSimpleName();
			if("INTEGER".equals(fortranType)) ref = int.class.getSimpleName();
			if("CHARACTER".equals(fortranType)) ref = char.class.getSimpleName();
			if("REAL".equals(fortranType)) ref = float.class.getSimpleName();
			if("DOUBLE".equals(fortranType)) ref = double.class.getSimpleName();
		}
		else
		{
			for(Class<?> clazz : new Class[]{CHARACTER.class,REAL.class,INTEGER.class,LOGICAL.class,DOUBLE.class})
				if(clazz.getSimpleName().equals(fortranType)) 
					ref = clazz.getSimpleName();
		}
		return ref+(isArray?"[]":"");
	}
}
