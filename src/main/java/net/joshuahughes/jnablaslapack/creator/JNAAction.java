package net.joshuahughes.jnablaslapack.creator;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import net.joshuahughes.jnablaslapack.pointer.CHARACTER;
import net.joshuahughes.jnablaslapack.pointer.COMPLEX;
import net.joshuahughes.jnablaslapack.pointer.COMPLEX16;
import net.joshuahughes.jnablaslapack.pointer.COMPLEX8;
import net.joshuahughes.jnablaslapack.pointer.DOUBLE;
import net.joshuahughes.jnablaslapack.pointer.DOUBLECOMPLEX;
import net.joshuahughes.jnablaslapack.pointer.INTEGER;
import net.joshuahughes.jnablaslapack.pointer.LOGICAL;
import net.joshuahughes.jnablaslapack.pointer.REAL;

import org.antlr.runtime.Token;
import org.antlr.v4.runtime.misc.Pair;

import fortran.ofp.parser.java.FortranParserActionNull;
import fortran.ofp.parser.java.IFortranParser;

public class JNAAction extends FortranParserActionNull {
	LinkedHashMap<String,Pair<Class<?>,Boolean>> map = new LinkedHashMap<>();
	boolean endOfArgs = false;
	Class<?> latestType;
	LinkedHashSet<String> argNameList = new LinkedHashSet<>();
	private String subprogramName;
	private Class<?> subprogramReturnType;
	boolean hasKindSelector = false;
	public JNAAction(String[] arg0, IFortranParser arg1, String arg2) {
		super(arg0, arg1, arg2);
	}
	public void intrinsic_type_spec(Token keyword1, Token keyword2,
			int type, boolean hasKindSelector)
	{
		latestType = convert(keyword1.getText());
		this.hasKindSelector = hasKindSelector;
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
			map.put(candidate, new Pair<>(latestType, hasArraySpec||hasKindSelector));
	}
	public void end() {
		try {
			JNAFotranClassCreator.writer.write("\tpublic "+returnByValue(subprogramReturnType).getSimpleName()+" "+subprogramName.toLowerCase()+"_(");
			
			boolean initial = true;
			for(String name : argNameList)
				if(map.containsKey(name)){
					Pair<Class<?>,Boolean> pair = map.get(name);
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
	public static Class<?> convert(String string)
	{
		for(Class<?> clazz : new Class<?>[]{LOGICAL.class,INTEGER.class,CHARACTER.class,REAL.class,DOUBLE.class,COMPLEX.class,DOUBLECOMPLEX.class})
			if(clazz.getSimpleName().equals(string))
				return clazz;
		return null;
	}
	public static Class<?> returnByValue(Class<?> fortranType)
	{
		if(LOGICAL.class.equals(fortranType)) return boolean.class;
		if(INTEGER.class.equals(fortranType)) return int.class;
		if(CHARACTER.class.equals(fortranType)) return char.class;
		if(REAL.class.equals(fortranType)) return float.class;
		if(DOUBLE.class.equals(fortranType)) return double.class;
		if(COMPLEX.class.equals(fortranType)) return COMPLEX8.class;
		if(DOUBLECOMPLEX.class.equals(fortranType)) return COMPLEX16.class;
		return void.class;
	}
	public static String passByReference(Class<?> fortranType, boolean isArray)
	{
		Class<?> ref = void.class;
		if(COMPLEX.class.equals(fortranType) || (REAL.class.equals(fortranType) && isArray))
		{
			isArray = true;
			ref = float.class;
		}
		else if(DOUBLECOMPLEX.class.equals(fortranType) || (DOUBLE.class.equals(fortranType) && isArray))
		{
			isArray = true;
			ref = double.class;
		}
		else if(isArray)
		{
			if(LOGICAL.class.equals(fortranType)) ref = boolean.class;
			if(INTEGER.class.equals(fortranType)) ref = int.class;
			if(CHARACTER.class.equals(fortranType)) ref = byte.class;
			if(REAL.class.equals(fortranType)) ref = float.class;
			if(DOUBLE.class.equals(fortranType)) ref = double.class;
		}
		else
		{
			for(Class<?> clazz : new Class[]{CHARACTER.class,REAL.class,INTEGER.class,LOGICAL.class,DOUBLE.class})
				if(clazz.equals(fortranType)) 
					ref = clazz;
		}
		return ref.getSimpleName()+(isArray?"[]":"");
	}
}
