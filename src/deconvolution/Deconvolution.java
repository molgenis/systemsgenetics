package deconvolution;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import org.renjin.eval.EvalException;
import org.renjin.sexp.Vector;

public class Deconvolution
{
  public static void main(String[] args)
    throws Exception
  {
    String expressionTablePath = args[0];
    String genotypeTablePath = args[1];
    String cellcountTablePath = args[2];
    
    ScriptEngineManager manager = new ScriptEngineManager();
    
    ScriptEngine engine = manager.getEngineByName("Renjin");
    if (engine == null) {
      throw new RuntimeException("Renjin Script Engine not found on the classpath.");
    }
    engine.eval("df <- data.frame(x=1:10, y=(1:10)+rnorm(n=10))");
    engine.eval("print(df)");
    engine.eval("print(lm(y ~ x, df))");
    List<List<String>> expressionTable = readTabDelimited(expressionTablePath, Boolean.valueOf(false));
    List<List<String>> genotypeTable = readTabDelimited(genotypeTablePath, Boolean.valueOf(false));
    List<List<String>> cellcountTable = readTabDelimited(cellcountTablePath, Boolean.valueOf(true));
    
    BufferedReader buf = new BufferedReader(new FileReader(cellcountTablePath));
    
    List<String> celltypes = new ArrayList(Arrays.asList(buf.readLine().replace("\n", "").split("\t")));
    celltypes.removeAll(Arrays.asList(new String[] { "", null }));
    buf.close();
    List<String> geneNames = (List)expressionTable.get(0);
    for (int i = 1; i < expressionTable.size(); i++)
    {
      List<String> expressionVectorString = (List)expressionTable.get(i);
      List<Double> expressionVector = StringListToDouble(expressionVectorString);
      List<String> genotypeVectorString = (List)genotypeTable.get(i);
      List<Double> genotypeVector = StringListToDouble(genotypeVectorString);
      try
      {
        engine.put("expressionVector", expressionVector);
        engine.put("genotypeVector", genotypeVector);
        engine.put("geneNames", geneNames);
        engine.put("cellcountTable", cellcountTable);
        engine.put("cellTypes", celltypes);
        engine.eval(new FileReader("deconvolution/src/deconvolution/deconFunctions.R"));
      }
      catch (EvalException e)
      {
        Vector condition = (Vector)e.getCondition();
        
        String msg = condition.getElementAsString(0);
        System.out.println(msg);
      }
    }
  }
  
  public static List<List<String>> readTabDelimited(String filepath, Boolean header)
  {
    List<List<String>> allColumns = new ArrayList();
    try
    {
      BufferedReader buf = new BufferedReader(new FileReader(filepath));
      String lineJustFetched = null;
      for (;;)
      {
        lineJustFetched = buf.readLine();
        if (lineJustFetched == null) {
          break;
        }
        lineJustFetched = lineJustFetched.replace("\n", "").replace("\r", "");
        for (int i = 0; i < lineJustFetched.split("\t").length; i++) {
          try
          {
            ((List)allColumns.get(i)).add(lineJustFetched.split("\t")[i]);
          }
          catch (IndexOutOfBoundsException e)
          {
            List<String> newColumn = new ArrayList();
            newColumn.add(lineJustFetched.split("\t")[i]);
            allColumns.add(newColumn);
          }
        }
      }
      buf.close();
    }
    catch (Exception e)
    {
      e.printStackTrace();
    }
    return allColumns;
  }
  
  public static List<Integer> range(int begin, final int end)
  {
    new AbstractList()
    {
      public Integer get(int index)
      {
        return Integer.valueOf(this.val$begin + index);
      }
      
      public int size()
      {
        return end - this.val$begin;
      }
    };
  }
  
  public static List<Double> StringListToDouble(List<String> ArrayToConvert)
  {
    List<Double> doubles = new ArrayList();
    for (String string : ArrayToConvert) {
      doubles.add(Double.valueOf(Double.parseDouble(string)));
    }
    return doubles;
  }
}
