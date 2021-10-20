package betaqtl.junk;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

public class DFTest {


    public static void main(String[] args) {
        double d = -3000000.646876876387;
//        d = 0.001;
//        d = 0.001;
        DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat dfDefault = new DecimalFormat("#.######", symbols);
        DecimalFormat dfPval = new DecimalFormat("#.####E0", symbols);
        System.out.println(dfDefault.format(d));
        System.out.println(d);
    }
}
