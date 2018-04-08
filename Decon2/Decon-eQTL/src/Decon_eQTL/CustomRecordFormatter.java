	package Decon_eQTL;
	
	import java.io.PrintWriter;
	import java.io.StringWriter;
	import java.util.logging.SimpleFormatter;
	import java.util.logging.LogRecord;
	
	class CustomRecordFormatter extends SimpleFormatter {
	    @Override
	    public String format(final LogRecord r) {
	        StringBuilder sb = new StringBuilder();
	        sb.append(formatMessage(r)).append(System.getProperty("line.separator"));
	        if (null != r.getThrown()) {
	            sb.append("Throwable occurred: "); //$NON-NLS-1$
	            Throwable t = r.getThrown();
	            PrintWriter pw = null;
	            try {
	                StringWriter sw = new StringWriter();
	                pw = new PrintWriter(sw);
	                t.printStackTrace(pw);
	                sb.append(sw.toString());
	            } finally {
	                if (pw != null) {
	                    try {
	                        pw.close();
	                    } catch (Exception e) {
	                        // ignore
	                    }
	                }
	            }
	        }
	        return sb.toString();
	    }
	}