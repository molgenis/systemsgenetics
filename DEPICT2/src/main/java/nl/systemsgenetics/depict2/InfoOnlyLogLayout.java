/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import org.apache.log4j.Layout;
import org.apache.log4j.Level;
import org.apache.log4j.spi.LoggingEvent;

/**
 *
 * @author patri
 */
public class InfoOnlyLogLayout extends Layout {
	
	StringBuffer sbuf = new StringBuffer(128);

	@Override
	public String format(LoggingEvent event) {
		if (event.getLevel() == Level.INFO) {
			sbuf.setLength(0);
			sbuf.append(event.getRenderedMessage());
			sbuf.append(LINE_SEP);
			return sbuf.toString();
		} else {
			return "";
		}
	}

	@Override
	public boolean ignoresThrowable() {
		return true;
	}

	@Override
	public void activateOptions() {
	}

}
