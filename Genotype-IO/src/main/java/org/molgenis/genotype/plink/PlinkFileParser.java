package org.molgenis.genotype.plink;

import java.io.Closeable;
import java.nio.charset.Charset;

public interface PlinkFileParser extends Closeable
{
	public static final Charset FILE_ENCODING = Charset.forName("UTF-8");
	public static final char DEFAULT_FIELD_SEPARATOR = ' ';
	public static final String DEFAULT_READ_FIELD_SEPARATORS = " \t";
	public static final String LINE_SEPARATOR = "\n";
}
