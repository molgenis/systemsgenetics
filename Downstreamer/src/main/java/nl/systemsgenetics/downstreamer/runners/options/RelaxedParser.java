package nl.systemsgenetics.downstreamer.runners.options;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import java.util.ArrayList;
import java.util.List;

/**
 * Parser that doesn't crash on unknown arguments.
 * Full credit to https://stackoverflow.com/questions/33874902/apache-commons-cli-1-3-1-how-to-ignore-unknown-arguments
 */
public class RelaxedParser extends PosixParser {

    @Override
    public CommandLine parse(final Options options, final String[] arguments, boolean stopAtNonOption) throws ParseException {
        final List<String> knownArgs = new ArrayList<>();
        for (int i = 0; i < arguments.length; i++) {
            if (options.hasOption(arguments[i])) {
                knownArgs.add(arguments[i]);
                if (i + 1 < arguments.length && options.getOption(arguments[i]).hasArg()) {
                    knownArgs.add(arguments[i + 1]);
                }
            }
        }
        return super.parse(options, knownArgs.toArray(new String[0]), false);
    }
}