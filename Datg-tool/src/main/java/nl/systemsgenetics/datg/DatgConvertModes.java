package nl.systemsgenetics.datg;

public enum DatgConvertModes {

    TXT_2_DATG("Convert tab separated .txt or .txt.gz files to a .datg file. The first row is expected to be header information and the first column should contain row names"),
    DATG_2_TXT("Convert a datg file to a txt"),
    DAT_2_DATG("Special mode for our old inhouse .dat files, not recommend for public use."),
    INSPECT("Get meta data of DATG file and show some elements"),
    ROW_CONCAT("Concatenate .txt matrices with identical columns into .datg. "),
    UPGRADE("This modes upgrades a .datg file to latest version."),
    TEST("Internal testing");

    private final String description;

    DatgConvertModes(String description) {
        this.description = description;
    }

    public String getDescription() {
        return description;
    }

    public static String getFullDescriptionString() {
        StringBuilder output = new StringBuilder();
        for (DatgConvertModes mode: DatgConvertModes.values()) {
            if(mode == TEST){
                continue;
            }
            output.append(mode.toString());
            output.append(" - ");
            output.append(mode.getDescription());
            output.append("\n");
        }

        return output.toString();
    }
}
