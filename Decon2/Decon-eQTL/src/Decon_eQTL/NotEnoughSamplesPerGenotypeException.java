package Decon_eQTL;

public class NotEnoughSamplesPerGenotypeException extends Exception {
	private static final long serialVersionUID = 7189762595398313333L;

	public NotEnoughSamplesPerGenotypeException() {
        super();
    }


    public NotEnoughSamplesPerGenotypeException(String message) {
        super(message);
    }

    public NotEnoughSamplesPerGenotypeException(String message, Throwable cause) {
        super(message, cause);
    }

    public NotEnoughSamplesPerGenotypeException(Throwable cause) {
        super(cause);
    }

     protected NotEnoughSamplesPerGenotypeException(String message, Throwable cause,
                        boolean enableSuppression,
                        boolean writableStackTrace) {
        super(message, cause, enableSuppression, writableStackTrace);
    }
}
