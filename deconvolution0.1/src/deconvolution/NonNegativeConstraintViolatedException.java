package deconvolution;

public class NonNegativeConstraintViolatedException extends Exception {
	private static final long serialVersionUID = 476078491106110717L;

	public NonNegativeConstraintViolatedException() {
        super();
    }


    public NonNegativeConstraintViolatedException(String message) {
        super(message);
    }

    public NonNegativeConstraintViolatedException(String message, Throwable cause) {
        super(message, cause);
    }

    public NonNegativeConstraintViolatedException(Throwable cause) {
        super(cause);
    }

     protected NonNegativeConstraintViolatedException(String message, Throwable cause,
                        boolean enableSuppression,
                        boolean writableStackTrace) {
        super(message, cause, enableSuppression, writableStackTrace);
    }
}
