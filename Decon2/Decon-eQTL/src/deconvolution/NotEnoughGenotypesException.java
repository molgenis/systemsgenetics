package deconvolution;

public class NotEnoughGenotypesException extends Exception {
	private static final long serialVersionUID = -2137984747579058008L;

	public NotEnoughGenotypesException() {
        super();
    }


    public NotEnoughGenotypesException(String message) {
        super(message);
    }

    public NotEnoughGenotypesException(String message, Throwable cause) {
        super(message, cause);
    }

    public NotEnoughGenotypesException(Throwable cause) {
        super(cause);
    }

     protected NotEnoughGenotypesException(String message, Throwable cause,
                        boolean enableSuppression,
                        boolean writableStackTrace) {
        super(message, cause, enableSuppression, writableStackTrace);
    }
}
