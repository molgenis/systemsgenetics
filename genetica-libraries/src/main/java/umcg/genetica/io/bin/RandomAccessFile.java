/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.bin;

/*
 * Copyright 1997, University Corporation for Atmospheric Research
 * See COPYRIGHT file for copying and redistribution conditions.
 *
 * RandomAccessFile.java.  By Russ Rew, based on
 * BufferedRandomAccessFile by Alex McManus, based on Sun's source code
 * for java.io.RandomAccessFile.  For Alex McManus version from which
 * this derives, see his <a href="http://www.aber.ac.uk/~agm/Java.html">
 * Freeware Java Classes</a>.  
 */
//DJ package ucar.netcdf;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.DataOutput;
import java.io.EOFException;
import java.io.File;
import java.io.FileDescriptor;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UTFDataFormatException;
import java.util.Date;
import java.util.Random;

/**
 * A buffered drop-in replacement for java.io.RandomAccessFile. Instances of
 * this class realise substantial speed increases over java.io.RandomAccessFile
 * through the use of buffering. This is a subclass of Object, as it was not
 * possible to subclass java.io.RandomAccessFile because many of the methods are
 * final. However, if it is necessary to use RandomAccessFile and
 * java.io.RandomAccessFile interchangeably, both classes implement the
 * DataInput and DataOutput interfaces.
 *
 * @author Alex McManus
 * @author Russ Rew
 * @version $Id: RandomAccessFile.java,v 1.2 1998/10/20 18:38:03 russ Exp $
 * @see DataInput
 * @see DataOutput
 * @see java.io.RandomAccessFile
 */
public class RandomAccessFile extends Object
        implements DataInput, DataOutput {

    /**
     * Read from the file. This is always implied.
     */
    public static final int READ = 1;
    /**
     * Write to the file.
     */
    public static final int WRITE = 2;
    /**
     * Create the file rather than overwriting it. This is ignored if mode is
     * not also WRITE.
     */
    public static final int CREATE = 4;
    /**
     * The default buffer size, in bytes.
     */
    protected static final int defaultBufferSize = 4096;
    /**
     * The underlying java.io.RandomAccessFile.
     */
    protected java.io.RandomAccessFile file;
    /**
     * The offset in bytes from the file start, of the next read or write
     * operation.
     */
    protected long filePosition;
    /**
     * The buffer used to load the data.
     */
    protected byte buffer[];
    /**
     * The offset in bytes of the start of the buffer, from the start of the
     * file.
     */
    protected long bufferStart;
    /**
     * The offset in bytes of the end of the data in the buffer, from the start
     * of the file. This can be calculated from
     * <code>bufferStart + dataSize</code>, but it is cached to speed up the
     * read( ) method.
     */
    protected long dataEnd;
    /**
     * The size of the data stored in the buffer, in bytes. This may be less
     * than the size of the buffer.
     */
    protected int dataSize;
    /**
     * True if we are at the end of the file.
     */
    protected boolean endOfFile;
    /**
     * The access mode of the file. This is a logical OR of READ, WRITE and
     * CREATE.
     */
    protected int mode;
    /**
     * True if the data in the buffer has been modified.
     */
    boolean bufferModified = false;

    /**
     * Create a new buffered random-access file with a default buffer size. Note
     * that the mode CREATE implies WRITE.
     *
     * @param filename the name of the file.
     * @param mode how the file is to be opened. This may be a combination
     * (logical OR) of CREATE, WRITE, and READ.
     * @exception IOException if an I/O error occurrs.
     * @exception SecurityException if a security manager exists, its checkRead
     * method is called with the name argument to see if the application is
     * allowed read access to the file. If the mode argument is WRITE, its
     * checkWrite method also is called with the name argument to see if the
     * application is allowed write access to the file. Either of these may
     * result in a security exception.
     */
    public RandomAccessFile(String filename, int mode)
            throws IOException {
        this(filename, mode, defaultBufferSize);
    }

    /**
     * Creates a random access file stream to read from, and optionally to write
     * to, a file with the specified name. <p> The mode argument must either be
     * equal to
     * <code>"r"</code> or
     * <code>"rw"</code>, indicating that the file is to be opened for input
     * only or for both input and output, respectively. If the mode is
     * <code>"rw"</code> and the file does not exist, then an attempt is made to
     * create it.
     *
     * @param name the system-dependent filename.
     * @param mode the access mode.
     * @exception IllegalArgumentException if the mode argument is not equal        to <code>"r"</code> or to <code>"rw"</code>.
     * @exception IOException if an I/O error occurs.
     * @exception SecurityException if a security manager exists, its
     * <code>checkRead</code> method is called with the name argument to see if
     * the application is allowed read access to the file. If the mode argument
     * is equal to <code>"rw"</code>, its <code>checkWrite</code> method also is
     * called with the name argument to see if the application is allowed write
     * access to the file. Either of these may result in a security exception.
     * @see java.lang.SecurityException
     * @see java.lang.SecurityManager#checkRead(java.lang.String)
     */
    public RandomAccessFile(String filename, String modeString)
            throws IOException {

        this(filename,
                modeString.equals("r") ? READ
                : (modeString.equals("rw") ? WRITE | READ : 0),
                defaultBufferSize);
    }

    /**
     * Creates a random access file stream to read from, and optionally to write
     * to, the file specified by the
     * <code>File</code> argument. A new {@link FileDescriptor} object is
     * created to represent this file connection. <p> The mode argument must
     * either be equal to
     * <code>"r"</code> or
     * <code>"rw"</code>, indicating that the file is to be opened for input
     * only or for both input and output, respectively. The write methods on
     * this object will always throw an
     * <code>IOException</code> if the file is opened with a mode of
     * <code>"r"</code>. If the mode is
     * <code>"rw"</code> and the file does not exist, then an attempt is made to
     * create it.
     *
     * @param file the file object.
     * @param mode the access mode.
     * @exception IllegalArgumentException if the mode argument is not equal        to <code>"r"</code> or to <code>"rw"</code>.
     * @exception IOException if an I/O error occurs.
     * @exception SecurityException if a security manager exists, its
     * <code>checkRead</code> method is called with the pathname of
     * the <code>File</code> argument to see if the application is allowed read
     * access to the file. If the mode argument is equal to <code>"rw"</code>,
     * its <code>checkWrite</code> method also is called with the pathname to
     * see if the application is allowed write access to the file.
     * @see java.io.File#getPath()
     * @see java.lang.SecurityManager#checkRead(java.lang.String)
     */
    public RandomAccessFile(File file, String modeString) throws IOException {
        this(file.getPath(), modeString);
    }

    /**
     * Create a new buffered random-access file with a specified buffer size.
     * Note that the mode CREATE implies WRITE, and the READ is always implied.
     *
     * @param filename the name of the file.
     * @param mode how the file is to be opened. This may be a combination
     * (logical OR) of CREATE, WRITE, and READ.
     * @param bufferSize the size of the temporary buffer, in bytes.
     * @exception IOException if an I/O error occurrs.
     * @exception SecurityException if a security manager exists, its checkRead
     * method is called with the name argument to see if the application is
     * allowed read access to the file. If the mode argument is WRITE, its
     * checkWrite method also is called with the name argument to see if the
     * application is allowed write access to the file. Either of these may
     * result in a security exception.
     */
    public RandomAccessFile(File file, String modeString, int bufferSize) throws IOException {
        this(file.getPath(),
        modeString.equals("r") ? READ
                : (modeString.equals("rw") ? WRITE | READ : 0),bufferSize);
    }
            
    
    /**
     * Create a new buffered random-access file with a specified buffer size.
     * Note that the mode CREATE implies WRITE, and the READ is always implied.
     *
     * @param filename the name of the file.
     * @param mode how the file is to be opened. This may be a combination
     * (logical OR) of CREATE, WRITE, and READ.
     * @param bufferSize the size of the temporary buffer, in bytes.
     * @exception IOException if an I/O error occurrs.
     * @exception SecurityException if a security manager exists, its checkRead
     * method is called with the name argument to see if the application is
     * allowed read access to the file. If the mode argument is WRITE, its
     * checkWrite method also is called with the name argument to see if the
     * application is allowed write access to the file. Either of these may
     * result in a security exception.
     */
    public RandomAccessFile(String filename, int mode, int bufferSize)
            throws IOException {
        this.mode = mode;

        // If we are CREATEing a file, we must also WRITE. READ is always
        // set.
        mode |= READ;
        if ((this.mode & CREATE) > 0) {
            this.mode |= WRITE;
        }

        // To match java.io.RandomAccessFile semantics, if we want to write 
        // a nonexistant file, create it first (even if CREATE not set)
        File checkfile = new File(filename);
        if ((this.mode & WRITE) > 0 && !checkfile.exists()) {
            mode |= CREATE;
        }

        // If a new file is to be created, delete any existing file with the same name.
        if ((this.mode & CREATE) > 0) {
            if (checkfile.exists()) {
                if (!checkfile.delete()) {
                    throw new IOException("Failed to delete " + filename);
                }
            }
        }

        // If only reading, check that the file exists.
        if (this.mode == READ && !(new File(filename)).exists()) {
            throw new FileNotFoundException(filename);
        }

        // Create the underlying file object.
        String modeString = ((this.mode & WRITE) > 0) ? "rw" : "r";
        file = new java.io.RandomAccessFile(filename, modeString);

        // Initialise the buffer;
        bufferStart = 0;
        dataEnd = 0;
        dataSize = 0;
        filePosition = 0;
        buffer = new byte[bufferSize];
        endOfFile = false;
    }

    /**
     * Close the file, and release any associated system resources.
     *
     * @exception IOException if an I/O error occurrs.
     */
    public void close()
            throws IOException {

        // If we are writing and the buffer has been modified, flush the contents
        // of the buffer.
        if ((mode | WRITE) > 0 && bufferModified) {
            file.seek(bufferStart);
            file.write(buffer, 0, (int) dataSize);
        }

        // Close the underlying file object.
        file.close();
    }

    /**
     * Set the position in the file for the next read or write.
     *
     * @param pos the offset (in bytes) from the start of the file.
     * @exception IOException if an I/O error occurrs.
     */
    public final void seek(long pos)
            throws IOException {

        // If the seek is into the buffer, just update the file pointer.
        if (pos >= bufferStart && pos < dataEnd) {
            filePosition = pos;
            return;
        }

        // If the current buffer is modified, write it to disk.
        if (bufferModified) {
            flush();
        }

        // Move to the position on the disk.
        file.seek(pos);
        filePosition = file.getFilePointer();
        bufferStart = filePosition;

        // Fill the buffer from the disk.
        dataSize = file.read(buffer);
        if (dataSize < 0) {
            dataSize = 0;
            endOfFile = true;
        } else {
            endOfFile = false;
        }

        // Cache the position of the buffer end.
        dataEnd = bufferStart + dataSize;
    }

    /**
     * Returns the current position in the file, where the next read or write
     * will occur.
     *
     * @return the offset from the start of the file in bytes.
     * @exception IOException if an I/O error occurrs.
     */
    public final long getFilePointer()
            throws IOException {
        return filePosition;
    }

    /**
     * Get the length of the file. The data in the buffer (which may not have
     * been written the disk yet) is taken into account.
     *
     * @return the length of the file in bytes.
     * @exception IOException if an I/O error occurrs.
     */
    public long length()
            throws IOException {
        long fileLength = file.length();
        if (fileLength < dataEnd) {
            return dataEnd;
        } else {
            return fileLength;
        }
    }

    /**
     * Returns the opaque file descriptor object associated with this file.
     *
     * @return the file descriptor object associated with this file.
     * @exception IOException if an I/O error occurs.
     */
    public final FileDescriptor getFD()
            throws IOException {
        return file.getFD();
    }

    /**
     * Copy the contents of the buffer to the disk.
     *
     * @exception IOException if an I/O error occurrs.
     */
    public void flush()
            throws IOException {
        file.seek(bufferStart);
        file.write(buffer, 0, dataSize);
        bufferModified = false;
    }

    //
    // Read primitives.
    //
    /**
     * Read a byte of data from the file, blocking until data is available.
     *
     * @return the next byte of data, or -1 if the end of the file is reached.
     * @exception IOException if an I/O error occurrs.
     */
    public final int read()
            throws IOException {

        // If the file position is within the data, return the byte...
        if (filePosition < dataEnd) {
            return (int) (buffer[(int) (filePosition++ - bufferStart)] & 0xff);

            // ...or should we indicate EOF...
        } else if (endOfFile) {
            return -1;

            // ...or seek to fill the buffer, and try again.
        } else {
            seek(filePosition);
            return read();
        }
    }

    /**
     * Read up to
     * <code>len</code> bytes into an array, at a specified offset. This will
     * block until at least one byte has been read.
     *
     * @param b the byte array to receive the bytes.
     * @param off the offset in the array where copying will start.
     * @param len the number of bytes to copy.
     * @return the actual number of bytes read, or -1 if there is not more data
     * due to the end of the file being reached.
     * @exception IOException if an I/O error occurrs.
     */
    private int readBytes(byte b[], int off, int len)
            throws IOException {

        // Check for end of file.
        if (endOfFile) {
            return -1;
        }

        // See how many bytes are available in the buffer - if none,
        // seek to the file position to update the buffer and try again.
        int bytesAvailable = (int) (dataEnd - filePosition);
        if (bytesAvailable < 1) {
            seek(filePosition);
            return readBytes(b, off, len);
        }

        // Copy as much as we can.
        int copyLength = (bytesAvailable >= len) ? len : bytesAvailable;
        System.arraycopy(buffer, (int) (filePosition - bufferStart),
                b, off, copyLength);
        filePosition += copyLength;

        // If there is more to copy...
        if (copyLength < len) {
            int extraCopy = len - copyLength;

            // If the amount remaining is more than a buffer's length, read it
            // directly from the file.
            if (extraCopy > buffer.length) {
                file.seek(filePosition);
                extraCopy = file.read(b, off + copyLength, len - copyLength);

                // ...or read a new buffer full, and copy as much as possible...
            } else {
                seek(filePosition);
                if (!endOfFile) {
                    extraCopy = (extraCopy > dataSize) ? dataSize : extraCopy;
                    System.arraycopy(buffer, 0, b, off + copyLength, extraCopy);
                } else {
                    extraCopy = -1;
                }
            }

            // If we did manage to copy any more, update the file position and
            // return the amount copied.
            if (extraCopy > 0) {
                filePosition += extraCopy;
                return copyLength + extraCopy;
            }
        }

        // Return the amount copied.
        return copyLength;
    }

    /**
     * Read up to
     * <code>len</code> bytes into an array, at a specified offset. This will
     * block until at least one byte has been read.
     *
     * @param b the byte array to receive the bytes.
     * @param off the offset in the array where copying will start.
     * @param len the number of bytes to copy.
     * @return the actual number of bytes read, or -1 if there is not more data
     * due to the end of the file being reached.
     * @exception IOException if an I/O error occurrs.
     */
    public int read(byte b[], int off, int len)
            throws IOException {
        return readBytes(b, off, len);
    }

    /**
     * Read up to
     * <code>b.length( )</code> bytes into an array. This will block until at
     * least one byte has been read.
     *
     * @param b the byte array to receive the bytes.
     * @return the actual number of bytes read, or -1 if there is not more data
     * due to the end of the file being reached.
     * @exception IOException if an I/O error occurrs.
     */
    public int read(byte b[])
            throws IOException {
        return readBytes(b, 0, b.length);
    }

    /**
     * Reads
     * <code>b.length</code> bytes from this file into the byte array. This
     * method reads repeatedly from the file until all the bytes are read. This
     * method blocks until all the bytes are read, the end of the stream is
     * detected, or an exception is thrown.
     *
     * @param b the buffer into which the data is read.
     * @exception EOFException if this file reaches the end before reading all
     * the bytes.
     * @exception IOException if an I/O error occurs.
     */
    public final void readFully(byte b[])
            throws IOException {
        readFully(b, 0, b.length);
    }

    /**
     * Reads exactly
     * <code>len</code> bytes from this file into the byte array. This method
     * reads repeatedly from the file until all the bytes are read. This method
     * blocks until all the bytes are read, the end of the stream is detected,
     * or an exception is thrown.
     *
     * @param b the buffer into which the data is read.
     * @param off the start offset of the data.
     * @param len the number of bytes to read.
     * @exception EOFException if this file reaches the end before reading all
     * the bytes.
     * @exception IOException if an I/O error occurs.
     */
    public final void readFully(byte b[], int off, int len)
            throws IOException {
        int n = 0;
        while (n < len) {
            int count = this.read(b, off + n, len - n);
            if (count < 0) {
                throw new EOFException();
            }
            n += count;
        }
    }

    /**
     * Skips exactly
     * <code>n</code> bytes of input. This method blocks until all the bytes are
     * skipped, the end of the stream is detected, or an exception is thrown.
     *
     * @param n the number of bytes to be skipped.
     * @return the number of bytes skipped, which is always <code>n</code>.
     * @exception EOFException if this file reaches the end before skipping all
     * the bytes.
     * @exception IOException if an I/O error occurs.
     */
    public int skipBytes(int n)
            throws IOException {
        seek(getFilePointer() + n);
        return n;
    }

    /**
     * Unread the last byte read. This method should not be used more than once
     * between reading operations, or strange things might happen.
     */
    public final void unread() {
        filePosition--;
    }

    //
    // Write primitives.
    //
    /**
     * Write a byte to the file. If the file has not been opened for writing, an
     * IOException will be raised only when an attempt is made to write the
     * buffer to the file. <p> Caveat: the effects of seek( )ing beyond the end
     * of the file are undefined.
     *
     * @exception IOException if an I/O error occurrs.
     */
    public final void write(int b)
            throws IOException {

        // If the file position is within the block of data...
        if (filePosition < dataEnd) {
            buffer[(int) (filePosition++ - bufferStart)] = (byte) b;
            bufferModified = true;

            // ...or (assuming that seek will not allow the file pointer
            // to move beyond the end of the file) get the correct block of
            // data... 
        } else {

            // If there is room in the buffer, expand it...
            if (dataSize != buffer.length) {
                buffer[(int) (filePosition++ - bufferStart)] = (byte) b;
                bufferModified = true;
                dataSize++;
                dataEnd++;

                // ...or do another seek to get a new buffer, and start again...
            } else {
                seek(filePosition);
                write(b);
            }
        }
    }

    /**
     * Write
     * <code>len</code> bytes from an array to the file.
     *
     * @param b the array containing the data.
     * @param off the offset in the array to the data.
     * @param len the length of the data.
     * @exception IOException if an I/O error occurrs.
     */
    public final void writeBytes(byte b[], int off, int len)
            throws IOException {

        // If the amount of data is small (less than a full buffer)...
        if (len < buffer.length) {

            // If any of the data fits within the buffer...
            int spaceInBuffer = 0;
            int copyLength = 0;
            if (filePosition >= bufferStart) {
                spaceInBuffer = (int) ((bufferStart + buffer.length) - filePosition);
            }
            if (spaceInBuffer > 0) {

                // Copy as much as possible to the buffer.
                copyLength = (spaceInBuffer > len) ? len : spaceInBuffer;
                System.arraycopy(b, off, buffer,
                        (int) (filePosition - bufferStart), copyLength);
                bufferModified = true;
                long myDataEnd = filePosition + copyLength;
                dataEnd = myDataEnd > dataEnd ? myDataEnd : dataEnd;
                dataSize = (int) (dataEnd - bufferStart);
                filePosition += copyLength;
            }

            // If there is any data remaining, move to the new position and copy to
            // the new buffer.
            if (copyLength < len) {
                seek(filePosition);
                System.arraycopy(b, off + copyLength, buffer,
                        (int) (filePosition - bufferStart),
                        len - copyLength);
                bufferModified = true;
                long myDataEnd = filePosition + (len - copyLength);
                dataEnd = myDataEnd > dataEnd ? myDataEnd : dataEnd;
                dataSize = (int) (dataEnd - bufferStart);
                filePosition += (len - copyLength);
            }

            // ...or write a lot of data...
        } else {

            // Flush the current buffer, and write this data to the file.
            if (bufferModified) {
                flush();
                bufferStart = dataEnd = dataSize = 0;
            }
            file.write(b, off, len);
            filePosition += len;
        }
    }

    /**
     * Writes
     * <code>b.length</code> bytes from the specified byte array starting at
     * offset
     * <code>off</code> to this file.
     *
     * @param b the data.
     * @exception IOException if an I/O error occurs.
     */
    public void write(byte b[]) throws IOException {
        writeBytes(b, 0, b.length);
    }

    /**
     * Writes
     * <code>len</code> bytes from the specified byte array starting at offset
     * <code>off</code> to this file.
     *
     * @param b the data.
     * @param off the start offset in the data.
     * @param len the number of bytes to write.
     * @exception IOException if an I/O error occurs.
     */
    public void write(byte b[], int off, int len) throws IOException {
        writeBytes(b, off, len);
    }

    //
    // DataInput methods.
    //
    /**
     * Reads a
     * <code>boolean</code> from this file. This method reads a single byte from
     * the file. A value of
     * <code>0</code> represents
     * <code>false</code>. Any other value represents
     * <code>true</code>. This method blocks until the byte is read, the end of
     * the stream is detected, or an exception is thrown.
     *
     * @return the <code>boolean</code> value read.
     * @exception EOFException if this file has reached the end.
     * @exception IOException if an I/O error occurs.
     */
    public final boolean readBoolean() throws IOException {
        int ch = this.read();
        if (ch < 0) {
            throw new EOFException();
        }
        return (ch != 0);
    }

    /**
     * Reads a signed 8-bit value from this file. This method reads a byte from
     * the file. If the byte read is
     * <code>b</code>, where
     * <code>0&nbsp;&lt;=&nbsp;b&nbsp;&lt;=&nbsp;255</code>, then the result is:      <ul><code>
    *     (byte)(b)
     * </code></ul> <p> This method blocks until the byte is read, the end of
     * the stream is detected, or an exception is thrown.
     *
     * @return the next byte of this file as a signed 8-bit <code>byte</code>.
     * @exception EOFException if this file has reached the end.
     * @exception IOException if an I/O error occurs.
     */
    public final byte readByte() throws IOException {
        int ch = this.read();
        if (ch < 0) {
            throw new EOFException();
        }
        return (byte) (ch);
    }

    /**
     * Reads an unsigned 8-bit number from this file. This method reads a byte
     * from this file and returns that byte. <p> This method blocks until the
     * byte is read, the end of the stream is detected, or an exception is
     * thrown.
     *
     * @return the next byte of this file, interpreted as an unsigned 8-bit
     * number.
     * @exception EOFException if this file has reached the end.
     * @exception IOException if an I/O error occurs.
     */
    public final int readUnsignedByte() throws IOException {
        int ch = this.read();
        if (ch < 0) {
            throw new EOFException();
        }
        return ch;
    }

    /**
     * Reads a signed 16-bit number from this file. The method reads 2 bytes
     * from this file. If the two bytes read, in order, are
     * <code>b1</code> and
     * <code>b2</code>, where each of the two values is between
     * <code>0</code> and
     * <code>255</code>, inclusive, then the result is equal to:      <ul><code>
    *     (short)((b1 &lt;&lt; 8) | b2)
     * </code></ul> <p> This method blocks until the two bytes are read, the end
     * of the stream is detected, or an exception is thrown.
     *
     * @return the next two bytes of this file, interpreted as a signed 16-bit
     * number.
     * @exception EOFException if this file reaches the end before reading two
     * bytes.
     * @exception IOException if an I/O error occurs.
     */
    public final short readShort() throws IOException {
        int ch1 = this.read();
        int ch2 = this.read();
        if ((ch1 | ch2) < 0) {
            throw new EOFException();
        }
        return (short) ((ch1 << 8) + (ch2 << 0));
    }

    /**
     * Reads an unsigned 16-bit number from this file. This method reads two
     * bytes from the file. If the bytes read, in order, are
     * <code>b1</code> and
     * <code>b2</code>, where
     * <code>0&nbsp;&lt;=&nbsp;b1, b2&nbsp;&lt;=&nbsp;255</code>, then the
     * result is equal to:      <ul><code>
    *     (b1 &lt;&lt; 8) | b2
     * </code></ul> <p> This method blocks until the two bytes are read, the end
     * of the stream is detected, or an exception is thrown.
     *
     * @return the next two bytes of this file, interpreted as an unsigned
     * 16-bit integer.
     * @exception EOFException if this file reaches the end before reading two
     * bytes.
     * @exception IOException if an I/O error occurs.
     */
    public final int readUnsignedShort() throws IOException {
        int ch1 = this.read();
        int ch2 = this.read();
        if ((ch1 | ch2) < 0) {
            throw new EOFException();
        }
        return (ch1 << 8) + (ch2 << 0);
    }

    /**
     * Reads a Unicode character from this file. This method reads two bytes
     * from the file. If the bytes read, in order, are
     * <code>b1</code> and
     * <code>b2</code>, where
     * <code>0&nbsp;&lt;=&nbsp;b1,&nbsp;b2&nbsp;&lt;=&nbsp;255</code>, then the
     * result is equal to:      <ul><code>
    *     (char)((b1 &lt;&lt; 8) | b2)
     * </code></ul> <p> This method blocks until the two bytes are read, the end
     * of the stream is detected, or an exception is thrown.
     *
     * @return the next two bytes of this file as a Unicode character.
     * @exception EOFException if this file reaches the end before reading two
     * bytes.
     * @exception IOException if an I/O error occurs.
     */
    public final char readChar() throws IOException {
        int ch1 = this.read();
        int ch2 = this.read();
        if ((ch1 | ch2) < 0) {
            throw new EOFException();
        }
        return (char) ((ch1 << 8) + (ch2 << 0));
    }

    /**
     * Reads a signed 32-bit integer from this file. This method reads 4 bytes
     * from the file. If the bytes read, in order, are
     * <code>b1</code>,
     * <code>b2</code>,
     * <code>b3</code>, and
     * <code>b4</code>, where
     * <code>0&nbsp;&lt;=&nbsp;b1, b2, b3, b4&nbsp;&lt;=&nbsp;255</code>, then
     * the result is equal to:      <ul><code>
    *     (b1 &lt;&lt; 24) | (b2 &lt;&lt; 16) + (b3 &lt;&lt; 8) + b4
     * </code></ul> <p> This method blocks until the four bytes are read, the
     * end of the stream is detected, or an exception is thrown.
     *
     * @return the next four bytes of this file, interpreted as an
     * <code>int</code>.
     * @exception EOFException if this file reaches the end before reading four
     * bytes.
     * @exception IOException if an I/O error occurs.
     */
    public final int readInt() throws IOException {
        int ch1 = this.read();
        int ch2 = this.read();
        int ch3 = this.read();
        int ch4 = this.read();
        if ((ch1 | ch2 | ch3 | ch4) < 0) {
            throw new EOFException();
        }
        return ((ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
    }

    /**
     * Reads a signed 64-bit integer from this file. This method reads eight
     * bytes from the file. If the bytes read, in order, are
     * <code>b1</code>,
     * <code>b2</code>,
     * <code>b3</code>,
     * <code>b4</code>,
     * <code>b5</code>,
     * <code>b6</code>,
     * <code>b7</code>, and
     * <code>b8,</code> where:      <ul><code>
    *     0 &lt;= b1, b2, b3, b4, b5, b6, b7, b8 &lt;=255,
     * </code></ul> <p> then the result is equal to:      <p><blockquote><pre>
    *     ((long)b1 &lt;&lt; 56) + ((long)b2 &lt;&lt; 48)
     *     + ((long)b3 &lt;&lt; 40) + ((long)b4 &lt;&lt; 32)
     *     + ((long)b5 &lt;&lt; 24) + ((long)b6 &lt;&lt; 16)
     *     + ((long)b7 &lt;&lt; 8) + b8
     * </pre></blockquote> <p> This method blocks until the eight bytes are
     * read, the end of the stream is detected, or an exception is thrown.
     *
     * @return the next eight bytes of this file, interpreted as a
     * <code>long</code>.
     * @exception EOFException if this file reaches the end before reading eight
     * bytes.
     * @exception IOException if an I/O error occurs.
     */
    public final long readLong() throws IOException {
        return ((long) (readInt()) << 32) + (readInt() & 0xFFFFFFFFL);
    }

    /**
     * Reads a
     * <code>float</code> from this file. This method reads an
     * <code>int</code> value as if by the
     * <code>readInt</code> method and then converts that
     * <code>int</code> to a
     * <code>float</code> using the
     * <code>intBitsToFloat</code> method in class
     * <code>Float</code>. <p> This method blocks until the four bytes are read,
     * the end of the stream is detected, or an exception is thrown.
     *
     * @return the next four bytes of this file, interpreted as a
     * <code>float</code>.
     * @exception EOFException if this file reaches the end before reading four
     * bytes.
     * @exception IOException if an I/O error occurs.
     * @see java.io.RandomAccessFile#readInt()
     * @see java.lang.Float#intBitsToFloat(int)
     */
    public final float readFloat() throws IOException {
        return Float.intBitsToFloat(readInt());
    }

    /**
     * Reads a
     * <code>double</code> from this file. This method reads a
     * <code>long</code> value as if by the
     * <code>readLong</code> method and then converts that
     * <code>long</code> to a
     * <code>double</code> using the
     * <code>longBitsToDouble</code> method in class
     * <code>Double</code>. <p> This method blocks until the eight bytes are
     * read, the end of the stream is detected, or an exception is thrown.
     *
     * @return the next eight bytes of this file, interpreted as a
     * <code>double</code>.
     * @exception EOFException if this file reaches the end before reading eight
     * bytes.
     * @exception IOException if an I/O error occurs.
     * @see java.io.RandomAccessFile#readLong()
     * @see java.lang.Double#longBitsToDouble(long)
     */
    public final double readDouble() throws IOException {
        return Double.longBitsToDouble(readLong());
    }

    /**
     * Reads the next line of text from this file. This method successively
     * reads bytes from the file until it reaches the end of a line of text. <p>
     * A line of text is terminated by a carriage-return character 
    * (
     * <code>'&#92;r'</code>), a newline character (
     * <code>'&#92;n'</code>), a carriage-return character immediately followed
     * by a newline character, or the end of the input stream. The
     * line-terminating character(s), if any, are included as part of the string
     * returned. <p> This method blocks until a newline character is read, a
     * carriage return and the byte following it are read (to see if it is a
     * newline), the end of the stream is detected, or an exception is thrown.
     *
     * @return the next line of text from this file.
     * @exception IOException if an I/O error occurs.
     */
    public final String readLine() throws IOException {
        StringBuffer input = new StringBuffer();
        int c;

        while (((c = read()) != -1) && (c != '\n')) {
            input.append((char) c);
        }
        if ((c == -1) && (input.length() == 0)) {
            return null;
        }
        return input.toString();
    }

    /**
     * Reads in a string from this file. The string has been encoded using a
     * modified UTF-8 format. <p> The first two bytes are read as if by
     * <code>readUnsignedShort</code>. This value gives the number of following
     * bytes that are in the encoded string, not the length of the resulting
     * string. The following bytes are then interpreted as bytes encoding
     * characters in the UTF-8 format and are converted into characters. <p>
     * This method blocks until all the bytes are read, the end of the stream is
     * detected, or an exception is thrown.
     *
     * @return a Unicode string.
     * @exception EOFException if this file reaches the end before reading all
     * the bytes.
     * @exception IOException if an I/O error occurs.
     * @exception UTFDataFormatException if the bytes do not represent valid
     * UTF-8 encoding of a Unicode string.
     * @see java.io.RandomAccessFile#readUnsignedShort()
     */
    public final String readUTF() throws IOException {
        return DataInputStream.readUTF(this);
    }

    //
    // DataOutput methods.
    //
    /**
     * Writes a
     * <code>boolean</code> to the file as a 1-byte value. The value
     * <code>true</code> is written out as the value
     * <code>(byte)1</code>; the value
     * <code>false</code> is written out as the value
     * <code>(byte)0</code>.
     *
     * @param v a <code>boolean</code> value to be written.
     * @exception IOException if an I/O error occurs.
     */
    public final void writeBoolean(boolean v) throws IOException {
        write(v ? 1 : 0);
    }

    /**
     * Writes a
     * <code>byte</code> to the file as a 1-byte value.
     *
     * @param v a <code>byte</code> value to be written.
     * @exception IOException if an I/O error occurs.
     */
    public final void writeByte(int v) throws IOException {
        write(v);
    }

    /**
     * Writes a
     * <code>short</code> to the file as two bytes, high byte first.
     *
     * @param v a <code>short</code> to be written.
     * @exception IOException if an I/O error occurs.
     */
    public final void writeShort(int v) throws IOException {
        write((v >>> 8) & 0xFF);
        write((v >>> 0) & 0xFF);
    }

    /**
     * Writes a
     * <code>char</code> to the file as a 2-byte value, high byte first.
     *
     * @param v a <code>char</code> value to be written.
     * @exception IOException if an I/O error occurs.
     */
    public final void writeChar(int v) throws IOException {
        write((v >>> 8) & 0xFF);
        write((v >>> 0) & 0xFF);
    }

    /**
     * Writes an
     * <code>int</code> to the file as four bytes, high byte first.
     *
     * @param v an <code>int</code> to be written.
     * @exception IOException if an I/O error occurs.
     */
    public final void writeInt(int v) throws IOException {
        write((v >>> 24) & 0xFF);
        write((v >>> 16) & 0xFF);
        write((v >>> 8) & 0xFF);
        write((v >>> 0) & 0xFF);
    }

    /**
     * Writes a
     * <code>long</code> to the file as eight bytes, high byte first.
     *
     * @param v a <code>long</code> to be written.
     * @exception IOException if an I/O error occurs.
     */
    public final void writeLong(long v) throws IOException {
        write((int) (v >>> 56) & 0xFF);
        write((int) (v >>> 48) & 0xFF);
        write((int) (v >>> 40) & 0xFF);
        write((int) (v >>> 32) & 0xFF);
        write((int) (v >>> 24) & 0xFF);
        write((int) (v >>> 16) & 0xFF);
        write((int) (v >>> 8) & 0xFF);
        write((int) (v >>> 0) & 0xFF);
    }

    /**
     * Converts the float argument to an
     * <code>int</code> using the
     * <code>floatToIntBits</code> method in class
     * <code>Float</code>, and then writes that
     * <code>int</code> value to the file as a 4-byte quantity, high byte first.
     *
     * @param v a <code>float</code> value to be written.
     * @exception IOException if an I/O error occurs.
     * @see java.lang.Float#floatToIntBits(float)
     */
    public final void writeFloat(float v) throws IOException {
        writeInt(Float.floatToIntBits(v));
    }

    /**
     * Converts the double argument to a
     * <code>long</code> using the
     * <code>doubleToLongBits</code> method in class
     * <code>Double</code>, and then writes that
     * <code>long</code> value to the file as an 8-byte quantity, high byte
     * first.
     *
     * @param v a <code>double</code> value to be written.
     * @exception IOException if an I/O error occurs.
     * @see java.lang.Double#doubleToLongBits(double)
     */
    public final void writeDouble(double v) throws IOException {
        writeLong(Double.doubleToLongBits(v));
    }

    /**
     * Writes the string to the file as a sequence of bytes. Each character in
     * the string is written out, in sequence, by discarding its high eight
     * bits.
     *
     * @param s a string of bytes to be written.
     * @exception IOException if an I/O error occurs.
     */
    public final void writeBytes(String s) throws IOException {
        int len = s.length();
        for (int i = 0; i < len; i++) {
            write((byte) s.charAt(i));
        }
    }

    /**
     * Writes the character array to the file as a sequence of bytes. Each
     * character in the string is written out, in sequence, by discarding its
     * high eight bits.
     *
     * @param b a character array of bytes to be written.
     * @param off the index of the first character to write.
     * @param len the number of characters to write.
     * @exception IOException if an I/O error occurs.
     */
    public final void writeBytes(char b[], int off, int len) throws IOException {
        for (int i = off; i < len; i++) {
            write((byte) b[i]);
        }
    }

    /**
     * Writes a string to the file as a sequence of characters. Each character
     * is written to the data output stream as if by the
     * <code>writeChar</code> method.
     *
     * @param s a <code>String</code> value to be written.
     * @exception IOException if an I/O error occurs.
     * @see java.io.RandomAccessFile#writeChar(int)
     */
    public final void writeChars(String s) throws IOException {
        int len = s.length();
        for (int i = 0; i < len; i++) {
            int v = s.charAt(i);
            write((v >>> 8) & 0xFF);
            write((v >>> 0) & 0xFF);
        }
    }

    /**
     * Writes a string to the file using UTF-8 encoding in a machine-independent
     * manner. <p> First, two bytes are written to the file as if by the
     * <code>writeShort</code> method giving the number of bytes to follow. This
     * value is the number of bytes actually written out, not the length of the
     * string. Following the length, each character of the string is output, in
     * sequence, using the UTF-8 encoding for each character.
     *
     * @param str a string to be written.
     * @exception IOException if an I/O error occurs.
     */
    public final void writeUTF(String str) throws IOException {
        int strlen = str.length();
        int utflen = 0;

        for (int i = 0; i < strlen; i++) {
            int c = str.charAt(i);
            if ((c >= 0x0001) && (c <= 0x007F)) {
                utflen++;
            } else if (c > 0x07FF) {
                utflen += 3;
            } else {
                utflen += 2;
            }
        }
        if (utflen > 65535) {
            throw new UTFDataFormatException();
        }

        write((utflen >>> 8) & 0xFF);
        write((utflen >>> 0) & 0xFF);
        for (int i = 0; i < strlen; i++) {
            int c = str.charAt(i);
            if ((c >= 0x0001) && (c <= 0x007F)) {
                write(c);
            } else if (c > 0x07FF) {
                write(0xE0 | ((c >> 12) & 0x0F));
                write(0x80 | ((c >> 6) & 0x3F));
                write(0x80 | ((c >> 0) & 0x3F));
            } else {
                write(0xC0 | ((c >> 6) & 0x1F));
                write(0x80 | ((c >> 0) & 0x3F));
            }
        }
    }

    /**
     * Create a string representation of this object.
     *
     * @return a string representation of the state of the object.
     */
    public String toString() {
        return "fp=" + filePosition + ", bs=" + bufferStart
                + ", de=" + dataEnd + ", ds=" + dataSize
                + ", bl=" + buffer.length + ", m=" + mode
                + ", bm=" + bufferModified;
    }

    /**
     * Test the byte operations of the RandomAccessFile class. These are the
     * methods that read/write on a byte-by-byte basis. The following checks are
     * made: <ul> <li>Writing random bytes to a file. <li>Checking the size of
     * the file is correct. <li>Checking that EOF is correctly raised.
     * <li>Reading the file back in and verifying its contents. </ul> The test
     * file is 4.5 times the size of the buffer, in order to test paging between
     * buffers, and using files that end in the middle of a buffer. A constant
     * seed value is used for the random number generator, to ensure any bugs
     * are reproduceable.
     *
     * @param filename the name of the test file to generate.
     * @param bufferSize the size of the buffer to use.
     */
    public static void testBytes(String filename, int bufferSize) {

        System.out.println("\nTesting byte operations...");
        int newFileSize = (int) (bufferSize * 4.5);

        try {

            // Create a test file.
            RandomAccessFile outFile = new RandomAccessFile(filename,
                    RandomAccessFile.WRITE
                    | RandomAccessFile.CREATE, bufferSize);
            try {
                Random random = new Random(0);
                byte b = 0;
                for (int i = 0; i < newFileSize; i++) {
                    b = (byte) (random.nextInt() % 256);
                    outFile.writeByte(b);
                }
            } finally {
                outFile.close();
            }

            // Check that the file length is correct.
            if ((new File(filename)).length() == newFileSize) {
                System.out.println(". File size correct (" + newFileSize + ").");
            } else {
                System.out.println("X New file size incorrect (should be " + newFileSize
                        + ", but is " + (new File(filename)).length() + ").");
            }

            // Read the file, verify and modify its contents.
            RandomAccessFile inoutFile = new RandomAccessFile(filename,
                    RandomAccessFile.READ
                    | RandomAccessFile.WRITE, bufferSize);

            boolean verified = true;
            int byteNo = 0;
            try {

                // Read each byte in the file.
                Random random = new Random(0);
                byte b = 0;
                for (byteNo = 0; byteNo < newFileSize; byteNo++) {
                    b = (byte) (random.nextInt() % 256);
                    byte currentByte = inoutFile.readByte();

                    // Check the value is correct.
                    if (currentByte != b) {
                        verified = false;
                    }

                    // Modify selected values.
                    if (currentByte >= 128) {
                        inoutFile.seek(inoutFile.getFilePointer() - 1);
                        inoutFile.writeByte(0);
                    }
                }

                // Check the EOF is correctly trapped.
                boolean foundEOF = false;
                try {
                    inoutFile.readByte();
                } catch (EOFException e) {
                    foundEOF = true;
                }
                if (foundEOF) {
                    System.err.println(". EOF found correctly");
                } else {
                    System.err.println("X No EOF found.");
                }

                // Trace a premature EOF.
            } catch (EOFException e) {
                e.printStackTrace();
                System.err.println("    At byte " + byteNo);
            } finally {
                inoutFile.close();
            }

            // Check that the read was verified.
            if (verified) {
                System.out.println(". Read/Write verified");
            } else {
                System.out.println("X Read/Write verification failed");
            }

            // Read the file and verify contents.
            RandomAccessFile inFile = new RandomAccessFile(filename,
                    RandomAccessFile.READ, bufferSize);

            verified = true;
            byteNo = 0;
            try {

                // Read each byte in the file.
                Random random = new Random(0);
                byte b = 0;
                for (byteNo = 0; byteNo < newFileSize; byteNo++) {
                    b = (byte) (random.nextInt() % 256);
                    byte currentByte = inFile.readByte();

                    // Account for the modification.
                    if (currentByte >= 128) {
                        currentByte = 0;
                    }

                    // Check the byte's value.
                    if (currentByte != b) {
                        verified = false;
                    }
                }

                // Trap a premature EOF.
            } catch (EOFException e) {
                e.printStackTrace();
                System.err.println("    At byte " + byteNo);
            } finally {
                inFile.close();
            }

            // Check that the read was verified.
            if (verified) {
                System.out.println(". Update verified");
            } else {
                System.out.println("X Update verification failed");
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Test the block operations of the RandomAccessFile class. These are the
     * methods that read/write blocks of data. The following checks are made:
     * <ul> <li>Writing blocks of data that are smaller than the buffer size.
     * <li>Writing blocks of data that are larger than the buffer size.
     * <li>Checking the size of the file is correct. <li>Reading small blocks of
     * the file back in and verifying its contents. <li>Reading large blocks of
     * the file back in and verifying its contents. </ul>
     *
     * @param filename the name of the test file to generate.
     */
    public static void testBlocks(String filename) {

        System.err.println("\nTesting block operations...");

        // Generate the data.
        int bufferSize = 10;
        byte data[] = new byte[256];
        for (int i = 0; i < data.length; i++) {
            data[i] = (byte) (i % 256);
        }

        try {

            // Write the data in small and large blocks.
            RandomAccessFile outFile = new RandomAccessFile(
                    filename, RandomAccessFile.WRITE
                    | RandomAccessFile.CREATE, bufferSize);
            for (int i = 0; i < data.length;) {
                int blockSize = (i < data.length / 2) ? 3
                        : 13;
                blockSize = (i + blockSize >= data.length) ? (data.length - i)
                        : blockSize;
                outFile.write(data, i, blockSize);
                i += blockSize;
            }

            outFile.close();

            // Check that the file length is correct.
            if ((new File(filename)).length() != data.length) {
                System.out.println("X New file size incorrect (should be " + data.length
                        + ", but is " + (new File(filename)).length() + ").");
            } else {
                System.out.println(". File size correct (" + data.length + ").");
            }

            // Reopen the file for reading.
            RandomAccessFile inFile = new RandomAccessFile(
                    filename, RandomAccessFile.READ, bufferSize);

            // Read and check random small blocks of data.
            boolean verified = true;
            int firstFailure = 256;
            Random random = new Random(0);
            byte block[] = new byte[(int) (bufferSize * 0.5)];
            for (int i = 0; i < 100; i++) {
                int index = Math.abs(random.nextInt()) % (data.length - block.length);
                inFile.seek(index);
                inFile.read(block);

                // Verify the block of data.
                for (int j = 0; j < block.length; j++) {
                    if (block[j] != data[index + j]) {
                        verified = false;
                        if (index + j < firstFailure) {
                            firstFailure = index + j;
                        }
                    }
                }
            }
            if (verified) {
                System.err.println(". Reading small blocks verified.");
            } else {
                System.err.println("X Reading small blocks failed (byte " + firstFailure + ").");
            }

            // Read and check random large (bigger than the bufferSize) blocks
            // of data.
            verified = true;
            random = new Random(0);
            block = new byte[(int) (bufferSize * 1.5)];
            for (int i = 0; i < 100; i++) {
                int index = Math.abs(random.nextInt()) % (data.length - block.length);
                inFile.seek(index);
                inFile.read(block);

                // Verify the block of data.
                for (int j = 0; j < block.length; j++) {
                    if (block[j] != data[j + index]) {
                        verified = false;
                    }
                }
            }
            if (verified) {
                System.err.println(". Reading large blocks verified.");
            } else {
                System.err.println("X Reading large blocks failed.");
            }

            // Close the input file.
            inFile.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    /**
     * Benchmark the performance of the new RandomAccessFile class. Its speed is
     * compared to that of a java.io.RandomAccessFile, based on reading and
     * writing a test file, byte by byte.
     *
     * @param filename the name of the test file.
     * @param bufferSize the buffer size to use.
     */
    public static void benchmark(String filename, int bufferSize) {
        System.out.println("\nBenchmarking...");

        // Start the clock, and open a file for reading and a file for writing.
        long time = (new Date()).getTime();
        try {
            RandomAccessFile inFile = new RandomAccessFile(filename,
                    RandomAccessFile.READ, bufferSize);
            RandomAccessFile outFile = new RandomAccessFile("temp.data",
                    RandomAccessFile.WRITE
                    | RandomAccessFile.CREATE, bufferSize);

            // Copy one file to the other.
            try {

                while (true) {
                    outFile.writeByte(inFile.readByte());
                }

            } catch (EOFException e) {
            } catch (IOException e) {
                e.printStackTrace();
            } finally {
                inFile.close();
                outFile.close();
            }
            System.out.println(". RandomAccessFile elapsed time="
                    + ((new Date()).getTime() - time));

            // Restart the clock, and open RandomAccessFiles for reading and writing.
            time = (new Date()).getTime();
            java.io.RandomAccessFile inFile2 = new java.io.RandomAccessFile(filename, "r");
            java.io.RandomAccessFile outFile2 = new java.io.RandomAccessFile("temp.data", "rw");

            // Copy one file to the other.
            try {

                while (true) {
                    outFile2.writeByte(inFile2.readByte());
                }

            } catch (EOFException e) {
            } catch (IOException e) {
                e.printStackTrace();
            } finally {
                inFile2.close();
                outFile2.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(". java.io.RandomAccessFile elapsed time=" + ((new Date()).getTime() - time));
    }

//    /**
//     * Test the RandomAccessFile class. This involves testing the byte methods,
//     * the block methods, and benchmarking the performance. By appending 'test'
//     * or 'benchmark' to the command-line, it can be limited to the tests or
//     * benchmarking alone. The test filename is only used for the benchmarking,
//     * the other tests create a file called "temp.data" in the current
//     * directory. Note that the size of the buffer determines the size of the
//     * test file (which is 4.5 times the size of the buffer).
//     *
//     * @param argv Usage: <testFilename> [bufferSize] [test | benchmark]
//     * @see testBytes
//     * @see testBlocks
//     * @see benchmark
//     */
//    public static void main(String argv[]) {
//
//        int defaultPageSize = 4096;
//
//        // Parse the command-line arguments.
//        String filename = null;
//        int bufferSize = 0;
//        boolean test = true;
//        boolean benchmark = true;
//        if (argv.length < 1) {
//            System.err.println("Usage: RandomAccessFile <filename> [buffer.length] [benchmark | test]");
//            System.exit(-1);
//        } else if (argv.length < 2) {
//            filename = argv[0];
//            bufferSize = defaultPageSize;
//        } else if (argv.length < 3) {
//            filename = argv[0];
//            bufferSize = Integer.parseInt(argv[1]);
//        } else {
//            filename = argv[0];
//            bufferSize = Integer.parseInt(argv[1]);
//            if (argv[2].equals("benchmark")) {
//                test = false;
//            } else if (argv[2].equals("test")) {
//                benchmark = false;
//            }
//        }
//
//        System.out.println("\nRandomAccessFile\n"
//                + "========================");
//        System.out.println("filename=" + filename
//                + ", bufferSize=" + bufferSize);
//        System.out.println("totalMemory="
//                + (Runtime.getRuntime().totalMemory() / 1000) + "k"
//                + " freeMemory=" + (Runtime.getRuntime().freeMemory() / 1000) + "k");
//
//        if (test) {
//            RandomAccessFile.testBytes("temp.data", bufferSize);
//            RandomAccessFile.testBlocks("temp.data");
//        }
//        if (benchmark) {
//            RandomAccessFile.benchmark(filename, bufferSize);
//        }
//
//        System.out.println("\nEND");
//    }
}
