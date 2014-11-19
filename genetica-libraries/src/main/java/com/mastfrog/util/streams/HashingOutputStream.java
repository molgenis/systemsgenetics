/* 
 * The MIT License
 *
 * Copyright 2013 Tim Boudreau.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package com.mastfrog.util.streams;

import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

/**
 * Output stream which wraps another output stream and computes a cryptographic
 * hash using some algorithm as the data arrives.
 *
 * @author Tim Boudreau
 */
public final class HashingOutputStream extends FilterOutputStream {
    private final MessageDigest digest;
    private volatile boolean closed;
    public HashingOutputStream(String algorithm, OutputStream out) throws NoSuchAlgorithmException {
        super (out);
        digest = HashingInputStream.createDigest(algorithm);
    }

    public static HashingOutputStream sha1(OutputStream out) throws NoSuchAlgorithmException {
        return new HashingOutputStream("SHA-1", out);
    }

    @Override
    public void close() throws IOException {
        closed = true;
        super.close();
    }

    public byte[] getDigest() {
        if (!closed) {
            throw new IllegalStateException ("Stream not closed");
        }
        return digest.digest();
    }

    @Override
    public void write(int b) throws IOException {
        if (closed) {
            throw new IOException ("Stream closed");
        }
        super.write(b);
        digest.update((byte) b);
    }

    public String getHashAsString() throws IOException {
        if (!closed) {
            close();
        }
        byte[] bytes = getDigest();
        return hashString(bytes);
    }
    
    public static String hashString(byte[] bytes) {
        // assumes hash is a multiple of 4
        IntBuffer ib = ByteBuffer.wrap(bytes).asIntBuffer();
        StringBuilder sb = new StringBuilder();
        while (ib.position() < ib.capacity()) {
            long val = ib.get();
            if (val < 0) {
                val = -val + Integer.MAX_VALUE;
            }
            sb.append(Long.toString(val, 36));
        }
        return sb.toString();
    }
    
}
