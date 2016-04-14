package deconvolution;


/*
 * From: https://github.com/tantaman/commons/blob/master/src/main/java/com/tantaman/commons/concurrent/NamedThreadFactory.java
 * 
 * Copyright 2010 Matt Crinklaw-Vogt
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


import java.util.concurrent.ThreadFactory;
import java.util.concurrent.atomic.AtomicLong;

public class NamedThreadFactory implements ThreadFactory {
	private static final AtomicLong THREAD_POOL_NUM = new AtomicLong(0);
	private final AtomicLong mThreadNum = new AtomicLong(0);
	private final String mPrefix;
	private final boolean mIsDaemon;
	private final long mPoolNum;
	
	public NamedThreadFactory(String pPrefix) {
		this(pPrefix, true);
	}
	
	public NamedThreadFactory(String pPrefix, boolean pIsDaemon) {
		mIsDaemon = pIsDaemon;
		mPrefix = pPrefix;
		mPoolNum = THREAD_POOL_NUM.incrementAndGet();
	}
	
	@Override
	public Thread newThread(Runnable r) {
		Thread t = new Thread(r, mPrefix + "-" + mPoolNum + "-Thread-" + mThreadNum.incrementAndGet());
		if (t.isDaemon() != mIsDaemon)
			t.setDaemon(mIsDaemon);
		
		return t;
	}
}