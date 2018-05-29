//******************************************************************************
//
// File:    NumericRuntimeException.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.NumericRuntimeException
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************

package main.java.decon_eQTL;

/**
 * Class NumericRuntimeException is the base class for an unchecked runtime
 * exception thrown by a method in a numerical class.
 *
 * @author  Alan Kaminsky
 * @version 06-Jul-2007
 */
public class NumericRuntimeException
	extends RuntimeException
	{

// Exported constructors.

	/**
	 * 
	 */
	private static final long serialVersionUID = 4282166190962135014L;

	/**
	 * Construct a new numeric runtime exception with no detail message and no
	 * cause.
	 */
	public NumericRuntimeException()
		{
		super();
		}

	/**
	 * Construct a new numeric runtime exception with the given detail message
	 * and no cause.
	 *
	 * @param  message  Detail message.
	 */
	public NumericRuntimeException
		(String message)
		{
		super (message);
		}

	/**
	 * Construct a new numeric runtime exception with the given cause and the
	 * default detail message.
	 *
	 * @param  cause  Cause.
	 */
	public NumericRuntimeException
		(Throwable cause)
		{
		super (cause);
		}

	/**
	 * Construct a new numeric runtime exception with the given detail message
	 * and the given cause.
	 * 
	 * @param message 	Detail message.
	 * @param  cause  Cause.
	 */
	public NumericRuntimeException
		(String message,
		 Throwable cause)
		{
		super (message, cause);
		}

	}
