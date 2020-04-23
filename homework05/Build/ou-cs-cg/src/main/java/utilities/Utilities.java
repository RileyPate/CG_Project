//******************************************************************************
// Copyright (C) 2016 University of Oklahoma Board of Trustees.
//******************************************************************************
// Last modified: Tue Feb 26 15:40:03 2019 by Chris Weaver
//******************************************************************************
// Major Modification History:
//
// 20160225 [weaver]:	Original file.
// 20190226 [weaver]:	Moved to utilities package, added coordinate mappings.
//
//******************************************************************************
// Notes:
//
//******************************************************************************

package utilities;

//import java.lang.*;
import java.awt.event.*;
import com.jogamp.opengl.*;
import com.jogamp.opengl.glu.GLU;
import com.jogamp.opengl.fixedfunc.GLMatrixFunc;

//******************************************************************************

/**
 * The <CODE>EventUtilities</CODE> class.<P>
 *
 * @author  Chris Weaver
 * @version %I%, %G%
 */
public final class Utilities
{
	//**********************************************************************
	// Private Methods (Coordinate Mapping)
	//**********************************************************************

	// Give this method integer x and y coordinates, e.g. of a mouse event.
	public static double[]	mapViewToScene(GL2 gl, double x, double y, double z)
	{
		GLU		glu = GLU.createGLU();
		double[]	mv = new double[16];		// modelview
		double[]	pr = new double[16];		// projection
		int[]		vp = new int[4];			// viewport

		gl.glGetDoublev(GLMatrixFunc.GL_MODELVIEW_MATRIX, mv, 0);
		gl.glGetDoublev(GLMatrixFunc.GL_PROJECTION_MATRIX, pr, 0);
		gl.glGetIntegerv(GL.GL_VIEWPORT, vp, 0);

		double[]	p = new double[3];			// Scene coordinates

		glu.gluUnProject(x, y, z, mv, 0, pr, 0, vp, 0, p, 0);

		return p;
	}

	// Receive from this method integer x and y pixel coordinates of scene item.
	public static double[]	mapSceneToView(GL2 gl, double x, double y, double z)
	{
		GLU		glu = GLU.createGLU();
		double[]	mv = new double[16];		// modelview
		double[]	pr = new double[16];		// projection
		int[]		vp = new int[4];			// viewport

		gl.glGetDoublev(GLMatrixFunc.GL_MODELVIEW_MATRIX, mv, 0);
		gl.glGetDoublev(GLMatrixFunc.GL_PROJECTION_MATRIX, pr, 0);
		gl.glGetIntegerv(GL.GL_VIEWPORT, vp, 0);

		double[]	p = new double[3];			// Screen coordinates

		glu.gluProject(x, y, z, mv, 0, pr, 0, vp, 0, p, 0);

		return p;
	}

	//**********************************************************************
	// Public Class Methods (Event Handling)
	//**********************************************************************

	public static boolean	isAltDown(InputEvent e)
	{
		return ((e.getModifiersEx() & InputEvent.ALT_DOWN_MASK) != 0);
	}

	public static boolean	isAltGraphDown(InputEvent e)
	{
		return ((e.getModifiersEx() & InputEvent.ALT_GRAPH_DOWN_MASK) != 0);
	}

	public static boolean	isControlDown(InputEvent e)
	{
		return ((e.getModifiersEx() & InputEvent.CTRL_DOWN_MASK) != 0);
	}

	public static boolean	isMetaDown(InputEvent e)
	{
		return ((e.getModifiersEx() & InputEvent.META_DOWN_MASK) != 0);
	}

	public static boolean	isShiftDown(InputEvent e)
	{
		return ((e.getModifiersEx() & InputEvent.SHIFT_DOWN_MASK) != 0);
	}

	public static boolean	isButton1Down(InputEvent e)
	{
		return (((e.getModifiersEx() & InputEvent.BUTTON1_DOWN_MASK) != 0) ||
				((e instanceof MouseEvent) &&
				 (((MouseEvent)e).getButton() == 1)));
	}

	public static boolean	isButton2Down(InputEvent e)
	{
		return (((e.getModifiersEx() & InputEvent.BUTTON2_DOWN_MASK) != 0) ||
				((e instanceof MouseEvent) &&
				 (((MouseEvent)e).getButton() == 2)));
	}

	public static boolean	isButton3Down(InputEvent e)
	{
		return (((e.getModifiersEx() & InputEvent.BUTTON3_DOWN_MASK) != 0) ||
				((e instanceof MouseEvent) &&
				 (((MouseEvent)e).getButton() == 3)));
	}
}

//******************************************************************************
