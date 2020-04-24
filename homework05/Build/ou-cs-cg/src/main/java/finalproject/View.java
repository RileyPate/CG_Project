//******************************************************************************
// Copyright (C) 2016-2019 University of Oklahoma Board of Trustees.
//******************************************************************************
// Last modified: Tue Mar 19 02:13:07 2019 by Chris Weaver
//******************************************************************************
// Major Modification History:
//
// 20160209 [weaver]:	Original file.
// 20190203 [weaver]:	Updated to JOGL 2.3.2 and cleaned up.
// 20190227 [weaver]:	Updated to use model and asynchronous event handling.
// 20190318 [weaver]:	Modified for homework04.
//
//******************************************************************************
// Notes:
//
//******************************************************************************

package finalproject;

//import java.lang.*;
import java.awt.*;
import java.awt.geom.Point2D;
import java.text.DecimalFormat;
import java.util.*;
import com.jogamp.opengl.*;
import com.jogamp.opengl.awt.GLJPanel;
import com.jogamp.opengl.glu.*;
import com.jogamp.opengl.util.FPSAnimator;
import com.jogamp.opengl.util.awt.TextRenderer;

import utilities.Utilities;

//******************************************************************************

/**
 * The <CODE>View</CODE> class.<P>
 *
 * @author  Chris Weaver
 * @version %I%, %G%
 */
public final class View
	implements GLEventListener
{
	//**********************************************************************
	// Private Class Members
	//**********************************************************************

	private static final int			DEFAULT_FRAMES_PER_SECOND = 60;
	private static final DecimalFormat	FORMAT = new DecimalFormat("0.000");

	//**********************************************************************
	// Public Class Members
	//**********************************************************************

	public static final int			MIN_SIDES = 3;
	public static final int			MAX_SIDES = 12;

	//**********************************************************************
	// Private Members
	//**********************************************************************

	// State (internal) variables
	private final GLJPanel				canvas;
	private int							w;			// Canvas width
	private int							h;			// Canvas height

	private TextRenderer				renderer;

	private final FPSAnimator			animator;
	private int							counter;	// Frame counter

	private final Model					model;

	private final KeyHandler			keyHandler;
	private final MouseHandler			mouseHandler;

	private final Deque<Point2D.Double>				special;
	private final ArrayList<Deque<Point2D.Double>>	regions;
	private final ArrayList<Point2D.Double>			holeLocations;

	private final Deque<Point2D.Double>				tracing;
	private final Deque<Point2D.Double>				bounces;

	//private double						dx = 1.0 / DEFAULT_FRAMES_PER_SECOND;
	//private double						dy = 1.0 / DEFAULT_FRAMES_PER_SECOND;
	private double						dx = 0.0;
	private double						dy = 0.0;
	private double						slowDownRate = .9;

	//**********************************************************************
	// Constructors and Finalizer
	//**********************************************************************

	public View(GLJPanel canvas)
	{
		this.canvas = canvas;

		// Initialize rendering
		counter = 0;
		canvas.addGLEventListener(this);

		// Initialize model (scene data and parameter manager)
		model = new Model(this);

		// Initialize container polygons
		special = createSpecialCourse();					// For N = 2
		regions = new ArrayList<Deque<Point2D.Double>>();	// For MIN to MAX

		for (int i=MIN_SIDES; i<=MAX_SIDES; i++)
			regions.add(createPolygon(i));

		
		//Defines all of the locations for the holes.
		holeLocations = new ArrayList<Point2D.Double>();
		holeLocations.add(new Point2D.Double(-.5, .5));
		holeLocations.add(new Point2D.Double(0, .4));
		holeLocations.add(new Point2D.Double(.2, -.4));
		holeLocations.add(new Point2D.Double(-.3, -.5));
		holeLocations.add(new Point2D.Double(0, -.5));
		holeLocations.add(new Point2D.Double(0, -.2));
		holeLocations.add(new Point2D.Double(.7, .5));
		holeLocations.add(new Point2D.Double(.5, .5));
		holeLocations.add(new Point2D.Double(.8, -.2));
		holeLocations.add(new Point2D.Double(.3, -.8));
		holeLocations.add(new Point2D.Double(.2, .8));
		

		tracing = new ArrayDeque<Point2D.Double>();
		bounces = new ArrayDeque<Point2D.Double>();

		// Initialize controller (interaction handlers)
		keyHandler = new KeyHandler(this, model);
		mouseHandler = new MouseHandler(this, model);

		// Initialize animation
		animator = new FPSAnimator(canvas, DEFAULT_FRAMES_PER_SECOND);
		animator.start();
	}

	//**********************************************************************
	// Getters and Setters
	//**********************************************************************

	public GLJPanel	getCanvas()
	{
		return canvas;
	}

	public int	getWidth()
	{
		return w;
	}

	public int	getHeight()
	{
		return h;
	}

	//**********************************************************************
	// Public Methods
	//**********************************************************************

	public void	clearAllTrace()
	{
		tracing.clear();
		bounces.clear();
	}

	//**********************************************************************
	// Override Methods (GLEventListener)
	//**********************************************************************

	public void	init(GLAutoDrawable drawable)
	{
		w = drawable.getSurfaceWidth();
		h = drawable.getSurfaceHeight();

		renderer = new TextRenderer(new Font("Monospaced", Font.PLAIN, 12),
									true, true);

		initPipeline(drawable);
	}

	public void	dispose(GLAutoDrawable drawable)
	{
		renderer = null;
	}

	public void	display(GLAutoDrawable drawable)
	{
		updatePipeline(drawable);

		update(drawable);
		render(drawable);
	}

	public void	reshape(GLAutoDrawable drawable, int x, int y, int w, int h)
	{
		this.w = w;
		this.h = h;
	}

	//**********************************************************************
	// Private Methods (Rendering)
	//**********************************************************************

	private void	update(GLAutoDrawable drawable)
	{
		counter++;									// Advance animation counter
		
		
		dx = model.getRealPowerX();
		dy = model.getRealPowerY();
		
		Deque<Point2D.Double>	polygon = getCurrentCourse();
		Point2D.Double			q = model.getObject();

		updatePointWithReflection(polygon, q);
		model.setObjectInSceneCoordinatesAlt(new Point2D.Double(q.x, q.y));

		while (tracing.size() > DEFAULT_FRAMES_PER_SECOND)
			tracing.removeFirst();

		while (bounces.size() > DEFAULT_FRAMES_PER_SECOND)
			bounces.removeFirst();
		
		dx*=slowDownRate;
		if(dx<.001 && dx>-.001) {
			model.setRealPowerX(0.0);
		}
		else{
			model.setRealPowerX(dx);
		}
		
		dy*=slowDownRate;
		if(dy<.001 && dy>-.001) {
			model.setRealPowerY(0.0);
		}
		else{
			model.setRealPowerY(dy);
		}
		
		
	}

	private void	render(GLAutoDrawable drawable)
	{
		GL2	gl = drawable.getGL().getGL2();

		gl.glClear(GL.GL_COLOR_BUFFER_BIT);		// Clear the buffer

		// Draw the scene
		drawMain(gl);								// Draw main content

		gl.glFlush();								// Finish and display
	}

	//**********************************************************************
	// Private Methods (Pipeline)
	//**********************************************************************

	private void	initPipeline(GLAutoDrawable drawable)
	{
		GL2	gl = drawable.getGL().getGL2();

		gl.glClearColor(0.0f, 0.0f, 0.0f, 0.0f);	// Black background
	}

	private void	updatePipeline(GLAutoDrawable drawable)
	{
		GL2			gl = drawable.getGL().getGL2();
		GLU			glu = GLU.createGLU();

		gl.glMatrixMode(GL2.GL_PROJECTION);		// Prepare for matrix xform
		gl.glLoadIdentity();						// Set to identity matrix
		glu.gluOrtho2D(-1.2, 1.2, -1.2, 1.2);		// 2D translate and scale
	}

	//**********************************************************************
	// Private Methods (Scene)
	//**********************************************************************

	private void	drawMain(GL2 gl)
	{
		drawGreen(gl);						// Container polygon
		drawHole(gl);
		drawBall(gl);						// The moving object
		drawCursor(gl);						// Cursor around the mouse point
		drawAim(gl);						//Draw the aiming cursor
	}

	
	
	
	
	
	
	// Fills and edges the golf green that is surrounding the moving object.
	private void	drawGreen(GL2 gl)
	{
		Deque<Point2D.Double>	polygon = getCurrentCourse();

		gl.glColor3f(0.0f, 0.75f, 0.0f);			// Very dark gray
		fillPolygon(gl, polygon);

		gl.glColor3f(1.0f, 1.0f, 1.0f);			// White
		edgePolygon(gl, polygon);
	}
	
	
	
	
	
	
	
	// Draws the hole for the ball to go into.
	private void	drawHole(GL2 gl)
	{

		
		gl.glColor3f(0.0f, 0.0f, 0.0f);			//Black
		gl.glBegin(GL2.GL_POLYGON);

		for (int i = 0; i < 32; i++)
		{
			double	theta = (2.0 * Math.PI) * (i / 32.0);

			gl.glVertex2d(getCurrentCourseHole().x + 0.02 * Math.cos(theta),
						  getCurrentCourseHole().y + 0.02 * Math.sin(theta));
			
		}
		gl.glEnd();

	}
	
	
		
		
		
		
		

	// If the cursor point is not null, draw something helpful around it.
	private void	drawCursor(GL2 gl)
	{
		Point2D.Double	cursor = model.getCursor();

		if (cursor == null)
			return;

		// Draw a circle of radius 0.025

		gl.glColor3f(0.5f, 0.5f, 0.5f);			// Medium gray

		gl.glBegin(GL.GL_LINE_LOOP);

		for (int i=0; i<32; i++)
		{
			double	theta = (2.0 * Math.PI) * (i / 32.0);

			gl.glVertex2d(cursor.x + 0.025 * Math.cos(theta),
						  cursor.y + 0.025 * Math.sin(theta));
		}
		gl.glEnd();
	}

	
	
	
	// Draws the golf ball.
	private void	drawBall(GL2 gl)
	{
		Point2D.Double	object = model.getObject();
		
		
		gl.glColor3f(1.0f, 1.0f, 1.0f);			//White
		gl.glBegin(GL2.GL_POLYGON);

		
		//Fills the golf ball
		for (int i = 0; i < 32; i++)
		{
			double	theta = (2.0 * Math.PI) * (i / 32.0);

			gl.glVertex2d(object.x + 0.02 * Math.cos(theta),
						  object.y + 0.02 * Math.sin(theta));
			
		}
		gl.glEnd();
		
		

		
		
		gl.glColor3f(0.0f, 0.0f, 0.0f);			// Black
		gl.glBegin(GL.GL_LINE_LOOP);

		//outlines the golf ball
		for (int i = 0; i < 32; i++)
		{
			double	theta = (2.0 * Math.PI) * (i / 32.0);

			gl.glVertex2d(object.x + 0.02 * Math.cos(theta),
						  object.y + 0.02 * Math.sin(theta));
		}
		gl.glEnd();	
	}
	
	//Draw the aim cursor
	private void drawAim(GL2 gl) {
		if((dx==0.0) && (dy==0.0)) {
			Point2D.Double			q = model.getObject();
			double radius = model.getPower();
			double angle = model.getAngle();
			double x = radius*Math.cos(angle) + q.x;
			double y = radius*Math.sin(angle) + q.y;
			gl.glColor3f(0.0f, 1.0f, 1.0f);
				
			gl.glBegin(GL2.GL_POLYGON);

				
			//Fills the aim dot
			gl.glVertex2d(x+.01, y+.01);
			gl.glVertex2d(x+.01, y-.01);
			gl.glVertex2d(x-.01, y-.01);
			gl.glVertex2d(x-.01, y+.01);
			gl.glVertex2d(x+.01, y+.01);
			gl.glEnd();
			
		}
		
	}




	//**********************************************************************
	// Private Methods (Polygons)
	//**********************************************************************

	// Custom polygon for the sides=2 case. Irregular but convex.
	private Deque<Point2D.Double>	createSpecialCourse()
	{
		Deque<Point2D.Double>	polygon = new ArrayDeque<Point2D.Double>(10);

		polygon.add(new Point2D.Double( 1.00, -0.86));
		polygon.add(new Point2D.Double( 1.00, -0.24));
		polygon.add(new Point2D.Double( 0.48,  0.90));
		polygon.add(new Point2D.Double( 0.05,  1.00));
		polygon.add(new Point2D.Double(-0.34,  0.87));

		polygon.add(new Point2D.Double(-0.86,  0.40));
		polygon.add(new Point2D.Double(-1.00,  0.04));
		polygon.add(new Point2D.Double(-0.93, -0.42));
		polygon.add(new Point2D.Double(-0.53, -0.84));
		polygon.add(new Point2D.Double( 0.71, -1.00));

		return polygon;
	}

	
	
	// Creates a regular N-gon with points stored in counterclockwise order.
	// The polygon is centered at the origin with first vertex at (1.0, 0.0).
	private Deque<Point2D.Double>	createPolygon(int sides)
	{
		Deque<Point2D.Double>	polygon = new ArrayDeque<Point2D.Double>(sides);
		double					inv = 2.0 * Math.PI / sides;

		for (int i=0; i<sides; i++)
		{
			double	theta = i * inv;

			polygon.add(new Point2D.Double(Math.cos(theta), Math.sin(theta)));
		}

		return polygon;
	}

	// Draws the sides of the specified polygon.
	private void	edgePolygon(GL2 gl, Deque<Point2D.Double> polygon)
	{
		gl.glLineWidth(10.0f);
		gl.glBegin(GL.GL_LINE_LOOP);

		
		
		for (Point2D.Double p : polygon)
			gl.glVertex2d(p.x, p.y);

		gl.glEnd();
		
		gl.glLineWidth(1.0f);
	}

	// Draws the interior of the specified polygon.
	private void	fillPolygon(GL2 gl, Deque<Point2D.Double> polygon)
	{
		gl.glBegin(GL2.GL_POLYGON);

		for (Point2D.Double p : polygon)
			gl.glVertex2d(p.x, p.y);

		gl.glEnd();
	}

	// Get the polygon that is currently containing the moving object.
	private Deque<Point2D.Double>	getCurrentCourse()
	{
		int	sides = model.getNumber();

		if (sides == 2)
		{
			return special;
		}
		else if ((MIN_SIDES <= sides) && (sides <= MAX_SIDES))
		{
			return regions.get(sides - MIN_SIDES);
		}
		else
			return null;
	}
	
	
	
	
	
	// Get the polygon that is currently containing the moving object.
		private Point2D.Double	getCurrentCourseHole()
		{
			int	sides = model.getNumber();

			return holeLocations.get(sides - MIN_SIDES + 1);
		}

		
		
		
		
		
		
	// Special method for privileged use by the Model class ONLY.
	public boolean	currentPolygonContains(Point2D.Double q)
	{
		return contains(getCurrentCourse(), q);
	}

	//**********************************************************************
	// Private Methods (Reflection)
	//**********************************************************************

	// Updates the x and y coordinates of point q. Adds a vector to the provided
	// point, reflecting as needed off the sides of the provided polygon to
	// determine the new coordinates. The new coordinates are "returned" in q.
	public void	updatePointWithReflection(Deque<Point2D.Double> polygon,
											  Point2D.Double q)
	{
		// Scale the reference vector by the current velocity factor.
		double		factor = model.getFactor();
		double		ddx = factor * dx;
		double		ddy = factor * dy;

		// Consume the scaled vector through reflections until nothing is left.
		while (true)
		{
			int				sides = polygon.size();

			// Crude workaround for very rare (in practice) "corner" cases.
			// Unfortunately, even this doesn't resolve all edge cases!
			for (int i=0; i<sides; i++)
			{
				// Test each corner/vertex.
				Point2D.Double	p1 = polygon.peekFirst();

				// If the point is *exactly* on a corner/vertex...
				if ((p1.x == q.x) && (p1.y == q.y))
				{
					// ...reverse instead of reflect the reference vector...
					dx = -dx;
					dy = -dy;
					model.setRealPowerX(dx);
					model.setRealPowerY(dy);

					// ...and the scaled vector too...
					ddx = -ddx;
					ddy = -ddy;

					// ...and remember to freak out. Just a little. :)
					System.out.println("***WARNING: EXACT CORNER!!!***");
				}

				polygon.offerLast(polygon.pollFirst());
			}

			// Calculate which side the point will reach first. These variables
			// store that side's vertices and the parametric time to hit it.
			Point2D.Double		pp1 = null;
			Point2D.Double		pp2 = null;
			double				tmin = Double.MAX_VALUE;

			// Go around the polygon counterclockwise, taking vertices pairwise.
			Point2D.Double		p1 = polygon.peekLast();	// First in pair

			for (int i=0; i<sides; i++)
			{
				Point2D.Double	p2 = polygon.peekFirst();	// Second in pair

				// Calculate the CCW (inward-pointing) perp vector for the pair.
				double			vdx = p2.x - p1.x;		// Calc side vector
				double			vdy = p2.y - p1.y;		// from p1 to p2
				double			ndx = -vdy;			// Calc perp vector:
				double			ndy = vdx;				// negate y and swap

				// See page 175 and the slide on "Intersection of a Line through
				// a Line". (Note: R and v on the slide are named q and w here.)
				double			wdx = p1.x - q.x;		// Calc test vector
				double			wdy = p1.y - q.y;		// from q to p1

				// Calculate the top part of the t_hit equation.
				double			dnw = dot(ndx, ndy, 0.0, wdx, wdy, 0.0);

				// Check if q is strictly on the inside of the polygon. The dot
				// product will be 0 if q is on a side, or slightly positive if
				// it is beyond it (which can happen due to roundoff error). See
				// Figure 4.37 and the dot products below it on page 176.
				if (dnw < 0.0)
				{
					// Calculate the bottom part of the t_hit equation.
					double	dnv = dot(ndx, ndy, 0.0, ddx, ddy, 0.0);

					// If the dot project is zero, the direction of motion is
					// parallel to the side. Disqualify it as a hit candidate
					// (even if the point is exactly ON the side).
					double	thit = ((dnv != 0.0) ? (dnw / dnv) : 0.0);

					// Remember the side with the smallest positive t_hit.
					// It's the side that the point will reach first.
					if ((0.0 < thit) && (thit < tmin))
					{
						pp1 = p1;
						pp2 = p2;
						tmin = thit;
					}
				}

				// Promote the second vertex to be the first, and advance CCW.
				p1 = p2;
				polygon.offerLast(polygon.pollFirst());
			}

			if (tmin > 1.0)	// If the smallest positive t_hit is over 1.0,
			{					// the point won't reach the closest side in
				q.x += ddx;	// this update. Simply add the velocity vector
				q.y += ddy;	// to translate the point to its new position.
				tracing.offerLast(new Point2D.Double(q.x, q.y));

				break;			// Now escape from this infinite-seeming loop!
			}
			else
			{
				// If it is under 1.0, there will be at least one reflection.
				// Translate the point to the reflection point along the side.
				q.x += ddx * tmin;
				q.y += ddy * tmin;
				tracing.offerLast(new Point2D.Double(q.x, q.y));
				bounces.offerLast(new Point2D.Double(q.x, q.y));

				// Calculate the CCW (inward-pointing) perp vector for the side.
				double		vdx = pp2.x - pp1.x;	// Calc side vector
				double		vdy = pp2.y - pp1.y;	// from p1 to p2
				double		ndx = -vdy;			// Calc perp vector:
				double		ndy = vdx;				// negate y and swap

				// Need a NORMALIZED perp vector for the reflection calculation.
				double		nn = Math.sqrt(ndx * ndx + ndy * ndy);

				ndx = ndx / nn;	// Divide each coordinate by the length to
				ndy = ndy / nn;	// make a UNIT vector normal to the side.

				// Calculate v_reflected. See pages 148-149 and the slide on
				// "Reflecting Trajectories". (Note: P and v on the slide are
				// named q and dd here.)
				double		dot = dot(ddx, ddy, 0.0, ndx, ndy, 0.0);
				double		vreflectedx = ddx - 2.0 * dot * ndx;
				double		vreflectedy = ddy - 2.0 * dot * ndy;

				// Reflect the update vector, and reduce it to compensate for
				// the distance the point moved to reach the side.
				ddx = vreflectedx * (1.0 - tmin);
				ddy = vreflectedy * (1.0 - tmin);

				// Also reflect the reference vector. It will change direction
				// but remain the same length.
				double		dot2 = dot(dx, dy, 0.0, ndx, ndy, 0.0);

				dx -= 2.0 * dot2 * ndx;
				dy -= 2.0 * dot2 * ndy;
				model.setRealPowerX(dx);
				model.setRealPowerY(dy);
			}
		}
	}

	//**********************************************************************
	// Private Methods (Vectors)
	//**********************************************************************

	// This might be a method to calculate a dot product. Sure seems like it.
	private double		dot(double vx, double vy, double vz,
							double wx, double wy, double wz)
	{
		return (vx * wx + vy * wy + vz * wz);
	}

	// Determines if point q is to the left of line p1->p2. If strict is false,
	// points exactly on the line are considered to be left of it.
	private boolean	isLeft(Point2D.Double p1, Point2D.Double p2,
							   Point2D.Double q, boolean strict)
	{
		// Calculate the CCW (inward-pointing) perp vector for the side.
		double		vdx = p2.x - p1.x;		// Calc side vector
		double		vdy = p2.y - p1.y;		// from p1 to p2
		double		ndx = -vdy;			// Calc perp vector:
		double		ndy = vdx;				// negate y and swap

		// See the dot product on the slide on "Testing Containment in 2D".
		double		wdx = q.x - p1.x;		// Calculate test vector
		double		wdy = q.y - p1.y;		// from p1 to q
		double		dot = dot(wdx, wdy, 0.0, ndx, ndy, 0.0);

		// If strict, disallow cases of points exactly along the side's line.
		return (strict ? (dot > 0.0) : (dot >= 0.0));
	}

	// Determines if point q is inside a polygon. The polygon must be convex
	// with points stored in counterclockwise order. Points exactly on any side
	// of the polygon are considered to be outside of it.
	private boolean	contains(Deque<Point2D.Double> polygon,
								 Point2D.Double q)
	{
		int			sides = polygon.size();
		boolean		contains = true;
		Point2D.Double	p1 = polygon.peekLast();

		for (int i=0; i<sides; i++)
		{
			Point2D.Double	p2 = polygon.peekFirst();

			if (!isLeft(p1, p2, q, true))
				contains = false;

			polygon.offerLast(polygon.pollFirst());	// Cycle around polygon
			p1 = p2;
		}

		return contains;
	}
}

//******************************************************************************
