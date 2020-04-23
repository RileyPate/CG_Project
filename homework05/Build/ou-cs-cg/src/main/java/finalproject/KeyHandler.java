//******************************************************************************
// Copyright (C) 2016-2019 University of Oklahoma Board of Trustees.
//******************************************************************************
// Last modified: Mon Mar 18 23:14:07 2019 by Chris Weaver
//******************************************************************************
// Major Modification History:
//
// 20160225 [weaver]:	Original file.
// 20190227 [weaver]:	Updated to use model and asynchronous event handling.
// 20190318 [weaver]:	Modified for homework04.
//
//******************************************************************************
// Notes:
//
//******************************************************************************

package finalproject;

//import java.lang.*;
import java.awt.Component;
import java.awt.event.*;
import java.awt.geom.Point2D;

import utilities.Utilities;

//******************************************************************************


public final class KeyHandler extends KeyAdapter
{
	//**********************************************************************
	// Private Members
	//**********************************************************************

	// State (internal) variables
	private final View		view;
	private final Model	model;

	//**********************************************************************
	// Constructors and Finalizer
	//**********************************************************************

	public KeyHandler(View view, Model model)
	{
		this.view = view;
		this.model = model;

		Component	component = view.getCanvas();

		component.addKeyListener(this);
	}

	//**********************************************************************
	// Override Methods (KeyListener)
	//**********************************************************************

	public void	keyPressed(KeyEvent e)
	{
		boolean	b = Utilities.isShiftDown(e);
		int		number = model.getNumber();
		double		factor = model.getFactor();

		switch (e.getKeyCode())
		{
			case KeyEvent.VK_C:	// Move object to center
				model.setObjectInSceneCoordinates(new Point2D.Double(0.0, 0.0));
				return;

			case KeyEvent.VK_Q:	// Sort through golf courses.
				model.setNumber(Math.max(number - 1, 2));
				break;

			case KeyEvent.VK_W:	// Sort through golf courses
				model.setNumber(Math.min(number + 1, View.MAX_SIDES));
				break;

		}
	}
}

//******************************************************************************
