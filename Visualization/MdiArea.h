// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

// DOC -> https://doc.qt.io/archives/qt-4.8/dnd.html

#include <QMdiArea>

#include <QDragEnterEvent>

namespace FEVV {

class SimpleWindow;

enum { bNone, bLeft, bRight };

/*! 
 * \class MdiArea
 * \brief MdiArea class.
 */
class MdiArea : public QMdiArea
{
	//Q_OBJECT

	public:
		/*!
		 * \brief Constructor.
		 *
		 * \param parent : parent window.
		 */
		MdiArea(QWidget *parent);
		/*!
		 * \fn setMainWindow(mainwindow* mw)
		 * \brief Set mainwindow pointer to mdiarea.
		 *
		 * \param mw mainwindow pointer.
		 */
		void setMainWindow(SimpleWindow* mw) { m_mw = mw; }

	protected:
		/*!
		 * \fn dragEnterEvent(QDragEnterEvent *event)
		 * \brief Accept a drop mesh file in this zone.
		 *
		 * \param event an event.
		 */
		void dragEnterEvent(QDragEnterEvent *event);
		/*!
		 * \fn dropEvent(QDropEvent *event)
		 * \brief Load drop mesh file.
		 *
		 * \param event an event.
		 */
		void dropEvent(QDropEvent *event);

	private:
		SimpleWindow* m_mw;	//!< mainwindow
		int bType;			//!< mouse bouton type (bNone, bLeft, bRight)
};

} // namespace FEVV

#ifndef Q_MOC_RUN
#include "Visualization/MdiArea.inl"
#endif // Q_MOC_RUN
