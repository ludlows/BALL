// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: pyWidget.C,v 1.44.6.25 2005/12/14 14:42:05 amoll Exp $
//

// This include has to be first in order to avoid collisions.
#include <Python.h>

#include <BALL/VIEW/WIDGETS/pyWidget.h>
#include <BALL/VIEW/KERNEL/mainControl.h>
#include <BALL/VIEW/DIALOGS/pythonSettings.h>
#include <BALL/PYTHON/pyInterpreter.h>
#include <BALL/VIEW/DIALOGS/preferences.h>
#include <BALL/FORMAT/lineBasedFile.h>

#include <qscrollbar.h>
#include <qfiledialog.h>
#include <qapplication.h>
#include <qdragobject.h>

// currently doesnt work right
//#undef BALL_QT_HAS_THREADS


namespace BALL
{
	namespace VIEW
	{


#ifdef BALL_QT_HAS_THREADS
		RunPythonThread::RunPythonThread()
			throw()
			: QThread(),
				input(),
				output()
		{}

		void RunPythonThread::run()
		{
			output = PyInterpreter::run(input, state);
		}
#endif


		const Hotkey& Hotkey::operator = (const Hotkey& hotkey)
			throw()
		{
			action 				= hotkey.action;
			key 	 				= hotkey.key;
			button_state 	= hotkey.button_state;
			return *this;
		}


		bool Hotkey::operator == (const Hotkey& hotkey) const
			throw()
		{
			return (action 				== hotkey.action &&
							key 	 				== hotkey.key    &&
							button_state 	== hotkey.button_state);
		}

		bool Hotkey::operator == (const QKeyEvent& e) const
			throw()
		{

			return (key 	 				== e.key() 		&&
							button_state 	== e.state()	);
		}

		bool Hotkey::set(const String& data) 
			throw()
		{
			vector<String> fields;
			Size nr = data.split(fields, String('#').c_str());
			if (nr < 3)
			{
				Log.error() << "Could not parse Hotkey " << data << std::endl;
				return false;
			}
			
			try
			{
				key = (Qt::Key) fields[0].toUnsignedInt();
				button_state = (Qt::ButtonState) fields[1].toUnsignedInt();
				action = fields[2];
			}
			catch(...)
			{
				Log.error() << "Could not parse Hotkey " << data << std::endl;
				return false;
			}

			return true;
		}


		void Hotkey::get(String& data) const 
			throw()
		{
			data = String(key) + "#" + String(button_state) + "#" + action;
		}

		// ==========================================================

		PyWidgetData::PyWidgetData(QWidget* parent, const char* name)
			: QTextEdit(parent, name),
				thread_(0),
				stop_script_(false)
		{
			setWrapPolicy(QTextEdit::Anywhere);
			#ifdef BALL_QT_HAS_THREADS
			  thread_ = new RunPythonThread();
			#endif
		}

		PyWidgetData::PyWidgetData(const PyWidgetData& /*widget*/)
			:	QTextEdit(),
				stop_script_(false)
		{
		}

		PyWidgetData::~PyWidgetData()
			throw()
		{	
#ifdef BALL_QT_HAS_THREADS
			if (thread_->running())
			{
				thread_->wait(500);
				thread_->terminate();
			}
			delete thread_;
#endif
		}


		void PyWidgetData::startInterpreter()
		{
			stop_script_ = false;
			// initialize the interpreter
			PyInterpreter::initialize();

			if (!PyInterpreter::isValid())
			{
				return;
			}

			// print the PyBALL version and clear
			// the widget's contents in case of a restart
			setText((String("BALL ") + VersionInfo::getVersion()).c_str());

			// print the first prompt 
			multi_line_mode_ = false;
			newPrompt_();
			history_position_ = history_.size() + 1;
		}


		void PyWidgetData::retrieveHistoryLine_(Position index)
		{
			if (index > history_.size()) 
			{
				history_position_ = history_.size();
				return;
			}

			int row, col;
			getCursorPosition(&row, &col);

			if (index == history_.size())
			{
				history_position_ = history_.size();
				removeParagraph(row);
				insertParagraph(getPrompt_(), row);
				setCursorPosition(paragraphs()-1, col);
				QScrollBar* sb = verticalScrollBar();
				if (sb != 0)
				{
					sb->setValue(sb->maxValue());
				}
				removeParagraph(row +1);
				return;
			}

			String line = getPrompt_()+ history_[index];

			// replace the line's contents
			removeParagraph(row);
			insertParagraph(line.c_str(), row);
			setCursorPosition(row, line.size()); 
			QScrollBar* sb = verticalScrollBar();
			if (sb != 0)
			{
				sb->setValue(sb->maxValue());
			}

			// update the history position
			history_position_ = index;
			removeParagraph(row +1);
		}


		bool PyWidgetData::returnPressed()
		{
			// check for an empty line (respect the prompt)
			int row, col;
			getCursorPosition(&row, &col);
			current_line_ = getCurrentLine_();
			QTextEdit::returnPressed();
			if (col < 5)
			{
				if (multi_line_mode_)
				{
					// in multi line mode: end of input - parse it!
					parseLine_();
				}
				else	
				{
					// return on an empty line is handled 
					// as in the interactive interpreter: do nothing and 
					// print another prompt
					newPrompt_();
				}
			} 
			else 
			{	
				// parse the line
				parseLine_();
			}

			return false;
		}


		bool PyWidgetData::parseLine_(String line, bool silent)
		{
			if (!Py_IsInitialized())
			{
				append("ERROR: no interpreter running!\n");
				return false;
			}

			// check for comments
			String temp(line);
			temp.trim();
			if (temp.hasPrefix("#")) return true;

			history_position_ = history_.size();

			line.trimRight();

			if (line == "quit" || line == "quit()")
			{
				stop_script_ = true;
				PyInterpreter::finalize();
				getMainControl()->quit();
				return true;
			}

			if (line.hasPrefix("run("))
			{
				// This code is probably the worst I've seen in a long, long time...
				// Needs to be replaced by something integrated with Py language concepts
				// in the future (e.g. calling run with a variable argument won't work!) -- OK -- 11/2006
				vector<String> tokens;
				Size nr = line.split(tokens, String("(\"')").c_str());
				if (nr < 2)		
				{
					return false;
				}
				PyInterpreter::runFile(tokens[1]);
				appendToHistory_(line);
				
				return true;
			}
			
			if (!multi_line_mode_)
			{
				if (line.isEmpty()) return true;
				if ((line.hasPrefix("for ") 		|| 
						 line.hasPrefix("def ") 		||
						 line.hasPrefix("class ") 	||
						 line.hasPrefix("while ")	  ||
						 line.hasPrefix("if "))
						&& line.hasSuffix(":"))
				{
						multi_line_mode_ = true;
						multi_line_text_ = line;
						multi_line_text_.append("\n");
						multi_lines_ = 1;
						appendToHistory_(line);
						if (!silent)
						{
							newPrompt_();
						}
						return true;
				}

				multi_lines_ = 0;
				appendToHistory_(line);
			}
			else // Multiline mode
			{
				multi_lines_ += 1;
				appendToHistory_(line);
				if (!line.isEmpty())
				{
					multi_line_text_ += line + "\n";
					if (!silent)
					{
						newPrompt_();	
					}
					return true;
				}

				line = multi_line_text_ + "\n";
			}

			bool state = true;
			bool ok = true;
		  String result;
#ifndef BALL_QT_HAS_THREADS
			result = PyInterpreter::run(line, state);
#else
			thread_->input = line;
			thread_->start();

			while (thread_->running()) 
			{
				qApp->wakeUpGuiThread();
				qApp->processEvents(500);
				thread_->wait(50);
			}

			result = thread_->output;
			state = thread_->state;
#endif

			if (result != "") 
			{
				if (result.hasSubstring("ERROR"))
				{
					setColor(red);
					ok = false;
				}
				else
				{
					setColor(blue);
				}
				append(result.c_str());
				setColor(black);
			}

				
			if (!multi_line_mode_)
			{
				results_.push_back(state);
			}
			else
			{
				for (Position p = 0; p <= multi_lines_ -1; p++)
				{
					results_.push_back(state);
				}
			}
			
			multi_line_mode_ = false;

			if (!silent)
			{
				newPrompt_();
			}
			return ok;
		}

		bool PyWidgetData::runString(String command)
		{
			if (!command.has('\n'))
			{
				return parseLine_(command);
			}

			vector<String> lines;
			Size nr = command.split(lines, String('\n').c_str());
			for (Position p = 0; p < nr; p++)
			{
				if (!parseLine_(lines[p])) return false;
			}

			return true;
		}
		

		void PyWidgetData::appendToHistory_(const String& line)
		{
			history_.push_back(line);
			history_position_ = history_.size();
		}

		const char* PyWidgetData::getPrompt_() const
		{
			return (multi_line_mode_ ? "... " : ">>> ");
		}

		void PyWidgetData::newPrompt_()
		{
			append(getPrompt_());
			setCursorPosition(paragraphs() - 1, 4);
		}

		bool PyWidgetData::parseLine_()
		{
			String line = current_line_.getSubstring(4);
			return parseLine_(line);
		}


		void PyWidgetData::keyPressEvent(QKeyEvent* e)
		{
			int row, col;
			getCursorPosition(&row, &col);

			if (row != paragraphs()-1 || col < 4)
			{
				setCursorPosition(paragraphs()-1, paragraphLength(paragraphs() -1) -1);
			}

			if (e->key() == Key_Left || e->key() == Key_Backspace)
			{
				if (col <= 4) return;
			}
			else if (e->key() == Key_Right)
			{
				setCursorPosition(paragraphs()-1, col+1);
				return;
			}
			else if (e->key() == Key_Up)
			{
				if (history_position_ != 0) retrieveHistoryLine_(history_position_ - 1);
				return;
			}
			else if (e->key() == Key_Down)
			{
				retrieveHistoryLine_(history_position_ + 1);
				return;
			}
			else if (e->key() == Key_Home)
			{
				setCursorPosition(row, 4);
				return;
			}
			else if (e->key() == Key_Return)
			{
				if (!returnPressed()) return;
			}
			else if (e->key() == Key_PageUp)
			{
				return;
			}
			else if (e->key() == Key_PageDown)
			{
				return;
			}

			QTextEdit::keyPressEvent(e);
		} 


		void PyWidgetData::dump(std::ostream& s, Size depth) const
			throw()
		{
			BALL_DUMP_STREAM_PREFIX(s);

			BALL_DUMP_HEADER(s, this, this);

			BALL_DUMP_DEPTH(s, depth);
			s << "multiline_mode : " << multi_line_mode_<< std::endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "multi_line_text  : " << multi_line_text_<< std::endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "history : "<< std::endl;
			
			for (Position i = 0; i < history_.size(); i++)
			{
				BALL_DUMP_DEPTH(s, depth);
				s << history_[i]<< std::endl;
			}

			BALL_DUMP_DEPTH(s, depth);
			s << "history_position : " << history_position_ << std::endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "current_line : " << current_line_ << std::endl;

			BALL_DUMP_STREAM_SUFFIX(s);
		}


		String PyWidgetData::getCurrentLine_()
		{
			int row, col;
			getCursorPosition(&row, &col);
			return String(text(row).ascii());
		}


		bool PyWidgetData::runFile(const String& filename)
		{
			stop_script_ = false;
			append(String("> executing script from " + filename + "\n").c_str());
			LineBasedFile file;
			try
			{
				file.open(filename);
			}
			catch	(...)
			{
				append(String("> Could not find file " + filename + "\n").c_str());
				newPrompt_();
				return false;
			}

			while (file.readLine())
			{
				// Call parse line with 'silent = true' to make sure we do not get prompts
				// for each line read.
				if (!parseLine_(file.getLine(), true))
				{
					String result_string = "> Error in line " + String(file.getLineNumber()) + " of file " + filename + "\n";
					append(result_string.c_str());
					newPrompt_();
					return false;
				}

				if (stop_script_) 
				{
					stop_script_ = false;
 					((PyWidget*)parent())->setStatusbarText("Aborted script");
					append("> aborted...");
					newPrompt_();
					return false;
				}
			}
			append("> Finished.");
			((PyWidget*)parent())->setStatusbarText("Finished script.");
			newPrompt_();
			return true;
		}

		void PyWidget::scriptDialog()
		{
			if (working_dir_ == "") working_dir_ = getWorkingDir();

			QString s = QFileDialog::getOpenFileName(
										working_dir_.c_str(),
										"Python Scripts(*.py)",
										getMainControl(),
										"Run Python Script",
										"Choose a Python script" );

		 	if (s == QString::null) return;
			setWorkingDirFromFilename_(s.ascii());
			working_dir_ = getWorkingDir();

			run(s.ascii());
		}

		void PyWidgetData::exportHistory()
		{
			PyWidget* p = (PyWidget*) parent();
			QString s = QFileDialog::getSaveFileName(
										p->getWorkingDir().c_str(),
										"",
										p->getMainControl(),
										"Export History",
										"");

		 	if (s == QString::null) return;
			String filename(s.ascii());
			p->setWorkingDirFromFilename_(filename);

			File file(filename, std::ios::out);
			if (!file.isOpen()) 
			{
				append(String("> Could not export history to file " + filename + "\n").c_str());
				newPrompt_();
				return;
			}
					
			for (Position p = 0; p < history_.size(); p++)
			{
				if (results_[p]) file << history_[p] << std::endl;
			}

			file.close();
		}

		void PyWidgetData::clear()
		{
			QTextEdit::clear();
			newPrompt_();
		}

		void PyWidgetData::paste()
		{
			int row, col;
			getCursorPosition(&row, &col);

			if (row != paragraphs()-1)
			{
				setCursorPosition(paragraphs()-1, paragraphLength(paragraphs() -1) -1);
			}

			QTextEdit::paste();
		}


		void PyWidgetData::contentsDragEnterEvent(QDragEnterEvent * e)
		{
			e->accept(QTextDrag::canDecode(e));
		}

		void PyWidgetData::contentsDropEvent(QDropEvent *e)
		{
			VIEW::processDropEvent(e);
		}

// ######################################################################################################

		PyWidget::PyWidget(QWidget *parent, const char *name)
			throw()
			: DockWidget(parent, name),
				text_edit_(new PyWidgetData(this)),
				working_dir_(""),
				valid_(false),
				started_startup_script_(false)
		{
		#ifdef BALL_VIEW_DEBUG
			Log.error() << "new PyWidget " << this << std::endl;
		#endif
			default_visible_ = false;
			setGuest(*text_edit_);
			registerWidget(this);
		}

		PyWidget::~PyWidget()
			throw()
		{}

		void PyWidget::initializeWidget(MainControl& main_control)
			throw()
		{
//   			insertMenuEntry(MainControl::TOOLS_PYTHON, "Restart Python", text_edit_, SLOT(startInterpreter()));

			DockWidget::initializeWidget(main_control);
			registerWidgetForHelpSystem(this, "pythonInterpreter.html");

			Index id1 = insertMenuEntry(MainControl::TOOLS_PYTHON, "Run Python Script", this , SLOT(scriptDialog()));
			Index id2 = insertMenuEntry(MainControl::TOOLS_PYTHON, "Abort Python Script", text_edit_, SLOT(abortScript()));
			Index id3 = insertMenuEntry(MainControl::TOOLS_PYTHON, "Export History", text_edit_, SLOT(exportHistory()));

			text_edit_->startInterpreter();

			valid_ = PyInterpreter::isValid();
	
			if (!valid_)
			{
				menuBar()->setItemEnabled(id1, false);
				menuBar()->setItemEnabled(id2, false);
				menuBar()->setItemEnabled(id3, false);
				text_edit_->setText("No Python support available:");
				text_edit_->runString("import BALL");
				text_edit_->setEnabled(false);
				setStatusbarText("No Python support available! (See PyWidget)", true);
				return;
			}
		}


		void PyWidget::finalizeWidget(MainControl& main_control)
			throw()
		{
			text_edit_->abortScript();
			PyInterpreter::finalize();

			DockWidget::finalizeWidget(main_control);
		}


		void PyWidget::initializePreferencesTab(Preferences &preferences)
			throw()
		{
			text_edit_->python_settings_= new PythonSettings();
			preferences.insertEntry(text_edit_->python_settings_);
		}

		void PyWidget::finalizePreferencesTab(Preferences &preferences)
			throw()
		{
			if (text_edit_->python_settings_ != 0)
			{
				preferences.removeEntry(text_edit_->python_settings_);
				text_edit_->python_settings_ = 0;
			}
		}

		void PyWidget::applyPreferences()
			throw()
		{
			DockWidget::applyPreferences();

			if (text_edit_->python_settings_ == 0)
			{
				return;	
			}

 			hotkeys_ = (text_edit_->python_settings_->getContent());

			if (started_startup_script_ || !isValid())
			{
				return;
			}

			started_startup_script_ = true;

			String startup = getDataPath() + "startup.py";
			if (!text_edit_->runFile(startup))
			{
				Log.error() << "Could not find startup script. Please set the correct path to the data path!" << std::endl;
				Log.error() << "To do so set the environment variable BALL_DATA_PATH or BALLVIEW_DATA_PATH." << std::endl;
			}
			
			String user_startup = text_edit_->python_settings_->getFilename();
			if (user_startup == "") 
			{
				return;	
			}

			text_edit_->runFile(user_startup);
		}

		void PyWidgetData::abortScript()
		{
			((PyWidget*)(parent()))->setStatusbarText("Aborting Python script");
			 stop_script_ = true;
		}

		PyWidget::PyWidget(const PyWidget& p)
			throw()
		 : DockWidget(p)
		{}

		bool PyWidget::toAbortScript() throw() 
		{
			return text_edit_->stop_script_;
		}

		void PyWidget::insertHotkey(const Hotkey& hotkey) 
			throw()
		{
			List<Hotkey>::iterator it = hotkeys_.begin();
			for (; it != hotkeys_.end(); it++)
			{
				if ((*it) == hotkey) return;
			}

			hotkeys_.push_back(hotkey);
		}

		void PyWidget::removeHotkey(const Hotkey& hotkey) 
			throw()
		{
			List<Hotkey>::iterator it = hotkeys_.begin();
			for (; it != hotkeys_.end(); it++)
			{
				if ((*it) == hotkey) return;
			}

			hotkeys_.erase(it);
		}

		void PyWidget::reactTo(const QKeyEvent& e) 
			throw() 
		{
			// doesnt work, no idea yet why:
			/*
			if (getMainControl()->compositesAreLocked() ||
					getMainControl()->getPrimitiveManager().updateRunning())
			{
				return;
			}
			*/

			List<Hotkey>::iterator it = hotkeys_.begin();
			for (; it != hotkeys_.end(); it++)
			{
				if ((*it) == e) 
				{
					text_edit_->runString((*it).action);
					return;
				}
			}
		}

		bool PyWidget::run(const String& filename) throw() 
		{
			last_script_ = filename;
			return text_edit_->runFile(filename);
		}
			
		bool PyWidget::runAgain()
		{
			return run(last_script_);
		}

	} // namespace VIEW
} // namespace BALL
