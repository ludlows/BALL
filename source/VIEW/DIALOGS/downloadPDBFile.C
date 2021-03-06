// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: downloadPDBFile.C,v 1.31.4.14 2005/12/20 12:46:21 amoll Exp $
//

#include <BALL/VIEW/DIALOGS/downloadPDBFile.h>

#include <BALL/SYSTEM/TCPTransfer.h>
#include <BALL/KERNEL/system.h>
#include <BALL/FORMAT/lineBasedFile.h>
#include <BALL/FORMAT/PDBFile.h>
#include <BALL/VIEW/KERNEL/mainControl.h>
#include <BALL/VIEW/KERNEL/message.h>
#include <BALL/VIEW/KERNEL/threads.h>

#include <qcombobox.h> 
#include <qlineedit.h> 
#include <qtextedit.h> 
#include <qfile.h>
#include <qurl.h>
#include <qradiobutton.h>
#include <qcheckbox.h>
#include <qimage.h>
#include <qpushbutton.h>
#include <qgroupbox.h>
#include <qapplication.h>

#ifdef BALL_HAS_SSTREAM
# include <sstream>
#else
# include <strstream>
#endif

namespace BALL
{
	namespace VIEW
	{

DownloadPDBFile::DownloadPDBFile(QWidget* parent, const char* name, bool modal, WFlags fl)
	throw()
	: DownloadPDBFileData(parent, name, modal, fl),
		ModularWidget(name),
		qb_(0),
		thread_(0),
		aborted_(false),
		error_(false)
{
#ifdef BALL_VIEW_DEBUG
	Log.error() << "new DownloadPDBFile" << this << std::endl;
#endif
	// register the widget with the MainControl
	registerWidget(this);
	hide();
//   	connect(results, SIGNAL(activated(const QString&)), this, SLOT(slotNewId(const QString&)));
	pdbId->setFocus();

#ifdef BALL_QT_HAS_THREADS
	thread_ = new FetchHTMLThread();
#endif
}

DownloadPDBFile::~DownloadPDBFile()
	throw()
{
#ifdef BALL_VIEW_DEBUG
	Log.info() << "Destructing object " << this << " of class DownloadPDBFile" << std::endl; 
#endif 

	if (qb_ != 0) delete qb_;

#ifdef BALL_QT_HAS_THREADS
	if (thread_ != 0) 
	{
		if (thread_->running())
		{
			thread_->terminate();
			thread_->wait();
		}
		delete thread_;
	}
#endif
}

void DownloadPDBFile::initializeWidget(MainControl&)
	throw()
{
	menu_id_ = insertMenuEntry(MainControl::FILE_OPEN, "Download PDB", this,
															 SLOT(show()), CTRL+Key_T, 1);
	setMenuHint("Download a PDB file from www.rcsb.org");
}

void DownloadPDBFile::slotSearch()
{
	/*
	results->clear();

	QString search_contents = searchField->text();
	QString filename;
	filename = "http://pdbbeta.rcsb.org/pdb/navbarsearch.do?newSearch=yes&isAuthorSearch=";
	if (authors->isChecked())
	{
		filename += "yes";
	}
	else
	{
		filename += "no";
	}

	filename += "&radioset=Structures&inputQuickSearch=";
	QUrl::encode(search_contents);
	filename += search_contents;
	filename += "&image.x=0&image.y=0";
	
	try 
	{
#ifndef BALL_QT_HAS_THREADS
		std::stringstream search_result;
		TCPTransfer tcp;
		tcp.set(search_result, filename.latin1());
		setProxyAndTransfer_(tcp);

		if (tcp.getStatusCode() != TCPTransfer::OK)
		{
			setStatusbarText(String("Failed to download file ") + filename.latin1() + 
											 ". ErrorCode " + String(tcp.getStatusCode()), true);
			return;
		}
#else   // =============================
		thread_->setFilename("");
		threadedDownload_(filename.latin1());
		downloadEnded_();
		if (error_) 
		{
			return;
		}

		std::stringstream& search_result = thread_->getStream();
#endif
		
		char buffer[10000];
		vector<String> result;
		while (!search_result.eof())
		{
			search_result.getline(buffer, 10000);
			String tmp(buffer);
			if (!tmp.hasSubstring("qrb_structid")) continue;

			Size pos = tmp.find(">");
			if (pos != string::npos)
			{
				result.push_back(tmp.substr(pos + 1, 4));
			}
		}

		for (Size i = 0; i < result.size(); i++)
		{
			results->insertItem(result[i].c_str());
		}

		if (result.size() != 0)
		{
			pdbId->setText(result[0].c_str());
		}
	}
	catch (...)
	{ 
		setStatusbarText("Could not download search result from PDB.org", true);
	}
	*/
}

bool DownloadPDBFile::threadedDownload_(const String& url)
{
	error_ = false;
#ifdef BALL_QT_HAS_THREADS
	downloadStarted_();
	thread_->setURL(url);
	thread_->start();
	Size last_bytes = 0;
	while (thread_->running())
	{
		qApp->processEvents();
		
		Size bytes = thread_->getTCPTransfer().getReceivedBytes();
		if (bytes != last_bytes)
		{
			setStatusbarText("received " + String(bytes) + " bytes", true);
			last_bytes = bytes;
		}
		
		thread_->wait(200);
	}

	if (aborted_) return false;

	if (thread_->getTCPTransfer().getReceivedBytes() == 0 ||
			thread_->getTCPTransfer().getStatusCode() != TCPTransfer::OK) 
	{
		if (thread_->getTCPTransfer().getReceivedBytes() == 0)
		{
			setStatusbarText("Could not download the given file. Maybe it does not exist on pdb.org?", true);
		}
		else
		{
			setStatusbarText(String("Failed to download file, ErrorCode ") + 
											 String(thread_->getTCPTransfer().getStatusCode()), true);
		}
		error_ = true;
		return false;
	}
#else
	// prevent compiler warnings
	url.isValid();
#endif

	return true;
}

void DownloadPDBFile::slotDownload()
{
	System *system = new System();

	try
	{
		String id = pdbId->text().latin1();
		String url = "http://pdbbeta.rcsb.org/pdb/displayFile.do?fileFormat=PDB&structureId=";
		url += id;
#ifndef BALL_QT_HAS_THREADS
		PDBFile pdb_file(url);
#else   // =============================
		String temp_filename = VIEW::createTemporaryFilename();
		thread_->setFilename(temp_filename);
		bool ok = threadedDownload_(url);
		downloadEnded_();

		if (!ok) 
		{
			delete system;
			return;
		}

		PDBFile pdb_file(temp_filename);
#endif

		pdb_file >> *system;
		pdb_file.close();

#ifdef BALL_QT_HAS_THREADS
 		removeFile_(temp_filename);
#endif

		if (system->countAtoms() == 0)
		{
			delete system;
			show();
			setStatusbarText("Could not fetch the given PDBFile", true);
			return;
		}
		else
		{
			setStatusbarText(String("read ") + String(system->countAtoms()) + " atoms for:  " + id, true);
		}

		if (system->getName() == "")
		{
			system->setName(pdbId->text().latin1());
		}

		system->setProperty("FROM_FILE", url);
		close();
		getMainControl()->insert(*system, pdbId->text().latin1());
		
		notify_(new CompositeMessage(*system, CompositeMessage::CENTER_CAMERA));
		
		download->setDefault(true);
		pdbId->setText("");
		pdbId->setFocus();
	}
	catch(...)
	{
		setStatusbarText("download PDB file failed", true);
		delete system;
	}
}

void DownloadPDBFile::slotShowDetail()
{
	if (pdbId->text() == "") return;

	setStatusbarText("Downloading information, please wait...", true);
	QString filename = "http://pdbbeta.rcsb.org/pdb/explore.do?structureId=";
	filename += pdbId->text();

	displayHTML(filename);
}

void DownloadPDBFile::slotNewId(const QString& new_id)
{
	pdbId->setText(new_id);
}

void DownloadPDBFile::displayHTML(const QString& url)
{
	if (url == "") return;

	try
	{
		QString filename;

		if (url.find("http://") == -1) filename = "http://www.rcsb.org/"+url;	
		else 													 filename = url;

		setStatusbarText(String("Reading ") + filename.ascii(), true);

#ifndef BALL_QT_HAS_THREADS
		std::stringstream search_result;
		TCPTransfer tcp;
		tcp.set(search_result, filename.latin1());
		setProxyAndTransfer_(tcp);

		if (tcp.getStatusCode() != TCPTransfer::OK)
		{
			setStatusbarText(String("Failed to download file ") + filename.latin1() + ". ErrorCode " + String(tcp.getStatusCode()), true);
			return;
		}
		Size size = tcp.getReceivedBytes();
#else
		thread_->setFilename("");
		threadedDownload_(filename.ascii());
		if (aborted_) return;

		std::stringstream& search_result = thread_->getStream();
		Size size = thread_->getTCPTransfer().getReceivedBytes();
#endif

		if (size == 0)
		{
			setStatusbarText(String("URL ") + filename.ascii() + " does not exist.", true);
			error_ = true;
			downloadEnded_();
			return;
		}

		setStatusbarText("Please wait, while loading images...", true);

		List<String> images;
		String result;
		char buffer[10000];

		while (!search_result.eof())
		{
			search_result.getline(buffer, 10000);
			String current_line(buffer);
			result += current_line + "\n";

			// find out all the images
			String tmp = current_line;
			tmp.toUpper();
			Size pos_1 = 0;
			Size pos_2 = 0;

//   			while ( (pos_1 = tmp.find("<IMG", pos_2)) != string::npos )
			while ( false)
			{
				pos_1 = tmp.find("SRC=\"", pos_1);
				pos_2 = current_line.find("\"", pos_1+9);

				String img_url = current_line.substr(pos_1+5, pos_2 - (pos_1+5));
				images.push_back(img_url);
				if (!unsupported_images_.has(img_url) && !image_cache_.has(img_url))
				{
					String tmp_filename = VIEW::createTemporaryFilename();
					File img(tmp_filename, std::ios::out);
					TCPTransfer tcp;
					tcp.set(img, "http://www.rcsb.org/" + img_url);
					setProxyAndTransfer_(tcp);
					img.close();

					QImage qi;
					qi.load(tmp_filename.c_str());
 					removeFile_(tmp_filename);

					if (qi.isNull())
					{
						unsupported_images_.insert(img_url);
					}
					else
					{
						image_cache_[img_url] = qi;
					}
				}
			}
		}

 		if (qb_ == 0) 
		{
			qb_ = new QTextBrowser();
			connect(qb_, SIGNAL(linkClicked(const QString&)), this, SLOT(displayHTML(const QString&)));
		}

		QImage empty;
		List<String>::Iterator it = images.begin();
		for (; it != images.end(); it++)
		{
			if (!unsupported_images_.has(*it))
			{
				qb_->mimeSourceFactory()->setImage((*it).c_str(), image_cache_[*it]);
			}
			else
			{
				qb_->mimeSourceFactory()->setImage((*it).c_str(), empty);
			}
		}

		qb_->setText(QString(result.c_str()));

		downloadEnded_();
		setStatusbarText("Finished download of HTML page", true);

		qb_->showMaximized();
		qb_->show();
	}
	catch (...)
	{ 
		downloadEnded_();
		setStatusbarText("Failed to download HTML page.", true);
	}

#ifdef BALL_QT_HAS_THREADS
	removeFile_(thread_->getFilename());
#endif
} 


void DownloadPDBFile::idChanged()
{
	download->setEnabled(pdbId->text() != "");
}

void DownloadPDBFile::abort()
{
#ifdef BALL_QT_HAS_THREADS

	if (thread_ == 0) return;
	aborted_ = true;

	thread_->abort();
	thread_->wait(5500);

	if (thread_->running())
	{
		thread_->terminate();
		thread_->wait();
	}
	removeFile_(thread_->getFilename());
	
	downloadEnded_();

#endif
}

void DownloadPDBFile::downloadStarted_()
	throw()
{
	aborted_ = false;
	error_   = false;
	setStatusbarText("Starting download ...", true);
	button_abort->setEnabled(true);
	download->setEnabled(false);
	pdbId->setEnabled(false);
//   	search_box->setEnabled(false);			
//   	showDetails->setEnabled(false);
	buttonClose->setEnabled(false);
}

void DownloadPDBFile::downloadEnded_()
	throw()
{
	if (!aborted_ && !error_)
	{
		setStatusbarText("Finished downloading, loading structure ...", true);
	}
	button_abort->setEnabled(false);
	download->setEnabled(true);
	pdbId->setEnabled(true);
//   	search_box->setEnabled(true);			
//   	showDetails->setEnabled(true);
	buttonClose->setEnabled(true);
	idChanged();
	qApp->processEvents();
	pdbId->setFocus();

#ifdef BALL_QT_HAS_THREADS
	if (error_)
	{
		removeFile_(thread_->getFilename());
	}
#endif
}

void DownloadPDBFile::removeFile_(const String& filename)
{
	try
	{
 		File::remove(filename);
	}
	catch(...) {}
}

void DownloadPDBFile::setProxyAndTransfer_(TCPTransfer& tcp)
{
	MainControl* mc = getMainControl();
 	tcp.setProxy(mc->getProxy(), mc->getProxyPort());
	tcp.transfer();
}

	}
}
