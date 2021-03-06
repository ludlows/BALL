// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: QTTimer.C,v 1.17 2003/08/29 15:37:56 amoll Exp $


#include <BALL/VIEW/KERNEL/QTTimer.h>
#include <BALL/COMMON/exception.h>
#include <BALL/DATATYPE/string.h>
#include <BALL/COMMON/rtti.h>
#include <stdio.h>

using namespace std;

namespace BALL
{

	namespace VIEW
	{

		QTTimer::NoValidInterval::NoValidInterval(const char* file, int line, const string& data)
				throw()
			:	Exception::GeneralException(file, line, string("NoValidInterval"), string("not a valid interval: ") + data)
		{
		}

		QTTimer::QTTimer(QObject* parent, const char* name)
				throw()
			:	QObject(parent, name),
			  timer_id_(0),
				interval_(100),
				timer_enabled_(false),
				timer_running_(false)
		{
		}

		QTTimer::QTTimer(const QTTimer &timer, QObject* parent, const char* name)
				throw()
			:	QObject(parent, name),
			  timer_id_(0),
				interval_(timer.interval_),
				timer_enabled_(false),
				timer_running_(false)
		{
		}

		QTTimer::~QTTimer()
				throw()
		{
#ifdef BALL_VIEW_DEBUG
				Log.info() << "Destructing object " << (void *)this 
					<< " of class " << RTTI::getName<QTTimer>() << endl;
#endif 

			clear();
		}

		void QTTimer::clear()
				throw()
		{
			stopTimer();

			timer_running_ = false;
			interval_ = 100;
		}

		void QTTimer::set(const QTTimer& timer)
				throw()
		{
			stopTimer();

			timer_running_ = false;
			interval_ = timer.interval_;
		}

		const QTTimer& QTTimer::operator =(const QTTimer& timer)
				throw()
		{
			set(timer);

			return *this;
		}

		void QTTimer::swap(QTTimer& timer)
				throw()
		{
			stopTimer();
			timer.stopTimer();

			int i = interval_;
			interval_ = timer.interval_;
			timer.interval_ = i;
		}

		void QTTimer::setInterval(int interval)
				throw(QTTimer::NoValidInterval)
		{
			if (interval <= 0)
			{
				char temp[20];

				sprintf(&temp[0], "%d", interval);
				throw NoValidInterval(__FILE__, __LINE__, temp);
			}

			interval_ = interval;
		}

		void QTTimer::timer()
		{
			return;
		}

		void QTTimer::timerEvent(QTimerEvent *timer_event)
		{
			if (timer_event->timerId() == timer_id_)
			{
				if (timer_running_ == false)
				{
					timer_running_ = true;
					timer();
					timer_running_ = false;
				}
			}
		}

		void QTTimer::dump(ostream& s, Size depth) const
				throw()
		{
			BALL_DUMP_STREAM_PREFIX(s);
			
			BALL_DUMP_DEPTH(s, depth);
			BALL_DUMP_HEADER(s, this, this);

			BALL_DUMP_DEPTH(s, depth);
			s << "interval: " << interval_ << endl;

			BALL_DUMP_DEPTH(s, depth);
			s << "timer running: " 
					 << ((timer_enabled_ == true) ? "yes" : "no") 
					 << endl;

			BALL_DUMP_STREAM_SUFFIX(s);
		}

#		ifdef BALL_NO_INLINE_FUNCTIONS
#			include <BALL/VIEW/KERNEL/QTTimer.iC>
#		endif 

	} // namespace VIEW

} // namespace BALL
