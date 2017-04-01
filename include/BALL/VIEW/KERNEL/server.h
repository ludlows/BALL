// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: server.h,v 1.11.6.2 2005/09/29 14:01:33 amoll Exp $
//

#ifndef BALL_VIEW_KERNEL_SERVER_H
#define BALL_VIEW_KERNEL_SERVER_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_DATATYPE_HASHMAP_H
#	include <BALL/DATATYPE/hashMap.h>
#endif

#ifndef BALL_VIEW_KERNEL_QTTIMER_H
#	include <BALL/VIEW/KERNEL/QTTimer.h>
#endif

#ifndef BALL_VIEW_KERNEL_MODULARWIDGET_H
#	include <BALL/VIEW/KERNEL/modularWidget.h>
#endif

#ifndef BALL_VIEW_KERNEL_OBJECTCREATOR_H
# include <BALL/VIEW/KERNEL/objectCreator.h>
#endif

class QLabel;

namespace BALL
{
	class Composite;
	class IOStreamSocket;
	class SockInetBuf;

	namespace VIEW
	{
		class ServerPreferences;
		class Preferences;

		/** Server class.
				The class Server handles all incoming PersistentObject objects,
				converts them into Composite objects (if possible) and sents
				them through the ConnectionObject tree with the message
				NewCompositeMessage. Also it stores all received Composite objects
				and replaces them if the same Composite object is received again.
				If a Composite object is replaced the message 
				RemovedCompositeMessage will be sent through the ConnectionObject
				tree and after that the the message NewCompositeMessage with the new
				received composite will be sent.
			\ingroup ViewKernelClient
		*/
		class BALL_VIEW_EXPORT Server
			: public QTTimer,
   			public ModularWidget
		{
			public:
			
			BALL_EMBEDDABLE(Server,ModularWidget)

			/**	@name	Constructors
			*/	
			//@{

			/** Default Constructor.
					The state of this server is:
					  - no object creator registered
						- server listening on <tt> VIEW_DEFAULT_PORT</tt> if activated
					\par
					\see         QTTimer
					\see         ModularWidget
			*/
			Server(QWidget* parent = 0, const char* name = 0)
				throw();

			//@}
			/** @name Destructors 
			*/
			//@{

			/** Destructor.
			*/
			virtual ~Server()
				throw();

			/** Explicit default initialization.
					Calls QTTimer::clear.
					Calls ConnectionObject::clear.
					\see QTTimer::clear
					\see ConnectionObject::clear
			*/
			virtual void clear()
				throw();
			//@}

			/**	@name	Exceptions
			*/
			//@{
			
			/** NotCompositeObject Exception class.
					This exeption will be thrown if a PersistentObject was received
					that was not a Composite object.
					\see         GeneralException			
			*/
			class BALL_VIEW_EXPORT NotCompositeObject:	public Exception::GeneralException
			{
				public:

				NotCompositeObject(const char* file, int line)
					throw();
			};

			//@}
			/**	@name	Accessors: inspectors and mutators 
			*/
			//@{

			/** Activates the server.
					Creates a new socket stream with the given port and enables the timer
					that will check
					every second whether an object will be available at the stream.
					After this method the <b> timer</b> method will be called every second.
					Must be called before other methods!
					Calls QTTimer::startTimer
					\see QTTimer::startTimer
					\see timer
			*/
			void activate()
				throw();

			/** Deactivates the server.
					If this server is already running this method stops the server
					and closes the socket stream.
					Calls QTTimer::stopTimer
					\see QTTimer::stopTimer
			*/
			void deactivate()
				throw();

			/**	Set the server port.
					Set port of this server. Must be called before activate
					to have any effect.
					\param  port the new port
			*/
			void setPort(const int port)
				throw();

			/**	Return the server port.
					Return the port of this server.
					\return int the port of this server
			*/
			int getPort() const
				throw();

			/** Register a ObjectCreator that is used for converting 
					PersistentObject objects into Composite objects.
					Every ObjectCreator, that is still registered, when a Server instance is destructed, will be deleted.
					\see ObjectCreator
			*/
			void registerObjectCreator(const ObjectCreator& s)
				throw();

			/** Reset the ObjectCreator.
					After calling this method PersistentObject objects will be converted
					using the default ObjectCreator.
					\see ObjectCreator
			*/
			void unregisterObjectCreator()
				throw();

			/**	Initialize the server widget.
					This method initializes the icon of this server and adds it
					to MainControl. This method will be called by show of 
					the class MainControl.
				  \see  ModularWidget
					\see  show
			*/
			virtual void initializeWidget(MainControl& main_control)
				throw();
			
			/**	Remove the server widget.
					This method deletes the icon of this server and removes it
					from MainControl.
					This method will be called by aboutToExit of 
					the class MainControl.
				  \see  ModularWidget
					\see  aboutToExit
			*/
			virtual void finalizeWidget(MainControl& main_control)
				throw();
			
			/** Initialize a preferences tab for the server.
					This method creates the preferences tab ServerPreferences for
					this server and inserts it into the Preferences dialog
					of the MainControl.
					This method is called automatically by the method show of 
					the class MainControl at the start of the application.
					See ModularWidget	for more information concerning preferences tabs.\par
					\param  preferences the Preferences dialog of the MainControl
					\see    show
					\see    ServerPreferences
					\see    Preferences
			*/
			virtual void initializePreferencesTab(Preferences &preferences)
				throw();
			
			/**	Remove the preferences tab.
					This method removes the ServerPreferences tab of this server
					from the Preferences dialog of the MainControl.
					This method is called automatically by the method aboutToExit
					method  of the class MainControl at the end of the application.
					See ModularWidget
					for more information concerning preferences tabs.\par
					\param  preferences the Preferences dialog of the MainControl
					\see    aboutToExit
					\see    ServerPreferences
					\see    Preferences
			*/
			virtual void finalizePreferencesTab(Preferences &preferences)
				throw();
			
			/** Apply the preferences of the specific tab.
					This method applies the preferences of the own tab ServerPreferences
					to this server.
					This method is called automatically by the method applyPreferencesTab of 
					the class MainControl.
					See ModularWidget	for more information concerning preferences tabs.\par
					\param  preferences the Preferences dialog of the MainControl
					\see    applyPreferencesTab
					\see    ServerPreferences
					\see    Preferences
			*/
			virtual void applyPreferences()
				throw();

			//@}
			/**	@name	debuggers and diagnostics
			*/
			//@{

			/** Internal state and consistency self-validation.
					Calls {ConnectionObject::isValid}.
					\return			bool <tt> true</tt> if the internal state of this server is correct
					\see        ConnectionObject::isValid
			*/
			virtual bool isValid() const
				throw();

			/** Dump the current state of this server to 
					the output ostream with a dumping depth.
					\param   s output stream where to output the state of this server
					\param   depth the dumping depth
					\see     ConnectionObject::dump
					\see     QTTimer::dump
			*/
			virtual void dump(std::ostream& s = std::cout, Size depth = 0) const
				throw();

			//@}	
			protected:

			/** @name Timer method.
			*/
			//@{
			
			/** Timer method.
					Virtually overridden method.
					This method handles the socket stream. Every second it checks whether
					a new object is available at the stream. If this is the case the stream
					will be accepted and the incoming object will be reveiced.
					At the moment only Composite objects will be accepted. If
					another object is received the exception NotCompositeObject
					will be thrown.				
					\see    QTTimer::timer
					\exception NotCompositeObject thrown if another object than Composite object is received
			*/
			virtual void timer();
						
			//@}

			private:

			/** private methodes used for reacting to client requests.
			*/
			void sendObject(IOStreamSocket& iostream_socket)
				throw(NotCompositeObject);

			void setCreatorValue(IOStreamSocket& iostream_socket);
			void getCreatorValue(IOStreamSocket& iostream_socket);
			void hasCreatorValue(IOStreamSocket& iostream_socket);


			/** private storage variables.
			*/
			ObjectCreator *object_creator_;

			Composite* received_composite_;

			typedef HashMap<unsigned long, Composite *> CompositeHashMap;

			CompositeHashMap composite_hashmap_;

			IOStreamSocket*	iostream_socket_;
			SockInetBuf*		sock_inet_buf_;
				
			// the port to bind to
			int							port_; 

			ServerPreferences  *server_preferences_;
			QLabel 					   *server_icon_;
			static const char  *mini_ray_xpm_[];
		};



#		ifndef BALL_NO_INLINE_FUNCTIONS
#			include <BALL/VIEW/KERNEL/server.iC>
#		endif
  
	}// namespace VIEW
}// namespace BALL

#endif // BALL_VIEW_KERNEL_SERVER_H
