// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: TCPTransfer.h,v 1.20.4.4 2005/07/29 12:38:14 amoll Exp $
//

#ifndef BALL_SYSTEM_TCPTRANSFER
#define BALL_SYSTEM_TCPTRANSFER

#ifndef BALL_COMMON_H
# include <BALL/common.h>
#endif

#ifndef BALL_DATATYPE_STRING_H
# include <BALL/DATATYPE/string.h>
#endif

#include <fstream>

#define BUFFER_SIZE 1024

namespace BALL
{
	/** Class for TCPTransfers.
	 * 	You can use http and ftp protocol to transfer files
	 * 	from the internet and store them in an ostream.
	 * 	Support for login, password and port is build in for
	 * 	FTP and HTTP. FTP transports are done by passive FTP, so
	 * 	no problems with local firewalls are encountered.
	 
	 \ingroup System
	 */
	class BALL_EXPORT TCPTransfer
	{
		public:

			/// typdef for Socket
			typedef int Socket;

			/** Enumeration of some status codes.
					In addition to these TCPTransfer can also return the
					status codes of HTTP and FTP.
			*/
			enum Status
			{
				OK  										= 0,
				GETHOSTBYNAME__ERROR 		,
				SOCKET__ERROR 					,
				CONNECT__ERROR 					,
				RECV__ERROR 						,
				OUTOFMEMORY__ERROR 			,
				BODY__ERROR 						,
				UNKNOWN__ERROR 					,
				ADDRESS__ERROR 					,
				UNINITIALIZED__ERROR 		,
				TRANSFER__ERROR					,
				SEND__ERROR							,
				PORT__ERROR							,
				UNKNOWN_PROTOCOL__ERROR ,
				LOGON__ERROR						,
				PROXY__ERROR 						,	
				FILENOTFOUND__ERROR 		= 404
			};

			/** Enumeration of the supported protocols.
			 */
			enum Protocol
			{
				///
				UNKNOWN_PROTOCOL 	= 0,
				
				///
				HTTP_PROTOCOL			= 1,
				
				///
				FTP_PROTOCOL			= 2
			};
			
			/** TransferFailed exception.
			 		This exception is thrown if a transfer fails.
			\ingroup System
			*/
			class BALL_EXPORT TransferFailed 
				: public Exception::GeneralException
			{
				 public:

					TransferFailed(const char* file, int line, Index error_code)
						throw();
				};
			
			/** Default constructor.
					The instance is set to UNINITIALIZED__ERROR.
			*/
			TCPTransfer()
				throw();
				
			/** Detailled constructor.
			 *	The file from the given address is downloaded and
			 *	stored in the ostream instance.
			 *	You can specify protocol, port and a login/pass by using
			 *	the syntax:  \par
			 *	http[ftp]://login:pass@server.com:port/fileaddress
			 */
			TCPTransfer(std::ostream& file, const String& address)
				throw(TransferFailed);
			
			/// Destructor.
			virtual ~TCPTransfer()
				throw();
		
			/** Clear method.
			 *  All attributes are set to default values and the status is set
			 *  to UNINITIALIZED__ERROR.
			 */
			void clear()
				throw();
		
			/** Set method.
			 * 	You can set a new file and address, but the transfer is not
			 * 	yet done. To do that, use transfer() afterwards.
			 */
			bool set(::std::ostream& file, const String& address)
				throw();
			
			
			/** Detailled set method.
			 *  @see set(ostream& file, const String& address)
			 */
			void set(::std::ostream& file, 
							 Protocol protocol, 
							 const String& host_address, 
							 const String& file_address,
							 const String& login,
							 const String& password,
							 Position 		 port = 80)
				throw();

			/** Return the host address.
			*/
			const String& getHostAddress() const
				throw()
			{
				return host_address_;
			}

			/// Return the file address.
			const String& getFileAddress() const
				throw()
			{
				return file_address_;
			}

			/// Return the port number.
			Position getPort() const
				throw()
			{
				return port_;
			}

			/** Return the status code.
			 * 	@see Status
			 */
			Status getStatusCode() const
				throw()
			{
				return status_;
			}
				
			/// Return the amount of received bytes.
			Position getReceivedBytes() const
				throw()
			{
				return received_bytes_;
			}
			
			/** Return the protocol.
			 * 	@see Protocol
			 */
			Protocol getProtocol() const
				throw()
			{
				return protocol_;
			}
			
			/// Return the login.
			const String& getLogin() const
				throw()
			{
				return login_;
			}

			/// Return the password.
			const String& getPassword() const
				throw()
			{
				return password_;
			}

			/// Get a pointer to the stream.
			const ::std::ostream* getStream() const
				throw()
			{
				return fstream_;
			}
					
			/** Return a pointer to the buffer.
			 		Remember that the buffer may contain NULL-chars.
					The length of the buffer is less than BUFFER_LENGTH.
					This method should only be used for debuging.
			*/
			const char* getBuffer() const
				throw()
			{
				return &(buffer_[0]);
			}

			/** Transfer method.
					If the address is specified, the file can be transfered
					with this method.
					@return Status the status of the transfer/instance
					@see Status
			*/
			Status transfer()
				throw();

			///
			void setProxy(const String proxy_address, Position port);

			///
			bool usingProxy() const;

			/// abort a running transfer
			void abort() { abort_ = true;}
			
			/** Dump the content of the buffer to an ostream.
					@param	s the stream to which we will dump
					@param	depth the indentation depth of the output
			*/
			void dump(std::ostream& s = std::cout, Size depth = 0) const
				throw();

			protected:
				
				String 			host_address_;
				String 			file_address_;
				Position 		port_;
				String			login_;
				String 			password_;
				Status			status_;
				int			 		received_bytes_;
				Protocol 		protocol_;
				char				buffer_[BUFFER_SIZE + 1];
				Socket			socket_;
				std::ostream*  fstream_;
				String 			proxy_address_;
				Position 		proxy_port_;
				bool 				abort_;
				
				/*_ Send data through the socket.
				 */
				Status sendData_(const String& query, Socket socket)
					throw();
				
				/*_ Logon to a server.
				 *  @param query string to send to the server as first contact
				 */
				Status	logon_(const String& query)
					throw();
				
				//_ Specified method for transfering per FTP-protocol
				Status	getFTP_()
					throw();
				
				//_ Specified method for transfering per HTTP-protocol
				Status	getHTTP_()
					throw();

				/*_ Read a complete status message form a FTP server
				    Return false if timeout of 20 seconds is exceeded or an other than the given message
						is received.
				*/						
				bool getFTPMessage_(Index status)
					throw();				
				
				//_ Compute the status of a ftp server from its response
				Status	getFTPStatus_()
					throw();

				//_ Compute the status of a http server from its response
				Status	getHTTPStatus_()
					throw();

				//_ Set a socket to blocking or nonblocking mode.
				Status	setBlock_(Socket socket, bool block = true)
					throw();

				//_ Wait a given time for output from the Socket.
				bool 		waitForOutput_(const String& key, Size seconds)
					throw();
				
				//_
				int getReceivedBytes_(Socket& socket);

			private:
				
				void operator = (TCPTransfer)
					throw();
	};
  
}

#endif // BALL_SYSTEM_TCPTRANSFER_H 
