// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: message.sip,v 1.9.4.1 2005/11/11 17:21:23 leonhardt Exp $
//

class Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End
  public:

	Message() throw();
	Message(const Message& /message/) throw();
//  virtual ~Message() throw();
	// void setSender(const void* object);
	// const void* getSender() const;
	// void setDeletable(bool flag = true);
	// bool isDeletable() const;
};


class CompositeMessage
 : Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:

	enum CompositeMessageType
	{
		UNDEFINED = -1,
		NEW_COMPOSITE,
		REMOVED_COMPOSITE,
		CHANGED_COMPOSITE,
		CHANGED_COMPOSITE_HIERARCHY,
		SELECTED_COMPOSITE,
		DESELECTED_COMPOSITE,
		CENTER_CAMERA,
		NEW_MOLECULE
	};

	CompositeMessage() throw();
	CompositeMessage(const Composite& , CompositeMessageType , bool = true) throw();
	CompositeMessage(const CompositeMessage& /message/) throw();
	//virtual ~CompositeMessage() throw();
	void setComposite(const Composite& /composite/) throw();
	Composite* getComposite() const throw();
	void setCompositeName(const String& /name/) throw();
	const String& getCompositeName() const throw();
	void setType(CompositeMessageType /type/) throw();
	CompositeMessageType getType() const throw();
	void setUpdateRepresentations(bool) throw(); 
	bool updateRepresentations() const throw();
	void setShowSelectionInfos(bool) throw();
	bool showSelectionInfos() throw();
};


class SceneMessage
	: Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:

	enum SceneMessageType
	{
		UNDEFINED = 0,
		REBUILD_DISPLAY_LISTS,
		REDRAW,
		UPDATE_CAMERA,
		REMOVE_COORDINATE_SYSTEM,
		EXPORT_PNG,
		EXPORT_POVRAY,
		ENTER_ROTATE_MODE,
		ENTER_PICKING_MODE,
		ENTER_MOVE_MODE
	};
	SceneMessage(SceneMessageType = UNDEFINED) throw();
	SceneMessage(const SceneMessage&) throw();
	//virtual ~SceneMessage() throw();
	
	void setType(SceneMessageType) throw();
	SceneMessageType getType() const throw();
	void setStage(Stage) throw();
	Stage& getStage() throw();
 	const Stage& getStage() const throw();
};


class GenericSelectionMessage
	: Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:

	GenericSelectionMessage() throw();
	GenericSelectionMessage(const GenericSelectionMessage& ) throw();
	//virtual ~GenericSelectionMessage() throw();
	void setSelection(const PyCompositeList&) throw();
	PyCompositeList getSelection() const throw();
	
};


class ControlSelectionMessage
	: GenericSelectionMessage
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:
	ControlSelectionMessage() throw();
};


class NewSelectionMessage
	: Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:
	NewSelectionMessage() throw();
	bool openItems() throw();
	void setOpenItems(bool) throw();

};


class GeometricObjectSelectionMessage
	: Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:
	GeometricObjectSelectionMessage() throw();
	//virtual ~GeometricObjectSelectionMessage() throw();
	void setSelection(const GeometricObjectList& ) throw();
	const GeometricObjectList& getSelection() const throw();
	void setSelected(bool ) throw();
	bool isSelected() const throw();
};


class RepresentationMessage
	: Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:

	enum RepresentationMessageType
	{
		UNDEFINED = -1,
		ADD = 0,
		REMOVE,
		UPDATE,
		SELECTED,
		STARTED_UPDATE,
		FINISHED_UPDATE,
		ADD_TO_GEOMETRIC_CONTROL
	};

	RepresentationMessage() throw();
	//virtual ~RepresentationMessage() throw();
	RepresentationMessage(Representation& , RepresentationMessageType ) throw();
	void setRepresentation(Representation& ) throw();
	Representation* getRepresentation() throw();
	void setType(RepresentationMessageType ) throw();
	RepresentationMessageType getType() const throw();
};


class MolecularTaskMessage
	: Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:

	enum MolecularTaskMessageType
	{
		UNDEFINED = -1,
		BUILD_BONDS,
		CHECK_RESIDUE,
		ADD_HYDROGENS,
		CREATE_DISTANCE_GRID
	};

	MolecularTaskMessage(MolecularTaskMessageType = UNDEFINED) throw();
	void setType(MolecularTaskMessageType) throw();
	MolecularTaskMessageType getType() const throw();
};


class ShowDisplayPropertiesMessage
	:Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:
		ShowDisplayPropertiesMessage();
};


class CreateRepresentationMessage
	:Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:
		CreateRepresentationMessage() throw();
		CreateRepresentationMessage(const PyCompositeList&, 
																ModelType , 
																ColoringMethod ) throw();   
		//const PyCompositeList& getComposites() const throw();
		ModelType getModelType() const throw();
		ColoringMethod getColoringMethod() const throw();
};
	

class NewTrajectoryMessage
	:CompositeMessage
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:
		NewTrajectoryMessage() throw();
		void setTrajectoryFile(TrajectoryFile&) throw();
		TrajectoryFile* getTrajectoryFile() throw();
};
	
class DeselectControlsMessage
	: Message
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:
		DeselectControlsMessage();
};


class RegularDataMessage
	: CompositeMessage
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:

		enum RegularDataMessageType
		{
			UNDEFINED = -1,
			NEW = 100,
			REMOVE,
			UPDATE,
			SELECTED
		};

		RegularDataMessage() throw();
};

class RegularData1DMessage
	: RegularDataMessage
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:

	RegularData1DMessage(RegularDataMessageType /*type = UNDEFINED*/) throw();
	void setData(RegularData1D& /*data*/) throw();
	RegularData1D* getData() throw();
};

class RegularData2DMessage
	: RegularDataMessage
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:

	RegularData2DMessage(RegularDataMessageType /*type = UNDEFINED*/) throw();
	void setData(RegularData2D& /*data*/) throw();
	RegularData2D* getData() throw();
};

class RegularData3DMessage
	: RegularDataMessage
{
%TypeHeaderCode
	#include <BALL/VIEW/KERNEL/message.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End

	public:

	RegularData3DMessage(RegularDataMessageType /*type = UNDEFINED*/) throw();
	void setData(RegularData3D& /*data*/) throw();
	RegularData3D* getData() throw();
};



