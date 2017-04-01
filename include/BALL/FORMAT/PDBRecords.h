// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: PDBRecords.h,v 1.1.4.1 2005/07/28 13:52:55 amoll Exp $
//

#ifndef BALL_FORMAT_PDBRECORDS_H
#define BALL_FORMAT_PDBRECORDS_H

#include <vector>
#include <stdexcept>

#ifndef BALL_DATATYPE_HASHMAP_H
#	include <BALL/DATATYPE/hashMap.h>
#endif

namespace BALL 
{
	
	/**	PDB record class.
			This class contains PDB records in an undigested
			format and provides some means of digesting this format.
			It is meant to capture all records not parsed by
			\link GenericPDBFile \endlink and \link PDBFile \endlink.
			\par
			The class fulfills the requirements of an STL container and
			behaves mostly like the vector of Strings it actually is.
	*/
	class BALL_EXPORT PDBRecords
	{
		public:
		
		/* STL compatibility typedefs */
		typedef std::vector<String>::iterator Iterator;
		typedef std::vector<String>::iterator iterator;
		typedef std::vector<String>::const_iterator ConstIterator;
		typedef std::vector<String>::const_iterator const_iterator;
		typedef std::vector<String>::reverse_iterator ReverseIterator;
		typedef std::vector<String>::reverse_iterator reverse_iterator;
		typedef std::vector<String>::const_reverse_iterator ConstReverseIterator;
		typedef std::vector<String>::const_reverse_iterator const_reverse_iterator;
		typedef String value_type;
		typedef String ValueType;
		typedef String& reference;
		typedef String& Reference;
		typedef const String& const_reference;
		typedef const String& ConstReference;
		typedef std::vector<String>::difference_type difference_type;
		typedef std::vector<String>::difference_type DifferenceType;
		typedef std::vector<String>::size_type size_type;
		typedef std::vector<String>::size_type SizeType;
		
		/**	@name Constructors and Destructor */
		//@{
		///
		PDBRecords() throw();
		///
		PDBRecords(const PDBRecords& pdbi) throw();
		///
		virtual ~PDBRecords() throw();
		//@}

		/**	Assignment
		*/
		//@{
		///
		PDBRecords& operator = (const PDBRecords& rhs);
		//@}

		/** STL container compatibility interface.
				These methods just wrap the corresponding methods
				of std::vector. Refer to STL documentation for details.
		*/
		//@{
		///
		ConstIterator begin() const { return records_.begin(); }
		///
		Iterator begin() { return records_.begin(); }
		///
		ConstIterator end() const { return records_.end(); }
		///
		Iterator end() { return records_.end(); }
		///
		ConstReverseIterator rbegin() const { return records_.rbegin(); }
		///
		ReverseIterator rbegin() { return records_.rbegin(); }
		///
		ConstReverseIterator rend() const { return records_.rend(); }
		///
		ReverseIterator rend() { return records_.rend(); }
		///
		SizeType size() const { return records_.size(); }
		///
		SizeType max_size() const { return records_.max_size(); }
		///
		SizeType capacity() const { return records_.capacity(); }
		///
		bool empty() const { return records_.empty(); }
		///
		void clear() { records_.clear(); }
		///
		void resize(SizeType sz, ValueType c = ValueType()) { records_.resize(sz, c); }
		///
		Reference front() { return records_.front(); }
		///
		ConstReference front() const { return records_.front(); }
		///
		Reference back() { return records_.back(); }
		///
		ConstReference back() const { return records_.back(); }
		///
		void push_back(ConstReference x) { records_.push_back(x); }
		///
		void insert(Iterator pos, ConstReference value) { records_.insert(pos, value); }
		///
		void insert(Iterator pos, SizeType n, ConstReference value) { records_.insert(pos, n, value); }
		///
		void pop_back() { records_.pop_back(); }
		///
		Iterator erase(Iterator pos) { return records_.erase(pos); }
		///
		Iterator erase(Iterator first, Iterator last) { return records_.erase(first, last); }
		///
		bool operator == (const PDBRecords& rhs) const { return records_ == rhs.records_; }
		///
		bool operator != (const PDBRecords& rhs) const { return records_ != rhs.records_; }
		///
		bool operator < (const PDBRecords& rhs) const { return records_ < rhs.records_; }
		///
		bool operator > (const PDBRecords& rhs) const { return records_ > rhs.records_; }
		///
		bool operator <= (const PDBRecords& rhs) const { return records_ <= rhs.records_; }
		///
		bool operator >= (const PDBRecords& rhs) const { return records_ >= rhs.records_; }
		///
		void swap(PDBRecords& rhs) { records_.swap(rhs.records_); }
		///
		ConstReference operator [] (SizeType n) const { return records_[n]; }
		///
		Reference operator [] (SizeType n) { return records_[n]; }
		///
		ConstReference at(SizeType n) const throw(std::out_of_range) { return records_.at(n); }
		///
		Reference at(SizeType n) throw(std::out_of_range) { return records_.at(n); }
		//@}

		protected:
		/// The PDB record buffer
		std::vector<String>	records_;
	};
	
} // namespace BALL

#endif // BALL_FORMAT_PDBRECORDS_H
