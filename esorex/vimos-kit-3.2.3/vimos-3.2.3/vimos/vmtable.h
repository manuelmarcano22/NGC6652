/* $Id: vmtable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2004 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-03-25 11:43:04 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VM_TABLE_H
#define VM_TABLE_H

#include <fitsio.h>

#include <pilmacros.h>

#include <vmtypes.h>


PIL_BEGIN_DECLS

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   defines of strings for tables

   Description:
   The following string are defined for use with VimosTables

   Values:
   VM_EMPTY_TABLE_STRING     string with name of empty table

   Updates:
   03 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
#define VM_EMPTY_TABLE_STRING  "Empty Table"

/* Length of the string containing the name of a VimosDescriptor */
#define VM_DESC_LENGTH  81


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   union VimosDescValue

   Description:
   VimosDescValue is the union that underlies the VimosDescriptors. Each
   VimosDescriptor has a pointer to a VimosDescValue that stores the value of
   the VimosDescriptor.

   Layout:
     int    i          use value as int
     VimosBool   b          use value as VimosBool
     float  f          use value as float
     double d          use value as double
     char   *s         use value as *char. Note that if a VimosDescriptor is
                       destructed using deleteDesctructor, and the type of the 
                       value stored in the VimosDescriptor is VM_STRING, the
                       memory pointed to by s is freed.
   
   Updates:
   03 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
typedef union _VIMOS_DESC_VAL_
{
  int   i;
  VimosBool   b;
  float  f;
  double d;
  char   *s;
  int *iar;
  float *far;
  double *dar;
} VimosDescValue;




/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   union VimosColumnValue

   Description:
   VimosColumnValue is the union that underlies the VimosColumn in the
   VimosTables. Each VimosColumn has a pointer to a VimosColumnValue that
   stores the array of the column. The destructor
   deleteColumnValue(VimosColumnValue *cVal) frees the memory used as the
   array in the ColumnValue.

   Layout:
     void *p               generic pointer
     int *iArray           use as array of ints
     float *fArray         use as array of floats
     double *dArray        use as array of doubles
     char   **sArray       use as array of Strings
   
   Updates:
   03 Nov 98: Created (TAO)
   04 Feb 99: String array added

-------------------------------------------------------------------------------- 
*/
typedef union _VIMOS_COLUMN_VALUE_
{
  void   *p;
  int    *iArray;
  float  *fArray;
  double *dArray;
  char   *cArray;
  char   **sArray;
} VimosColumnValue;



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   struct VimosDescriptor

   Description:
   Name-value-comment triplet to be used to store descriptors from a table 
   or an image. VimosDescriptor is built as a linked list.

   Layout: 
     VimosVarType   descType      Type of value stored in the descriptor. Only
                                  VM_INT, VM_BOOL, VM_FLOAT, VM_DOUBLE or
                                  VM_STRING are allowed.
     char           *descName     Name of the descriptor
     VimosDescValue *descValue    Pointer to the VimosDescValue containing the 
                                  value
     char           *descComment  Comment associated to the descriptor
     VimosDescriptor *prev        Link for linked list
     VimosDescriptor *next        Link for linked list
   
   Updates:
   03 Nov 98: Created (TAO)
   07 Jan 00: Added the comment field (MS)

-------------------------------------------------------------------------------- 
*/
typedef struct _VIMOS_DESCRIPTOR_
{
  VimosVarType    descType;
  char            *descName;
  int             len;
  VimosDescValue  *descValue;
  char            *descComment;
  struct _VIMOS_DESCRIPTOR_ *prev;
  struct _VIMOS_DESCRIPTOR_ *next;
} VimosDescriptor;




typedef struct _VIMOS_COLUMN_
{
  VimosVarType   colType;
  char           *colName;
  int            len;
  VimosColumnValue  *colValue;
  struct _VIMOS_COLUMN_ *prev;
  struct _VIMOS_COLUMN_ *next;
} VimosColumn;

typedef struct _VIMOS_TABLE_ 
{
  char            name[VM_DESC_LENGTH];
  VimosDescriptor *descs;
  int             numColumns;
  VimosColumn     *cols;
  fitsfile *fptr;
} VimosTable;

/* Constructor of VimosDescValue */
VimosDescValue *newDescValue();

/* Desctructor of VimosDescValue */
void deleteDescValue(VimosDescValue *dValue);


/* Constructor of VimosDescriptor */
VimosDescriptor *newDescriptor();

/* Desctructor of VimosDescriptor 

   Deletes a single descriptor (does NOT traverse the linked list)

   Updates:
   11 Nov 98: added handling of VM_STRING for descriptor value (TAO)
*/
void deleteDescriptor(VimosDescriptor *vDesc);

/* Desctructor of VimosDescriptor 

   Deletes a single descriptors in a linked list 
*/
void deleteAllDescriptors(VimosDescriptor *vDesc);

/*
 * The next one is an extension of deleteDescriptor (in vimosTable.c).
 * Tipically a descriptor belongs to a linked list, and just deleting
 * it would leave the list broken: the gap should be filled, and the
 * previous descriptor should point directly to the next (and the
 * next to the previous). The extracted descriptor is then deleted
 * using deleteDescriptor. To simplify the operation, a findDescriptor
 * is also called, and vDesc may just point to any point of the list
 * coming before the descriptor to be removed, that can therefore be
 * specified by name. Returns the number of descs removed.
 */

int removeDescriptor(VimosDescriptor **vDesc, const char *dscName);

/*   
 *   Descriptors containing the given substring are searched in the 
 *   list. When a match is found the descriptor is removed from the 
 *   list. The character "*" (with its usual meaning) must be used at 
 *   the beginning and/or at the end of the substring, if it is not
 *   used the substring is interpreted as an entire descriptor name.
 *   Returns the number of descriptors deleted.
 */

int deleteSetOfDescriptors (VimosDescriptor **desc, char *substring);
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosDescriptor *newStringDescriptor(const char *name, const char *value,
                         conts char *comment)

  Description:
  Allocates a new VimosDescriptor of type VM_STRING with name *name
  and sets the value to *value and the comment to *comment.

  Input:
  const char *name
  String with the name of the new Descriptor

  const char *value
  String with the value of the new Descriptor

  const char *comment
  String with the comment of the new Descriptor

  Return Value (success): 
  Pointer to the new Descriptor

  Return Value (error)
  NULL

  Updates:
  18 Nov 98: Created (TAO)
  07 Jan 00: Added comment field (MS)
  13 Jun 00: Return NULL on error instead of exiting (Maura)

--------------------------------------------------------------------------------
*/

VimosDescriptor *newStringDescriptor(const char *name, const char *value, 
				     const char *comment);
VimosDescriptor *newIntDescriptor(const char *name, int value, 
				  const char *comment);
VimosDescriptor *newFloatDescriptor(const char *name, float value, 
				    const char *comment);
VimosDescriptor *newBoolDescriptor(const char *name, VimosBool value, 
				   const char *comment);
VimosDescriptor *newDoubleDescriptor(const char *name, double value, 
				     const char *comment);

VimosDescriptor *newIntArrayDescriptor(const char *name, int *value, 
				       const char *comment, int len);
VimosDescriptor *newFloatArrayDescriptor(const char *name, float *value, 
					 const char *comment, int len);
VimosDescriptor *newDoubleArrayDescriptor(const char *name, double *value, 
                                          const char *comment, int len);


/* Constructor of VimosColumnValue */
VimosColumnValue *newColumnValue();

/* Desctructor of VimosColumnValue */
void deleteColumnValue(VimosColumnValue *vColVal);


/* Constructor of VimosColumn */
VimosColumn *newColumn();
VimosColumn *newIntColumn(int size, const char *name);
VimosColumn *newFloatColumn(int size, const char *name);
VimosColumn *newDoubleColumn(int size, const char *name);
VimosColumn *newCharacterColumn(int size, const char *name);
VimosColumn *newStringColumn(int size, const char *name);

/* Desctructor of VimosColumn */
void deleteColumn(VimosColumn *vCol);


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosTable *newTable()

  Description:
  Allocates a new VimosTable. The pointers in the table are set to NULL, the
  field numColumns is set to 0, and the name of the table is set to
  CM_EMPTY_TABLE_STRING. 

  Return Value (success): 
  Pointer to the new table

  Return Value (error)
  NULL

  Updates:
  04 Nov 98: Created (TAO)
  09 Nov 98: Added field colLength (TAO)
  13 Jun 00: Return NULL on error instead of exiting (Maura)
  21 Aug 00: Change field colLength to numColumns (CI)
--------------------------------------------------------------------------------
*/
VimosTable *newTable();



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void deleteAllColumns(VimosTable *table)


  Description:
  Removes all columns from a VimosTable

  Return Value:
  void

  Updates:
  09 Nov 98: Created (TAO)
--------------------------------------------------------------------------------
*/
void deleteAllColumns(VimosTable *table);



void deleteTable(VimosTable *table);

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosDescriptor *copyDescriptor(VimosDescriptor *desc)

  Description:
  Allocates a new VimosDescriptor and copies the value of the input descriptor 
  into the new descriptor. The links prev and next of the new copy are set to
  NULL. 

  Return Value (success):
  Pointer to newly created copy of input descriptor

  Return Value (error):
  Null

  Updates:
  10 Nov 98: Created (TAO)
  13 Jun 00: Return NULL on error instead of exiting (Maura)

--------------------------------------------------------------------------------
*/
VimosDescriptor *copyOfDescriptor(VimosDescriptor *desc);


VimosBool copyTableDescriptors(VimosTable *inTable, VimosTable *outTable);

VimosBool copyAllDescriptors(VimosDescriptor *inDesc, 
                             VimosDescriptor **outDesc);

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosBool addDesc2Table(VimosTable *table, VimosDescriptor *desc)

  Description:
  Add a descriptor to the linked list of descriptors of the Table

  Input:
  VimosTable *table
  Table to add descriptor to

  VimosDescriptor *desc
  descriptor to be added

  Return Value (succes):
  VM_TRUE

  Return Value (error):
  VM_FALSE

  Updates: 
  18 Nov 98: Created (TAO)
  13 Jun 00: Return VimosBool instead of void (Maura)
--------------------------------------------------------------------------------
*/
VimosBool addDesc2Table(VimosDescriptor *desc, VimosTable *table);

VimosBool addDesc2Desc(VimosDescriptor *inDesc, VimosDescriptor **outDesc);


/* Get a pointer to the descriptor with name name */
VimosDescriptor *findDescriptor(VimosDescriptor *vDesc, const char *name);



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  int getDescriptorLength(VimosDescriptor *desc, const char *name)
  
  Description:
  get the length of the data of the Descriptor

  Input:
  VimosDescriptor *desc
  Descriptor to search for descriptor

  const char *name
  Name of descriptor to search for

  Return Value:
  void

  Updates: 
  29 Apr 99: Created (TAO)
--------------------------------------------------------------------------------
*/
int getDescriptorLength(VimosDescriptor *desc, const char *name);


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosBool readStringDescriptor(VimosDescriptor *desc, const char *name,
                            char *sVal, char *comment)

  Description:
  Searches *desc for a decriptor with name name and returns the value of that 
  descriptor, if it exists and is of the right type, and the associated 
  comment.

  Input:
  VimosDescriptor *desc
  Descriptor to be searched

  const char *name
  Name of the descriptor to search for

  int len (for array versions):
  maximum length of array to read

  Output:
  VimosBool   *bVal
  int    *iVal
  float  *fVal
  double *dVal
  char   *sVal
  Pointer to the variable where the value should be written in 
  For Array versions:
  pointer to existing array of length len
  
  char *comment

  Return Value (success):
  VM_TRUE

  Return Value (error):
  VM_FALSE

  Updates: 
  24 Nov 98: Created (TAO)
  07 Jan 00: Added comment field (MS)

--------------------------------------------------------------------------------
*/
VimosBool readBoolDescriptor(VimosDescriptor *desc, const char *name, 
			     VimosBool *bVal, char *comment);
VimosBool readIntDescriptor(VimosDescriptor *desc, const char *name, 
			    int *iVal, char *comment);
VimosBool readFloatDescriptor(VimosDescriptor *desc, const char *name, 
                         float *fVal, char *comment);
VimosBool readDoubleDescriptor(VimosDescriptor *desc, const char *name, 
                          double *dVal, char *comment);
VimosBool readStringDescriptor(VimosDescriptor *desc, const char *name, 
                          char *sVal, char *comment);
VimosBool readIntArrayDescriptor(VimosDescriptor *desc, const char *name, 
                                 int *iVal, char *comment, int len);
VimosBool readFloatArrayDescriptor(VimosDescriptor *desc, const char *name, 
                                   float *fVal, char *comment, int len);
VimosBool readDoubleArrayDescriptor(VimosDescriptor *desc, const char *name, 
                                    double *dVal, char *comment, int len);


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosBool readStringDescFromTable(VimosTable *table, const char *name
                               char *sVal, char *comment)

  Description:
  Searches *table for a decriptor with name name and returns the value of that 
  descriptor, if it exists and is of the right type.

  Input:
  VimosTable *table
  Table to be searched

  const char *name
  Name of the descriptor to search for

  Output:
  VimosBool   *bVal
  int    *iVal
  float  *fVal
  double *dVal
  char   *sVal
  Pointer to the variable where the value should be written in 

  char *comment

  Return Value (success):
  VM_TRUE

  Return Value (error):
  VM_FALSE

  Updates: 
  11 Nov 98: Created (TAO)
  24 Nov 98: Implemented in terms of read<type>Descriptor (TAO)
  07 Jan 00: Added comment field (MS)

--------------------------------------------------------------------------------
*/
VimosBool readBoolDescFromTable(VimosTable *table, const char *name, 
				VimosBool *bVal, char *comment);
VimosBool readIntDescFromTable(VimosTable *table, const char *name, 
			       int *iVal, char *comment);
VimosBool readFloatDescFromTable(VimosTable *table, const char *name, 
                                 float *fVal, char *comment);
VimosBool readDoubleDescFromTable(VimosTable *table, const char *name, 
                                  double *dVal, char *comment);
VimosBool readStringDescFromTable(VimosTable *table, const char *name, 
                                  char *sVal, char *comment);
VimosBool readIntArrayDescFromTable(VimosTable *table, const char *name, 
                                    int *iVal, char *comment, int len);
VimosBool readFloatArrayDescFromTable(VimosTable *table, const char *name, 
                                      float *fVal, char *comment, int len);
VimosBool readDoubleArrayDescFromTable(VimosTable *table, const char *name, 
                                       double *dVal, char *comment, int len);



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosBool writeIntDescriptor(VimosDescriptor **desc, const char *name
                               int value, const char *comment)
  VimosBool writeFloatDescriptor(VimosDescriptor **desc, const char *name
                                 float value, const char *comment)
  VimosBool writeDoubleDescriptor(VimosDescriptor **desc, const char *name
                                  double value, const char *comment)
  VimosBool writeStringDescriptor(VimosDescriptor **desc, const char *name
                                  const char *value, const char *comment)

  Description:
  Writes a descriptor with name name, value value, and comment comment
  to the list pointed to by *desc. If a descriptor of name name already 
  exists, the value of this one is overwritten, if it does not exist, or it 
  exists and has the wrong type, a new descriptor is created and appended 
  to the list. 

  Input:
  VimosDescriptor **desc
  Descriptor to be searched

  const char *name
  Name of the descriptor to search for

  Input:
  int    value
  float  value
  double value
  const char   *value
  Value of descriptor

  const char *comment
  comment field of the descriptor

  Return Value (succes):
  VM_TRUE

  Return Value (error):
  VM_FALSE

  Updates:
  24 Nov 98: Created (TAO)
  07 Jan 00: Added comment field (MS)
  13 Jun 00: Return VimosBool instead of void (Maura)

--------------------------------------------------------------------------------
*/

VimosBool writeIntDescriptor(VimosDescriptor **desc, const char *name, 
                             int value, const char *comment);

VimosBool writeFloatDescriptor(VimosDescriptor **desc, const char *name, 
			       float value, const char *comment);

VimosBool writeDoubleDescriptor(VimosDescriptor **desc, const char *name, 
			        double value, const char *comment);

VimosBool writeStringDescriptor(VimosDescriptor **desc, const char *name,
                                const char *value, const char *comment);

VimosBool copyFromHeaderToHeader(VimosDescriptor *fromDesc, const char *inName,
                     VimosDescriptor **toDesc, const char *outName);

VimosBool insertDescriptor(VimosDescriptor **desc, const char *name,
                           VimosDescriptor *newDesc, int before);

VimosBool insertIntDescriptor(VimosDescriptor **desc, const char *name,
     int value, const char *comment, const char *refName, int before);

VimosBool insertFloatDescriptor(VimosDescriptor **desc, const char *name,
     float value, const char *comment, const char *refName, int before);

VimosBool insertDoubleDescriptor(VimosDescriptor **desc, const char *name,
     double value, const char *comment, const char *refName, int before);

VimosBool insertStringDescriptor(VimosDescriptor **desc, const char *name,
     const char *value, const char *comment, const char *refName, int before);


/* Get a pointer to the Column with name name 
 */
VimosColumn *findColumn(VimosColumn *vCol, const char *name);

/* Get a pointer to the  descrpitor with name name in a VimosTable 
 */
VimosDescriptor *findDescInTab(VimosTable *table, const char *name);

/* Get a pointer to the column with name name in a VimosTable 
 */
VimosColumn *findColInTab(VimosTable *table, const char *name);

VimosBool insertHistoryDescriptor(VimosDescriptor **desc, const char *name,
				  const char *value, const char *comment);


VimosDescriptor *vimosDscFind(VimosDescriptor *, const char *);
int vimosDscCopy(VimosDescriptor **, VimosDescriptor *, const char *,
                 const char *);
int vimosDscErase(VimosDescriptor **, const char *);

int colGetSize(VimosColumn *);
int colSetName(VimosColumn *, const char *);

char **colGetStringData(VimosColumn *);
int *colGetIntData(VimosColumn *);
float *colGetFloatData(VimosColumn *);
double *colGetDoubleData(VimosColumn *);

int tblGetSize(VimosTable *table, const char *name);

char **tblGetStringData(VimosTable *, const char *);
int *tblGetIntData(VimosTable *, const char *);
float *tblGetFloatData(VimosTable *, const char *);
double *tblGetDoubleData(VimosTable *, const char *);

int tblSetStringValue(VimosTable *, const char *, unsigned int,
		      const char *);
int tblSetIntValue(VimosTable *, const char *, unsigned int, int);
int tblSetFloatValue(VimosTable *, const char *, unsigned int, float);
int tblSetDoubleValue(VimosTable *, const char *, unsigned int, double);

VimosColumn *tblCopyColumn(VimosTable *, const char *);
VimosColumn *tblRemoveColumn(VimosTable *, const char *);
int tblAppendColumn(VimosTable *, VimosColumn *);

/* FITS file should be open, and correct extension already selected,
   before using this function!!! */

VimosBool readDescsFromFitsTable(VimosDescriptor **desc, fitsfile *tableID);

/* FITS should exist and be open, and correct extension already
   selected, before using this function!!! */

VimosBool writeDescsToFitsTable(VimosDescriptor *desc, fitsfile *tableID);

/*
 * Create a disk table FITS file from structure VimosTable
 */

VimosBool createFitsTable(char *filename, VimosTable *, const char *tableType);

/* 
 * Open a Fits Table  
 */

VimosTable *openOldFitsTable(const char *tableName, int flag);

/*
 * read a Fits Table 
 */

VimosBool readFitsTable(VimosTable *, fitsfile *tableID);

VimosBool writeFitsTable(VimosTable *table, fitsfile *fptr);

/*
 * Close a Fits table
 */

VimosBool closeFitsTable (VimosTable *table, int flag);

int copyToPrimary(char *fileName, const char *keyName);

PIL_END_DECLS

#endif /* VM_TABLE_H */
