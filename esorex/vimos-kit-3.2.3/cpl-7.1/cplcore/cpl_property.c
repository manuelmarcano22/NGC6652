/*
 * This file is part of the ESO Common Pipeline Library
 * Copyright (C) 2001-2017 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string.h>
#include <complex.h>

#include <cxmemory.h>
#include <cxmessages.h>
#include <cxstrutils.h>

#include "cpl_error_impl.h"
#include "cpl_property.h"


/**
 * @defgroup cpl_property Properties
 *
 * This module implements the property type. The type @c cpl_property is
 * basically a variable container which consists of a name, a type identifier
 * and a specific value of that type. The type identifier always determines
 * the type of the associated value. A property is similar to an ordinary
 * variable and its current value can be set or retrieved through its name. In
 * addition a property may have a descriptive comment associated.
 * A property can be created for the basic types char, bool (int), int,
 * long, float and double. Also C strings are supported. Support for
 * arrays in general is currently not available.
 *
 * @par Synopsis:
 * @code
 *   #include <cpl_property.h>
 * @endcode
 */

/**@{*/

/*
 * The property value and property types
 */

typedef struct _cpl_property_value_ cpl_property_value;

struct _cpl_property_value_ {
    cpl_type type;
    cxsize size;
    cxptr data;
};


struct _cpl_property_ {
    cxchar *name;
    cxchar *comment;
    cpl_property_value *value;
};


/*
 * Private methods
 */

inline static cpl_property_value *
_cpl_property_value_new(cpl_type type, cxsize nelements, cxptr data)
{

    cpl_property_value *value;

    cxsize size = (cxsize)cpl_type_get_sizeof(type);


    if (size == 0 || nelements == 0)
        return NULL;

    if (!(type & CPL_TYPE_FLAG_ARRAY))
        nelements = 1;

    value = (cpl_property_value *) cx_malloc_clear(sizeof(cpl_property_value));

    value->type = type;
    value->size = nelements;

    value->data = cx_malloc_clear(nelements * size);

    if (data) {
        memcpy(value->data, data, nelements * size);
    } else {
        memset(value->data, 0, nelements * size);
    }

    return value;

}


inline static void
_cpl_property_value_delete(cpl_property_value *value)
{

    cx_assert(value != NULL);

    if (value->data != NULL)
        cx_free(value->data);

    cx_free(value);
    return;

}


inline static cpl_property_value *
_cpl_property_value_resize(cpl_property_value *value, cpl_type type,
                           cxsize nelements)
{

    cxsize sz = 0;

    cxptr tmp;


    cx_assert(value != NULL);
    cx_assert(value->type == type);
    cx_assert(nelements > 0);

    sz = (cxsize)cpl_type_get_sizeof(value->type);

    if (value->size != nelements) {


        sz *= nelements;
        tmp = cx_malloc_clear(sz);

        if (value->data != NULL)
            cx_free(value->data);

        value->data = tmp;
        value->size = nelements;
    }

    return value;

}


inline static int
_cpl_property_value_set(cpl_property_value *value, cpl_type type,
                        cxsize nelements, cxcptr data)
{

    cxsize sz = 0;


    cx_assert(value != NULL);
    cx_assert(nelements > 0);

    sz = (cxsize)cpl_type_get_sizeof(value->type);

    if (value->type != type) {
        return 1;
    }


    /*
     * Resize the value if necessary.
     */

    value = _cpl_property_value_resize(value, type, nelements);

    memcpy(value->data, data, sz * nelements);

    return 0;

}


inline static cxptr
_cpl_property_value_get(cpl_property_value *value, cpl_type type,
                        const cxchar *name)
{

    cx_assert(value != NULL);

    if (value->type != type) {

        cx_assert(name != NULL);

        cpl_error_set(name, CPL_ERROR_TYPE_MISMATCH);

        return NULL;
    }

    return value->data;

}


/*
 * Public methods
 */

/**
 * @brief
 *    Create an empty property of a given type.
 *
 * @param name     Property name.
 * @param type     Property type flag.
 *
 * @return
 *   The newly created property, or @c NULL if it could not be created. In the
 *   latter case an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>name</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_INVALID_TYPE</td>
 *       <td class="ecr">
 *         The requested property type <i>type</i> is not supported.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function allocates memory for a property of type @em type and assigns
 * the identifier string @em name to the newly created property.
 *
 * The returned property must be destroyed using the property destructor
 * @b cpl_property_delete().
 *
 * @see cpl_property_delete()
 */

cpl_property *
cpl_property_new(const char *name, cpl_type type)
{

    cxsize sz;

    cpl_property *property;


    if (name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }


    /* Type size is 0 for invalid types */

    sz = (cxsize)cpl_type_get_sizeof(type);

    if (sz == 0) {
        cpl_error_set_(CPL_ERROR_INVALID_TYPE);
        return NULL;
    }

    property = cx_malloc(sizeof(cpl_property));

    property->name = cx_strdup(name);
    property->value = _cpl_property_value_new(type, 1, NULL);
    property->comment = NULL;

    return property;

}


/**
 * @brief
 *    Create an empty property of a given type and size.
 *
 * @param name     Property name.
 * @param type     Property type flag.
 * @param size     Size of the property value.
 *
 * @return
 *   The newly created property, or @c NULL if it could not be created. in
 *   the latter case an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>name</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function allocates memory for a property of type @em type and assigns
 * the identifier string @em name to the newly created property. The
 * property value is created such that @em size elements of type @em type
 * can be stored.
 *
 * The returned property must be destroyed using the property destructor
 * @b cpl_property_delete().
 *
 * @see cpl_property_delete()
 */

cpl_property *
cpl_property_new_array(const char *name, cpl_type type, cpl_size size)
{

    cxsize sz;

    cpl_property *property;


    if (name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }


    /* Type size is 0 for invalid types */

    sz = (cxsize)cpl_type_get_sizeof(type);

    if (sz == 0) {
        return NULL;
    }

    property = cx_malloc(sizeof(cpl_property));

    property->name = cx_strdup(name);
    property->value = _cpl_property_value_new(type, (cxsize)size, NULL);
    property->comment = NULL;

    return property;

}


/**
 * @brief
 *    Create a copy of a property.
 *
 * @param other  The property to copy.
 *
 * @return
 *   The copy of the given property, or @c NULL in case of an error.
 *   In the latter case an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function returns a copy of the property @em self. The copy is a
 * deep copy, i.e. all property members are copied.
 */

cpl_property *
cpl_property_duplicate(const cpl_property *other)
{

    cpl_property *copy;


    if (other == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    copy = cx_malloc(sizeof(cpl_property));
    copy->name = cx_strdup(other->name);

    if (other->value) {
        copy->value = _cpl_property_value_new(other->value->type,
                                              other->value->size,
                                              other->value->data);
    }

    if (other->comment) {
        copy->comment = cx_strdup(other->comment);
    }
    else {
        copy->comment = NULL;
    }

    return copy;

}


/**
 * @brief
 *    Destroy a property.
 *
 * @param self  The property.
 *
 * @return Nothing.
 *
 * The function destroys a property of any kind. All property members
 * including their values are properly deallocated. If the property @em self
 * is @c NULL, nothing is done and no error is set.
 */

void
cpl_property_delete(cpl_property *self)
{

    if (self) {

        if (self->comment) {
            cx_free(self->comment);
        }

        if (self->value) {
            _cpl_property_value_delete(self->value);
        }

        if (self->name) {
            cx_free(self->name);
        }

        cx_free(self);

    }

    return;

}


/**
 * @brief
 *   Get the current number of elements a property contains.
 *
 * @param self  The property.
 *
 * @return
 *   The current number of elements or -1 in case of an error. If an
 *   error occurred an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function returns the current number of elements the property
 * @em self contains. This is the array size of the property's
 * value and in particular, for a property of the type
 * @c CPL_TYPE_STRING, it is the length of the string as
 * given by the @b strlen() function plus 1.
 */

cpl_size
cpl_property_get_size(const cpl_property *self)
{

    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return -1;
    }

    return (cpl_size)self->value->size;

}


/**
 * @brief
 *   Get the type of a property.
 *
 * @param self  The property.
 *
 * @return
 *   The type code of this property. In case of an error the returned
 *   type code is @c CPL_TYPE_INVALID and an appropriate error code is
 *   set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * This function returns the type of the data value stored in the
 * property @em self.
 */

cpl_type
cpl_property_get_type(const cpl_property *self)
{

    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return CPL_TYPE_INVALID;
    }

    return self->value->type;

}


/**
 * @brief
 *   Modify the name of a property.
 *
 * @param self  The property.
 * @param name  New name.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current name of @em self with a copy
 * of the string @em name. The function returns an error if @em name is
 * @c NULL.
 */

cpl_error_code
cpl_property_set_name(cpl_property *self, const char *name)
{

    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    if (!name) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    if (self->name) {
        cx_free(self->name);
    }

    self->name = cx_strdup(name);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Modify a property's comment.
 *
 * @param self     The property.
 * @param comment  New comment.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>comment</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current comment field of @em self with a
 * copy of the string @em comment. The new comment may be @c NULL. In this
 * case the function effectively deletes the current comment.
 */

cpl_error_code
cpl_property_set_comment(cpl_property *self, const char *comment)
{

    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    if (self->comment) {
        cx_free(self->comment);
    }

    if (comment) {
        self->comment = cx_strdup(comment);
    }
    else {
        self->comment = NULL;
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of a character property.
 *
 * @param self   The property.
 * @param value  New value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_CHAR</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current character value of @em self with a
 * copy of the value @em value.
 */

cpl_error_code
cpl_property_set_char(cpl_property *self, char value)
{

    cxint status = 0;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_assert(self->value != NULL);

    status = _cpl_property_value_set(self->value, CPL_TYPE_CHAR,
                                     1, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of a boolean property.
 *
 * @param self   The property.
 * @param value  New value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_BOOL</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current boolean value of @em self with a
 * a 0, if @em value is equal to 0. If @em value is different from 0 any
 * previous value of @em self is replaced by a 1.
 */

cpl_error_code
cpl_property_set_bool(cpl_property *self, int value)
{

    cxint status = 0;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_assert(self->value != NULL);

    value = value ? TRUE : FALSE;

    status = _cpl_property_value_set(self->value, CPL_TYPE_BOOL,
                                     1, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of an integer property.
 *
 * @param self   The property.
 * @param value  New value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_INT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current integer value of @em self with a
 * copy of the value @em value.
 */

cpl_error_code
cpl_property_set_int(cpl_property *self, int value)
{

    cxint status = 0;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_assert(self->value != NULL);

    status = _cpl_property_value_set(self->value, CPL_TYPE_INT,
                                     1, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of a long property.
 *
 * @param self   The property.
 * @param value  New value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_LONG</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current long value of @em self with a
 * copy of the value @em value.
 */

cpl_error_code
cpl_property_set_long(cpl_property *self, long value)
{

    cxint status = 0;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_assert(self->value != NULL);

    status = _cpl_property_value_set(self->value, CPL_TYPE_LONG,
                                     1, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of a long long property.
 *
 * @param self   The property.
 * @param value  New value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_LONG</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current long long value of @em self with a
 * copy of the value @em value.
 */

cpl_error_code
cpl_property_set_long_long(cpl_property *self, long long value)
{

    cxint status = 0;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_assert(self->value != NULL);

    status = _cpl_property_value_set(self->value, CPL_TYPE_LONG_LONG,
                                     1, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of a float property.
 *
 * @param self   The property.
 * @param value  New value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_FLOAT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current float value of @em self with a
 * copy of the value @em value.
 */

cpl_error_code
cpl_property_set_float(cpl_property *self, float value)
{

    cxint status = 0;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_assert(self->value != NULL);

    status = _cpl_property_value_set(self->value, CPL_TYPE_FLOAT,
                                     1, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of a double property.
 *
 * @param self   The property.
 * @param value  New value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_DOUBLE</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current double value of @em self with a
 * copy of the value @em value.
 */

cpl_error_code
cpl_property_set_double(cpl_property *self, double value)
{


    cxint status = 0;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_assert(self->value != NULL);

    status = _cpl_property_value_set(self->value, CPL_TYPE_DOUBLE,
                                     1, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of a string property.
 *
 * @param self   The property.
 * @param value  New value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>value</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_STRING</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current string value of @em self with a
 * copy of the value @em value.
 */

cpl_error_code
cpl_property_set_string(cpl_property *self, const char *value)
{


    cxint status = 0;

    cxsize len;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    if (value == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_assert(self->value != NULL);


    /*
     * Determine the length of the string value and add 1 for the
     * trailing 0.
     */

    len = (cxsize) strlen(value) + 1;

    status = _cpl_property_value_set(self->value, CPL_TYPE_STRING,
                                     len, (cxcptr)value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of a complex float property.
 *
 * @param self   The property.
 * @param value  New value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type
 *            <tt>CPL_TYPE_COMPLEX_FLOAT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current complex float value of @em self with a
 * copy of the value @em value.
 */

cpl_error_code
cpl_property_set_float_complex(cpl_property *self, float complex value)
{

    cxint status = 0;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_assert(self->value != NULL);

    status = _cpl_property_value_set(self->value, CPL_TYPE_FLOAT_COMPLEX,
                                     1, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of a double complex property.
 *
 * @param self   The property.
 * @param value  New value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success and an appropriate
 *   error code if an error occurred. In the latter case the error code is
 *   also set in the error handling system.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type
 *            <tt>CPL_TYPE_DOUBLE_COMPLEX</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function replaces the current double value of @em self with a
 * copy of the value @em value.
 */

cpl_error_code
cpl_property_set_double_complex(cpl_property *self, double complex value)
{

    cxint status = 0;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_assert(self->value != NULL);

    status = _cpl_property_value_set(self->value, CPL_TYPE_DOUBLE_COMPLEX,
                                     1, (cxcptr)&value);

    if (status)
    {
        return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Get the property name.
 *
 * @param self  The property.
 *
 * @return
 *   The name of the property or @c NULL if an error occurred. In the
 *   latter case an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function returns a handle for the read-only identifier string of @em self.
 */

const char *
cpl_property_get_name(const cpl_property *self)
{

    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    return self->name;

}


/**
 * @brief
 *   Get the property comment.
 *
 * @param self  The property.
 *
 * @return
 *   The comment of the property if it is present or @c NULL. If an error
 *   occurrs the function returns @c NULL and sets an appropriate error code.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function returns a handle for the read-only comment string of @em self.
 */

const char *
cpl_property_get_comment(const cpl_property *self)
{

    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (!self->comment) {
        return NULL;
    }

    return self->comment;

}


/**
 * @brief
 *   Get the value of a character property.
 *
 * @param self  The property.
 *
 * @return
 *   The current property value. If an error occurs the function returns
 *   '\\0' and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_CHAR</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function retrieves the character value currently stored in the
 * property @em self.
 */

char
cpl_property_get_char(const cpl_property *self)
{

    cxptr result;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return '\0';
    }

    cx_assert(self->value != NULL);

    result = _cpl_property_value_get(self->value, CPL_TYPE_CHAR, cpl_func);

    if (!result) {
        return '\0';
    }

    return *((char *)result);

}


/**
 * @brief
 *   Get the value of a boolean property.
 *
 * @param self  The property.
 *
 * @return
 *   The current property value. If an error occurs the function returns
 *   0 and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_BOOL</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function retrieves the boolean value currently stored in the
 * property @em self.
 */

int
cpl_property_get_bool(const cpl_property *self)
{



    cxptr result;
    cxbool val;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0;
    }

    cx_assert(self->value != NULL);

    result = _cpl_property_value_get(self->value, CPL_TYPE_BOOL, cpl_func);

    if (!result) {
        return 0;
    }

    val = *((cxbool *)result);

    return val == TRUE ? 1 : 0;

}


/**
 * @brief
 *   Get the value of an integer property.
 *
 * @param self  The property.
 *
 * @return
 *   The current property value. If an error occurs the function returns
 *   0 and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_INT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function retrieves the integer value currently stored in the
 * property @em self.
 */

int
cpl_property_get_int(const cpl_property *self)
{


    cxptr result;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0;
    }

    cx_assert(self->value != NULL);

    result = _cpl_property_value_get(self->value, CPL_TYPE_INT, cpl_func);

    if (!result) {
        return 0;
    }

    return *((int *)result);

}


/**
 * @brief
 *   Get the value of a long property.
 *
 * @param self  The property.
 *
 * @return
 *   The current property value. If an error occurs the function returns
 *   0 and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_LONG</tt>
 *         or <tt>CPL_TYPE_INT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function retrieves the value currently stored in the property @em self
 * as a long integer. The function accepts properties of integer type
 * with a rank less or equal than the function's return type.
 */

long
cpl_property_get_long(const cpl_property *self)
{

    cxlong value = 0;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return value;
    }

    cx_assert(self->value != NULL);


    /*
     * If the property self is of the same kind, but of a smaller rank,
     * promote the property's value to the target type.
     */

    switch (self->value->type) {

        case CPL_TYPE_INT:
        {
            cxint *result = _cpl_property_value_get(self->value,
                                                    CPL_TYPE_INT, cpl_func);

            value = *result;
            break;
        }

        case CPL_TYPE_LONG:
        {
            cxlong *result = _cpl_property_value_get(self->value,
                                                     CPL_TYPE_LONG, cpl_func);

            value = *result;
            break;
        }

        default:
        {
            cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
            break;
        }

    }

    return value;

}


/**
 * @brief
 *   Get the value of a long long property.
 *
 * @param self  The property.
 *
 * @return
 *   The current property value. If an error occurs the function returns
 *   0 and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_LONG_LONG</tt>,
 *         <tt>CPL_TYPE_LONG</tt>, or <tt>CPL_TYPE_INT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function retrieves the value currently stored in the property @em self
 * as a long long integer. The function accepts properties of integer type
 * with a rank less or equal than the function's return type.
 */

long long
cpl_property_get_long_long(const cpl_property *self)
{


    long long value = 0;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return value;
    }

    cx_assert(self->value != NULL);


    /*
     * If the property self is of the same kind, but of a smaller rank,
     * promote the property's value to the target type.
     */

    switch (self->value->type) {

        case CPL_TYPE_INT:
        {
            cxint *result = _cpl_property_value_get(self->value,
                                                    CPL_TYPE_INT,
                                                    cpl_func);

            value = *result;
            break;
        }

        case CPL_TYPE_LONG:
        {
            cxlong *result = _cpl_property_value_get(self->value,
                                                     CPL_TYPE_LONG,
                                                     cpl_func);

            value = *result;
            break;
        }

        case CPL_TYPE_LONG_LONG:
        {
            long long *result = _cpl_property_value_get(self->value,
                                                        CPL_TYPE_LONG_LONG,
                                                        cpl_func);

            value = *result;
            break;
        }

        default:
        {
            cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
            break;
        }

    }

    return value;

}


/**
 * @brief
 *   Get the value of a float property.
 *
 * @param self  The property.
 *
 * @return
 *   The current property value. If an error occurs the function returns
 *   0. and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_FLOAT</tt>
 *         or <tt>CPL_TYPE_DOUBLE</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function retrieves the value currently stored in the property @em self
 * as a float. The function also accepts properties of type double and returns
 * the property value as float.
 *
 * @note
 *  Calling the function for a property of type double may lead to
 *  truncation errors!
 */

float
cpl_property_get_float(const cpl_property *self)
{

    cxfloat value = 0.;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return value;
    }

    cx_assert(self->value != NULL);


    /*
     * If the property self is of type double cast it to float.
     * This may lead to truncation errors!
     */

    switch (self->value->type) {

        case CPL_TYPE_FLOAT:
        {
            cxfloat *result = _cpl_property_value_get(self->value,
                                                      CPL_TYPE_FLOAT,
                                                      cpl_func);

            value = *result;
            break;
        }

        case CPL_TYPE_DOUBLE:
        {
            cxdouble *result = _cpl_property_value_get(self->value,
                                                       CPL_TYPE_DOUBLE,
                                                       cpl_func);

            value = (cxfloat)(*result);
            break;
        }

        default:
        {
            cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
            break;
        }

    }

    return value;

}


/**
 * @brief
 *   Get the value of a double property.
 *
 * @param self  The property.
 *
 * @return
 *   The current property value. If an error occurs the function returns
 *   0. and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_DOUBLE</tt>
 *         or <tt>CPL_TYPE_FLOAT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function retrieves the value currently stored in the property @em self
 * as a double. The function accepts properties of a floating-point type
 * with a rank less or equal than the function's return type.
 */

double
cpl_property_get_double(const cpl_property *self)
{

    cxdouble value = 0.;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return value;
    }

    cx_assert(self->value != NULL);


    /*
     * If the property self is of the same kind, but of a smaller rank,
     * promote the property's value to the target type.
     */

    switch (self->value->type) {

        case CPL_TYPE_FLOAT:
        {
            cxfloat *result = _cpl_property_value_get(self->value,
                                                      CPL_TYPE_FLOAT, cpl_func);

            value = *result;
            break;
        }

        case CPL_TYPE_DOUBLE:
        {
            cxdouble *result = _cpl_property_value_get(self->value,
                                                       CPL_TYPE_DOUBLE,
                                                       cpl_func);

            value = *result;
            break;
        }

        default:
        {
            cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
            break;
        }

    }

    return value;

}


/**
 * @brief
 *   Get the value of a string property.
 *
 * @param self  The property.
 *
 * @return
 *   The current property value. If an error occurs the function returns
 *   @c NULL and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type <tt>CPL_TYPE_STRING</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function retrieves the string value currently stored in the
 * property @em self.
 */

const char *
cpl_property_get_string(const cpl_property *self)
{

    cxptr result;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    cx_assert(self->value != NULL);

    result = _cpl_property_value_get(self->value, CPL_TYPE_STRING, cpl_func);

    if (!result) {
        return NULL;
    }

    return (const char *)result;

}

/**
 * @brief
 *   Get the value of a float complex property.
 *
 * @param self  The property.
 *
 * @return
 *   The current property value. If an error occurs the function returns
 *   0. and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type
 *            <tt>CPL_TYPE_FLOAT_COMPLEX</tt> or
 *            <tt>CPL_TYPE_DOUBLE_COMPLEX</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function retrieves the value currently stored in the property @em self
 * as a float complex. The function also accepts properties of type
 * double complex and returns the property value as float complex.
 *
 * @note
 *  Calling the function for a property of type double complex may lead to
 *  truncation errors!
 */

float complex
cpl_property_get_float_complex(const cpl_property *self)
{

    float complex value = 0.;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return value;
    }

    cx_assert(self->value != NULL);


    /*
     * If the property self is of type double cast it to float.
     * This may lead to truncation errors!
     */

    switch (self->value->type) {

        case CPL_TYPE_FLOAT_COMPLEX:
        {
            float complex *result =
                    _cpl_property_value_get(self->value,
                                            CPL_TYPE_FLOAT_COMPLEX,
                                            cpl_func);

            value = *result;
            break;
        }

        case CPL_TYPE_DOUBLE_COMPLEX:
        {
            double complex *result =
                    _cpl_property_value_get(self->value,
                                            CPL_TYPE_DOUBLE_COMPLEX,
                                            cpl_func);

            value = (float complex)(*result);
            break;
        }

        default:
        {
            cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
            break;
        }

    }

    return value;

}


/**
 * @brief
 *   Get the value of a double complex property.
 *
 * @param self  The property.
 *
 * @return
 *   The current property value. If an error occurs the function returns
 *   0. and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property <i>self</i> is not of type
 *            <tt>CPL_TYPE_DOUBLE_COMPLEX</tt> or
 *            <tt>CPL_TYPE_FLOAT_COMPLEX</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function retrieves the value currently stored in the property @em self
 * as a double complex. The function accepts properties of a complex type
 * with a rank less or equal than the function's return type.
 */

double complex
cpl_property_get_double_complex(const cpl_property *self)
{

    double complex value = 0.;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return value;
    }

    cx_assert(self->value != NULL);


    /*
     * If the property self is of the same kind, but of a smaller rank,
     * promote the property's value to the target type.
     */

    switch (self->value->type) {

        case CPL_TYPE_FLOAT_COMPLEX:
        {
            float complex *result =
                    _cpl_property_value_get(self->value,
                                            CPL_TYPE_FLOAT_COMPLEX,
                                            cpl_func);

            value = *result;
            break;
        }

        case CPL_TYPE_DOUBLE_COMPLEX:
        {
            double complex *result =
                    _cpl_property_value_get(self->value,
                                            CPL_TYPE_DOUBLE_COMPLEX,
                                            cpl_func);

            value = *result;
            break;
        }

        default:
        {
            cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
            break;
        }

    }

    return value;

}
/**@}*/
