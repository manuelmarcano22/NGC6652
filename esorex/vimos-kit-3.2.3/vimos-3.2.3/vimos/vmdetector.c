/* $Id: vmdetector.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string.h>
#include <math.h>

#include <pilmemory.h>
#include <pilstrutils.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include "vmdetector.h"
#include "cpl.h"


#define MAX_COMMENT_LENGTH (80)

/**
 * @name vimosDetector
 *
 * @doc
 *   The module vimosDetector collects low/medium level functions
 *   related to detector properties.
 */

/**@{*/

/**
 * @memo
 *   Creates a new window.
 *
 * @return Pointer to the newly allocated window.
 *
 * @doc
 *   The function allocates and initializes a new list of windows,
 *   containing one empty window.
 *
 * @author C. Izzo
 */

VimosWindow *newWindow(void)
{
  VimosWindow *window;

  window = (VimosWindow *) cpl_malloc(sizeof(VimosWindow));
  if (window) {
    window->startX = 0;
    window->startY = 0;
    window->nX = 0;
    window->nY = 0;
    window->prev = window->next = NULL;
  }
  return(window);
}


/**
 * @memo
 *   Destroys a window.
 *
 * @return Nothing.
 *
 * @param window Pointer to the window to destroy.
 *
 * @doc
 *   The function removes the window from a list, and destroys it
 *   leaving the list unbroken.
 *
 * @author C. Izzo
 */

void deleteWindow(VimosWindow *window)
{
  if (window) {

   /*
    * Detach this element from its linked list
    */

    if (window->prev) window->prev->next = window->next;
    if (window->next) window->next->prev = window->prev;

   /*
    * Destroy it
    */

    cpl_free(window);
    window = NULL;
  }
  return;
}


/**
 * @memo
 *   Destroys a list of windows.
 *
 * @return Nothing.
 *
 * @param window Pointer to the first window of the list to destroy.
 *
 * @doc
 *   The function destroys a list of windows starting from any given
 *   element on, and leaving the rest of the list intact.
 *
 * @author C. Izzo
 */

void deleteWindowList(VimosWindow *window)
{
  if (window) {
    deleteWindowList(window->next);
    deleteWindow(window);
  }
  return;
}


/**
 * @memo
 *   Creates a new port.
 *
 * @return Pointer to the newly allocated port.
 *
 * @doc
 *   The function allocates and initializes a new list of ports,
 *   containing one empty port.
 *   A port consists of four windows: the window defining the port
 *   itself, the readout window (i.e. the port portion that is actually
 *   read), the prescan region, and the overscan region. The windows
 *   coordinates are given in terms of the pixel coordinates of a
 *   related image (starting from (0,0)). Adding the given shifts 
 *   in X and Y the coordinates are converted into chip pixels.
 *
 * @author C. Izzo
 */

VimosPort *newPort(void)
{
  VimosPort *port;

  port = (VimosPort *) cpl_malloc(sizeof(VimosPort));
  if (port) {
    port->window        = newWindow();
    port->prScan        = newWindow();
    port->ovScan        = newWindow();
    port->readOutWindow = newWindow();
    port->shiftX        = 0;
    port->shiftY        = 0;
    port->prev = port->next = NULL;
  }
  return(port);
}


/**
 * @memo
 *   Destroys a port.
 *
 * @return Nothing.
 *
 * @param port Pointer to the port to destroy.
 *
 * @doc
 *   The function removes the port from a list, and destroys it
 *   leaving the list unbroken.
 *
 * @author C. Izzo
 */

void deletePort(VimosPort *port)
{
  if (port) {

   /*
    * Detach this element from its linked list
    */

    if (port->prev) port->prev->next = port->next;
    if (port->next) port->next->prev = port->prev;

   /*
    * Destroy it
    */

    deleteWindow(port->window);
    deleteWindow(port->prScan);
    deleteWindow(port->ovScan);
    deleteWindow(port->readOutWindow);

    cpl_free(port);
  }
  return;
}


/**
 * @memo
 *   Destroys a list of ports.
 *
 * @return Nothing.
 *
 * @param port Pointer to the first port of the list to destroy.
 *
 * @doc
 *   The function destroys a list of ports starting from any given
 *   element on, and leaving the rest of the list intact.
 *
 * @author C. Izzo
 */

void deletePortList(VimosPort *port)
{
  if (port) {
    deletePortList(port->next);
    deletePort(port);
  }
  return;
}

/**
 * @memo
 *   Generate a list of ports from the content of a data image header.
 *
 * @return A list of ports.
 *
 * @param image  Image to be read the port structure from.
 * @param nports Returned number of found ports.
 *
 * @doc
 *   This function returns correct results, as long as the following
 *   assumptions are correct:
 *   \begin{itemize}
 *     \item ESO.DET.OUTPUTS contains the number of readout ports.

 *     \item ESO.DET.OUTi.X, Y contain the starting position of the 
 *           i-th port in chip coordinates (if there is just one 
 *           readout port, the index i is omitted). The starting 
 *           position may correspond to any of the corners of the 
 *           port, and the corner may differ from port to port 
 *           (as in the FORS case), therefore it is internally 
 *           converted always to the lower left corner of the port.
 *
 *     \item ESO.DET.OUTi.NX, NY contain the number of pixel in the 
 *           X and Y directions of the i-th port (if there is just 
 *           one readout port, the index i is omitted).
 *   \end{itemize}
 *
 * Readout windows can be defined on the chip surface independently
 * on the subdivision into ports. In other words, a given readout
 * window may extend over different ports, and it includes prescan 
 * and overscans if present. 
 *
 *   \begin{itemize}
 *     \item ESO.DET.WINDOWS contains the number of readout windows.
 *           It is assumed that just one readout window can be defined
 *           for a given data image (Stefan Bogun, private communication).
 *
 *     \item ESO.DET.WINi.STRX, STRY, NX, NY contain the lower left corner
 *           and the size of the readout window in detector coordinates. 
 *           Detector coordinates are not necessarily coincident with 
 *           chip coordinates, as adding a prescan would define a larger 
 *           "virtual detector". Specifically, if a readout window is 
 *           entirely contained in the chip area, then its coordinates 
 *           coincide with chip coordinates. But if a prescan region 
 *           exists and the readout window extending beyond the chip, 
 *           then the window coordinates will have a shift with respect 
 *           to chip coordinates, equal to the extension of the prescan.
 *
 *     \item ESO.DET.OUT.PRSCX, PRSCY, OVSCX, OVSCY supply the number of
 *           pixels in the prescan and overscan regions (beyond the chip).
 *   \end{itemize}
 *
 * This routine returns a list of ports, and for each port its 
 * own portion of the whole readout window with the corresponding 
 * prescan/overscan window (if a readout window is extending into 
 * different ports, then it is split into more (rectangular) regions, 
 * one for each involved port).
 *
 * All the start coordinates for each region are given in image 
 * pixels (NOT chip pixels), where the lower left pixel is assigned 
 * coordinates (0,0).
 *
 * @author C. Izzo
 */

VimosPort *getPorts(VimosImage *image, int *nports)
{
  char          modName[]            = "getPorts";

  VimosPort    *lastPort             = NULL;
  VimosPort    *currentPort          = NULL;
  VimosPort    *headPort             = NULL;
  VimosWindow  *lastWindow           = NULL;
  VimosWindow  *currentWindow        = NULL;
  VimosWindow  *headWindow           = NULL;
  VimosWindow  *currentReadOutWindow = NULL;
  VimosWindow  *newReadOutWindow     = NULL;
  const char   *dscName              = NULL;
  int           swap[2]              = {0,0};
  int           tx                   = 99999999;
  int           ty                   = 99999999;
  int           noIndex              = 0;
  int           vertical             = 0;
  int           i, j, k, nwin;
  int           status;
  int           npix;
  int           x1, x2, nx1, nx2;

 /*
  * Determine number of readout ports
  */

  *nports = 0;
  if (readIntDescriptor(image->descs, pilTrnGetKeyword("NumberOfPorts"), 
                        nports, image->descs->descComment) == VM_TRUE) {
    cpl_msg_debug(modName, "Number of ports found is %d.", *nports);
    if (*nports < 1 || *nports > MAX_RDPORTS) {
      cpl_msg_error(modName, "Max. number allowed is %d.", MAX_RDPORTS);
      return NULL;
    }
  }
  else {
    cpl_msg_error(modName, "Keyword %s not found.", 
                pilTrnGetKeyword("NumberOfPorts"));
    return NULL;
  }

 /*
  * Determine all the port regions
  */

  status = 0;

  if (findDescriptor(image->descs, pilTrnGetKeyword("PortStartX"))) {
    noIndex = 1;
  }
  for (i = 0; i < *nports; i++) {
    currentPort = newPort();
    if (i == 0) {
      headPort = currentPort;   /* Keep track of list's head */
    }
    else {
      lastPort->next = currentPort;   /* Append a new port to list */
      currentPort->prev = lastPort;
    }
    lastPort = currentPort;

    if (noIndex) {

     /*
      * The current block is just necessary to avoid indexing if there
      * is just one port.
      */

     /*
      * Get start X,Y coords of current port
      */

      if (readIntDescriptor(image->descs, pilTrnGetKeyword("PortStartX"),
                            &(currentPort->window->startX), 
                            image->descs->descComment) == VM_FALSE) {
        status++;
        cpl_msg_error(modName, "Keyword %s not found.", 
                    pilTrnGetKeyword("PortStartX"));
      }

      if (readIntDescriptor(image->descs, pilTrnGetKeyword("PortStartY"),
                             &(currentPort->window->startY), 
                             image->descs->descComment) == VM_FALSE) {
        status++;
        cpl_msg_error(modName, "Keyword %s not found.", 
                    pilTrnGetKeyword("PortStartY"));
      }

     /*
      * Get X and Y size of current port
      */

      if (readIntDescriptor(image->descs, pilTrnGetKeyword("PortSizeX"),
                            &(currentPort->window->nX), 
                            image->descs->descComment) == VM_FALSE) {
        status++;
        cpl_msg_error(modName, "Keyword %s not found.", 
                    pilTrnGetKeyword("PortSizeX"));
      }

      if (readIntDescriptor(image->descs, pilTrnGetKeyword("PortSizeY"),
                            &(currentPort->window->nY), 
                            image->descs->descComment) == VM_FALSE) {
        status++;
        cpl_msg_error(modName, "Keyword %s not found.", 
                    pilTrnGetKeyword("PortSizeY"));
      }

     /*
      * Get extension of X or Y prescan. Note: one of the two should be zero.
      */

      if (readIntDescriptor(image->descs, pilTrnGetKeyword("PortPrscX"),
                            &npix, image->descs->descComment) == VM_TRUE) {
        if (npix > 0) {
          currentPort->prScan->nX     = npix;
          currentPort->prScan->startX = currentPort->window->startX
                                      - currentPort->prScan->nX;
          currentPort->prScan->startY = currentPort->window->startY;
          currentPort->prScan->nY     = currentPort->window->nY;
        }
      }

      if (readIntDescriptor(image->descs, pilTrnGetKeyword("PortPrscY"),
                            &npix, image->descs->descComment) == VM_TRUE) {
        if (npix > 0) {
          if (currentPort->prScan->nX) { 
            cpl_msg_error(modName, "Illegal prescan region");
            return(NULL);
          }
          currentPort->prScan->nY     = npix;
          currentPort->prScan->startY = currentPort->window->startY
                                      - currentPort->prScan->nY;
          currentPort->prScan->startX = currentPort->window->startX;
          currentPort->prScan->nX     = currentPort->window->nX;
        }
      }

     /*
      * Get extension of X or Y overscan. Note: one of the two
      * must be zero.
      */

      if (readIntDescriptor(image->descs, pilTrnGetKeyword("PortOvscX"),
                            &npix, image->descs->descComment) == VM_TRUE) {
        if (npix > 0) {
          currentPort->ovScan->nX     = npix;
          currentPort->ovScan->startX = currentPort->window->startX
                                      + currentPort->window->nX;
          currentPort->ovScan->startY = currentPort->window->startY;
          currentPort->ovScan->nY     = currentPort->window->nY;
        }
      }

      if (readIntDescriptor(image->descs, pilTrnGetKeyword("PortOvscY"),
                            &npix, image->descs->descComment) == VM_TRUE) {
        if (npix > 0) {
          if (currentPort->ovScan->nX) { 
            cpl_msg_error(modName, "Illegal overscan region");
            return(NULL);
          }
          currentPort->ovScan->nY     = npix;
          currentPort->ovScan->startY = currentPort->window->startY
                                      + currentPort->window->nY;
          currentPort->ovScan->startX = currentPort->window->startX;
          currentPort->ovScan->nX     = currentPort->window->nX;
        }
      }
    }
    else {

     /*
      * Here is where ports are in total more than one (or just one, 
      * but indexed).
      */

     /*
      * Get start X,Y coords of current port
      */

      if (readIntDescriptor(image->descs, 
                            pilTrnGetKeyword("SeqPortStartX", i + 1),
                            &(currentPort->window->startX), 
                            image->descs->descComment) == VM_FALSE) {
        status++;
        cpl_msg_error(modName, "Keyword %s not found.", 
                    pilTrnGetKeyword("SeqPortStartX", i + 1));
      }

      if (readIntDescriptor(image->descs, 
                            pilTrnGetKeyword("SeqPortStartY", i + 1),
                            &(currentPort->window->startY), 
                            image->descs->descComment) == VM_FALSE) {
        status++;
        cpl_msg_error(modName, "Keyword %s not found.", 
                    pilTrnGetKeyword("SeqPortStartY", i + 1));
      }

     /*
      * Get X and Y size of current port
      */

      if (readIntDescriptor(image->descs, 
                            pilTrnGetKeyword("SeqPortSizeX", i + 1),
                            &(currentPort->window->nX), 
                            image->descs->descComment) == VM_FALSE) {
        status++;
        cpl_msg_error(modName, "Keyword %s not found.", 
                    pilTrnGetKeyword("SeqPortSizeX", i + 1));
      }

      if (readIntDescriptor(image->descs, 
                            pilTrnGetKeyword("SeqPortSizeY", i + 1),
                            &(currentPort->window->nY), 
                            image->descs->descComment) == VM_FALSE) {
        status++;
        cpl_msg_error(modName, "Keyword %s not found.", 
                    pilTrnGetKeyword("SeqPortSizeY", i + 1));
      }

     /*
      * Here we standardize the start coordinates of the port:
      * we want them to be always at the lower left corner of 
      * the port, and this is not always the case. What follows
      * is valid just under the assumption that at least the first 
      * port has its start coordinates at its lower left corner, as
      * is hopefully always the case: we MUST make this assumption,
      * as this information is NOT currently present in the data
      * header! In this way, we can always assume that for each
      * new port found, the previous one was already assigned the
      * starting point at the lower left corner.
      */

      if (i > 0) {  /* Not necessary to enter this block for the first port */

        swap[0] = swap[1] = 0;

       /*
        * swap[] is meant to record the modifications done to
        * standardize the port starting X and Y coordinates to 
        * its lower left corner, and therefore understand whether 
        * a prescan should be seen as an overscan, and viceversa.
        */

        for (j = 0; j < 2; j++) {
          if (j) {                                     /* Y coordinate */
            x1  = currentPort->prev->window->startY;
            x2  = currentPort->window->startY;
            nx1 = currentPort->prev->window->nY;
            nx2 = currentPort->window->nY;
          } 
          else {                                       /* X coordinate */
            x1  = currentPort->prev->window->startX;
            x2  = currentPort->window->startX;
            nx1 = currentPort->prev->window->nX;
            nx2 = currentPort->window->nX;
          }
          if (x2 - x1 == nx1 + nx2 - 1) {
            x2 -= (nx2 - 1);
            swap[j] = 1;
          }
          else if (x2 - x1 == nx1 + nx2) {

           /*
            * PLEASE NOTE: this part was added here just to patch the 
            * wrong start Y coordinates for ports 3 and 4 in FORS, set
            * to 2049 instead of 2048. This block should in principle
            * be removed the day that that mistake is solved, but it
            * may be even wiser to leave it here, as it fixes easily
            * this kind of common mistakes.
            */

            x2 -= nx2;
            swap[j] = 1;
          }
          else if (x2 - x1 == nx1 - 1 || x2 - x1 == nx1) {

           /*                        ^^^^^^^^^^^^^^^^^^
            * PLEASE NOTE: the underscored part was added just to patch 
            * wrong start Y coordinates for ports 3 and 4 in FORS, set
            * to 2049 instead of 2048. This part should in principle
            * be removed the day that that mistake is solved, but it
            * may be even wiser to leave it here, as it fixes easily 
            * this kind of common mistakes.
            */

            x2 = x1;
            swap[j] = 1;
          }
          else if (x1 - x2 == 1) {
            x2 -= (nx2 - 1);
            swap[j] = 1;
          }
          if (j) {
            currentPort->window->startY = x2;
          } 
          else {
            currentPort->window->startX = x2;
          }
        }
      }

     /*
      * Get extension of X or Y prescan. Note: one of the two
      * must be zero.
      */

      if (swap[0]) {

       /* 
        * If the current X corner coordinate was swapped, read
        * the overscan as a prescan.
        */

        dscName = pilTrnGetKeyword("SeqPortOvscX", i + 1);
      }
      else {
        dscName = pilTrnGetKeyword("SeqPortPrscX", i + 1);
      }

      if (readIntDescriptor(image->descs, dscName,
                            &npix, image->descs->descComment) == VM_TRUE) {
        if (npix > 0) {
          vertical = 0;
          currentPort->prScan->nX     = npix;
          currentPort->prScan->startX = currentPort->window->startX
                                      - currentPort->prScan->nX;
          currentPort->prScan->startY = currentPort->window->startY;
          currentPort->prScan->nY     = currentPort->window->nY;
        }
      }

      if (swap[1]) {

       /* 
        * If the current Y corner coordinate was swapped, read
        * the overscan as a prescan.
        */

        dscName = pilTrnGetKeyword("SeqPortOvscY", i + 1);
      }
      else {
        dscName = pilTrnGetKeyword("SeqPortPrscY", i + 1);
      }

      if (readIntDescriptor(image->descs, dscName,
                            &npix, image->descs->descComment) == VM_TRUE) {
        if (npix > 0) {
          vertical = 1;
          if (currentPort->prScan->nX) {
            cpl_msg_error(modName, "Illegal prescan region");
            return NULL;
          }
          currentPort->prScan->nY     = npix;
          currentPort->prScan->startY = currentPort->window->startY
                                      - currentPort->prScan->nY;
          currentPort->prScan->startX = currentPort->window->startX;
          currentPort->prScan->nX     = currentPort->window->nX;
        }
      }

     /*
      * Get extension of X or Y overscan. Note: one of the two
      * must be zero.
      */

      if (swap[0]) {

       /* 
        * If the current X corner coordinate was swapped, read
        * the prescan as an overscan.
        */

        dscName = pilTrnGetKeyword("SeqPortPrscX", i + 1);
      }
      else {
        dscName = pilTrnGetKeyword("SeqPortOvscX", i + 1);
      }

      if (readIntDescriptor(image->descs, dscName,
                            &npix, image->descs->descComment) == VM_TRUE) {
        if (npix > 0) {
          if (vertical) {
             cpl_msg_error(modName, 
                         "Overscan region inconsistent with prescan region");
             return NULL;
          }
          currentPort->ovScan->nX     = npix;
          currentPort->ovScan->startX = currentPort->window->startX
                                      + currentPort->window->nX;
          currentPort->ovScan->startY = currentPort->window->startY;
          currentPort->ovScan->nY     = currentPort->window->nY;
        }
      }
      if (swap[1]) {

       /* 
        * If the current Y corner coordinate was swapped, read
        * the prescan as an overscan.
        */

        dscName = pilTrnGetKeyword("SeqPortPrscY", i + 1);
      }
      else {
        dscName = pilTrnGetKeyword("SeqPortOvscY", i + 1);
      }

      if (readIntDescriptor(image->descs, dscName,
                            &npix, image->descs->descComment) == VM_TRUE) {
        if (npix > 0) {
          if (currentPort->ovScan->nX) {
            cpl_msg_error(modName, "Illegal overscan region");
            return NULL;
          }
          if (!vertical) {
             cpl_msg_error(modName, 
                         "Overscan region inconsistent with prescan region");
             return NULL;
          }
          currentPort->ovScan->nY     = npix;
          currentPort->ovScan->startY = currentPort->window->startY
                                      + currentPort->window->nY;
          currentPort->ovScan->startX = currentPort->window->startX;
          currentPort->ovScan->nX     = currentPort->window->nX;
        }
      }
    }

   /*
    * tx and ty should contain the lowest coordinate in image:
    */

    if (currentPort->prScan->nX) {  /* i.e., if ANY (X or Y) prescan exists! */
      if (vertical) {
        if (currentPort->prScan->startY < ty) ty = currentPort->prScan->startY;
      }
      else {
        if (currentPort->prScan->startX < tx) tx = currentPort->prScan->startX;
      }
    }
  }

 /*
  * Patch (an attempt to cope with the twisted logic of the keyword header):
  * size of unic port is always the size of the chip. This is solving a
  * compatibility problem with FORS data.
  */

 /* FIXME: it must be disabled to let vmDet work. It's actually vmDet that
  *        should be probably modified

  if (*nports == 1) {
    if (getChipSize(image, &(headPort->window->nX), &(headPort->window->nY)) 
                     == EXIT_FAILURE) {
      status++;
      cpl_msg_error(modName, "Cannot read chip size.");
    }
  }

  * End Of FIXME */

 /*
  * Determine all the readout regions
  */

  nwin = 0;
  if (readIntDescriptor(image->descs, pilTrnGetKeyword("NumberOfWindows"), 
                        &nwin, image->descs->descComment) == VM_TRUE) {
    cpl_msg_debug(modName, "Number of readout windows found is %d", nwin);
    if (nwin != 1) {
      cpl_msg_error(modName, 
                  "Just ONE read out window is allowed - %d found", nwin);
      return NULL;
    }
  }
  else {
    cpl_msg_error(modName, "Keyword %s not found", 
                pilTrnGetKeyword("NumberOfWindows"));
    return NULL;
  }

 /*
  * Using loop for easier (future) generalization - just one readout 
  * window is currently expected.
  */

  for (i = 0; i < nwin; i++) {
    currentWindow = newWindow();
    if (i == 0) {
      headWindow = currentWindow;         /* Keep track of list's head   */
    }
    else {
      lastWindow->next = currentWindow;   /* Append a new window to list */
      currentWindow->prev = lastWindow;
    }
    lastWindow = currentWindow;

    if (readIntDescriptor(image->descs, 
                          pilTrnGetKeyword("SeqWindowStartX", i + 1),
                          &(currentWindow->startX), 
                          image->descs->descComment) == VM_FALSE) {
      status++;
      cpl_msg_error(modName, "Keyword %s not found.", 
                  pilTrnGetKeyword("SeqWindowStartX", i + 1));
    }

    if (readIntDescriptor(image->descs, 
                          pilTrnGetKeyword("SeqWindowStartY", i + 1),
                          &(currentWindow->startY), 
                          image->descs->descComment) == VM_FALSE) {
      status++;
      cpl_msg_error(modName, "Keyword %s not found.", 
                  pilTrnGetKeyword("SeqWindowStartY", i + 1));
    }

    if (readIntDescriptor(image->descs, 
                          pilTrnGetKeyword("SeqWindowSizeX", i + 1),
                          &(currentWindow->nX), 
                          image->descs->descComment) == VM_FALSE) {
      status++;
      cpl_msg_error(modName, "Keyword %s not found.", 
                  pilTrnGetKeyword("SeqWindowSizeX", i + 1));
    }

    if (readIntDescriptor(image->descs, 
                          pilTrnGetKeyword("SeqWindowSizeY", i + 1),
                          &(currentWindow->nY), 
                          image->descs->descComment) == VM_FALSE) {
      status++;
      cpl_msg_error(modName, "Keyword %s not found.", 
                  pilTrnGetKeyword("SeqWindowSizeY", i + 1));
    }

   /*
    * tx and ty should contain the lowest coordinate in image
    */

    if (currentWindow->startX < tx) tx = currentWindow->startX;
    if (currentWindow->startY < ty) ty = currentWindow->startY;
  }

 /*
  * If any of the above entries were not found, abort the program.
  */

  if (status > 0) 
    return NULL;

 /*
  * Now find the intersections of the readout window with each port
  */

  cpl_msg_debug(modName, "The following regions are found:");
  cpl_msg_indent_more();
  cpl_msg_indent_more();
  for (currentPort = headPort, i = 1, k = 1; currentPort;
                                        currentPort=currentPort->next, i++) {
    currentReadOutWindow = currentPort->readOutWindow;
    cpl_msg_indent_less();
    for (currentWindow = headWindow; currentWindow;
                                       currentWindow = currentWindow->next) {
      currentReadOutWindow->startX =
        MAX(currentPort->window->startX, currentWindow->startX);
      currentReadOutWindow->startY =
        MAX(currentPort->window->startY, currentWindow->startY);
      currentReadOutWindow->nX =
        MIN(currentPort->window->startX + currentPort->window->nX,
            currentWindow->startX + currentWindow->nX)
            - currentPort->readOutWindow->startX;
      currentReadOutWindow->nY =
        MIN(currentPort->window->startY + currentPort->window->nY,
            currentWindow->startY + currentWindow->nY)
            - currentPort->readOutWindow->startY;

      if (currentReadOutWindow->nX > 0 && currentReadOutWindow->nY > 0) {

       /*
        * If the intersection is not empty, convert all start
        * coordinates to image coordinates (starting from (0,0) ).
        */

        currentReadOutWindow->startX -= tx;
        currentReadOutWindow->startY -= ty;

       /*
        * Monitoring (debug):
        */

        cpl_msg_debug(modName,
                "Window %d             : start = %4d,%-4d  npix = %4d,%-4d\n",
                 k,
                 currentReadOutWindow->startX,
                 currentReadOutWindow->startY,
                 currentReadOutWindow->nX,
                 currentReadOutWindow->nY);
        k++;

       /*
        * Create, and position to, next readout window of current port.
        */

        newReadOutWindow = newWindow();
        currentReadOutWindow->next = newReadOutWindow;
        newReadOutWindow->prev = currentReadOutWindow;
        currentReadOutWindow = newReadOutWindow;
      }
      else {
        cpl_msg_debug(modName,
                "Readout window %d,%d,%d,%d not including port %d,%d,%d,%d",
                 currentWindow->startX,
                 currentWindow->startY,
                 currentWindow->nX,
                 currentWindow->nY,
                 currentPort->window->startX,
                 currentPort->window->startY,
                 currentPort->window->nX,
                 currentPort->window->nY);
        currentReadOutWindow->next = NULL;
      }
    }
    cpl_msg_indent_more();
    if (currentReadOutWindow->prev) {
      currentReadOutWindow->prev->next = NULL;
      deleteWindow(currentReadOutWindow);
    }

   /*
    * Now convert all start coordinates to image coordinates
    * (starting from (0,0) ).
    */

    currentPort->window->startX -= tx;
    currentPort->window->startY -= ty;

    if (currentPort->prScan->nX > 0) {
      currentPort->prScan->startX -= tx;
      currentPort->prScan->startY -= ty;
    }
    else {
      currentPort->prScan->startX = 0;
      currentPort->prScan->startY = 0;
    }

    if (currentPort->ovScan->nX > 0) {
      currentPort->ovScan->startX -= tx;
      currentPort->ovScan->startY -= ty;
    }
    else {
      currentPort->ovScan->startX = 0;
      currentPort->ovScan->startY = 0;
    }

    currentPort->shiftX = tx;
    currentPort->shiftY = ty;

    cpl_msg_debug(modName,
    "in Port %d          : start = %4d,%-4d  npix = %4d,%-4d",i,
    currentPort->window->startX,currentPort->window->startY,
    currentPort->window->nX,currentPort->window->nY);
    cpl_msg_debug(modName,
    "with prescan region: start = %4d,%-4d  npix = %4d,%-4d",
    currentPort->prScan->startX,currentPort->prScan->startY,
    currentPort->prScan->nX,currentPort->prScan->nY);
    cpl_msg_debug(modName,
    "and overscan region: start = %4d,%-4d  npix = %4d,%-4d\n",
    currentPort->ovScan->startX,currentPort->ovScan->startY,
    currentPort->ovScan->nX,currentPort->ovScan->nY);
  }
  cpl_msg_indent_less();
  cpl_msg_indent_less();

  return(headPort);
}

/**
 * @memo
 *   Get total size of read out window through all ports.
 *
 * @return Total number of pixels contained in read out window.
 *
 * @param ports         List of ports.
 * @param startWindowX  (returned) Start X position of total read out window.
 * @param startWindowY  (returned) Start Y position of total read out window.
 * @param sizeWindowX   (returned) X size of total read out window.
 * @param sizeWindowY   (returned) Y size of total read out window.
 *
 * @doc
 *   Given a ports list, position and size of the whole read out
 *   window (sum of the read out windows within each port) are
 *   returned.
 *
 * @author C. Izzo
 */

int getTotalReadoutWindow(VimosPort *ports, int *startWindowX,
    int *startWindowY, int *sizeWindowX, int *sizeWindowY)
{

  VimosPort *currPort;
  int        startPortsX, startPortsY, sizePortsX, sizePortsY;
  int        endPortsX, endPortsY;

  if (!ports)
    return 0;

  startPortsX = ports->readOutWindow->startX;
  startPortsY = ports->readOutWindow->startY;
  endPortsX = ports->readOutWindow->startX + ports->readOutWindow->nX;
  endPortsY = ports->readOutWindow->startY + ports->readOutWindow->nY;

  for (currPort = ports->next; currPort; currPort = currPort->next) {

    if (startPortsX > currPort->readOutWindow->startX)
      startPortsX = currPort->readOutWindow->startX;

    if (startPortsY > currPort->readOutWindow->startY)
      startPortsY = currPort->readOutWindow->startY;

    if (endPortsX < currPort->readOutWindow->startX
                  + currPort->readOutWindow->nX)
      endPortsX = currPort->readOutWindow->startX
                + currPort->readOutWindow->nX;

    if (endPortsY < currPort->readOutWindow->startY
                  + currPort->readOutWindow->nY)
      endPortsY = currPort->readOutWindow->startY
                + currPort->readOutWindow->nY;

  }

  sizePortsX = endPortsX - startPortsX;
  sizePortsY = endPortsY - startPortsY;

  *startWindowX = startPortsX;
  *startWindowY = startPortsY;
  *sizeWindowX = sizePortsX;
  *sizeWindowY = sizePortsY;

  return sizePortsX * sizePortsY;

}


/**
 * @memo
 *   The bias level is estimated on the found overscan regions.
 *
 * @return Array of bias levels, one for each port.
 *
 * @param image     Image where the bias levels must be estimated.
 * @param port      Image port structure.
 *
 * @doc
 *   The function computes for each port the average bias from all 
 *   its available overscan regions. If no overscan region is found
 *   for a port, a null pointer is returned.
 *
 * @author C. Izzo
 */

VimosFloatArray *estimateImageBias(VimosImage *image, VimosPort *port)
{
  const char       modName[] = "estimateImageBias";
  VimosPort       *currentPort;
  VimosBool        success;
  VimosFloatArray *biasLevel;
  float           *bias;
  float            prscBiasLevel, ovscBiasLevel;
  int              sizePrsc;
  int              sizeOvsc;
  int              nPorts = 0;
  int              i = 0;

  if (!(image && port)) {
    cpl_msg_debug(modName, "NULL input(s)");
    return NULL;
  }

  for (currentPort = port; currentPort; currentPort = currentPort->next) 
    nPorts++;

  if (!(biasLevel = newFloatArray(nPorts))) {
    cpl_msg_debug(modName, "Cannot allocate output");
    return NULL;
  }

  for (currentPort = port; currentPort; currentPort = currentPort->next) {
    success = VM_FALSE;

    sizePrsc = 0;
    sizeOvsc = 0;

    if (currentPort->prScan->nX > 0) {

     /*
      * Prescan region exists - extract it and get its average
      */

      bias = extractFloatImage(image->data, image->xlen, image->ylen,
             currentPort->prScan->startX, currentPort->prScan->startY,
             currentPort->prScan->nX, currentPort->prScan->nY);

      if (!bias) {
        cpl_msg_debug(modName, "Memory allocation failure");
        return NULL;
      }

      sizePrsc = currentPort->prScan->nX * currentPort->prScan->nY;
      prscBiasLevel = computeAverageFloat(bias, sizePrsc);
      cpl_free(bias);
      success = VM_TRUE;
    }
    else {
      prscBiasLevel = 0.;
    }

    if (currentPort->ovScan->nX > 0) {

     /*
      * Overscan region exists - extract it and get its average
      */

      bias = extractFloatImage(image->data, image->xlen, image->ylen,
             currentPort->ovScan->startX, currentPort->ovScan->startY,
             currentPort->ovScan->nX, currentPort->ovScan->nY);

      if (!bias) {
        cpl_msg_debug(modName, "Memory allocation failure");
        return NULL;
      }

      sizeOvsc = currentPort->ovScan->nX * currentPort->ovScan->nY;
      ovscBiasLevel = computeAverageFloat(bias, sizeOvsc);
      cpl_free(bias);
      success = VM_TRUE;
    }
    else {
      ovscBiasLevel = 0.;
    }

    if (success) {
      biasLevel->data[i] = 
                   (prscBiasLevel * sizePrsc + ovscBiasLevel * sizeOvsc)
                   / (sizePrsc + sizeOvsc);
      i++;
    }

  }

  if (i != nPorts) {
    deleteFloatArray(biasLevel);
    return(NULL);
  }

  return(biasLevel);

}


/**
 * @memo
 *   Subtract an estimated bias level from input image.
 *
 * @return VM_TRUE on success, VM_FALSE otherwise.
 *
 * @param a    Image data.
 * @param nx   Image x size in pixel.
 * @param ny   Image y size in pixel.
 * @param port Image port structure.
 *
 * @doc
 *   The average bias level is found for each overscan of each port.
 *   Each value is subtracted from the corresponding overscan region 
 *   it was derived from, while the average bias level computed from 
 *   all overscan regions of a given port is subtracted from the port
 *   window.
 *
 * @author C. Izzo
 */

VimosBool subtractOverscan(float a[], int nx, int ny, VimosPort *port)
{
  VimosPort *currentPort;
  VimosBool  success = VM_FALSE;
  float     *bias;
  float     *portWindow;
  float      prscBiasLevel, ovscBiasLevel, biasLevel;
  int        sizePrsc;
  int        sizeOvsc;
  int        i,n;

  for (currentPort = port; currentPort; currentPort = currentPort->next) {
    sizePrsc = 0;
    sizeOvsc = 0;
    if (currentPort->prScan->nX > 0) {

     /*
      * Prescan region exists - extract it and get its average
      */

      bias = extractFloatImage(a, nx, ny,
             currentPort->prScan->startX,
             currentPort->prScan->startY,
             currentPort->prScan->nX,
             currentPort->prScan->nY);

      sizePrsc = currentPort->prScan->nX * currentPort->prScan->nY;
      prscBiasLevel = computeAverageFloat(bias, sizePrsc);

      for (i = 0; i < sizePrsc; i++) bias[i] -= prscBiasLevel;

      insertFloatImage(a, nx, ny,
             currentPort->prScan->startX,
             currentPort->prScan->startY,
             currentPort->prScan->nX,
             currentPort->prScan->nY, bias);

      cpl_free(bias);
      success = VM_TRUE;
    }
    else {
      prscBiasLevel = 0.;
    }
    if (currentPort->ovScan->nX > 0) {

     /*
      * Overscan region exists - extract it and get its average
      */

      bias = extractFloatImage(a, nx, ny,
             currentPort->ovScan->startX,
             currentPort->ovScan->startY,
             currentPort->ovScan->nX,
             currentPort->ovScan->nY);

      sizeOvsc = currentPort->ovScan->nX * currentPort->ovScan->nY;
      ovscBiasLevel = computeAverageFloat(bias, sizeOvsc);

      for (i = 0; i < sizeOvsc; i++) bias[i] -= ovscBiasLevel;

      insertFloatImage(a, nx, ny,
             currentPort->ovScan->startX,
             currentPort->ovScan->startY,
             currentPort->ovScan->nX,
             currentPort->ovScan->nY, bias);

      cpl_free(bias);
      success = VM_TRUE;
    }
    else {
      ovscBiasLevel = 0.;
    }

    if (success) {
      biasLevel = (prscBiasLevel * sizePrsc + ovscBiasLevel * sizeOvsc)
                / (sizePrsc + sizeOvsc);

      portWindow = extractFloatImage(a, nx, ny,
             currentPort->window->startX,
             currentPort->window->startY,
             currentPort->window->nX,
             currentPort->window->nY);

      n = currentPort->window->nX * currentPort->window->nY;
      for (i = 0; i < n; i++) portWindow[i] -= biasLevel;

      insertFloatImage(a, nx, ny,
             currentPort->window->startX,
             currentPort->window->startY,
             currentPort->window->nX,
             currentPort->window->nY, portWindow);

      cpl_free(portWindow);
    }
    else break;
  }
  return(success);
}


/**
 * @memo
 *   Get scan direction of a given port.
 *
 * @return 1 = the scan direction is vertical, 0 = the scan direction
 *         is horizontal.
 *
 * @param port  Input port.
 *
 * @doc
 *   The function derive the scan direction for a given port
 *   (not for a list of ports!). This is derived from the position
 *   of the prescan/overscan regions, under the assumption that
 *   they are always extensions of the readout window along the 
 *   readout direction. If neither a prescan nor an overscan region
 *   is present, the scan direction is always returned as vertical.
 *
 * @author C. Izzo
 */

int getReadoutDirection(VimosPort *port)
{
  int vertical = 1;

  if (port->prScan->nX > 0 
                           && port->window->startY == port->prScan->startY) {
    vertical = 0;
  }
  else if (port->ovScan->nX > 0 
                           && port->window->startY == port->ovScan->startY) {
    vertical = 0;
  }

  return(vertical);
}


/**
 * @memo
 *   Get the chip extension along the x and y directions.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param image     Input image
 * @param chipSizeX X chip size to be returned.
 * @param chipSizeY Y chip size to be returned.
 *
 * @doc
 *   The function simply reads the information from the image header.
 *   With an instrument consisting of more than one chip, the appropriate
 *   chip index is searched and found.
 *
 * @author C. Izzo
 */

int getChipSize(VimosImage *image, int *chipSizeX, int *chipSizeY) 
{

  const char modName[] = "getChipSize";

  char      *dscNameChipX = NULL;
  char      *dscNameChipY = NULL;
  char       comment[MAX_COMMENT_LENGTH];
  int        status = EXIT_FAILURE;
  int        i, nChips;

  dscNameChipX = cpl_strdup(pilTrnGetKeyword("CHIP1.NX"));
  if (findDescriptor(image->descs, dscNameChipX)) {
    dscNameChipY = cpl_strdup(pilTrnGetKeyword("CHIP1.NY"));
  }
  else {
    cpl_free((char *) dscNameChipX);

   /*
    * With more than one chip around, look for the sequence
    * number of the current chip.
    */

    if (readIntDescriptor(image->descs, pilTrnGetKeyword("NCHIPS"), 
                          &nChips, comment) == VM_TRUE) {
      for (i = 0; i < nChips; i++) {
        dscNameChipX = cpl_strdup(pilTrnGetKeyword("CHIPi.NX", i + 1));
        if (findDescriptor(image->descs, dscNameChipX)) {
          dscNameChipY = cpl_strdup(pilTrnGetKeyword("CHIPi.NY", i + 1));
          break;     /* Related descriptors are found */
        }
        else {
          cpl_free(dscNameChipX);
          dscNameChipX = NULL;
        }
      } 
    }
    else {
      cpl_msg_debug(modName, "Unable to read keyword %s",
                  pilTrnGetKeyword("NCHIPS"));
    }
  }

  if (readIntDescriptor(image->descs, dscNameChipX, chipSizeX, comment) 
      == VM_TRUE) {
    if (readIntDescriptor(image->descs, dscNameChipY, chipSizeY, comment)
        == VM_TRUE) {
      status = EXIT_SUCCESS;
    }
  }

  cpl_free(dscNameChipX);
  cpl_free(dscNameChipY);

  return status;

}


/**
 * @memo
 *   The RON is estimated on the found overscan regions.
 *
 * @return Array of RONs, one for each port, measured in ADU.
 *
 * @param image     Image where the RON must be estimated.
 * @param port      Image port structure.
 *
 * @doc
 *   The function computes for each port the variance from all
 *   its available overscan regions. For a given port, all the 
 *   found variances are averaged, and the square root of this
 *   average is written to an array. If no overscan region is 
 *   found for a port, a null pointer is returned.
 *
 * @author C. Izzo
 */

VimosFloatArray *estimateImageRon(VimosImage *image, VimosPort *port)
{
  const char       modName[] = "estimateImageRon";
  VimosPort       *currentPort;
  VimosBool        success;
  VimosFloatArray *ron;
  float           *bias;
  float            prscRon, ovscRon;
  int              sizePrsc;
  int              sizeOvsc;
  int              nPorts = 0;
  int              i = 0;

  if (!(image && port)) {
    cpl_msg_debug(modName, "NULL input(s)");
    return NULL;
  }

  for (currentPort = port; currentPort; currentPort = currentPort->next) 
    nPorts++;

  if (!(ron = newFloatArray(nPorts))) {
    cpl_msg_debug(modName, "Cannot allocate output");
    return NULL;
  }

  for (currentPort = port; currentPort; currentPort = currentPort->next) {
    success = VM_FALSE;

    sizePrsc = 0;
    sizeOvsc = 0;
    if (currentPort->prScan->nX > 0) {

     /*
      * Prescan region exists - extract it and get its variance
      */

      bias = extractFloatImage(image->data, image->xlen, image->ylen,
             currentPort->prScan->startX, currentPort->prScan->startY,
             currentPort->prScan->nX, currentPort->prScan->nY);

      if (!bias) {
        cpl_msg_debug(modName, "Memory allocation failure");
        return NULL;
      }

      sizePrsc = currentPort->prScan->nX * currentPort->prScan->nY;
      prscRon = computeVarianceFloat2D(bias, currentPort->prScan->nX, 
                                       currentPort->prScan->nY);

      cpl_free(bias);
      success = VM_TRUE;
    }
    else {
      prscRon = 0.;
    }

    if (currentPort->ovScan->nX > 0) {

     /*
      * Overscan region exists - extract it and get its average
      */

      bias = extractFloatImage(image->data, image->xlen, image->ylen,
             currentPort->ovScan->startX, currentPort->ovScan->startY,
             currentPort->ovScan->nX, currentPort->ovScan->nY);

      if (!bias) {
        cpl_msg_debug(modName, "Memory allocation failure");
        return NULL;
      }

      sizeOvsc = currentPort->ovScan->nX * currentPort->ovScan->nY;
      ovscRon = computeVarianceFloat2D(bias, currentPort->ovScan->startX, 
                                       currentPort->ovScan->startY);

      cpl_free(bias);
      success = VM_TRUE;
    }
    else {
      ovscRon = 0.;
    }

    if (success) {
      ron->data[i] = sqrt((prscRon * sizePrsc + ovscRon * sizeOvsc)
                   / (sizePrsc + sizeOvsc));
      i++;
    }

  }

  if (i != nPorts) {
    deleteFloatArray(ron);
    return NULL;
  }

  return ron;

}

VimosFloatArray *getImageRon(VimosImage *image)
{
  const char       modName[] = "getImageRon";
  char             comment[MAX_COMMENT_LENGTH];
  VimosFloatArray *ron;
  double           dValue;
  int              nPorts = 0;
  int              i = 0;

  if (!image) {
    cpl_msg_debug(modName, "NULL input");
    return NULL;
  }

  if ((readIntDescriptor(image->descs, pilTrnGetKeyword("NumberOfPorts"),
                         &nPorts, comment) == VM_FALSE)) {
    return NULL;
  }

  if (!(ron = newFloatArray(nPorts))) {
    cpl_msg_debug(modName, "Cannot allocate output");
    return NULL;
  }

  for (i = 0; i < nPorts; i++) {
    if (readDoubleDescriptor(image->descs, 
                             pilTrnGetKeyword("SeqReadNoise", i + 1),
                             &dValue, comment) == VM_FALSE) {
      deleteFloatArray(ron);
      return NULL;
    }
    ron->data[i] = (float)dValue;
  }

  return ron;
}
/**@}*/
