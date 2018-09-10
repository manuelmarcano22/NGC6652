 /*
 				define.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN, IAP & Leiden observatory
*
*	Contents:	global definitions.
*
*	Last modify:	13/10/99
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*------------------------ what, who, when and where ------------------------*/

#define         BANNER          "PSFEx"
#define         VERSION         "1.8.1 (October 13, 1999)"
#define         COPYRIGHT       "Emmanuel BERTIN (bertin@iap.fr)"
#define         INSTITUTE       "IAP"

/*----------------------------- Internal constants --------------------------*/
#define		OUTPUT		stderr		/* where all msgs are sent */
#define		BIG		1e+30		/* a huge number */
#define		MAXCHAR		256		/* max. number of characters */
#define		MAXCHECK	16		/* max. # of CHECKimages */
#define		MAXCONTEXT	8		/* max. # of context keys */

/*------------ Set defines according to machine's specificities -------------*/

#ifdef DEC_ALPHA
#define BSWAP
#define INT64
#define	UNIX
#endif

#ifdef HP_UX
#define UNIX
#endif

#ifdef PC_LINUX
#define BSWAP
#define UNIX
#endif

#ifdef SUN_ULTRASPARC
#define INT64
#define	UNIX
#endif

#if 0
#define	NO_ENVVAR
#endif

/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)

/*--------------------- in case of missing constants ------------------------*/

#ifndef PI
#define PI		3.1415926535898	/* never met before? */
#endif

#ifndef         SEEK_SET
#define         SEEK_SET        0
#endif
#ifndef         SEEK_CUR
#define         SEEK_CUR        1
#endif

#ifndef EXIT_SUCCESS
#define 	EXIT_SUCCESS	0
#endif
#ifndef EXIT_FAILURE
#define		EXIT_FAILURE	-1
#endif

/*------------------- a few definitions to read FITS parameters ------------*/

#define	FITSTOF(k, def)	((point = fitsnfind(buf, k, n))? \
					 atof(strncpy(st, &point[10], 70)) \
					:(def))
#define	FITSTOI(k, def)	((point = fitsnfind(buf, k, n))? \
					 atoi(strncpy(st, &point[10], 70)) \
					:(def))

#define	FITSTOS(k, str, def) \
                { if (fitsread(buf,k,str,H_STRING,T_STRING)!= RETURN_OK) \
                    strcpy(str, (def)); \
                }

/*------------------------------- Other Macros -----------------------------*/

#define	DEXP(x)	exp(2.30258509299*(x))	/* 10^x */

#define QFREAD(ptr, size, afile, fname) \
		if (fread(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while reading ", fname)

#define QFWRITE(ptr, size, afile, fname) \
		if (fwrite(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while writing ", fname)

#define	QFSEEK(afile, offset, pos, fname) \
		if (fseek(afile, (offset), pos)) \
		  error(EXIT_FAILURE,"*Error*: file positioning failed in ", \
			fname)

#define	QFTELL(pos, afile, fname) \
		if ((pos=ftell(afile))==-1) \
		  error(EXIT_FAILURE,"*Error*: file position unknown in ", \
			fname)

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ)))) \
		   error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
                  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
                    error(EXIT_FAILURE, "Not enough memory for ", \
                        #ptrout " (" #nel " elements) !"); \
                   memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ));};;}

#define	RINT(x)	(int)(floor(x+0.5))

#define	PIX(x, y)	field->strip[(((int)y)%field->stripheight) \
				*field->width +(int)x]

#define	NPRINTF		if (prefs.verbose_type == NORM) fprintf

#define	NFPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[1M> %s\n\33[1A",x); \
			else if (prefs.verbose_type == FULL) \
				fprintf(w, "%s.\n", x);}

#define	FPRINTF		if (prefs.verbose_type == FULL)	fprintf
