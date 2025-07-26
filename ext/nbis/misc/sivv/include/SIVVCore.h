/*******************************************************************************

License: 
This software and/or related materials was developed at the National Institute
of Standards and Technology (NIST) by employees of the Federal Government
in the course of their official duties. Pursuant to title 17 Section 105
of the United States Code, this software is not subject to copyright
protection and is in the public domain. 

This software and/or related materials have been determined to be not subject
to the EAR (see Part 734.3 of the EAR for exact details) because it is
a publicly available technology and software, and is freely distributed
to any interested party with no licensing requirements.  Therefore, it is 
permissible to distribute this software as a free download from the internet.

Disclaimer: 
This software and/or related materials was developed to promote biometric
standards and biometric technology testing for the Federal Government
in accordance with the USA PATRIOT Act and the Enhanced Border Security
and Visa Entry Reform Act. Specific hardware and software products identified
in this software were used in order to perform the software development.
In no case does such identification imply recommendation or endorsement
by the National Institute of Standards and Technology, nor does it imply that
the products and equipment identified are necessarily the best available
for the purpose.

This software and/or related materials are provided "AS-IS" without warranty
of any kind including NO WARRANTY OF PERFORMANCE, MERCHANTABILITY,
NO WARRANTY OF NON-INFRINGEMENT OF ANY 3RD PARTY INTELLECTUAL PROPERTY
or FITNESS FOR A PARTICULAR PURPOSE or for any purpose whatsoever, for the
licensed product, however used. In no event shall NIST be liable for any
damages and/or costs, including but not limited to incidental or consequential
damages of any kind, including economic damage or injury to property and lost
profits, regardless of whether NIST shall be advised, have reason to know,
or in fact shall know of the possibility.

By using this software, you agree to bear all risk relating to quality,
use and performance of the software and/or related materials.  You agree
to hold the Government harmless from any claim arising from your use
of the software.

*******************************************************************************/

/*******************************************************************************
	LIBRARY:	SIVVCore - Spectral Image Validation/Verification (Core)

	FILE:		SIVVCORE.H

	AUTHORS:	Joseph C. Konczal 
				John D. Grantham
				Modified by Jiong Zhang (Seventh Sense AI) for C++ compatibility
	
	DATE:		01/05/2009 (JCK)
	UPDATED:	01/19/2009 (JDG)
				09/10/2009 (JDG)
				09/03/2010 (JDG)
				07/07/2011 (JDG)

	DESCRIPTION:

		Contains the core set of functions necessary for the processing images
		using the SIVV method (as described in NISTIR 7599). 

********************************************************************************/

#ifndef BNIS_SIVV_CORE_H
#define BNIS_SIVV_CORE_H

#include <vector>
#include <string>
#include <opencv2/core/types_c.h>

using std::vector;
using std::string;

/* Peak-Valley-Pair (extrenum) data structure */
struct extrenum {
	double min, max;
	int min_loc, max_loc;
};

/*******************************************************************************
FUNCTION NAME:	dump_image()

AUTHOR:			Joseph C. Konczal

DESCRIPTION:	Print image characteristics on stdout and create an
				image window displaying the image.

	INPUT:
		name			- A unique name used to label the image information
						printed on stdout and the window displaying the image
		img				- A pointer to the image to display
	    autosize		- A flag for turning autosizing of window on and off
        verbose			- A flag for turning printf's on and off

	OUTPUT:
		name			- An image window containing the given image (img)

*******************************************************************************/
// Modern C++ version using cv::Mat
/*******************************************************************************
FUNCTION NAME:	sivv()

AUTHOR:			John D. Grantham

DESCRIPTION:	Performs the entire standard SIVV process, described in NISTIR
				7599, on a given image, using either default (see overload 
				below) or given parameters


	INPUT:
		img				- An image to be processed by SIVV
		window			- A flag for turning on/off the windowing step in the
						SIVV process
		smoothscale		- A flag for setting the number of points to be used
						in the signal smoothing algorithm (see the
						smooth_sums() function definition above)
		verbose			- A flag for turning on/off "verbose" mode, useful for
						debugging the SIVV process
		textonly		- A flag for turning on/off the "textonly mode", which 
						determines whether the output is displayed graphically
						or as text only

	OUTPUT:
		return			- A string containing the results of the SIVV process,
						separated by commas in the following order:
							(1) Number of peak-valley-pairs (PVP's) found
							(2) Ordinal number of largest PVP found
							(3) Power difference between the valley and peak
							(4) Frequency difference between the valley and peak
							(5) Slope between the valley and peak
							(6) Frequency of the midpoint between the valley
								and peak


*******************************************************************************/

/* SIVV Function Overload for quick use with default values */
string sivv(const cv::Mat &src);

/*******************************************************************************
FUNCTION NAME:	caseInsensitiveStringCompare

AUTHOR:			John D. Grantham

DESCRIPTION:	Provides a case-insensitive comparison of strings


	INPUT:
		str1			- The first string to be compared
		str2			- The second string to be compared

	OUTPUT:
		return			- A boolean value indicating whether or not the strings
						are equal (returns true if strings are equal)

*******************************************************************************/
bool caseInsensitiveStringCompare(const std::string& str1, const std::string& str2);

cv::Point2i find_fingerprint_center_morph(const cv::Mat& src, int* xbound_min, int* xbound_max, int* ybound_min, int* ybound_max);

#endif /* !BNIS_SIVV_CORE_H */
