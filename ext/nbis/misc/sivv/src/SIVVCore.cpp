/*******************************************************************************

License: 
This software was developed at the National Institute of Standards and 
Technology (NIST) by employees of the Federal Government in the course 
of their official duties. Pursuant to title 17 Section 105 of the 
United States Code, this software is not subject to copyright protection 
and is in the public domain. NIST assumes no responsibility  whatsoever for 
its use by other parties, and makes no guarantees, expressed or implied, 
about its quality, reliability, or any other characteristic. 

This software has been determined to be outside the scope of the EAR
(see Part 734.3 of the EAR for exact details) as it has been created solely
by employees of the U.S. Government; it is freely distributed with no
licensing requirements; and it is considered public domain.� Therefore,
it is permissible to distribute this software as a free download from the
internet.

Disclaimer: 
This software was developed to promote biometric standards and biometric
technology testing for the Federal Government in accordance with the USA
PATRIOT Act and the Enhanced Border Security and Visa Entry Reform Act.
Specific hardware and software products identified in this software were used
in order to perform the software development.  In no case does such
identification imply recommendation or endorsement by the National Institute
of Standards and Technology, nor does it imply that the products and equipment
identified are necessarily the best available for the purpose.  

*******************************************************************************/

/*******************************************************************************
	LIBRARY:	SIVVCore - Spectral Image Validation/Verification (Core)

	FILE:		SIVVCORE.CPP

	AUTHORS:	Joseph C. Konczal 
				John D. Grantham

				Modified by Jiong Zhang (Seventh Sense AI) for C++ compatibility
	
	DATE:		01/05/2009 (JCK)
	UPDATED:	01/19/2009 (JDG)
				09/10/2009 (JDG)
				01/25/2010 (JDG)
				09/03/2010 (JDG)
				07/07/2011 (JDG)
				07/25/2025 (JZ)

	DESCRIPTION:

		Contains the core set of functions necessary for the processing images
		using the SIVV method (as described in NISTIR 7599). 
********************************************************************************/

/* Win32 includes */
#ifdef WIN32
/* Intentionally blank -- a placeholder for any Win32-specific includes */

/* Linux includes */
#else 
/* Intentionally blank -- a placeholder for any Linux-specific includes */

#endif

/* C++ Includes */
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;

/* OpenCV includes */
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/opencv.hpp>

/* SIVV Includes */
#include "SIVVCore.h"


/*******************************************************************************
 * FUNCTION NAME:	init_blackman_1d_filter_row()
 * @param arr the Mat object to be initialized (row)
 * @param alpha the alpha parameter for the Blackman filter
 */
void init_blackman_1d_filter_row(cv::Mat& arr, const double alpha)
{
    const double a0 = (1.0 - alpha) / 2.0;
    const double a1 = 1.0/2.0;
    const double a2 = alpha / 2.0;
    const double f = 2.0 * CV_PI / (arr.cols - 1.0);
    
    if (arr.rows != 1) {
        fprintf(stderr, "error: height != 1\n");
        exit(EXIT_FAILURE);
    }
    
    auto* dptr = arr.ptr<double>(0);
    for (int n = 0; n < arr.cols; n++) {
        dptr[n] = a0 - a1 * cos(f * (double)n) + a2 * cos(2.0 * f * (double)n);
    }
}

/*******************************************************************************
 * FUNCTION NAME:	init_blackman_1d_filter_col()
 * @param arr the Mat object to be initialized (column)
 * @param alpha the alpha parameter for the Blackman filter
 */
void init_blackman_1d_filter_col(cv::Mat& arr, const double alpha)
{
    const double a0 = (1.0 - alpha) / 2.0;
    const double a1 = 1.0/2.0;
    const double a2 = alpha / 2.0;
    const double f = 2.0 * CV_PI / (arr.rows - 1.0);

    if (arr.cols != 1) {
        fprintf(stderr, "error: width != 1\n");
        exit(EXIT_FAILURE);
    }

    for (int n = 0; n < arr.rows; n++) {
        auto* dptr = arr.ptr<double>(n);
        dptr[0] = a0 - a1 * cos(f * (double)n) + a2 * cos(2.0 * f * (double)n);
    }
}

/*******************************************************************************
 * FUNCTION NAME:	init_tukey_1d_filter_row()
 * @param arr the Mat object to be initialized (row)
 * @param alpha the alpha parameter for the Tukey filter
 */
void init_tukey_1d_filter_row(cv::Mat& arr, const double alpha)
{
    const double a1 = 0.5;
    const double a2 = (2.0 * CV_PI / alpha);
    const double a3 = (alpha / 2.0);
    
    if (arr.rows != 1) {
        fprintf(stderr, "error: height != 1\n");
        exit(EXIT_FAILURE);
    }
    
    auto* dptr = arr.ptr<double>(0);
    for (int n = 0; n < arr.cols; n++) {
        double x = (double)n / (double)arr.cols;
        
        if (x < a3) {
            dptr[n] = a1 * (1.0 + cos(a2 * (x - a3)));
        }
        else if (x < (1.0 - a3)) {
            dptr[n] = 1.0;
        }
        else if (x <= 1.0) {
            dptr[n] = a1 * (1.0 + cos(a2 * (x - 1.0 + a3)));
        }
    }
}

/*******************************************************************************
 * FUNCTION NAME:	init_tukey_1d_filter_col()
 * @param arr the Mat object to be initialized (col)
 * @param alpha the alpha parameter for the Tukey filter
 */
void init_tukey_1d_filter_col(cv::Mat& arr, const double alpha)
{
    const double a1 = 0.5;
    const double a2 = (2.0 * CV_PI / alpha);
    const double a3 = (alpha / 2.0);

    if (arr.cols != 1) {
        fprintf(stderr, "error: width != 1\n");
        exit(EXIT_FAILURE);
    }

    for (int n = 0; n < arr.rows; n++) {
        auto* dptr = arr.ptr<double>(n);
        double x = (double)n / (double)arr.rows;

        if (x < a3) {
            dptr[0] = a1 * (1.0 + cos(a2 * (x - a3)));
        }
        else if (x < (1.0 - a3)) {
            dptr[0] = 1.0;
        }
        else if (x <= 1.0) {
            dptr[0] = a1 * (1.0 + cos(a2 * (x - 1.0 + a3)));
        }
    }
}

/*******************************************************************************
FUNCTION NAME:	apply_tukey_window()

AUTHOR:			John D. Grantham. Modified for C++ by Jiong Zhang (Seventh Sense AI)

DESCRIPTION:	Create an appropriately sized 2D Tukey window, and apply 
				it to the image.

	INPUT:
		src				- Points to an allocated image structure containing the
						input image.

	OUTPUT:
		dst				- Points to allocated image structure of the same size 
						and type as the input (src) to receive the filtered 
						result. Src and dst may be equal.

*******************************************************************************/
void apply_tukey_window(const cv::Mat& src, cv::Mat& dst)
{
    cv::Mat hor = cv::Mat::zeros(1, src.cols, CV_64F);
    cv::Mat ver = cv::Mat::zeros(src.rows, 1, CV_64F);

    init_tukey_1d_filter_row(hor, 0.25);
    init_tukey_1d_filter_col(ver, 0.25);

    cv::Mat tka;
    cv::gemm(ver, hor, 1.0, cv::Mat(), 0.0, tka);

    cv::multiply(src, tka, dst);
}

/*******************************************************************************
FUNCTION NAME:	apply_blackman_window()

AUTHOR:			Joseph C. Konczal. Modified for C++ by Jiong Zhang (Seventh Sense AI)

DESCRIPTION:	Create an appropriately sized 2D Blackman window, and apply 
				it to the image.

	INPUT:
		src				- Points to an allocated image structure containing the
						input image.

	OUTPUT:
		dst				- Points to allocated image structure of the same size 
						and type as the input (src) to receive the filtered 
						result. Src and dst may be equal.

*******************************************************************************/
void apply_blackman_window(const cv::Mat& src, cv::Mat& dst)
{
    cv::Mat hor = cv::Mat::zeros(1, src.cols, CV_64F);
    cv::Mat ver = cv::Mat::zeros(src.rows, 1, CV_64F);

    init_blackman_1d_filter_row(hor, 0.16);
    init_blackman_1d_filter_col(ver, 0.16);

    cv::Mat bwa;
    cv::gemm(ver, hor, 1.0, cv::Mat(), 0.0, bwa);

    cv::multiply(src, bwa, dst);
}

/*******************************************************************************
FUNCTION NAME:	polar_transform()

AUTHOR:			John D. Grantham. Modified for C++ by Jiong Zhang (Seventh Sense AI)

DESCRIPTION:	Performs a polar transformation to the input array (src) and
				writes the result into the output array (dst). Both arrays must 
				be two dimensional and the values must be 32-bit floating point.

	INPUT:
		src				- Points to an array structure containing the 2D, 
						32-bit floating point input data.
		flags			- OpenCV Interpolation flag (ex: CV_INTER_CUBIC), 
						determines the type of interpolation used in remapping 

	OUTPUT:
		dst				- Points to allocated array structure of the same size
						and type as the input (src) to receive the transformed 
						result. Src and dst may not be equal.

*******************************************************************************/
void polar_transform(const cv::Mat& src, cv::Mat& dst, int flags)
{
    if (&dst == &src) {
        printf("ERROR: polar_transform input arrays cannot be equal!\n");
        exit(EXIT_FAILURE);      
    }
    
    cv::Mat mapx = cv::Mat::zeros(src.size(), CV_32F);
    cv::Mat mapy = cv::Mat::zeros(src.size(), CV_32F);
    
    int width = src.cols;
    int height = src.rows;
    
    int x1 = 0, x2 = width;
    int y1 = 0, y2 = height;
    int xdiff = x2 - x1;
    int ydiff = y2 - y1;
    
    double xm = xdiff / 2.0;
    double ym = ydiff / 2.0;
    double angl = 0;
    long circle = 100;
    
    for (int row = y1; row < y2; row++) {
        auto* mapx_ptr = mapx.ptr<float>(row);
        auto* mapy_ptr = mapy.ptr<float>(row);
        
        for (int col = x1; col < x2; col++) {
            double phi = (2 * CV_PI) * (col - x1) / xdiff;
            phi = fmod(phi + angl, 2 * CV_PI);
            
            double phi2;
            if (phi >= 1.5 * CV_PI)
                phi2 = 2 * CV_PI - phi;
            else if (phi >= CV_PI)
                phi2 = phi - CV_PI;
            else if (phi >= 0.5 * CV_PI)
                phi2 = CV_PI - phi;
            else
                phi2 = phi;
            
            double xx = tan(phi2);
            double m = (xx != 0) ? (1.0 / xx) : 0;
            
            double xmax, ymax;
            if (m <= ((double)(ydiff) / (double)(xdiff))) {
                if (phi2 == 0) {
                    xmax = 0;
                    ymax = ym - y1;
                } else {
                    xmax = xm - x1;
                    ymax = m * xmax;
                }
            } else {
                ymax = ym - y1;
                xmax = ymax / m;
            }
            
            double rmax = sqrt((xmax*xmax) + (ymax*ymax));
            double t = ((ym - y1) < (xm - x1)) ? (ym - y1) : (xm - x1);
            rmax = (rmax - t) / 100.0 * (100 - circle) + t;
            
            double r = rmax * (double)((y2 - row) / (double)(ydiff));
            xx = r * sin(phi2);
            double yy = r * cos(phi2);
            
            double x, y;
            if (phi >= 1.5 * CV_PI) {
                x = xm - xx;
                y = ym - yy;
            } else if (phi >= CV_PI) {
                x = xm - xx;
                y = ym + yy;
            } else if (phi >= 0.5 * CV_PI) {
                x = xm + xx;
                y = ym + yy;
            } else {
                x = xm + xx;
                y = ym - yy;
            }
            
            mapx_ptr[col] = static_cast<float>(x);
            mapy_ptr[col] = static_cast<float>(y);
        }
    }
    
    cv::remap(src, dst, mapx, mapy, flags, cv::BORDER_CONSTANT, cv::Scalar::all(0));
}

/*******************************************************************************
FUNCTION NAME:	diagonal_shuffle_quadrants()

AUTHOR:			Joseph C. Konczal. Modified for C++ by Jiong Zhang (Seventh Sense AI)

DESCRIPTION:	Rearrange the quadrants of a 2D array by swapping them 
				diagonally.

	INPUT:
		src				- Points to an array structure containing the input data

	OUTPUT:
		dst				- Points to allocated structure of the same size and
						type as the input (src) to receive the shuffled result.
						Src and dst may be equal.

*******************************************************************************/
void diagonal_shuffle_quadrants(const cv::Mat& src, cv::Mat& dst)
{
    int width = src.cols, height = src.rows;
    int wx = width/2, hy = height/2;
    
    if (dst.cols != width || dst.rows != height) {
        fprintf(stderr, "Unmatched sizes: src %d x %d != dst %d x %d\n", 
                src.cols, src.rows, dst.cols, dst.rows);
        exit(EXIT_FAILURE);
    }
    
    cv::Mat srcq[4], dstq[4];
    
    // Define quadrants
    srcq[0] = src(cv::Rect(0, 0, wx, hy));      // Top-left
    srcq[1] = src(cv::Rect(wx, 0, wx, hy));     // Top-right
    srcq[2] = src(cv::Rect(wx, hy, wx, hy));    // Bottom-right
    srcq[3] = src(cv::Rect(0, hy, wx, hy));     // Bottom-left
    
    if (&src == &dst) {
        cv::Mat tmp = cv::Mat::zeros(hy, wx, src.type());
        
        // Swap quadrants diagonally
        srcq[0].copyTo(tmp);
        srcq[2].copyTo(srcq[0]);
        tmp.copyTo(srcq[2]);
        
        srcq[1].copyTo(tmp);
        srcq[3].copyTo(srcq[1]);
        tmp.copyTo(srcq[3]);
    } else {
        dstq[0] = dst(cv::Rect(0, 0, wx, hy));
        dstq[1] = dst(cv::Rect(wx, 0, wx, hy));
        dstq[2] = dst(cv::Rect(wx, hy, wx, hy));
        dstq[3] = dst(cv::Rect(0, hy, wx, hy));
        
        srcq[0].copyTo(dstq[2]);
        srcq[1].copyTo(dstq[3]);
        srcq[2].copyTo(dstq[0]);
        srcq[3].copyTo(dstq[1]);
    }
}

/*******************************************************************************
FUNCTION NAME:	log_power_spectrum()

AUTHORS:		Joseph C. Konczal and John D. Grantham,  Modified for C++ by Jiong Zhang (Seventh Sense AI)

DESCRIPTION:	Compute the log power spectrum of the DFT image. 

	INPUT:
		src				- Points to an array structure containing the input data

	OUTPUT:
		dst		        - Points to the output
*******************************************************************************/
void log_power_spectrum(const cv::Mat& src, cv::Mat& dst)
{
    // Compute log polar power spectrum
    cv::Mat dft_comb = cv::Mat::zeros(src.size(), CV_64FC2);
    cv::Mat dft_real = cv::Mat::zeros(src.size(), CV_64F);
    cv::Mat dft_dpy = cv::Mat::zeros(src.size(), CV_64FC3);
    cv::Mat dft_imgy = cv::Mat::zeros(src.size(), CV_64F);
    
    // Prepare complex image for DFT
    vector<cv::Mat> planes = {src, dft_imgy};
    cv::merge(planes, dft_comb);
    
    // Perform DFT
    cv::dft(dft_comb, dft_comb, cv::DFT_COMPLEX_OUTPUT);
    
    // Split real and imaginary parts
    cv::split(dft_comb, planes);
    dft_real = planes[0];
    dft_imgy = planes[1];
    
    // Create display image
    cv::Mat dft_zero = cv::Mat::zeros(src.size(), CV_64F);
    vector<cv::Mat> display_planes = {dft_zero, dft_real, dft_imgy};
    cv::merge(display_planes, dft_dpy);
    
    diagonal_shuffle_quadrants(dft_dpy, dft_dpy);
    
    cv::split(dft_dpy, display_planes);
    dft_real = display_planes[1];
    dft_imgy = display_planes[2];
    
    // Compute magnitude spectrum
    cv::pow(dft_real, 2.0, dft_real);
    cv::pow(dft_imgy, 2.0, dft_imgy);
    cv::add(dft_real, dft_imgy, dft_real);
    cv::sqrt(dft_real, dft_real);
    
    // Compute log power spectrum
    cv::add(dft_real, cv::Scalar::all(1.0), dft_imgy);
    cv::log(dft_imgy, dft_imgy);
    
    // Convert to output
    dft_imgy.convertTo(dst, dst.type());
}

/*******************************************************************************
FUNCTION NAME:	findmax()

AUTHOR:			John D. Grantham, Modified for C++ by Jiong Zhang (Seventh Sense AI)

DESCRIPTION:	Finds and returns the maximum value of an image. 

	INPUT:
		src				- Points to an array structure containing the input data

	OUTPUT:
		return  		- A double containing the maximum value found in the 
						input image.

*******************************************************************************/
double findmax(const cv::Mat& src)
{
    double min_val, max_val;
    cv::minMaxLoc(src, &min_val, &max_val);
    return max_val;
}

/*******************************************************************************
FUNCTION NAME:	sum_rows()

AUTHOR:			John D. Grantham, Modified for C++ by Jiong Zhang (Seventh Sense AI)

DESCRIPTION:	Sums each of the rows in a given image and stores the sums in a
				given vector of doubles. 

	INPUT:
		src				- Points to an array structure containing the input data

	OUTPUT:
		rowsums  		- A vector of doubles containing the sum of each row of
						the input image, stored in top-to-bottom order

*******************************************************************************/
void sum_rows(const cv::Mat& src, vector<double>& rowsums)
{
    rowsums.resize(src.rows);
    
    for (int i = 0; i < src.rows; i++) {
        cv::Scalar sum = cv::sum(src.row(i));
        rowsums[i] = sum[0];
    }
}

/*******************************************************************************
FUNCTION NAME:	normalize_sums()

AUTHOR:			John D. Grantham

DESCRIPTION:	Normalizes the sums stored in a vector to the 0th term 

	INPUT:
		sums			- A vector of doubles containing sums to be normalized

	OUTPUT:
		sums  			- A vector of doubles, normalized to the 0th value

*******************************************************************************/
void normalize_sums(vector<double>& sums)
{
    if (sums.empty()) return;
    
    double normterm = sums[0];
    for (auto& sum : sums) {
        sum = sum / normterm;
    }
}

/*******************************************************************************
FUNCTION NAME:	smooth_sums()

AUTHOR:			John D. Grantham,  Modified for C++ by Jiong Zhang (Seventh Sense AI)

DESCRIPTION:	Smooths the signal using an n-point moving average filter, where 
				n is an odd number greater than 1 (if given an even number, the
				function will choose an odd number 1 lower than the given value).
				The filtering algorithm performs zero-phase digital filtering by
				processing the input data in both forward and reverse directions
				as described in NISTIR 7599.


	INPUT:
		rowsums			- A vector of doubles containing sums which represent 
						the spectrum/signal to be smoothed
		numpoints		- The number of points over which to the signal will be
						smoothed (higher values mean more smoothing)

	OUTPUT:
		rowsums  			- A vector of doubles, smoothed by n points

*******************************************************************************/
void smooth_sums(vector<double> &rowsums, int num_points)	
{
    int num_rows = static_cast<int>(rowsums.size());
    vector<double> firstpass(num_rows), firstreverse(num_rows), secondpass(num_rows);
    int distance = (num_points % 2 == 0) ? num_points / 2 : (num_points - 1) / 2;
    int maxindex = num_rows - 1;

    // First pass
    for (int i = 0; i <= maxindex; i++) {
        double smoothsum = 0;
        
        if (i < distance) {
            for (int j = i+1; j < (i + distance); j++)
                smoothsum += rowsums[j]; 
            for (int j = 0; j < i; j++)
                smoothsum += rowsums[j];
            for (int j = i+1; j <= distance; j++)
                smoothsum += rowsums[j];
        } else if (i > (maxindex - distance)) {
            for (int j = i-1; j > (i - distance); j--)
                smoothsum += rowsums[j]; 
            for (int j = maxindex; j > i; j--)
                smoothsum += rowsums[j];
            for (int j = i-1; j >= maxindex - distance; j--)
                smoothsum += rowsums[j];
        } else {
            for (int j = (i - distance); j < (i + distance); j++) {
                if (j != i)
                    smoothsum += rowsums[j];
            }
        }
        
        double smoothavg = smoothsum / (distance * 2);
        firstpass[i] = smoothavg;
    }

    // Reverse signal for second pass
    for (int i = 0; i <= maxindex; i++) {
        firstreverse[i] = firstpass[maxindex-i];
    }

    // Second pass (similar logic)
    for (int i = 0; i <= maxindex; i++) {
        double smoothsum = 0;
        
        if (i < distance) {
            for (int j = i+1; j < (i + distance); j++)
                smoothsum += firstreverse[j]; 
            for (int j = 0; j < i; j++)
                smoothsum += firstreverse[j];
            for (int j = i+1; j <= distance; j++)
                smoothsum += firstreverse[j];
        } else if (i > (maxindex - distance)) {
            for (int j = i-1; j > (i - distance); j--)
                smoothsum += firstreverse[j]; 
            for (int j = maxindex; j > i; j--)
                smoothsum += firstreverse[j];
            for (int j = i-1; j >= maxindex - distance; j--)
                smoothsum += firstreverse[j];
        } else {
            for (int j = (i - distance); j < (i + distance); j++) {
                if (j != i)
                    smoothsum += firstreverse[j];
            }
        }
        
        double smoothavg = smoothsum / (distance * 2);
        secondpass[i] = smoothavg;
    }

    // Write results back (undo reversal)
    for (int i = 0; i <= maxindex; i++) {
        rowsums[i] = secondpass[maxindex-i];
    }
}

/*******************************************************************************
FUNCTION NAME:	find_global_minmax()

AUTHOR:			John D. Grantham

DESCRIPTION:	Finds the global minimum and maximum of a given vector of 
				doubles, and stores the values (along with their locations) in 
				a given extrenum structure.


	INPUT:
		rowsums			- A vector of doubles containing sums 
		global_minmax	- The extrenum structure in which to store the global
						minimum and maximum points found in the vector of sums

	OUTPUT:
		global_minmax	- An extrenum structure containing the global minimum
						and maximum values, along with their locations (points)

*******************************************************************************/
void find_global_minmax(const vector<double>& rowsums, extrenum& global_minmax)
{
    if (rowsums.empty()) return;
    
    global_minmax.min = rowsums[0];
    global_minmax.max = rowsums[0];
    global_minmax.min_loc = 0;
    global_minmax.max_loc = 0;

    for (int i = 0; i < (int)rowsums.size(); i++) {
        if (rowsums[i] <= global_minmax.min) {
            global_minmax.min = rowsums[i];
            global_minmax.min_loc = i;
        } else if (rowsums[i] >= global_minmax.max) {
            global_minmax.max = rowsums[i];
            global_minmax.max_loc = i;
        }
    }
}

/*******************************************************************************
FUNCTION NAME:	cvp_distance()

AUTHOR:			John D. Grantham

DESCRIPTION:	Calculates and returns the distance between two points


	INPUT:
		a				- A CvPoint (origin)
		b				- Another CvPoint (destination)

	OUTPUT:
		return			- The distance between the two given CvPoints (from
						origin to destination)

*******************************************************************************/
double cvp_distance(const CvPoint a, const CvPoint b)
{
	return sqrt((double)((b.x - a.x)*(b.x - a.x)) + ((b.y - a.y)*(b.y - a.y)));
}

/*******************************************************************************
FUNCTION NAME:	peak_finder()

AUTHOR:			John D. Grantham,  Modified for C++ by Jiong Zhang (Seventh Sense AI)

DESCRIPTION:	Finds peak-valley-pairs (PVP's) within a given signal above the
				given peak size threshold and within a specified distance
				between the	valley and peak


	INPUT:
		peaks			- A vector in which to store any PVP's found
		rowsums			- A vector containing the signal to be processed
		global_minmax	- An extrenum structure containing the global minimum
						and maximum values of the signal
		peak_threshold	- The given peak size threshold (the minimum relative
						size of peaks to be found in the signal)
		step			- The maximum distance allowable between a peak and
						its preceeding valley

	OUTPUT:
		peaks			- A vector containing any PVP's found (within step and
						above peak_threshold)

*******************************************************************************/
void peak_finder(vector<extrenum>& peaks, const vector<double>& rowsums, const extrenum& global_minmax, double peak_threshold, int step)
{
    double min = global_minmax.min;
    double max = global_minmax.max;
    int num_rows = static_cast<int>(rowsums.size());

    int lmin_loc = 0;
    int lmax_loc = 0;
    extrenum pvp{};

    for (int i = 0; i < (num_rows - step); i++) {
        double lmin = max;
        double lmax = min;

        if (rowsums[i] < lmin && i > lmax_loc)
        {
            lmin = rowsums[i];
            lmin_loc = i;

            for (int j = i; j < (i + step); j++)
            {
                if (rowsums[j] > lmax)
                {
                    lmax = rowsums[j];
                    lmax_loc = j;
                }
            }

            if (lmax != min)
            {
                for (int j = i; j < (i+step); j++)
                {
                    if (rowsums[j] < lmin && j < lmax_loc)
                    {
                        lmin = rowsums[j];
                        lmin_loc = j;
                    }
                }
            }
        }

        if ((lmin != max && lmax != min) && (lmax - lmin) > peak_threshold)
        {
            /* Add peak to vector of PVP's */
            pvp.min_loc = lmin_loc;
            pvp.max_loc = lmax_loc;
            pvp.min = lmin;
            pvp.max = lmax;
            peaks.push_back(pvp);
        }
    }
}

/*******************************************************************************
FUNCTION NAME:	sivv()

AUTHOR:			John D. Grantham, Modified for C++ by Jiong Zhang (Seventh Sense AI)

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
string sivv(const cv::Mat& src, int smoothscale, int verbose, int textonly, vector<double>* signal, string window, const string& graphfile, int* fail)
{
    cv::Mat img, img_dfp, img_polar, polar_trans;
    double grey_scale_factor;
    string results;

    // Window selection
    if (caseInsensitiveStringCompare(window, "blackman")) {
        window = "blackman";
    } else if (caseInsensitiveStringCompare(window, "tukey")) {
        window = "tukey";
    } else {
        window = "";
    }

    if (smoothscale < 1) smoothscale = 1;
    if (verbose != 0) verbose = 1;
    if (textonly != 0) textonly = 1;

    // Convert to grayscale if needed
    if (src.channels() > 1) {
        cv::cvtColor(src, img, cv::COLOR_BGR2GRAY);
    } else {
        img = src.clone();
    }

    // Convert to double precision floating point [0.0, 1.0]
    grey_scale_factor = 1.0 / 255.0;
    img.convertTo(img_dfp, CV_64F, grey_scale_factor);

    // Apply window function
    if (window == "blackman") {
        apply_blackman_window(img_dfp, img_dfp);
    } else if (window == "tukey") {
        apply_tukey_window(img_dfp, img_dfp);
    }

    // Compute log polar power spectrum
    cv::Mat img_lps = cv::Mat::zeros(img_dfp.size(), CV_64F);

    log_power_spectrum(img_dfp, img_lps);

	if (textonly == 0)
	{
		//dump_image("Log Power Spectrum", img_lps, 0, verbose);
	}

    img_lps.convertTo(img_polar, CV_32F, (1.0 / findmax(img_lps)));

    // Polar transform using bicubic interpolation
    polar_trans = cv::Mat::zeros(img_polar.size(), CV_32F);
    polar_transform(img_polar, polar_trans, cv::INTER_CUBIC);
	if (textonly == 0)
	{
		//dump_image("Polar Transform", polar_trans, false, verbose);
	}

    // Reduce LogPolar to angles 0 - 180
    cv::Mat polar_trans_half = polar_trans(cv::Rect(0, 0, polar_trans.cols / 2, polar_trans.rows));

    // Sum rows of polar transform (0 - 180)
    vector<double> rowsums;
    cv::Mat flipped;
    cv::flip(polar_trans_half, flipped, 0);

    sum_rows(flipped, rowsums);

    // Smooth sums
    if (smoothscale > 1)
        smooth_sums(rowsums, smoothscale); 

    // Normalize sums to DC (0th) term
    normalize_sums(rowsums);

    // Find global min and max
    extrenum global_minmax{};
    find_global_minmax(rowsums, global_minmax);

    if (verbose != 0)
        printf("Global maximum: %0.10f, Global minimum: %0.10f\nGlobal max location: %d, Global min location: %d\n", 
               global_minmax.max, global_minmax.min, global_minmax.max_loc, global_minmax.min_loc);

    // Find peaks
    vector<extrenum> peaks;
    int num_rows = static_cast<int>(rowsums.size());
    peak_finder(peaks, rowsums, global_minmax, 0.005, (num_rows / 5));

    double cpp_scale = 1 / (2 * (double)num_rows);

    // Process peaks
    double maxdiff = 0.0;
    extrenum largestpvp{};
    int lpvp_peaknum = 0;

    for (int i = 0; i < (int)peaks.size(); i++) {
        double diff = peaks[i].max - peaks[i].min;
        if (diff > maxdiff) {
            maxdiff = diff;
            largestpvp = peaks[i];
            lpvp_peaknum = i+1;
        }
    }

    if (verbose != 0)
        printf("Peak value: %10f, Valley value: %10f\nPeak location (freq): %f, Valley location (freq): %f\n", 
               largestpvp.max, largestpvp.min, (largestpvp.max_loc * cpp_scale), (largestpvp.min_loc * cpp_scale));

    // Calculate peak statistics
    double dx = 0, dy = 0, slope = 0, center_freq = 0, peak_freq = 0;
    int num_peaks = 0;

    if (!peaks.empty()) {
        num_peaks = static_cast<int>(peaks.size());
        dx = (largestpvp.max_loc - largestpvp.min_loc) * cpp_scale;
        dy = (largestpvp.max - largestpvp.min);
        center_freq = ((largestpvp.min_loc * cpp_scale) + (dx / 2));
        peak_freq = largestpvp.max_loc * cpp_scale;
        
        slope = (dy == 0.0 || dx == 0.0) ? 0.0 : dy / dx;
    }

    // Set fail flag
    if (fail != nullptr) {
        *fail = (dy == 0) ? 1 : 0;
    }

    stringstream out;
    out << lpvp_peaknum << ", " << num_peaks << ", " << dy << ", " << dx << ", " 
        << slope << ", " << center_freq << ", " << peak_freq;
    results = out.str();

    if (signal != nullptr) {
        *signal = rowsums;
    }

    return results;
}

// Overloaded function for quick use with default values
string sivv(const cv::Mat& src)
{
    int fail = 0;
//    string results = sivv(src, 7, 0, 0, nullptr, "blackman", "black_man_graph.png", &fail);
    string results = sivv(src, 7, 0, 1, nullptr, "blackman", "", &fail);

    if (fail == 1) {
        results = sivv(src, 7, 0, 1, nullptr, "tukey", "", &fail);
    }

    return results;
}

/*******************************************************************************
FUNCTION NAME:	find_fingerprint_center_morph()

AUTHOR:			John D. Grantham, Modified for C++ by Jiong Zhang (Seventh Sense AI)
				(with credit to Bruce Bandini for his morphology code)

DESCRIPTION:	Finds the center of the area of the center of a potential 
				fingerprint, based on several morphology operations


	INPUT:
		src				- An image to be processed
		xbound_min		- A pointer to a location in which to store the
						minimum x-value bounding the area of highest edge
						density (if found)
		xbound_max		- A pointer to a location in which to store the
						maximum x-value bounding the area of highest edge
						density (if found)
		ybound_min		- A pointer to a location in which to store the
						minimum y-value bounding the area of highest edge
						density (if found)
		ybound_max		- A pointer to a location in which to store the
						maximum y-value bounding the area of highest edge
						density (if found)

	OUTPUT:
		xbound_min		- The minimum x-value bounding the area of highest 
						edge density (if found)
		xbound_max		- The maximum x-value bounding the area of highest 
						edge density (if found)
		ybound_min		- The minimum y-value bounding the area of highest 
						edge density (if found)
		ybound_max		- The maximum y-value bounding the area of highest 
						edge density (if found)

*******************************************************************************/
cv::Point2i find_fingerprint_center_morph(const cv::Mat& src, int* xbound_min, int* xbound_max, int* ybound_min, int* ybound_max)
{
    cv::Mat inputImg, topHatImg, blackWhiteImg, closedOutputImg;

    // Convert to grayscale if needed
    if (src.channels() > 1) {
        cv::cvtColor(src, inputImg, cv::COLOR_BGR2GRAY);
    } else {
        inputImg = src.clone();
    }

    // Create structuring element
    cv::Size sz(14, 14);
    cv::Mat strEl = cv::getStructuringElement(cv::MORPH_ELLIPSE, sz);

    // Top Hat Filter
    cv::morphologyEx(inputImg, topHatImg, cv::MORPH_TOPHAT, strEl);

    // Convert to binary using OTSU threshold
    cv::threshold(topHatImg, blackWhiteImg, 128, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);

    // Morphological closing
    cv::morphologyEx(blackWhiteImg, closedOutputImg, cv::MORPH_CLOSE, strEl);

    // Apply mask and invert image values
    for (int v = 0; v < closedOutputImg.total(); v++) {
        auto* closed_ptr = closedOutputImg.ptr<uchar>();
        auto* bw_ptr = blackWhiteImg.ptr<uchar>();

        if (closed_ptr[v] != 0) {
            bw_ptr[v] = abs(bw_ptr[v] - 255);
        } else {
            bw_ptr[v] = 0;
        }
    }

    // Calculate moments to find center
    cv::Moments moments = cv::moments(closedOutputImg, true);

    int posX = static_cast<int>(moments.m10 / moments.m00);
    int posY = static_cast<int>(moments.m01 / moments.m00);

    // Set boundaries
    *xbound_min = std::max(0, posX - 200);
    *xbound_max = std::min(src.cols, posX + 200);
    *ybound_min = std::max(0, posY - 200);
    *ybound_max = std::min(src.rows, posY + 200);

    return {posX, posY};
}

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
bool caseInsensitiveStringCompare(const string& str1, const string& str2)
{
    string str1Cpy(str1);
    string str2Cpy(str2);

    transform(str1Cpy.begin(), str1Cpy.end(), str1Cpy.begin(), ::tolower);
    transform(str2Cpy.begin(), str2Cpy.end(), str2Cpy.begin(), ::tolower);

    return (str1Cpy == str2Cpy);
}
