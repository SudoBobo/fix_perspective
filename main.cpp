#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
//#include <tclDecls.h>

using namespace cv;
//using namespace std;

Mat source; Mat src_gray;
int thresh = 100;
int max_thresh = 255;
RNG rng(12345);

// Finds the intersection of two lines, or returns false.
// The lines are defined by (o1, p1) and (o2, p2).
Point2f* intersection(Point2f o1, Point2f p1, Point2f o2, Point2f p2)
{
	Point2f x = o2 - o1;
	Point2f d1 = p1 - o1;
	Point2f d2 = p2 - o2;

	float cross = d1.x*d2.y - d1.y*d2.x;
	if (abs(cross) < /*EPS*/1e-8)
		return NULL;

	double t1 = (x.x * d2.y - x.y * d2.x)/cross;
	Point2f r = o1 + d1 * t1;
	return new Point2f(r);
}

double distanse(Point2f* p1, Point2f* p2) {
	return std::sqrt(std::pow((p2->x - p1->x), 2) +
			 std::pow((p2->y - p1->y), 2));
}

bool same(double a, double b)
{
	return fabs(a - b) < 0.0001;
}

/** @function main */
int main( int argc, char** argv )
{
	/// Load source image and convert it to gray
	source = imread( argv[1], 1 );



	/// Create Window
	char* source_window = "Source";
	namedWindow( source_window, WINDOW_AUTOSIZE );
	imshow( source_window, source );

	Point center = Point(source.cols/2, source.rows/2);
	double angle = 180.0;
	double scale = 1;
	Mat rot_mat = getRotationMatrix2D( center, angle, scale );
	warpAffine(source, source, rot_mat, source.size() );

	source_window = "Source1";
	namedWindow( source_window, WINDOW_AUTOSIZE );
	imshow( source_window, source );

	/// Convert image to gray and blur it
	cvtColor( source, src_gray, COLOR_BGR2GRAY );
	blur( src_gray, src_gray, Size(3,3) );

	Mat canny_output;
	std::vector<std::vector<Point>> contours;
	std::vector<Vec4i> hierarchy;

	/// Detect edges using canny
	Canny( src_gray, canny_output, thresh, thresh*2, 3 );
	/// Find contours
	findContours( canny_output, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE, Point(0, 0) );

	int contour_id = 0;
	double max_area = 0;
	for (int i = 0; i < contours.size(); i++) {
		double contour_area = contourArea(contours[i]);
		if (contour_area > max_area) {
			max_area = contour_area;
			contour_id = i;
		}
	}

	Mat contour_drawing = Mat::zeros( canny_output.size(), CV_8UC3);
	Scalar color = Scalar( 255, 255, 255);
	drawContours(contour_drawing, contours, contour_id, color, 2, 8,
		     hierarchy, 0, Point() );

	Mat dst;
	Canny(contour_drawing, dst, 50, 200, 3);

	// get lines using Haff
	std::vector<Vec2f> lines; // will hold the results of the detection
	std::vector<std::pair<Point2f, Point2f>> vertical_lines;
	std::vector<std::pair<Point2f, Point2f>> horizontal_lines;
	HoughLines(dst, lines, 1, CV_PI/180, 150, 0, 0 ); // runs the actual detection
	// Draw the lines
	Mat cdst = Mat::zeros( canny_output.size(), CV_8UC3 );

	for( size_t i = 0; i < lines.size(); i++ )
	{
		float rho = lines[i][0], theta = lines[i][1];
		Point2f pt1, pt2;
		double a = cos(theta), b = sin(theta);
		double x0 = a*rho, y0 = b*rho;
		pt1.x = cvRound(x0 + 1000*(-b));
		pt1.y = cvRound(y0 + 1000*(a));
		pt2.x = cvRound(x0 - 1000*(-b));
		pt2.y = cvRound(y0 - 1000*(a));

		std::pair<Point2f, Point2f> new_pair(pt1, pt2);
		if (abs(pt1.x - pt2.x) < abs(pt1.y - pt2.y))
			vertical_lines.push_back(new_pair);
		else
			horizontal_lines.push_back(new_pair);
//		line( cdst, pt1, pt2, Scalar(255,255,255), 1, CV_AA);
	}

	// put all horizon
	// find corner dotes
	std::vector<Point2f> intersection_points;
	Point2f* new_point;
	Point2f mass_center(0, 0);
	// todo we assume all coordinates to be positive
	for (int i = 0; i < vertical_lines.size(); i++) {
		for (int j = 0; j < horizontal_lines.size(); j++) {
			Point2f vertical_1 = vertical_lines[i].first;
			Point2f vertical_2 = vertical_lines[i].second;

			Point2f horizontal_1 = horizontal_lines[j].first;
			Point2f horizontal_2 = horizontal_lines[j].second;

			new_point = intersection(vertical_1, vertical_2,
						 horizontal_1, horizontal_2);

//			circle(cdst, *new_point, 2, color);

			if (new_point != NULL) {
				intersection_points.push_back(*new_point);
				mass_center += *new_point;
			}
		}
	}

	mass_center /= (double) intersection_points.size();

	// todo check that points operations work as expected

	// todo before first launch perform an amount of visual checks

	// reduce corner point to amount of four
	Point2f left_top(0, 0);
	Point2f left_bottom(0, 0);
	Point2f right_top(0, 0);
	Point2f right_bottom(0, 0);

	int lt_size = 0;
	int lb_size = 0;
	int rt_size = 0;
	int rb_size = 0;

	double d[4];
	for (int i = 0; i < intersection_points.size(); i++) {
		Point2f* p = &intersection_points[i];
		if ((p->x > mass_center.x) && (p->y > mass_center.y)) {
			rt_size++;
			right_top += *p;
		}

		if ((p->x > mass_center.x) && (p->y < mass_center.y)) {
			rb_size++;
			right_bottom += *p;
		}

		if ((p->x < mass_center.x) && (p->y > mass_center.y)) {
			lt_size++;
			left_top += *p;
		}

		if ((p->x < mass_center.x) && (p->y < mass_center.y)) {
			lb_size++;
			left_bottom += *p;
		}
	}

	left_top /= lt_size;
	left_bottom /= lb_size;
	right_top /= rt_size;
	right_bottom /= rb_size;

	std::cout << "lt " << left_top << std::endl;
	std::cout << "lb " << left_bottom << std::endl;
	std::cout << "rt " << right_top << std::endl;
	std::cout << "rb " << right_bottom << std::endl;



//	Mat schema = Mat::zeros( canny_output.size(), CV_8UC3 );
	Mat schema = source;
	line( schema, left_top, left_bottom, Scalar(255,255,255), 1, CV_AA);
	line( schema, left_bottom, right_bottom, Scalar(255,255,255), 1, CV_AA);
	line( schema, right_bottom, right_top, Scalar(255,255,255), 1, CV_AA);
	line( schema, right_top, left_top, Scalar(255,255,255), 1, CV_AA);

	Scalar blue (255, 0, 0);

	// diags intersection
	Point2f* diags_intersection = intersection(left_top, right_bottom,
						   right_top, left_bottom);

	circle(schema, *diags_intersection, 2, blue);

//	line( schema, right_bottom, left_top, red, 1, CV_AA);
//	line( schema, right_top, left_bottom, red, 1, CV_AA);

	Point2f corner_points_mass_center =
		(left_top + left_bottom + right_bottom + right_top) / 4;

	circle(schema, corner_points_mass_center, 2, Scalar(255, 255, 255));


	// todo check what should be taken dX = mass - diag or dX = diag - mass
	double top_width = right_top.x - left_top.x;
	double bottom_width = right_bottom.x - left_bottom.x;
//	double dX = corner_points_mass_center.x - diags_intersection->x;
	double dX = diags_intersection->x - corner_points_mass_center.x;
//
	double fixed_width = (top_width + bottom_width) / 2 + 2 * dX;

	double left_height = left_top.y - left_bottom.y;
	double right_height = right_top.y - right_bottom.y;
//	double dY = corner_points_mass_center.y - diags_intersection->y;
	double dY = diags_intersection->y - corner_points_mass_center.y;

	double fixed_higth = (left_height + right_height) /2 + 2 * dY;


	// new lt, new rt, new central
	Mat warp_dst = Mat::zeros(source.rows, source.cols, source.type());
	Point2f srcTri[3];
	Point2f dstTri[3];

	// todo working
//	Point2f fixed_rt (right_bottom.x, right_bottom.y + fixed_higth);
//	Point2f fixed_lt (left_bottom.x, left_bottom.y + fixed_higth);
//	srcTri[0] = right_top;
//	srcTri[1] = left_top;
//	srcTri[2] = *diags_intersection;
//
//	dstTri[0] = fixed_rt;
//	dstTri[1] = fixed_lt;
//	dstTri[2] = corner_points_mass_center;

	// todo not working but fun

	//	dstTri[0] = Point2f(right_bottom.x, right_top.y);
//	dstTri[1] = Point2f(left_bottom.x, left_bottom.y);
//	dstTri[2] = corner_points_mass_center;

	// todo best but not working
//	srcTri[0] = right_top;
//	srcTri[1] = left_bottom;
//	srcTri[2] = *diags_intersection;
//
//
//	dstTri[0] = Point2f(
//		(corner_points_mass_center.x + fixed_width / 2),
//		(corner_points_mass_center.y + fixed_higth / 2));
//	dstTri[1] = Point2f(
//		(corner_points_mass_center.x - fixed_width / 2),
//		(corner_points_mass_center.y - fixed_higth / 2));
//	dstTri[2] = corner_points_mass_center;

	srcTri[0] = right_bottom;
	srcTri[1] = right_top;
	srcTri[2] = left_top;
	srcTri[3] = left_bottom;

	dstTri[0] = Point2f(
		(corner_points_mass_center.x + fixed_width / 2),
		(corner_points_mass_center.y - fixed_higth / 2));
	dstTri[1] = Point2f(
		(corner_points_mass_center.x + fixed_width / 2),
		(corner_points_mass_center.y + fixed_higth / 2));
	dstTri[2] = Point2f(
		(corner_points_mass_center.x - fixed_width / 2),
		(corner_points_mass_center.y + fixed_higth / 2));

	dstTri[3] = Point2f(
		(corner_points_mass_center.x - fixed_width / 2),
		(corner_points_mass_center.y - fixed_higth / 2));


	Scalar red (0, 0, 255);
	Scalar green(0, 255, 0);

	circle(schema, srcTri[0], 2, red);
	circle(schema, srcTri[1], 2, green);


//	srcTri[0] = left_bottom;
//	srcTri[1] = right_bottom;
//	srcTri[2] = *diags_intersection;
//
//	dstTri[0] = Point2f(left_top.x, right_bottom.y);
//	dstTri[1] = Point2f(right_top.x, right_bottom.y);
//	dstTri[2] = corner_points_mass_center;
//


//	srcTri[0] = right_top;
//	srcTri[1] = left_top;
//	srcTri[0] = left_bottom;
//	srcTri[1] = right_bottom;
//	srcTri[0] = right_top;
//	srcTri[1] = left_bottom;
//	srcTri[2] = *diags_intersection;
//
////	dstTri[0] = fixed_rt;
////	dstTri[1] = fixed_lt;
//	dstTri[0] = Point2f(corner_points_mass_center.x + fixed_width / 2,
//			    corner_points_mass_center.y + fixed_higth / 2);
//	dstTri[1] = Point2f(corner_points_mass_center.x - fixed_width / 2,
//			    corner_points_mass_center.y - fixed_higth / 2);
//	dstTri[2] = corner_points_mass_center;

//	circle(schema, fixed_rt, 10, red);
//	circle(schema, fixed_lt, 10, red);

//	line( schema, fixed_lt, left_bottom, red, 1, CV_AA);
//	line( schema, dstTri[0], dstTri[1], red, 1, CV_AA);
//
//	line( schema, left_bottom, right_bottom, red, 1, CV_AA);
//	line( schema, right_bottom, fixed_rt, red, 1, CV_AA);
//	line( schema, fixed_rt, left_top, red, 1, CV_AA);
//
//	line( schema, right_bottom, fixed_lt, red, 1, CV_AA);
//	line( schema, fixed_rt, left_bottom, red, 1, CV_AA);


	imshow("Detected Lines (in red) - Standard Hough Line Transform", schema);


//	Mat warp_mat = getAffineTransform( srcTri, dstTri );
	Mat warp_p_mat = getPerspectiveTransform(srcTri, dstTri);

//	warp_mat = getAffineTransform( srcTri, dstTri );
//	warpAffine( source, warp_dst, warp_mat, warp_dst.size() );
//	warpAffine( source, warp_dst, warp_mat, source.size() );

	warpPerspective(source, warp_dst, warp_p_mat, warp_dst.size());


	circle(warp_dst, dstTri[0], 10, red);
	circle(warp_dst, dstTri[1], 10, green);
	circle(warp_dst, dstTri[2], 10, blue);

//	Point2f rt2, rb2, lt2, lb2;
//	rt2 = dstTri[0];
//	lb2 = dstTri[1];
//	rb2 = Point2f(rt2.x, rt2.y - fixed_higth);
//	lt2 = Point2f(lb2.x, lb2.y + fixed_higth);
//
//	line( warp_dst, lb2, rb2, red, 1, CV_AA);
//	line( warp_dst, rb2, rt2, red, 1, CV_AA);
//	line( warp_dst, rt2, lt2, red, 1, CV_AA);
//	line( warp_dst, lt2, lb2, red, 1, CV_AA);

//	srcTri[0] = right_bottom;
//	srcTri[1] = right_top;
//	srcTri[2] = *diags_intersection;
//
//	dstTri[0] = Point2f(
//		(corner_points_mass_center.x + fixed_width / 2),
//		(corner_points_mass_center.y - fixed_higth / 2));
//	dstTri[1] = Point2f(
//		(corner_points_mass_center.x + fixed_width / 2),
//		(corner_points_mass_center.y + fixed_higth / 2));
//	dstTri[2] = corner_points_mass_center;




	const char* warp_window = "Warp";
	namedWindow( warp_window, WINDOW_AUTOSIZE );
	imshow( warp_window, warp_dst );

	center = Point(warp_dst.cols/2, warp_dst.rows/2);
	angle = 180.0;
	scale = 1;
	rot_mat = getRotationMatrix2D( center, angle, scale );
	warpAffine(warp_dst, warp_dst, rot_mat, warp_dst.size() );

	std::string fixed_filename =  argv[1];
	fixed_filename = fixed_filename.substr (0, fixed_filename.size() - 4);

	fixed_filename.append("_fixed.jpg");
	std::cout << fixed_filename << std::endl;
	imwrite(fixed_filename, warp_dst );

	waitKey(0);

	return(0);
}
