#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"

using namespace cv;

/*
 * Finds the intersection of two lines.
 * Lines are set with two dots (o, p).
 * function returns NULL if lines do not intersect.
 */
Point2f *intersection(const Point2f &o1, const Point2f &p1, const Point2f &o2, const Point2f &p2) {
	Point2f x = o2 - o1;
	Point2f d1 = p1 - o1;
	Point2f d2 = p2 - o2;

	float cross = d1.x * d2.y - d1.y * d2.x;
	if (abs(cross) < 1e-8)
		return NULL;

	float t1 = (x.x * d2.y - x.y * d2.x) / cross;
	Point2f r = o1 + d1 * t1;
	return new Point2f(r);
}

/*
 * Conceptually procedure consists of the following steps:
 * 1) Find document contour as a set of points.
 * 2) Find horizontal and vertical borders of the document, using
 *    Hough algorithm.
 * 3) Find corner points of document (as intersection points of
 *    vertical and horizontal Hough lines)
 * 4) Reduce corner points (find center of mass of each corner
 *    set of points. For example find center of mass of all
 *    points on the top-right, then use this center of mass as
 *    top-righ corner point.
 * 5) Find intersection of contours diagonals and center of mass
 *    of reduced corner points.
 * 6) Find parameters for perspective transformation which would
 *    сombine two points mentioned in step 5.
 * 7) Perform transformation.
 */
int main(int argc, char **argv) {
	Mat source = imread(argv[1], 1);

	// OpenCV uses coordinate system with Y axis directed
	// down, so we rotate the image before and
	// after processing and assume that Y axis  directed up.
	Point center = Point(source.cols / 2, source.rows / 2);
	float angle = 180.0;
	float scale = 1;
	Mat rot_mat = getRotationMatrix2D(center, angle, scale);
	warpAffine(source, source, rot_mat, source.size());

	Mat source_gray;
	cvtColor(source, source_gray, COLOR_BGR2GRAY);

	Mat canny_output;
	std::vector<std::vector<Point>> contours;
	std::vector<Vec4i> hierarchy;
	int thresh = 100;
	Canny(source_gray, canny_output, thresh, thresh * 2, 3);

	findContours(canny_output, contours, hierarchy, RETR_TREE,
		     CHAIN_APPROX_SIMPLE);

	// Find contour with maximum area.
	int contour_id = 0;
	double max_area = 0;
	for (int i = 0; i < contours.size(); i++) {
		double contour_area = contourArea(contours[i]);
		if (contour_area > max_area) {
			max_area = contour_area;
			contour_id = i;
		}
	}

	Mat contour_drawing = Mat::zeros(canny_output.size(), CV_8UC3);
	Scalar white = Scalar(255, 255, 255);
	drawContours(contour_drawing, contours, contour_id, white, 2, 8,
		     hierarchy, 0, Point());

	Mat dst;
	Canny(contour_drawing, dst, 50, 200, 3);

	// Get document border lines using Hough algorithm.
	std::vector<Vec2f> lines;
	std::vector<std::pair<Point2f, Point2f>> vertical_lines;
	std::vector<std::pair<Point2f, Point2f>> horizontal_lines;
	HoughLines(dst, lines, 1, CV_PI / 180, 150, 0, 0);

	for (size_t i = 0; i < lines.size(); i++) {
		float rho = lines[i][0];
		float theta = lines[i][1];
		Point2f pt1(0, 0);
		Point2f pt2(0, 0);
		double a = cos(theta);
		double b = sin(theta);
		double x0 = a * rho;
		double y0 = b * rho;
		pt1.x = cvRound(x0 + 1000 * (-b));
		pt1.y = cvRound(y0 + 1000 * (a));
		pt2.x = cvRound(x0 - 1000 * (-b));
		pt2.y = cvRound(y0 - 1000 * (a));

		std::pair<Point2f, Point2f> new_pair(pt1, pt2);
		if (abs(pt1.x - pt2.x) < abs(pt1.y - pt2.y))
			vertical_lines.push_back(new_pair);
		else
			horizontal_lines.push_back(new_pair);
	}

	std::vector<Point2f> intersection_points;
	Point2f *new_point;
	Point2f mass_center(0, 0);
	for (int i = 0; i < vertical_lines.size(); i++) {
		for (int j = 0; j < horizontal_lines.size(); j++) {
			Point2f vertical_1 = vertical_lines[i].first;
			Point2f vertical_2 = vertical_lines[i].second;

			Point2f horizontal_1 = horizontal_lines[j].first;
			Point2f horizontal_2 = horizontal_lines[j].second;

			new_point = intersection(vertical_1, vertical_2,
						 horizontal_1, horizontal_2);

			if (new_point != NULL) {
				intersection_points.push_back(*new_point);
				mass_center += *new_point;
			}
		}
	}

	mass_center /= (double) intersection_points.size();
	Point2f left_top(0, 0);
	Point2f left_bottom(0, 0);
	Point2f right_top(0, 0);
	Point2f right_bottom(0, 0);

	int lt_size = 0;
	int lb_size = 0;
	int rt_size = 0;
	int rb_size = 0;

	for (int i = 0; i < intersection_points.size(); i++) {
		Point2f *p = &intersection_points[i];
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

	Point2f *diags_intersection = intersection(left_top, right_bottom,
						   right_top, left_bottom);

	Point2f corner_points_mass_center =
		(left_top + left_bottom + right_bottom + right_top) / 4;

	double top_width = right_top.x - left_top.x;
	double bottom_width = right_bottom.x - left_bottom.x;
	double d_x = diags_intersection->x - corner_points_mass_center.x;
	double fixed_width = (top_width + bottom_width) / 2 + 2 * d_x;

	double left_height = left_top.y - left_bottom.y;
	double right_height = right_top.y - right_bottom.y;
	double d_y = diags_intersection->y - corner_points_mass_center.y;
	double fixed_higth = (left_height + right_height) / 2 + 2 * d_y;

	Mat warp_dst = Mat::zeros(source.rows, source.cols, source.type());
	Point2f orig_points[4];
	Point2f fixed_points[4];

	orig_points[0] = right_bottom;
	orig_points[1] = right_top;
	orig_points[2] = left_top;
	orig_points[3] = left_bottom;

	fixed_points[0] = Point2f(
		(corner_points_mass_center.x + fixed_width / 2),
		(corner_points_mass_center.y - fixed_higth / 2));
	fixed_points[1] = Point2f(
		(corner_points_mass_center.x + fixed_width / 2),
		(corner_points_mass_center.y + fixed_higth / 2));
	fixed_points[2] = Point2f(
		(corner_points_mass_center.x - fixed_width / 2),
		(corner_points_mass_center.y + fixed_higth / 2));

	fixed_points[3] = Point2f(
		(corner_points_mass_center.x - fixed_width / 2),
		(corner_points_mass_center.y - fixed_higth / 2));


	Mat warp_p_mat = getPerspectiveTransform(orig_points, fixed_points);
	warpPerspective(source, warp_dst, warp_p_mat, warp_dst.size());

	center = Point(warp_dst.cols / 2, warp_dst.rows / 2);
	angle = 180.0;
	scale = 1;
	rot_mat = getRotationMatrix2D(center, angle, scale);
	warpAffine(warp_dst, warp_dst, rot_mat, warp_dst.size());

	std::string fixed_filename = argv[1];
	fixed_filename = fixed_filename.substr(0, fixed_filename.size() - 4);

	fixed_filename.append("_fixed.jpg");
	std::cout << fixed_filename << std::endl;
	imwrite(fixed_filename, warp_dst);

	return (0);
}
