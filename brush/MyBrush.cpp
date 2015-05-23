#include "csc321.h"
#include "MyBrush.h"
#include "BrushInterface.h"
#include <cmath>
#include <iostream>
using namespace std;

void MyBrush::changedBrush() {
	// this is called anytime the brush type or brush radius changes
	// it should recompute the brush mask appropriately
	const int radius = brushUI->getRadius();
	//float mask[radius + 1];
	//mask.resize(radius*2+1, vector<float>(radius*2+1));
	mask.resize(radius*2+1);
	for (int i=0; i< mask.size(); i++) {
		mask[i].resize(radius*2+1);
	}



	switch(brushUI->getBrushType()) {
	case BRUSH_CONSTANT:
		for (int i = 0; i < radius*2+1; i++) {
			for (int j = 0; j < radius*2+1; j++) {
				mask[i][j] = ((i-radius)*(i-radius)+(j-radius)*(j-radius) < radius*radius) ? 1.0 : 0.0;
			}
		}
		break;
	case BRUSH_LINEAR:
		for (int i = 0; i < radius*2+1; i++) {
			for (int j = 0; j < radius*2+1; j++) {
				mask[i][j] = ((i-radius)*(i-radius)+(j-radius)*(j-radius) < radius*radius) ? 1 - sqrt(static_cast<float>((i-radius)*(i-radius) + (j-radius)*(j-radius)))/radius: 0.0;
			}
		}
		break;
	case BRUSH_QUADRATIC:
		for (int i = 0; i < radius*2+1; i++) {
			for (int j = 0; j < radius*2+1; j++) {
				mask[i][j] = ((i-radius)*(i-radius)+(j-radius)*(j-radius) < radius*radius) ? pow(1 - sqrt(static_cast<float>((i-radius)*(i-radius) + (j-radius)*(j-radius)))/radius, 2): 0.0;
			}
		}
		break;
	case BRUSH_GAUSSIAN:
		for (int i = 0; i < radius*2+1; i++) {
			for (int j = 0; j < radius*2+1; j++) {
				mask[i][j] = ((i-radius)*(i-radius)+(j-radius)*(j-radius) < radius*radius) ? exp( -0.5 / 100 * pow(sqrt(static_cast<float>((i-radius)*(i-radius) + (j-radius)*(j-radius))), 2)) : 0.0;
			}
		}
		break;
	case BRUSH_SPECIAL:  //special mixed brush
		float x = mouseDrag[0];
		float y = mouseDrag[1];
		for (int i = 0; i < radius*2+1; i++) {
			for (int j = 0; j < radius*2+1; j++) {
				mask[i][j] = ((i-radius)*(i-radius)+(j-radius)*(j-radius) < radius*radius) ? pow(1 - sqrt((i-x)*(i-x) + (j-y)*(j-y))/radius, 2): 0.0;
			}
		}
		break;
	}
}



void MyBrush::drawBrush( ) {
	// apply the current brush mask to image location (x,y)
	// the mouse location is in mouseDrag

	const int radius = brushUI->getRadius();
	const float pixelFlow = brushUI->getFlow();
	const Color colBrush = brushUI->getColor();

	//check boundary and assign indices
	int xmin = (mouseDrag[0] - radius < 0) ? 0 : mouseDrag[0] - radius;
	int xmax = (mouseDrag[0] + radius > imageWidth) ? imageWidth : mouseDrag[0] + radius;
	int ymin = (mouseDrag[1] - radius < 0) ? 0 : mouseDrag[1] - radius;
	int ymax = (mouseDrag[1] + radius > imageHeight) ? imageHeight : mouseDrag[1] + radius;

	//apply the mask
	for (int i = xmin; i < xmax; i++) {
		for (int j = ymin; j < ymax; j++) {
			putPixel(i, j, mask[i-(mouseDrag[0]-radius)][j-(mouseDrag[1]-radius)]*pixelFlow*colBrush + (1-mask[i-(mouseDrag[0]-radius)][j-(mouseDrag[1]-radius)]*pixelFlow)*getPixel(i, j));
		}
	}


}

void MyBrush::drawLine( ) {
	// draw a thick line from mouseDown to mouseDrag
	// the width of the line is given by the current brush radius
	const int radius = brushUI->getRadius();
	const Color colBrush = brushUI->getColor();
	//start point (x0,y0), end point (x1,y1)
	int x0 = mouseDown[0], y0 = mouseDown[1], x1 = mouseDrag[0], y1 = mouseDrag[1]; 

	if (x1 == x0){// draw a vertical line with a width of radius
		int x;
		for(int w = -radius/2; w <= radius/2; w++){
			x= x0 + w; 
			for (int y = min(y0, y1); y < max(y0, y1); y++)
				if(x >= 0 && x < screenWidth && y >= 0 && y < screenHeight) // check boundary
					putPixel(x,y,colBrush);
		}
	}
	else if(y1 == y0){// draw a horizontal line with a width of radius
		int y;
		for(int w = -radius/2; w <= radius/2; w++){
			y= y0 + w;
			for(int x = min(x0,x1); x < max(x0,x1); x++)
				if(x >= 0 && x < screenWidth && y >= 0 && y < screenHeight) // check boundary
					putPixel(x,y,colBrush);
		}
	}
	else{
		float s1 = (y1-y0)/float(x1-x0); // slope of the two border lines in a rectangle
		float s2 = -1/s1; // slope of the other two border lines in a rectangle

		float xtemp = (radius/(float)(2*sqrt(s2*s2 + 1))) + x0;
		float ytemp = s2*(xtemp-x0) + y0;
		int x00 =(int) xtemp;
		int y00 = (int) (ytemp);
		xtemp = -(radius/(float)(2*sqrt(s2*s2+1)))+x0;
		ytemp = s2*(xtemp-x0)+y0;
		int x01 = (int) xtemp,y01 = (int) ytemp;
		xtemp = (radius/float(2*sqrt(s2*s2+1)))+x1; 
		ytemp = s2*(xtemp-x1)+y1;
		int x10 = (int) xtemp,y10 = (int) ytemp;
		xtemp = -(radius/float(2*sqrt(s2*s2+1)))+x1; 
		ytemp = s2*(xtemp-x1)+y1;
		int x11 = (int) xtemp,y11 = (int) ytemp;

		//draw the border lines of a rectangle
		if(abs(s1) < 1){
			drawFlatLine(colBrush,x00,y00,x10,y10);
			drawFlatLine(colBrush,x01,y01,x11,y11);  
		}
		else {
			drawSteepLine(colBrush,x00,y00,x10,y10);
			drawSteepLine(colBrush,x01,y01,x11,y11); 
		}

		if(abs(s2) < 1){
			drawFlatLine(colBrush,x00,y00,x01,y01);
			drawFlatLine(colBrush,x10,y10,x11,y11);  
		}
		else {
			drawSteepLine(colBrush,x00,y00,x01,y01);
			drawSteepLine(colBrush,x00,y00,x01,y01); 
		}


		// fill the polygon betweeen the lines
		fillLine(colBrush,x00,y00,x01,y01,x10,y10,x11,y11);
	}	

}



//draw a line with 1px width when -1<= slope <=1
void MyBrush::drawFlatLine(Color col, int x0, int y0, int x1, int y1){
	if (x0 > x1){ // swap start point and end point if the start point of the line is on the right to the end point
		int temp = x1;
		x1 = x0;
		x0 = temp;
		temp = y1;
		y1 = y0;
		y0 = temp;
	}
	int dx = x1 - x0, dy = y1 - y0;
	int y_incr = 1; // increment of y in each decision
	if (dy < 0){ // if the start point is higher than the end point, negate
		y_incr = -1;
		dy = -dy;
	}

	int incrE = 2*dy;
	int incrNE = 2*(dy - dx);;
	int d = 2*dy - dx;
	int y = y0;

	for (int x = x0; x <= x1; x++){
		if(x >= 0 && x < screenWidth && y >= 0 && y < screenHeight) // check boundary
			putPixel(x,y,col);
		if (d <= 0)
			d += incrE;
		else{
			d += incrNE;
			y += y_incr;
		}
	}
}

//draw a line with 1px width when slope > 1 or slope < -1
void MyBrush::drawSteepLine(Color col, int x0, int y0, int x1, int y1){
	// swap x,y value for start and end points individually so E/NE algorithm can be used
	int temp = x1;
	x1 = y1;
	y1 = temp;
	temp = x0;
	x0 = y0;
	y0 = temp;

	if (x0 > x1){// swap start point and end point if the start point of the line is on the right to the end point 
		temp = x1;
		x1 = x0;
		x0 = temp;
		temp = y1;
		y1 = y0;
		y0 = temp;
	}

	int dx = x1 - x0, dy = y1 - y0;
	int y_incr = 1; // increment of y in each decision
	if (dy < 0){ // if the start point is higher than the end point, negate
		y_incr = -1;
		dy = -dy;
	}

	int incrE = 2*dy;
	int incrNE = 2*(dy - dx);;
	int d = 2*dy - dx;
	int y = y0;

	for (int x = x0; x <= x1; x++){
		if(y >= 0 && y < screenWidth && x >= 0 && x < screenHeight) // check boundary
			putPixel(y,x,col);
		if (d <= 0)
			d += incrE;
		else{
			d += incrNE;
			y += y_incr;
		}
	}
}

// fill the polygon between four border lines
void MyBrush::fillLine(Color col, int x00, int y00, int x01, int y01, int x10, int y10, int x11, int y11){
	float m1 = (y01-y00)/float(x01-x00);
	float b1 = y00-m1*x00;

	float m2 = (y10-y00)/float(x10-x00);
	float b2 = y00-m2*x00;

	float m3 = (y11-y10)/float(x11-x10);
	float b3 = y10-m3*x10;

	float m4 = (y11-y01)/float(x11-x01);
	float b4 = y01-m4*x01;

	int y_min = min(min(y00,y01),min(y10,y11));
	int y_max = max(max(y00,y01),max(y10,y11));

	vector<int> fillInters; 
	int numfillInters;


	for(int y = y_min; y <= y_max; y++){
		if(y >= min(y00,y01) && y <= max(y00,y01)){
			fillInters.push_back((y-b1)/m1);
		}
		if(y >= min(y00,y10) && y <= max(y00,y10)){
			fillInters.push_back((y-b2)/m2);
		}
		if(y >= min(y10,y11) && y <= max(y10,y11)){
			fillInters.push_back((y-b3)/m3);
		}
		if(y >= min(y11,y01) && y <= max(y11,y01)){
			fillInters.push_back((y-b4)/m4);
		}
		numfillInters = fillInters.size();
		if(numfillInters == 2){
			int x_min = min(fillInters[0],fillInters[1]),x_max = max(fillInters[0],fillInters[1]);
			for(int x=x_min; x<=x_max; x++)
				if(x >= 0 && x < screenWidth && y >= 0 && y < screenHeight)
					putPixel(x,y,col);
		}
		else if(numfillInters == 1){
			int x = fillInters[0];
			if(x >= 0 && x < screenWidth && y >= 0 && y < screenHeight)
				putPixel(x,y,col);
		}

		else if(numfillInters == 3){
			int x1 = fillInters[0], x2 = fillInters[1], x3 = fillInters[2];
			int x_min = min(min(x1,x2),x3),x_max = max(max(x1,x2),x3);
			for(int x = x_min; x<= x_max; x++)
				if(x >= 0 && x < screenWidth && y >= 0 && y < screenHeight)
					putPixel(x,y,col);
		}		
		fillInters.resize(0);
	}
}





void MyBrush::drawCircle() {
	// draw a thick circle at mouseDown with radius r
	// the width of the circle is given by the current brush radius
	float x0, y0, x1, y1;
	int R, W; // radius of the circle, width of the circle
	const Color colBrush = brushUI->getColor();
	const int width = brushUI->getRadius();

	x0 = mouseDown[0];
	y0 = mouseDown[1];
	x1 = mouseDrag[0];
	y1 = mouseDrag[1];
	R = sqrt(pow(x0 - x1,2) + pow(y0 - y1, 2));
	if (R <= width) {
		mask.resize(width*2+1);
		for (int i=0; i< mask.size(); i++) {
			mask[i].resize(width*2+1);
		}
		for (int i = 0; i < width*2+1; i++) {
			for (int j = 0; j < width*2+1; j++) {
				if ((i-width)*(i-width)+(j-width)*(j-width) < width*width) {
					mask[i][j] = 1.0;
				} else mask[i][j] = 0.0;
			}
		} 

		//check boundary and assign indices
		int xmin = max(x0 - width, 0);
		int xmax = min(x0 + width, imageWidth);
		int ymin = max(y0 - width, 0);
		int ymax = min(y0 + width, imageHeight);

		//apply the mask
		for (int i = xmin; i < xmax; i++) {
			for (int j = ymin; j < ymax; j++) {
				putPixel(i, j, mask[i-(x0-width)][j-(y0-width)]*colBrush + (1-mask[i-(x0-width)][j-(y0-width)])*getPixel(i,j));
			}
		}

	}
	else {
		mask.resize(R*2+1);
		for (int i=0; i< mask.size(); i++) {
			mask[i].resize(R*2+1);
		}
		for (int i = 0; i < R*2+1; i++) {
			for (int j = 0; j < R*2+1; j++) {
				if ((i-R)*(i-R)+(j-R)*(j-R) < R*R && (i-R)*(i-R)+(j-R)*(j-R) > (R-width)*(R-width)) {
					mask[i][j] = 1.0;
				} else mask[i][j] = 0.0;
			}
		}
		//check boundary and assign indices
		int xmin = max(x0 - R, 0);
		int xmax = min(x0 + R, imageWidth);
		int ymin = max(y0 - R, 0);
		int ymax = min(y0 + R, imageHeight);

		//apply the mask
		for (int i = xmin; i < xmax; i++) {
			for (int j = ymin; j < ymax; j++) {
				putPixel(i, j, mask[i-(x0-R)][j-(y0-R)]*colBrush + (1-mask[i-(x0-R)][j-(y0-R)])*getPixel(i,j));
			}
		}
	}

}


void MyBrush::drawPolygon() {
	// draw a polygon with numVertices whos coordinates are stored in the
	// polygon array: {x0, y0, x1, y1, ...., xn-1, yn-1}
	const float pixelFlow = brushUI->getFlow();
	const Color colBrush = brushUI->getColor();

}

void MyBrush::filterRegion( ) {
	// apply the filter indicated by filterType to the square
	// defined by the two corner points mouseDown and mouseDrag
	// these corners are not guarenteed to be in any order
	// The filter width is given by the brush radius

	//int radius = brushUI->getRadius();
	//if (radius % 2 == 0) {
	//	radius++;	
	//}

	//mask.resize(radius);
	//for (int i=0; i< mask.size(); i++) {
	//	mask[i].resize(radius);
	//}

	//for (int i = 0; i < radius; i++) {
	//	for (int j = 0; j < radius; j++) {
	//		mask[i][j] = exp( -0.5 / 100 * pow(sqrt(static_cast<float>((i-radius/2)*(i-radius/2) + (j-radius/2)*(j-radius/2))), 2));
	//	}
	//}


	////Set up x and y
	//int x0, y0, x1, y1, xStart, yStart, xEnd, yEnd, iSrc, jSrc;
	//long k;
	//Color colSum;
	//vector<Color> colorRegion;


	//x0 = mouseDown[0];
	//y0 = mouseDown[1];
	//x1 = mouseDrag[0];
	//y1 = mouseDrag[1];
	//xStart = min(x0, x1);
	//yStart = min(y0, y1);
	//xEnd = max(x0, x1);
	//yEnd = max(y0, y1);
	//for ( int x = xStart; x <= xEnd; x++ ) {
	//	for ( int y = yStart; y <= yEnd; y++ ) {
	//		float r_sum = 0.0;
	//		float g_sum = 0.0;
	//		float b_sum = 0.0;
	//		//colSum = (0,0,0);
	//		for ( int i = -radius/2; i <= radius/2; i++ ) {
	//			iSrc = x + i;
	//			if ( iSrc < 0 ) iSrc = 0;
	//			if ( iSrc >= imageWidth ) iSrc = imageWidth - 1;
	//			for ( int j = -radius/2; j <= radius/2; j++ ) {
	//				jSrc = y + j;
	//				if ( jSrc < 0 ) jSrc = 0;
	//				if ( jSrc >= imageHeight ) jSrc = imageHeight - 1;
	//				r_sum += getPixel(iSrc, jSrc)[0] * mask[i+radius/2][j+radius/2];
	//				g_sum += getPixel(iSrc, jSrc)[1] * mask[i+radius/2][j+radius/2];
	//				b_sum += getPixel(iSrc, jSrc)[2] * mask[i+radius/2][j+radius/2];
	//			}
	//		}
	//		Color temp(r_sum, g_sum, b_sum);
	//		//colorRegion[k] = temp; // copy the calculated pixel
	//		k++;
	//	}
	//
	//}

	//// change pixels of the filtered region
	//for ( int x = xStart; x <= xEnd; x++ ) {
	//	for ( int y = yStart; y <= yEnd; y++ ) {
	//		putPixel(x, y, colorRegion[k]);
	//		k++;
	//	}
	//}
}
