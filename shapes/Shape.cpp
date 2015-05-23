#include "../csc321.h"
#include "Shape.h"


Shape::~Shape() {
}

void Shape::draw() {
	glBegin(GL_TRIANGLES);
	for (int i=0; i < vertices.size(); i +=3) {
		glNormal3dv(&normals[i][0]);
		glVertex3dv(&vertices[i][0]);
		glNormal3dv(&normals[i+1][0]);
		glVertex3dv(&vertices[i+1][0]);
		glNormal3dv(&normals[i+2][0]);
		glVertex3dv(&vertices[i+2][0]);
	}
	glEnd();
}

Shape::Shape() : vertices(), normals() {
}

HitRecord Shape::intersect(const Point3& o, const Vector3& dir) {
	return HitRecord();
}

void Shape::addTriangle(const Point3& p1, const Point3& p2, const Point3& p3,
	const Vector3& n1, const Vector3& n2, const Vector3& n3) {
		vertices.push_back(p1);
		vertices.push_back(p2);
		vertices.push_back(p3);
		normals.push_back(n1);
		normals.push_back(n2);
		normals.push_back(n3);
}

Cylinder::Cylinder(int n, int m) : Shape() {
	double ca = 0.5;
	double sa = 0.0;
	for (int i = 1; i <= n; i++) {
		double a = (2 * M_PI * i) / n;
		double ca1 = 0.5 * cos(a);
		double sa1 = 0.5 * sin(a);

		//cap
		addTriangle(Point3(0.0, 0.5, 0.0), Point3(ca1, 0.5, sa1), Point3(ca, 0.5, sa), 
			Vector3(0.0, 1.0, 0.0), Vector3(0.0, 1.0, 0.0), Vector3(0.0, 1.0, 0.0));
		//bottom
		addTriangle(Point3(0.0, -0.5, 0.0), Point3(ca, -0.5, sa), Point3(ca1, -0.5, sa1), 
			Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0));

		//barrel
		double cb = m;
		for (int j = 1; j <= m; j++) {
			double b = 1.0 / m;
			double cb1 = 0.5 - j * b;
			double cb2 = 0.5 - (j-1) * b;

			addTriangle(Point3(ca1, cb2, sa1), Point3(ca, cb1, sa), Point3(ca, cb2, sa), 
				Vector3(ca1, 0.0, sa1), Vector3(ca, 0.0, sa), Vector3(ca, 0.0, sa));
			addTriangle(Point3(ca1, cb2, sa1), Point3(ca1, cb1, sa1), Point3(ca, cb1, sa),
				Vector3(ca1, 0.0, sa1),  Vector3(ca1, 0.0, sa1), Vector3(ca, 0.0, sa));

		}

		ca = ca1;
		sa = sa1;
	}
}

// cylinder intersect method
HitRecord Cylinder::intersect(const Point3& o, const Vector3& d) {
	HitRecord hr;

	// body
	// solve a*t^2 + b*t + c = 0
	double a = d[0]*d[0] + d[2]*d[2];
	double b = 2.0*( d[0]*o[0] + d[2]*o[2] );
	double c = o[0]*o[0] + o[2]*o[2] - 0.25;

	double dis = b*b - 4.0*a*c;
	if (dis > 0) {
		dis = sqrt(dis);
		double t1 = ( dis - b) / (2.0*a);
		double t2 = (-dis - b) / (2.0*a);

		Point3 p;
		Vector3 n;
		if (t1 > 0.0) {
			p = o + t1*d;
			if ((p[1] < 0.5) && (p[1] > -0.5)) {
				n = Vector3(p[0], 0, p[2]);
				n.normalize();
				hr.addHit(t1, 0.0, 0.0, p, n);
			}
		}
		if (t2 > 0.0) {  // else
			p = o + t2*d;
			if ((p[1] < 0.5) && (p[1] > -0.5)) {
				n = Vector3(p[0], 0, p[2]);
				n.normalize();
				hr.addHit(t2, 0.0, 0.0, p, n);
			}
		}
	}

	// cap
	double tc = (0.5 - o[1])/d[1];
	Point3 p2;
	Vector3 n2;
	p2 = o + tc*d;
	if ( tc > 0  && (p2[0]*p2[0] + p2[2]*p2[2] < 0.25)) {
		n2 = Vector3(0, 0.5, 0);
		n2.normalize();
		hr.addHit(tc, 0.0, 0.0, p2, n2);
	}
	// bottom
	double tb = (0.5 + o[1])/(-d[1]);
	Point3 p3;
	Vector3 n3;
	p3 = o + tb*d;
	if ( tb > 0 && (p3[0]*p3[0] + p3[2]*p3[2] < 0.25)) {
		n3 = Vector3(0, -0.5, 0);
		n3.normalize();
		hr.addHit(tb, 0.0, 0.0, p3, n3);
	}
	hr.sortHits();
	return hr;
}


Cone::Cone(int n, int m) : Shape() {
	double ca = 0.5;
	double sa = 0.0;
	if (n < 3) {
		n = 3;
	}
	for (int i = 1; i <= n; i++) {
		double a = (2 * M_PI * i) / n;
		double ca1 = 0.5 * cos(a);
		double sa1 = 0.5 * sin(a);

		//cap
		addTriangle(Point3(0.0, -0.5, 0.0), Point3(ca, -0.5, sa), Point3(ca1, -0.5, sa1), 
			Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0));

		double ca2 = ca;
		double sa2 = sa;
		double ca3 = ca1;
		double sa3 = sa1;
		//sphere
		for (int j = 1; j <= m; j++) {
			double d = j * 1.0 / m;		
			double ca4 = (1.0-d) / 1.0 * ca; // upper x
			double sa4 = (1.0-d) / 1.0 * sa; // upper y
			double ca5 = (1.0-d) / 1.0 * ca1; // upper x
			double sa5 = (1.0-d) / 1.0 * sa1; // upper y
			double h1 = -0.5 + (j - 1) *1.0/ (float) m;
			double h2 = -0.5 + d;

			addTriangle(Point3(ca3, h1, sa3), Point3(ca2, h1, sa2), Point3(ca5, h2, sa5),
				Vector3(ca3, 0.5*sqrt(ca3*ca3 + sa3*sa3), sa3), Vector3(ca2, 0.5*sqrt(ca2*ca2 + sa2*sa2), sa2), Vector3(ca3, 0.5*sqrt(ca3*ca3 + sa3*sa3), sa3));

			if (j != m) {

				addTriangle(Point3(ca5, h2, sa5), Point3(ca2, h1, sa2), Point3(ca4, h2, sa4), 
					Vector3(ca3, 0.5*sqrt(ca3*ca3 + sa3*sa3), sa3), Vector3(ca2, 0.5*sqrt(ca2*ca2 + sa2*sa2), sa2), Vector3(ca2, 0.5*sqrt(ca2*ca2 + sa2*sa2), sa2));
			}
			ca2 = ca4;
			sa2 = sa4;
			ca3 = ca5;
			sa3 = sa5;

		}


		ca = ca1;
		sa = sa1;

	}
}

// cone intersection method
HitRecord Cone::intersect(const Point3& o, const Vector3& d) {
	HitRecord hr;

	// body
	// solve a*t^2 + b*t + c = 0
	double a = d[0]*d[0] + d[2]*d[2] - 0.25*d[1]*d[1];
	double b = 2.0*( d[0]*o[0] + d[2]*o[2] ) + 0.5*d[1]*(0.5-o[1]);
	double c = o[0]*o[0] + o[2]*o[2] - 0.25*(0.5-o[1])*(0.5-o[1]);

	double dis = b*b - 4.0*a*c;
	if (dis > 0) {
		dis = sqrt(dis);
		double t1 = ( dis - b) / (2.0*a);
		double t2 = (-dis - b) / (2.0*a);

		Point3 p;
		Vector3 n;
		if (t1 > 0.0) {
			p = o + t1*d;
			if ((p[1] < 0.5) && (p[1] > -0.5)) {
				n = Vector3(p[0], 0.5*sqrt(p[0]*p[0]+p[2]*p[2]), p[2]);
				n.normalize();
				hr.addHit(t1, 0.0, 0.0, p, n);
			}
		}
		if (t2 > 0.0) {  // else
			p = o + t2*d;
			if ((p[1] < 0.5) && (p[1] > -0.5)) {
				n = Vector3(p[0], 0.5*sqrt(p[0]*p[0]+p[2]*p[2]), p[2]);
				n.normalize();
				hr.addHit(t2, 0.0, 0.0, p, n);
			}
		}
	}


	// bottom
	double tb = (0.5 + o[1])/(-d[1]);
	Point3 p3;
	Vector3 n3;
	p3 = o + tb*d;
	if (tb > 0 && (p3[0]*p3[0] + p3[2]*p3[2] < 0.25)) {
		n3 = Vector3(0, p3[1], 0);
		n3.normalize();
		hr.addHit(tb, 0.0, 0.0, p3, n3);
	}

	hr.sortHits();
	return hr;
}



Sphere::Sphere(int n, int m) : Shape() {
	if (n < 3) {
		n = 3;
	}
	for (int i = 0; i < n; i++) {
		double a = (2 * M_PI * i) / n;
		double a2 = (2 * M_PI * (i+1)) / n;


		for (int j = 0; j < m; j++) {
			double b = (M_PI/2.0 * j) / m;
			double b2 = (M_PI/2.0 * (j+1)) / m;	

			addTriangle(Point3(0.5*cos(b)*cos(a2), 0.5*sin(b), 0.5*cos(b)*sin(a2)), Point3(0.5*cos(b)*cos(a), 0.5*sin(b), 0.5*cos(b)*sin(a)), Point3(0.5*cos(b2)*cos(a), 0.5*sin(b2), 0.5*cos(b2)*sin(a)), 
				Vector3(0.5*cos(b)*cos(a2), 0.5*sin(b), 0.5*cos(b)*sin(a2)), Vector3(0.5*cos(b)*cos(a), 0.5*sin(b), 0.5*cos(b)*sin(a)), Vector3(0.5*cos(b2)*cos(a), 0.5*sin(b2), 0.5*cos(b2)*sin(a)));

			addTriangle(Point3(0.5*cos(b2)*cos(a), 0.5*sin(b2), 0.5*cos(b2)*sin(a)), Point3(0.5*cos(b2)*cos(a2), 0.5*sin(b2), 0.5*cos(b2)*sin(a2)), Point3(0.5*cos(b)*cos(a2), 0.5*sin(b), 0.5*cos(b)*sin(a2)), 
				Vector3(0.5*cos(b2)*cos(a), 0.5*sin(b2), 0.5*cos(b2)*sin(a)), Vector3(0.5*cos(b2)*cos(a2), 0.5*sin(b2), 0.5*cos(b2)*sin(a2)), Vector3(0.5*cos(b)*cos(a2), 0.5*sin(b), 0.5*cos(b)*sin(a2)));

			addTriangle(Point3(0.5*cos(b)*cos(a), -0.5*sin(b), 0.5*cos(b)*sin(a)), Point3(0.5*cos(b)*cos(a2), -0.5*sin(b), 0.5*cos(b)*sin(a2)), Point3(0.5*cos(b2)*cos(a), -0.5*sin(b2), 0.5*cos(b2)*sin(a)), 
				Vector3(0.5*cos(b)*cos(a), -0.5*sin(b), 0.5*cos(b)*sin(a)), Vector3(0.5*cos(b)*cos(a2), -0.5*sin(b), 0.5*cos(b)*sin(a2)), Vector3(0.5*cos(b2)*cos(a), -0.5*sin(b2), 0.5*cos(b2)*sin(a)));

			addTriangle(Point3(0.5*cos(b2)*cos(a2), -0.5*sin(b2), 0.5*cos(b2)*sin(a2)), Point3(0.5*cos(b2)*cos(a), -0.5*sin(b2), 0.5*cos(b2)*sin(a)), Point3(0.5*cos(b)*cos(a2), -0.5*sin(b), 0.5*cos(b)*sin(a2)), 
				Vector3(0.5*cos(b2)*cos(a2), -0.5*sin(b2), 0.5*cos(b2)*sin(a2)), Vector3(0.5*cos(b2)*cos(a), -0.5*sin(b2), 0.5*cos(b2)*sin(a)), Vector3(0.5*cos(b)*cos(a2), -0.5*sin(b), 0.5*cos(b)*sin(a2)));

		}	
	}

}

// sphere intersection method
HitRecord Sphere::intersect(const Point3& o, const Vector3& d) {
	HitRecord hr;

	// solve a*t^2 + b*t + c = 0
	double a = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
	double b = 2.0*( d[0]*o[0] + d[1]*o[1] + d[2]*o[2] );
	double c = o[0]*o[0] + o[1]*o[1] + o[2]*o[2] - 0.25;

	double dis = b*b - 4.0*a*c;
	if (dis >= 0) {
		dis = sqrt(dis);
		double t1 = ( dis - b) / (2.0*a);
		double t2 = (-dis - b) / (2.0*a);

		Point3 p;
		Vector3 n;
		if (t1 > 0.0) {
			p = o + t1*d;
			if ((p[1] <= 0.5) && (p[1] >= -0.5)) {
				n = Vector3(p[0], p[1], p[2]);
				n.normalize();
				hr.addHit(t1, 0.0, 0.0, p, n);
			}
		}
		if (t2 > 0.0) {  // else
			p = o + t2*d;
			if ((p[1] <= 0.5) && (p[1] >= -0.5)) {
				n = Vector3(p[0], p[1], p[2]);
				n.normalize();
				hr.addHit(t2, 0.0, 0.0, p, n);
			}
		}
	}

	hr.sortHits();
	return hr;
}



Cube::Cube(int n) : Shape() {
	double z1 = -0.5; // left
	double x1 = 0.5; // close
	for (int i = 0; i < n+1; i++) {
		x1 = 0.5;
		double z2 = -0.5 + (double) i* 1.0 / (double) n; // right
		for (int j = 0; j < n+1; j++) {
			double x2 = 0.5 - (double) j * 1.0 / (double) n; // far
			// front face
			addTriangle(Point3(x1, -0.5, z1), Point3(x1, -0.5, z2), Point3(x2, -0.5, z1), Vector3(0.0, -1, 0.0), Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0));
			addTriangle(Point3(x1, -0.5, z2), Point3(x2, -0.5, z2), Point3(x2, -0.5, z1), Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0));
			// back face
			addTriangle(Point3(x1, 0.5, z1), Point3(x2, 0.5, z1), Point3(x1, 0.5, z2), Vector3(0.0, 1.0, 0.0), Vector3(0.0, 1.0, 0.0), Vector3(0.0, 1.0, 0.0));
			addTriangle(Point3(x1, 0.5, z2), Point3(x2, 0.5, z1), Point3(x2, 0.5, z2), Vector3(0.0, 1.0, 0.0), Vector3(0.0, 1.0, 0.0), Vector3(0.0, 1.0, 0.0));
			// left face
			addTriangle(Point3(x1, z2, -0.5), Point3(x1, z1, -0.5), Point3(x2, z1, -0.5), Vector3(0.0, 0.0, -1.0), Vector3(0.0, 0.0, -1.0), Vector3(0.0, 0.0, -1.0));
			addTriangle(Point3(x2, z2, -0.5), Point3(x1, z2, -0.5), Point3(x2, z1, -0.5), Vector3(0.0, 0.0, -1.0), Vector3(0.0, 0.0, -1.0), Vector3(0.0, 0.0, -1.0));
			// right face
			addTriangle(Point3(x1, z1, 0.5), Point3(x1, z2, 0.5), Point3(x2, z1, 0.5), Vector3(0.0, 0.0, 1.0), Vector3(0.0, 0.0, 1.0), Vector3(0.0, 0.0, 1.0));
			addTriangle(Point3(x1, z2, 0.5), Point3(x2, z2, 0.5), Point3(x2, z1, 0.5), Vector3(0.0, 0.0, 1.0), Vector3(0.0, 0.0, 1.0), Vector3(0.0, 0.0, 1.0));
			// top face
			addTriangle(Point3(-0.5, x1, z2), Point3(-0.5, x1, z1), Point3(-0.5, x2, z1), Vector3(-1.0, 0.0, 0.0), Vector3(-1.0, 0.0, 0.0), Vector3(-1.0, 0.0, 0.0));
			addTriangle(Point3(-0.5, x2, z2), Point3(-0.5, x1, z2), Point3(-0.5, x2, z1), Vector3(-1.0, 0.0, 0.0), Vector3(-1.0, 0.0, 0.0),  Vector3(-1.0, 0.0, 0.0));
			// bottom face
			addTriangle(Point3(0.5, x1, z1), Point3(0.5, x1, z2), Point3(0.5, x2, z1), Vector3(1.0, 0.0, 0.0), Vector3(1.0, 0.0, 0.0), Vector3(1.0, 0.0, 0.0));
			addTriangle(Point3(0.5, x1, z2), Point3(0.5, x2, z2), Point3(0.5, x2, z1), Vector3(1.0, 0.0, 0.0), Vector3(1.0, 0.0, 0.0),  Vector3(1.0, 0.0, 0.0));
			x1 = x2;
		}
		z1 = z2; // move left
	}


}

// cube intersection method
HitRecord Cube::intersect(const Point3& o, const Vector3& d) {
	HitRecord hr;
	//n = (1, 0, 0)
	double t1 = (0.5-o[0])/d[0];
	//n = (0, 1, 0)
	double t2 = (0.5-o[1])/d[1];
	//n = (0, 0, 1)
	double t3 = (0.5-o[2])/d[2];
	//n = (-1, 0, 0)
	double t4 = (-0.5-o[0])/d[0];
	//n = (0, -1, 0)
	double t5 = (-0.5-o[1])/d[1];
	//n = (0, 0, -1)
	double t6 = (-0.5-o[2])/d[2];

	Point3 p;
	Vector3 n;
	p = o + t1*d;
	if (d[0] != 0.0 && t1 > 0.0 && p[1] < 0.5 && p[1] > -0.5 && p[2] < 0.5 && p[2] > -0.5){

		n = Vector3(1, 0, 0);
		n.normalize();
		hr.addHit(t1, 0.0, 0.0, p, n);

	}

	p = o + t4*d;
	if (d[0] != 0.0 && t4 > 0.0 && p[1] < 0.5 && p[1] > -0.5 && p[2] < 0.5 && p[2] > -0.5){

		n = Vector3(-1, 0, 0);
		n.normalize();
		hr.addHit(t4, 0.0, 0.0, p, n);

	}

	Point3 p2;
	Vector3 n2;
	p2 = o + t2*d;
	if (d[1] != 0.0 && t2 > 0.0 && p2[0] < 0.5 && p2[0] > -0.5 && p2[2] < 0.5 && p2[2] > -0.5){

		n2 = Vector3(0, 1, 0);
		n2.normalize();
		hr.addHit(t2, 0.0, 0.0, p2, n2);

	}

	p2 = o + t5*d;
	if (d[1] != 0.0 && t5 > 0.0 && p2[0] < 0.5 && p2[0] > -0.5 && p2[2] < 0.5 && p2[2] > -0.5){

		n2 = Vector3(0, -1, 0);
		n2.normalize();
		hr.addHit(t5, 0.0, 0.0, p2, n2);

	}

	Point3 p3;
	Vector3 n3;
	p3 = o + t3*d;
	if (d[2] != 0.0 && t3 > 0.0 && p3[1] < 0.5 && p3[1] > -0.5 && p3[0] < 0.5 && p3[0] > -0.5){

		n3 = Vector3(0, 0, 1);
		n3.normalize();
		hr.addHit(t3, 0.0, 0.0, p3, n3);

	}

	p3 = o + t6*d;
	if (d[2] != 0.0 && t6 > 0.0 && p3[1] < 0.5 && p3[1] > -0.5 && p3[0] < 0.5 && p3[0] > -0.5){

		n3 = Vector3(0, 0, -1);
		n3.normalize();
		hr.addHit(t6, 0.0, 0.0, p3, n3);

	}


	hr.sortHits();
	return hr;
}


Special1::Special1(int n, int m) :Shape() {	
	double ca = 0.5;
	double sa = 0.0;
	if (n < 3) {
		n = 3;
	}
	for (int i = 1; i <= n; i++) {
		double a = (2 * M_PI * i) / n;
		double ca1 = 0.5 * cos(a);
		double sa1 = 0.5 * sin(a);


		addTriangle(Point3(0.0, -0.5, 0.0), Point3(ca, -0.5, sa), Point3(ca1, -0.5, sa1), 
			Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0));

		double ca2 = ca;
		double sa2 = sa;
		double ca3 = ca1;
		double sa3 = sa1;

		for (int j = 1; j <= m; j++) {
			double d = j * 1.0 / m;		
			double ca4 = (1.0-d) / 1.0 * ca2; // upper x
			double sa4 = (1.0-d) / 1.0 * sa2; // upper y
			double ca5 = (1.0-d) / 1.0 * ca3; // upper x
			double sa5 = (1.0-d) / 1.0 * sa3; // upper y
			double h1 = -0.5 + (j - 1) *1.0/m;
			double h2 = -0.5 + d;

			addTriangle(Point3(ca3, h1, sa3), Point3(ca2, h1, sa2), Point3(ca5, h2, sa5),
				Vector3(ca3, 0.5*sqrt(ca3*ca3 + sa3*sa3), sa3), Vector3(ca2, 0.5*sqrt(ca2*ca2 + sa2*sa2), sa2), Vector3(ca3, 0.5*sqrt(ca3*ca3 + sa3*sa3), sa3));

			addTriangle(Point3(ca4, h2, sa4), Point3(ca5, h2, sa5), Point3(ca2, h1, sa2),
				Vector3(ca4, 0.5*sqrt(ca4*ca4 + sa4*sa4), sa4), Vector3(ca5, 0.5*sqrt(ca5*ca5 + sa5*sa5), sa5), Vector3(ca2, 0.5*sqrt(ca2*ca2 + sa2*sa2), sa2));
			ca2 = ca4;
			sa2 = sa4;
			ca3 = ca5;
			sa3 = sa5;

		}

		ca = ca1;
		sa = sa1;

	}
}


Special2::Special2(int n, int m) :Shape() {	
	double ca = 0.1;
	double sa = 0.0;
	for (int i = 1; i <= n; i++) {
		double a = (2 * M_PI * i) / n;
		double ca1 = 0.1 * cos(a);
		double sa1 = 0.1 * sin(a);

		//top
		addTriangle(Point3(0.0, 1.0, 0.0), Point3(ca1, 0.5, sa1), Point3(ca, 0.5, sa), 
			Vector3(0.5*ca, 0.5*(abs(ca)+abs(sa)), 0.5*sa), Vector3(ca, 0.5*(abs(ca)+abs(sa)), sa), Vector3(1.5*ca, 0.5*(abs(ca)+abs(sa)), 1.5*sa));
		//bottom
		addTriangle(Point3(0.0, -0.5, 0.0), Point3(ca, -0.5, sa), Point3(ca1, -0.5, sa1), 
			Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0), Vector3(0.0, -1.0, 0.0));


		//body
		double cb = m;
		for (int j = 1; j <= m; j++) {
			double b = 1.0 / m;
			double cb1 = 0.5 - j * b;
			double cb2 = 0.5 - (j-1) * b;

			addTriangle(Point3(ca1, cb2, sa1), Point3(ca, cb1, sa), Point3(ca, cb2, sa), 
				Vector3(ca1, 0.0, sa1), Vector3(ca, 0.0, sa), Vector3(ca, 0.0, sa));
			addTriangle(Point3(ca1, cb2, sa1), Point3(ca1, cb1, sa1), Point3(ca, cb1, sa),
				Vector3(ca1, 0.0, sa1),  Vector3(ca1, 0.0, sa1), Vector3(ca, 0.0, sa));

		}

		ca = ca1;
		sa = sa1;
	}
}