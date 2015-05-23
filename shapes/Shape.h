#ifndef _SHAPE_H_
#define _SHAPE_H_

#include "../csc321.h"
#include "../intersection/HitRecord.h"
#include <vector>
#include "Point3.h"
#include "Vector3.h"

class Shape {
public:
	virtual ~Shape();
	virtual void draw();
	virtual HitRecord intersect(const Point3&, const Vector3&);
protected:
	Shape();
	void addTriangle(const Point3& p1, const Point3& p2, const Point3& p3,
					 const Vector3& n1, const Vector3& n2, const Vector3& n3);
	std::vector<Point3> vertices;
	std::vector<Vector3> normals;
};

class Cylinder : public Shape {
public:
	Cylinder(int, int);
	HitRecord intersect(const Point3&, const Vector3&);
};

class Cone : public Shape {
public:
	Cone(int, int);
	HitRecord intersect(const Point3&, const Vector3&);
};

class Sphere : public Shape {
public:
	Sphere(int, int);
	HitRecord intersect(const Point3&, const Vector3&);
};

class Cube : public Shape {
public:
	Cube(int);
	HitRecord intersect(const Point3&, const Vector3&);
};

class Special1 : public Shape {
public:
	Special1(int, int);
};

class Special2 : public Shape {
public:
	Special2(int, int);
};


#endif /* _SHAPE_H_ */