#include "../csc321.h"
#include "../sceneview/MyScene.h"
#include "RenderingInterface.h"
#include <FL/gl.h>
#include <cfloat>

const int RECURSIVE_LIMIT = 5;

void MyScene::render(int type, int width, int height, unsigned char* pixels) {
	if (!isLoaded) {
		progress = 1.0;
		return;
	}
	bDoRender = true;
	// Add your rendering code here.
	// Keep track of your progress as a value between 0 and 1
	// so the progress bar can update as the rendering progresses
	progress = 0.0;
	switch (type) {
	case RenderingUI::RENDER_SCANLINE:  scanline(width, height, pixels); break;
	case RenderingUI::RENDER_RAY_TRACING:  raytrace(width, height, pixels); break;
	case RenderingUI::RENDER_PATH_TRACING:  break;
	default: break;
	}
	progress = 1.0;
}

void MyScene::stopRender()
{
	// Because this is threaded code, this function
	// can be called in the middle of your rendering code.
	// You should then stop at the next scanline
	bDoRender = false;
}

double MyScene::getRenderProgress() {
	// return the current progress as a value between 0 and 1
	return progress;
}

// add extra methods here
void MyScene::scanline(int w, int h, unsigned char* pixels) {
	resize(w, h);

	glPushAttrib( GL_ALL_ATTRIB_BITS );
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);

	glClearColor( background[0], background[1], background[2], 1.0f);

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd( &getCamera().getProjection()(0,0) );
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd( &getCamera().getWorldToCamera()(0,0) );

	glEnable(GL_LIGHTING);
	glPolygonMode(GL_FRONT, GL_FILL);

	draw();

	glPopAttrib();
	glReadPixels( 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &pixels[0] );
	progress = 1.0;
}

void MyScene::raytrace(int w, int h, unsigned char* pixels) {
	resize(w, h);

	for (int y=0; y < h; y++) {
		if (bDoRender == false)
			break;

		for (int x = 0; x < w; x++) {
			// determine the color of the pixel (x,y) by raytracing

			// form the ray
			Point3 pixel (x, y, -1.0);
			pixel[0] = -1.0 + 2.0 * pixel[0] / (camera.getWidth() - 1);
			pixel[1] = -1.0 + 2.0 * pixel[1] / (camera.getHeight() - 1);
			pixel = camera.getCameraToWorld() * pixel;
			Vector3 dir = pixel - camera.getEye();
			Point3 o = camera.getEye();

			// trace the ray
			Color c = traceRay(o, dir, RECURSIVE_LIMIT);

			// clamp and store the color value
			c[0] = (c[0] > 0.0) ? ((c[0] < 1.0) ? c[0] : 1.0) : 0.0;
			c[1] = (c[1] > 0.0) ? ((c[1] < 1.0) ? c[1] : 1.0) : 0.0;
			c[2] = (c[2] > 0.0) ? ((c[2] < 1.0) ? c[2] : 1.0) : 0.0;
			*pixels++ = (unsigned char) (c[0] * 255.0);
			*pixels++ = (unsigned char) (c[1] * 255.0);
			*pixels++ = (unsigned char) (c[2] * 255.0);
		}

		progress = (double) y / (double) h;
		Fl::check();
	}

	progress = 1.0;
}

Color MyScene::traceRay(Point3& o, Vector3& dir, int depth) {
	FirstHitRecord fhr = masters->get("root")->intersect(o, dir);
	if (!fhr.hit()) {
		return background;
	}

	Vector3 L, V, R; //
	Color c = fhr.node->object->ambient;// * ambientLight; 
	fhr.n.normalize();

	//reflect
	V = (-dir*fhr.n) * fhr.n;
	R = 2*V + dir;
	R.normalize();
	if (depth > 0 && fhr.node->object->reflect.getMax() > 0.0) {
		 c += fhr.node->object->reflect * traceRay(fhr.p+R*0.0001, R, depth-1);	
	} 

	//transparent
	if (depth > 1 && fhr.node->object->transparent.getMax() > 0.0) {
		c += fhr.node->object->transparent * traceRay(fhr.p+R*0.0001, dir, depth-1);
	} 

	for (int i = 0; i < lights.size(); i++) {
		Color cli;  // color contribution of light i
		Vector3 l, v, r;  // light vector pointing from point to source
		double ray;
		l = lights[i].getPos() - fhr.p;
		double lightDistance = l.length();
		//l.normalize();

		// shadow
		FirstHitRecord fhr2 = masters->get("root")->intersect(fhr.p+ 0.00001*l, l);
		if (!fhr2.hit()|| (fhr.p-fhr2.p).lengthSquared() > l.lengthSquared()) {    
			// diffuse 
			// check if in front
			l.normalize();
			double d = fhr.n * l;
			if ( d <= 0) {
				continue; //continue;
			}
			cli += fhr.node->object->diffuse * d;


			//specular and shine

			v = (-l*fhr.n) * fhr.n;
			r = 2*v +l;
			r.normalize();
			ray = dir * r;
			if (ray <=0) {
				ray = 0;
			}
			cli += fhr.node->object->specular * pow(ray, fhr.node->object->shine);


			//attenuation
			double fatt;
			Point3 falloff = lights[i].getFalloff();
			fatt = 1.0/(falloff[0]+falloff[1]*lightDistance+falloff[2]*lightDistance*lightDistance);

			// include this light's contribution in total color
			c += cli * lights[i].getColor() * fatt; 
		}
	}

	return c;

}