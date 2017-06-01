/*========================================================================
* COSC 363  Computer Graphics (2017)
* Ray tracer 
* See Lab07.pdf for details.
*=========================================================================
*/
#include <iostream>
#include "stdio.h"
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include "Sphere.h"
#include "SceneObject.h"
#include "Ray.h"
#include "Plane.h"
#include "TextureBMP.h"
#include <GL/glut.h>

using namespace std;

TextureBMP tex1;
TextureBMP tex2;
const float WIDTH = 20.0;  
const float HEIGHT = 20.0;
const float EDIST = 40.0;
const int NUMDIV = 500;
const int MAX_STEPS = 5;
const float XMIN = -WIDTH * 0.5;
const float XMAX =  WIDTH * 0.5;
const float YMIN = -HEIGHT * 0.5;
const float YMAX =  HEIGHT * 0.5;


vector<SceneObject*> sceneObjects;  //A global list containing pointers to objects in the scene

glm::vec3 spherePattern(Ray ray){
	if ((int) (ray.xpt.y) % 2 == 0){
		return glm::vec3(0, 0, 1);
	}else {
		return sceneObjects[ray.xindex]->getColor();
	}
}

glm::vec3 planePattern(Ray ray){
	if ((int) (ray.xpt.z) % 2 == 0 && (int) (ray.xpt.x) % 2 == 0){
		return glm::vec3(1, 1, 1);
	
	}else {
		return sceneObjects[ray.xindex]->getColor();
	}
}


//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{
	glm::vec3 backgroundCol(0);
	glm::vec3 light(80, 200, 200);
	glm::vec3 light2(-80, 200, 20);
	float ambientTerm = 0.2;
    ray.closestPt(sceneObjects);		//Compute the closest point of intersetion of objects with the ray

    if(ray.xindex == -1) return backgroundCol;      //If there is no intersection return background colour

    glm::vec3 col = sceneObjects[ray.xindex]->getColor(); //else return object's colour
    glm::vec3 normalVector = sceneObjects[ray.xindex]->normal(ray.xpt);
    glm::vec3 lightVector = light - ray.xpt;
    glm::vec3 lightVector2 = light2 - ray.xpt;
    glm::vec3 lightNormal = glm::normalize(lightVector);
    glm::vec3 lightNormal2 = glm::normalize(lightVector2);
    float lDotn = glm::dot(normalVector, lightNormal);
    float lDotn2 = glm::dot(normalVector, lightNormal2);

	if (ray.xindex == 9)
	{
		col = spherePattern(ray);
	}
	

	if (ray.xindex == 7)
	{
		col = planePattern(ray);
	}
	
	
	glm::vec3 reflVector = glm::reflect(-lightNormal, normalVector);
	glm::vec3 reflVector2 = glm::reflect(-lightNormal2, normalVector);
	
	

	glm::vec3 spec(dot(reflVector, -ray.dir)); 
	glm::vec3 spec2(dot(reflVector2, -ray.dir)); 
	

	reflVector = glm::normalize(reflVector);
	
	if (ray.xindex == 1){
		float ucoord = asin(normalVector.x)/M_PI + 0.5;
		float vcoord = asin(normalVector.y)/M_PI + 0.5;
		col = tex1.getColorAt(ucoord, vcoord);
	}
	if(ray.xindex == 8)
	{
		float u = (ray.xpt.x + 400)/(800);
		float v = (ray.xpt.y + 250)/(500);
		col = tex2.getColorAt(u, v);
	}
	
	//~ float c = glm::dot(reflVector, normalVector);
	Ray shadow(ray.xpt, lightNormal);
	Ray shadow2(ray.xpt, lightNormal2);
	shadow.closestPt(sceneObjects);
	shadow2.closestPt(sceneObjects);
	
	float d = glm::length(light);
	float d2 = glm::length(light2);
	glm::vec3 colorSum;
	
	if ((lDotn < 0 || (shadow.xindex>-1 && shadow.xdist<d)) || (lDotn2 < 0 || (shadow2.xindex>-1 && shadow2.xdist<d2))){
		if ((shadow.xindex == 2 && shadow.xdist<d) || (shadow2.xindex == 2 && shadow.xdist<d2))
		{
			ambientTerm = 0.6;
		}
		if ((shadow.xindex == 10 && shadow.xdist<d) || (shadow2.xindex == 10 && shadow.xdist<d2))
		{
			ambientTerm = 0.6;
		}
		colorSum = ambientTerm * col;
	} else {
		float spec1 = 100; 
		float specTerm = glm::pow(dot(reflVector,-ray.dir), spec1);
		float specTerm2 = glm::pow(dot(reflVector2,-ray.dir), spec1);
		colorSum = ambientTerm*col + lDotn*col + specTerm + lDotn2*col + specTerm2;
	}
	
	if(ray.xindex == 0 && step < MAX_STEPS) {
		//~ glm::vec3 reflectedDir = glm::refract(reflVector, -ray.xpt, 0.99f);
		glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
		Ray reflectedRay(ray.xpt, reflectedDir);
		glm::vec3 reflectedCol = trace(reflectedRay, step+1); //Recursion!
		colorSum = colorSum + (0.8f*reflectedCol);
	 }
	 //shiettty transparency
	if(ray.xindex == 2 && step < MAX_STEPS && ray.xindex > -1) {
		//~ reflVector = glm::refract(ray.xpt, normalVector, 0.99f);
		//~ reflVector2 = glm::refract(ray.xpt, normalVector, 0.9f);
		float eta = 1.0f;
		glm::vec3 g = glm::refract(ray.dir, normalVector, eta);
		//~ glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
		Ray refractedRay(ray.xpt, g);
		refractedRay.closestPt(sceneObjects);
		if (refractedRay.xindex != -1) {
		glm::vec3 m = sceneObjects[refractedRay.xindex]->normal(refractedRay.xpt);
		glm::vec3 h = glm::refract(g, -m, 1.0f/eta);
		Ray reflected(refractedRay.xpt, h);
		colorSum = (0.2f*colorSum) + trace(reflected, step+1); //Recursion!
		}
	}
	 //shiettty refraction
	if(ray.xindex == 10 && step < MAX_STEPS) {
		//~ reflVector = glm::refract(ray.xpt, normalVector, 0.99f);
		//~ reflVector2 = glm::refract(ray.xpt, normalVector, 0.9f);
		float eta = 1.0f/1.5f;
		glm::vec3 g = glm::refract(ray.dir, normalVector, eta);
		//~ glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
		Ray refractedRay(ray.xpt, g);
		refractedRay.closestPt(sceneObjects);
		
		glm::vec3 m = sceneObjects[refractedRay.xindex]->normal(refractedRay.xpt);
		glm::vec3 h = glm::refract(g, -m, 1.0f/eta);
		Ray reflected(refractedRay.xpt, h);
		return  colorSum + trace(reflected, step+1); //Recursion!
		
	 }
	 
	return colorSum;
}


//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display()
{
	float xp, yp;  //grid point
	float cellX = (XMAX-XMIN)/NUMDIV;  //cell width
	float cellY = (YMAX-YMIN)/NUMDIV;  //cell height

	glm::vec3 eye(0., 0., 0.);  //The eye position (source of primary rays) is the origin

	glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	glBegin(GL_QUADS);  //Each cell is a quad.
	


	for(int i = 0; i < NUMDIV; i++)  	//For each grid point xp, yp
	{
		xp = XMIN + i*cellX;
		for(int j = 0; j < NUMDIV; j++)
		{
			yp = YMIN + j*cellY;

		    //~ glm::vec3 dir(xp+0.5*cellX, yp+0.5*cellY, -EDIST);	//direction of the primary ray
		    glm::vec3 dir1(xp+0.25*cellX, yp+0.25*cellY, -EDIST);	//direction of the primary ray
		    glm::vec3 dir2(xp+0.75*cellX, yp+0.25*cellY, -EDIST);	//direction of the primary ray
		    glm::vec3 dir3(xp+0.75*cellX, yp+0.75*cellY, -EDIST);	//direction of the primary ray
		    glm::vec3 dir4(xp+0.25*cellX, yp+0.75*cellY, -EDIST);	//direction of the primary ray

		    //~ Ray ray = Ray(eye, dir);		//Create a ray originating from the camera in the direction 'dir'
		    Ray ray1 = Ray(eye, dir1);		//Create a ray originating from the camera in the direction 'dir'
		    Ray ray2 = Ray(eye, dir2);		//Create a ray originating from the camera in the direction 'dir'
		    Ray ray3 = Ray(eye, dir3);		//Create a ray originating from the camera in the direction 'dir'
		    Ray ray4 = Ray(eye, dir4);		//Create a ray originating from the camera in the direction 'dir'
			//~ ray.normalize();				//Normalize the direction of the ray to a unit vector
			ray1.normalize();				//Normalize the direction of the ray to a unit vector
			ray2.normalize();				//Normalize the direction of the ray to a unit vector
			ray3.normalize();				//Normalize the direction of the ray to a unit vector
			ray4.normalize();				//Normalize the direction of the ray to a unit vector
		    //~ glm::vec3 col = trace (ray, 1); //Trace the primary ray and get the colour value
		    glm::vec3 col1 = trace (ray1, 1); //Trace the primary ray and get the colour value
		    glm::vec3 col2 = trace (ray2, 1); //Trace the primary ray and get the colour value
		    glm::vec3 col3 = trace (ray3, 1); //Trace the primary ray and get the colour value
		    glm::vec3 col4 = trace (ray4, 1); //Trace the primary ray and get the colour value
		    
			
			//~ glColor3f(col.r, col.g, col.b);
			glColor3f((col1.r + col2.r + col3.r + col4.r)/4, (col1.g + col2.g + col3.g + col4.g)/4, (col1.b + col2.b + col3.b + col4.b)/4);
			glVertex2f(xp, yp);				//Draw each cell with its color value
			glVertex2f(xp+cellX, yp);
			glVertex2f(xp+cellX, yp+cellY);
			glVertex2f(xp, yp+cellY);
			
        }
    }

    glEnd();
    glFlush();
}


//---This function initializes the scene ------------------------------------------- 
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize()
{

    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
    glClearColor(0, 0, 0, 1);

//------------------------Spheres---------------------------------------
    Sphere *sphere1 = new Sphere(glm::vec3(-5.0, -5.0, -90.0), 15.0, glm::vec3(0, 0, 1));
	Sphere *patternedSpehere = new Sphere(glm::vec3(-14.0, -14.0, -75.0), 3.5, glm::vec3(75.0f/255.0f, 0.0f/255.0f, 130.0f/255.0f));
    Sphere *sphere2 = new Sphere(glm::vec3(1.0, 1.0, -60.0), 2.0, glm::vec3(0, 1, 0));
    Sphere *sphere3 = new Sphere(glm::vec3(3.0, -10.0, -58.0), 3.0, glm::vec3(1, 0.960, 0.101));
    Sphere *sphere4 = new Sphere(glm::vec3(7.0, 0.0, -58.0), 3.0, glm::vec3(0, 0, 0));
//------------------------Spheres------------------------------------------


//------------------------Plane-----------------------------------------
	Plane *plane = new Plane (glm::vec3(-400., -25, 400), //Point A
	 glm::vec3(400., -25, 400), //Point B
	 glm::vec3(400., -25, -2500), //Point C
	 glm::vec3(-400., -25, -2500), //Point D
	 glm::vec3(0.0, 0.0, 0)); //Colour
	 
	Plane *background = new Plane (glm::vec3(-400., -250, -1000), //Point A
	 glm::vec3(400., -250, -1000), //Point B
	 glm::vec3(400., 250, -1000), //Point C
	 glm::vec3(-400., 250, -1000), //Point D
	 glm::vec3(0.5, 0.5, 0.5)); //Colour
//------------------------Plane-----------------------------------------
	//https://static.pexels.com/photos/269888/pexels-photo-269888.jpeg
    tex1 = TextureBMP("earth_.bmp");
    tex2 = TextureBMP("background.bmp");

//------------------------Tetra------------------------------------------
	glm::vec3 a =  glm::vec3((3.)+10, -25, (-30)-90); // Front top left 
	glm::vec3 b = glm::vec3((13.)+10, -25, (-20)-90); // Front top right 
	glm::vec3 c = glm::vec3((8.)+10, -25, (-15)-90); // Front bottom left 
	glm::vec3 d = glm::vec3((23.0f/3.0f)+10, -15, (-65.0f/3.0f)-90); // Front bottom right 
	glm::vec3 ab = glm::vec3((5.)+10, -25, (-25)-90); // Back top left 
	glm::vec3 cd = glm::vec3((15.0f/2.0f)+10, -20.0f , (-18.2f)-90); // Back bottom left 
	a = a * 0.2f;
	b = b * 0.2f;
	c = c * 0.2f;
	d = d * 0.2f;
	ab = ab * 0.2f;
	cd = cd * 0.2f;
	
	Plane *triBot = new Plane (a, //Point A
	 ab, //Point B
	 b, //Point C
	 c, //Point D
	 glm::vec3(0.952, 0.223, 0.286)); //Colour
	 
	Plane *triSide1 = new Plane (a, //Point A
	 ab, //Point B
	 b, //Point C
	 d, //Point D
	 glm::vec3(0.952, 0.223, 0.286)); //Colour
	 
	Plane *triSide2 = new Plane (a, //Point A
	 d, //Point B
	 cd, //Point C
	 c, //Point D
	 glm::vec3(0.952, 0.223, 0.286)); //Colour
	 
	Plane *triSide3 = new Plane (b, //Point A
	 d, //Point B
	 cd, //Point C
	 c, //Point D
	 glm::vec3(0.952, 0.223, 0.286)); //Colour
	 
//------------------------Tetra------------------------------------------

//------------------------Cube------------------------------------------
	glm::vec3 frontTopLeft =  glm::vec3(8., 5, -55); // Front top left 
	glm::vec3 frontTopRight = glm::vec3(13., 5, -55); // Front top right 
	glm::vec3 frontBotLeft = glm::vec3(8., 10, -55); // Front bottom left 
	glm::vec3 frontBotRight = glm::vec3(13., 10, -55); // Front bottom right 
	glm::vec3 backTopLeft = glm::vec3(8., 5, -65); // Back top left 
	glm::vec3 backTopRight = glm::vec3(13., 5, -65); // Back top right 
	glm::vec3 backBottomLeft = glm::vec3(8., 10, -65); // Back bottom left 
	glm::vec3 backBottomRight = glm::vec3(13., 10, -65); // Back bottom right 
	 
	Plane *front = new Plane (frontTopLeft, //Point A
	 frontTopRight, //Point B
	 frontBotRight, //Point C
	 frontBotLeft, //Point D
	 glm::vec3(0.419, 0.8, 0.380)); //Colour
	 
	Plane *top = new Plane (frontTopLeft, //Point A
	 frontTopRight, //Point B
	 backTopRight, //Point C
	 backTopLeft, //Point D
	 glm::vec3(0.419, 0.8, 0.380)); //Colour
	 
	Plane *bottom = new Plane (frontBotLeft, //Point A
	 backTopRight, //Point B
	 backBottomRight, //Point C
	 backBottomLeft, //Point D
	 glm::vec3(0.419, 0.8, 0.380)); //Colour
	 
	Plane *back = new Plane (backTopLeft, //Point A
	 backTopRight, //Point B
	 backBottomRight, //Point C
	 backBottomLeft, //Point D
	 glm::vec3(0.419, 0.8, 0.380)); //Colour
	 
	Plane *right = new Plane (frontTopRight, //Point A
	 backTopRight, //Point B
	 backBottomRight, //Point C
	 frontBotRight, //Point D
	 glm::vec3(0.419, 0.8, 0.380)); //Colour

	Plane *left = new Plane (frontTopLeft, //Point A
	 backTopLeft, //Point B
	 backBottomLeft, //Point C
	 frontBotLeft, //Point D
	 glm::vec3(0.419, 0.8, 0.380)); //Colour
//------------------------Cube------------------------------------------



	//--Add the above to the list of scene objects.
	sceneObjects.push_back(sphere1); //0
	sceneObjects.push_back(sphere2); //1
	sceneObjects.push_back(sphere3); //2
	sceneObjects.push_back(front); //3
	sceneObjects.push_back(top); //4
	sceneObjects.push_back(right); //5
	sceneObjects.push_back(left); //6
	sceneObjects.push_back(plane); //7
	sceneObjects.push_back(background); //8
	sceneObjects.push_back(patternedSpehere); //9
	sceneObjects.push_back(sphere4); //10
	sceneObjects.push_back(bottom); //11
	sceneObjects.push_back(back); //12
	sceneObjects.push_back(triBot); //12
	sceneObjects.push_back(triSide1); //12
	sceneObjects.push_back(triSide2); //12
	sceneObjects.push_back(triSide3); //12
	
}



int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Raytracer");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}
