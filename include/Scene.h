#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Matrix4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<std::vector<double> > depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	explicit Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName);
	void forwardRenderingPipeline(Camera *camera);

    Matrix4 calculateModelingTransformationMatrix(Mesh* mesh);
    static Matrix4 calculateCameraTransformationMatrix(Camera *camera);
    static Matrix4 calculateOrthographicTransformationMatrix(Camera *camera);
    static Matrix4 calculatePerspectiveTransformationMatrix(Camera *camera);
    static Matrix4 calculateViewportTransformationMatrix(Camera *camera);

    void rasterizeSolid(Camera* camera, const Vec4& pointA, const Vec4& pointB, const Vec4& pointC);

    void rasterizeWireframe(Camera *camera, const Vec4 &pointA, const Vec4 &pointB);
};

#endif
