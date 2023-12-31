#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>

#include "../include/tinyxml2.h"
#include "../include/Triangle.h"
#include "../include/Helpers.h"
#include "../include/Scene.h"

// Outcodes for Cohen-Sutherland algorithm
#define INSIDE 0 // 0000
#define LEFT 1   // 0001
#define RIGHT 2  // 0010
#define BOTTOM 4 // 0100
#define TOP 8    // 1000

using namespace tinyxml2;
using namespace std;

// Parses XML file
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

// Initializes image with background color
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

// If given value is less than 0, converts value to 0.
// If given value is more than 255, converts value to 255.
// Otherwise, returns value itself.
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

// Writes contents of image (Color**) into a PPM file.
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

// Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	command = "convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}

#pragma region Transformation
Matrix4 Scene::calculateModelingTransformationMatrix(Mesh* mesh){
    auto modelingTransformation = getIdentityMatrix();

    for (int i = 0; i < mesh->numberOfTransformations; ++i) {

        // Translation
        if (mesh->transformationTypes[i] == 't')
        {
            auto translation = translations[mesh->transformationIds[i] - 1];
            auto x = translation->tx;
            auto y = translation->ty;
            auto z = translation->tz;
            double translationMatrix[4][4] =
                                 {{1,0,0,x},
                                  {0,1,0,y},
                                  {0,0,1,z},
                                  {0,0,0,1}};
            modelingTransformation = multiplyMatrixWithMatrix(translationMatrix, modelingTransformation);
        }

        // Scaling
        if (mesh->transformationTypes[i] == 's')
        {
            auto scaling = scalings[mesh->transformationIds[i] - 1];
            auto x = scaling->sx;
            auto y = scaling->sy;
            auto z = scaling->sz;
            double scalingMatrix[4][4] =
                    {{x,0,0,0},
                     {0,y,0,0},
                     {0,0,z,0},
                     {0,0,0,1}};
            modelingTransformation = multiplyMatrixWithMatrix(scalingMatrix, modelingTransformation);
        }

        // Rotation
        if (mesh->transformationTypes[i] == 'r')
        {
            auto rotation = rotations[mesh->transformationIds[i] - 1];

            auto u = Vec3(rotation->ux, rotation->uy, rotation->uz, -1);
            auto m = std::min(std::min(rotation->ux, rotation->uy), abs(rotation->uz));
            auto v = normalizeVec3(((m == abs(rotation->ux)) ? Vec3(0, -1 * rotation->uz, rotation->uy, -1) :
                                    (m == abs(rotation->uy)) ? Vec3(-1 * rotation->uz, 0, rotation->ux, -1) :
                                    Vec3(-1 * rotation->uy, rotation->ux, 0, -1)));
            auto w = normalizeVec3(crossProductVec3(u, v));
            double mMatrix[4][4] =
                    {{u.x,  u.y,    u.z,    0},
                     {v.x,  v.y,    v.z,    0},
                     {w.x,  w.y,    w.z,    0},
                     {0,    0,      0,      1}};
            double mMatrix_inverse[4][4] =
                    {{u.x,  v.x,    w.x,    0},
                     {u.y,  v.y,    w.y,    0},
                     {u.z,  v.z,    w.z,    0},
                     {0,    0,      0,      1}};

            auto ra = rotation->angle * M_PI/180;
            double rotationMatrix[4][4] =
                    {{1,      0,          0,          0},
                     {0,      cos(ra),    -sin(ra),   0},
                     {0,      sin(ra),    cos(ra),    0},
                     {0,      0,          0,          1}};

            auto r1 = multiplyMatrixWithMatrix(rotationMatrix, mMatrix);
            auto r2 = multiplyMatrixWithMatrix(mMatrix_inverse, r1);
            modelingTransformation = multiplyMatrixWithMatrix(r2, modelingTransformation);
        }
    }

    return modelingTransformation;
}

Matrix4 Scene::calculateCameraTransformationMatrix(Camera* camera){
    auto x = camera->position.x;
    auto y = camera->position.y;
    auto z = camera->position.z;
    auto u = camera->u;
    auto v = camera->v;
    auto w = camera->w;

    double translate[4][4] = {
            {1, 0, 0, -x},
            {0, 1, 0, -y},
            {0, 0, 1, -z},
            {0, 0, 0, 1}
    };

    double rotate[4][4] = {
            {u.x,   u.y,    u.z,    0},
            {v.x,   v.y,    v.z,    0},
            {w.x,   w.y,    w.z,    0},
            {0,     0,      0,      1}
    };

    return multiplyMatrixWithMatrix(rotate, translate);
}

Matrix4 Scene::calculateOrthographicTransformationMatrix(Camera* camera){
    auto l = camera->left;
    auto r = camera->right;
    auto t = camera->top;
    auto b = camera->bottom;
    auto n = camera->near;
    auto f = camera->far;

    double orthographic[4][4] = {
            {2/(r - l),             0,                  0,                  -((r + l) / (r - l))},
            {0,                     2/(t - b),          0,                  -((t + b) / (t - b))},
            {0,                     0,                  -(2/(f - n)),       -((f + n) / (f - n))},
            {0,                     0,                  0,                  1                   }
    };

    return orthographic;
}

Matrix4 Scene::calculatePerspectiveTransformationMatrix(Camera* camera){
    auto l = camera->left;
    auto r = camera->right;
    auto t = camera->top;
    auto b = camera->bottom;
    auto n = camera->near;
    auto f = camera->far;

    double perspective[4][4] = {
            {(2*n) / (r - l),   0,                  (r + l) / (r - l),      0                       },
            {0,                 (2*n) / (t - b),    (t + b) / (t - b),      0                       },
            {0,                 0,                  -((f + n) / (f - n)),   -((2 * f * n) / (f - n))},
            {0,                 0,                  -1,                     0                       }
    };

    return perspective;
}

Matrix4 Scene::calculateViewportTransformationMatrix(Camera* camera){
    auto h = camera->horRes;
    auto v = camera->verRes;

    double viewport[4][4] = {
            {h/2.0,     0,          0,          (h-1)/2.0},
            {0,         v/2.0,      0,          (v-1)/2.0},
            {0,         0,          0.5,        0.5      },
            {0,         0,          0,          1        }
    };

    return viewport;
}
#pragma endregion

#pragma region Rasterization

// Compute the outcode for a point (x, y, z)
int computeOutcode(double x, double y, double z) {
    int code = INSIDE; // Initialize as inside

    double x_min = -1, y_min = -1, z_min = -1;

    double x_max = 1, y_max = 1, z_max = 1;

    if (x <= x_min)      // to the left of clip window
        code |= LEFT;
    else if (x >= x_max) // to the right of clip window
        code |= RIGHT;

    if (y <= y_min)      // below the clip window
        code |= BOTTOM;
    else if (y >= y_max) // above the clip window
        code |= TOP;

    if (z <= z_min)      // behind the near clipping plane
        code |= 16;    // Additional bit for z
    else if (z >= z_max) // beyond the far clipping plane
        code |= 32;    // Additional bit for z

    return code;
}

// Clip a line using Cohen-Sutherland algorithm
void clipLine(Vec4& line1, Color* c1, Vec4& line2, Color* c2) {
    // Compute outcodes for the two endpoints of the line
    int outcode1 = computeOutcode(line1.x, line1.y, line1.z);
    int outcode2 = computeOutcode(line2.x, line2.y, line2.z);

    // Initially assume both endpoints are inside the clip window
    bool accept = false;

    double x_min = -1, y_min = -1, z_min = -1;

    double x_max = 1, y_max = 1, z_max = 1;

    while (true) {
        // If both endpoints are inside the clip window, accept the line
        if ((outcode1 == INSIDE) && (outcode2 == INSIDE)) {
            accept = true;
            break;
        }
        else if (outcode1 & outcode2) {
            // If the logical AND is not 0, the line is completely outside the clip window
            break;
        }
        else {
            // Calculate intersection point
            double x, y, z;

            // Pick the endpoint outside the clip window
            int outcodeOut = outcode1 ? outcode1 : outcode2;

            // Find the intersection point
            if (outcodeOut & TOP) {
                x = line1.x + (line2.x - line1.x) * (y_max - line1.y) / (line2.y - line1.y);
                y = y_max;
                z = line1.z + (line2.z - line1.z) * (y_max - line1.y) / (line2.y - line1.y);
            }
            else if (outcodeOut & BOTTOM) {
                x = line1.x + (line2.x - line1.x) * (y_min - line1.y) / (line2.y - line1.y);
                y = y_min;
                z = line1.z + (line2.z - line1.z) * (y_min - line1.y) / (line2.y - line1.y);
            }
            else if (outcodeOut & RIGHT) {
                y = line1.y + (line2.y - line1.y) * (x_max - line1.x) / (line2.x - line1.x);
                x = x_max;
                z = line1.z + (line2.z - line1.z) * (x_max - line1.x) / (line2.x - line1.x);
            }
            else if (outcodeOut & LEFT) {
                y = line1.y + (line2.y - line1.y) * (x_min - line1.x) / (line2.x - line1.x);
                x = x_min;
                z = line1.z + (line2.z - line1.z) * (x_min - line1.x) / (line2.x - line1.x);
            }
            else if (outcodeOut & 16) {
                z = line1.z + (line2.z - line1.z) * (z_min - line1.y) / (line2.y - line1.y);
                x = line1.x + (line2.x - line1.x) * (z_min - line1.z) / (line2.z - line1.z);
                y = line1.y + (line2.y - line1.y) * (z_min - line1.z) / (line2.z - line1.z);
            }
            else if (outcodeOut & 32) {
                z = line1.z + (line2.z - line1.z) * (z_max - line1.y) / (line2.y - line1.y);
                x = line1.x + (line2.x - line1.x) * (z_max - line1.z) / (line2.z - line1.z);
                y = line1.y + (line2.y - line1.y) * (z_max - line1.z) / (line2.z - line1.z);
            }

            // Update the endpoint outside the clip window
            if (outcodeOut == outcode1) {
                line1.x = x;
                line1.y = y;
                line1.z = z;
                outcode1 = computeOutcode(line1.x, line1.y, line1.z);
            }
            else {
                line2.x = x;
                line2.y = y;
                line2.z = z;
                outcode2 = computeOutcode(line2.x, line2.y, line2.z);
            }
        }
    }
}

void Scene::rasterizeWireframe(Camera* camera, const Vec4& pointA, const Vec4& pointB) {
    // Calculate bounds for efficient rasterization
    auto x0 = (int)pointA.x;
    auto y0 = (int)pointA.y;
    auto x1 = (int)pointB.x;
    auto y1 = (int)pointB.y;

    // Determine differences and directions
    auto dx = abs(x1 - x0);
    auto dy = abs(y1 - y0);
    auto sx = x0 < x1 ? 1 : -1;
    auto sy = y0 < y1 ? 1 : -1;
    auto err = (dx > dy ? dx : -dy) / 2;
    auto e2 = err;

    // Interpolate colors between pointA and pointB
    auto c0 = colorsOfVertices[pointA.colorId - 1];
    auto c1 = colorsOfVertices[pointB.colorId - 1];

    // Iterate over the line pixels and interpolate color
    while (true) {
        // Interpolate color
        double t = sqrt((x0 - pointA.x) * (x0 - pointA.x) + (y0 - pointA.y) * (y0 - pointA.y)) /
                   sqrt((pointB.x - pointA.x) * (pointB.x - pointA.x) + (pointB.y - pointA.y) * (pointB.y - pointA.y));
        auto color = Color(
                (1 - t) * c0->r + t * c1->r,
                (1 - t) * c0->g + t * c1->g,
                (1 - t) * c0->b + t * c1->b
        );

        // Round color
        color = color.round();

        // Draw pixel with the interpolated color
        if (x0 >= 0 && y0 >= 0 && x0 < camera->horRes && y0 < camera->verRes)
            image[x0][y0] = color;

        if (x0 == x1 && y0 == y1)
            break;

        e2 = err;
        if (e2 > -dx) {
            err -= dy;
            x0 += sx;
        }
        if (e2 < dy) {
            err += dx;
            y0 += sy;
        }
    }
}


void Scene::rasterizeSolid(Camera* camera, const Vec4& pointA, const Vec4& pointB, const Vec4& pointC){
    // Calculate bounds for efficient rasterization
    auto xMin = (int)max(0.0, min(pointA.x, min(pointB.x, pointC.x)));
    auto yMin = (int)max(0.0, min(pointA.y, min(pointB.y, pointC.y)));

    auto xMax = (int)max(0.0, min((double)camera->horRes - 1, max(pointA.x, max(pointB.x, pointC.x))));
    auto yMax = (int)max(0.0, min((double)camera->verRes - 1, max(pointA.y, max(pointB.y, pointC.y))));

    // Derive line equations
    auto x0 = pointA.x;
    auto x1 = pointB.x;
    auto x2 = pointC.x;
    auto y0 = pointA.y;
    auto y1 = pointB.y;
    auto y2 = pointC.y;

#pragma region Rasterization Helpers
    auto f01 = [&](double x, double y) {
        return (x * (y0 - y1)) + (y * (x1 - x0)) + (x0 * y1) - (y0 * x1);
    };
    auto f12 = [&](double x, double y) {
        return (x * (y1 - y2)) + (y * (x2 - x1)) + (x1 * y2) - (y1 * x2);
    };
    auto f20 = [&](double x, double y) {
        return (x * (y2 - y0)) + (y * (x0 - x2)) + (x2 * y0) - (y2 * x0);
    };
    auto calculateDepth = [&](double alpha, double beta, double gamma) {
        double depthV0 = pointA.z;
        double depthV1 = pointB.z;
        double depthV2 = pointC.z;

        double interpolatedDepth = depthV0 * alpha +
                                   depthV1 * beta +
                                   depthV2 * gamma;

        return interpolatedDepth;
    };
    auto calculateColor = [&](double alpha, double beta, double gamma) {
        auto c0 = colorsOfVertices[pointA.colorId - 1];
        auto c1 = colorsOfVertices[pointB.colorId - 1];
        auto c2 = colorsOfVertices[pointC.colorId - 1];

        auto color = Color(
                (alpha * (c0->r)) + (beta * (c1->r)) + (gamma * (c2->r)),
                (alpha * (c0->g)) + (beta * (c1->g)) + (gamma * (c2->g)),
                (alpha * (c0->b)) + (beta * (c1->b)) + (gamma * (c2->b))
        );

        return color;
    };
#pragma endregion

    // Iterate over the pixels in the bounds
    for (auto y = yMin; y <= yMax; ++y)
    {
        for (auto x = xMin; x <= xMax; ++x)
        {
            // Calculate barycentric coordinates
            auto alpha = f12(x, y) / f12(x0, y0);
            auto beta = f20(x, y) / f20(x1, y1);
            auto gamma = f01(x, y) / f01(x2, y2);

            // Check whether the pixel is inside the triangle
            if (alpha >= 0 && beta >= 0 && gamma >= 0){

                // Test depth
                auto calculatedDepth = calculateDepth(alpha, beta, gamma);
                if (calculatedDepth < depth[x][y]){
                    depth[x][y] = calculatedDepth;

                    // Calculate color
                    auto color = calculateColor(alpha, beta, gamma);

                    // Round color
                    color = color.round();

                    // Draw pixel
                    image[x][y] = color;
                }
            }
        }
    }
}
#pragma endregion

// Transformations, clipping, culling, rasterization are done here.
void Scene::forwardRenderingPipeline(Camera *camera)
{
    // Calculate camera transformation
    auto camTr = calculateCameraTransformationMatrix(camera);

    // Calculate projection transformation
    auto projTr = (camera->projectionType == 0) ?
                  calculateOrthographicTransformationMatrix(camera) :
                  calculatePerspectiveTransformationMatrix(camera);

    // Calculate viewport transformation
    auto viewTr = calculateViewportTransformationMatrix(camera);

    for (auto & mesh : meshes)
    {
        // Calculate modeling transformation
        auto modelTr = calculateModelingTransformationMatrix(mesh);

        // Precalculate combined transformation
        auto combined = multiplyMatrixWithMatrix(projTr, multiplyMatrixWithMatrix(camTr, modelTr));

        for (auto & triangle : mesh->triangles)
        {
            // Apply transformations
            Vec4 applied[3];
            for (int i = 0; i < 3; ++i)
            {
                auto vertexId = triangle.vertexIds[i];
                auto vertex = vertices[vertexId - 1];
                applied[i] = multiplyMatrixWithVec4(combined, Vec4(vertex->x, vertex->y, vertex->z, 1, vertex->colorId));
            }

            // Apply culling
            if (cullingEnabled)
            {
                auto v0 = Vec3(applied[0].x, applied[0].y, applied[0].z);
                auto v1 = Vec3(applied[1].x, applied[1].y, applied[1].z);
                auto v2 = Vec3(applied[2].x, applied[2].y, applied[2].z);
                auto normal = normalizeVec3(crossProductVec3(subtractVec3(v1, v0), subtractVec3(v2, v0)));
                if (dotProductVec3(normal, v0) < 0){
                    continue;
                }
            }

            // Wireframe
            if (mesh->type == 0){
                // Apply clipping
                Color *c0 = this->colorsOfVertices[applied[0].colorId - 1];
                Color *c1 = this->colorsOfVertices[applied[1].colorId - 1];
                Color *c2 = this->colorsOfVertices[applied[2].colorId - 1];

                clipLine(applied[0], c0, applied[1], c1);
                clipLine(applied[1], c1, applied[2], c2);
                clipLine(applied[2], c2, applied[0], c0);

                // Apply perspective division
                for (auto & point : applied)
                {
                    auto perspective = point.t;
                    point = divideVec4ByScalar(point, perspective);
                }

                // Apply viewport transformation
                for (auto & point : applied)
                {
                    point = multiplyMatrixWithVec4(viewTr, point);
                }

                // Apply rasterization
                rasterizeWireframe(camera, applied[0], applied[1]);
                rasterizeWireframe(camera, applied[1], applied[2]);
                rasterizeWireframe(camera, applied[2], applied[0]);

            }

            // Solid
            if (mesh->type == 1){
                // Apply perspective division
                for (auto & point : applied)
                {
                    auto perspective = point.t;
                    point = divideVec4ByScalar(point, perspective);
                }

                // Apply viewport transformation
                for (auto & point : applied)
                {
                    point = multiplyMatrixWithVec4(viewTr, point);
                }

                // Apply rasterization
                rasterizeSolid(camera, applied[0], applied[1], applied[2]);
            }
        }
    }
}