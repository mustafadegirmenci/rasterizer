#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "../include/tinyxml2.h"
#include "../include/Triangle.h"
#include "../include/Helpers.h"
#include "../include/Scene.h"

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

	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}

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
                auto vertex = vertices[vertexId];
                applied[i] = multiplyMatrixWithVec4(combined, Vec4(vertex->x, vertex->y, vertex->z, vertex->colorId));
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
                // Generate lines
                // TODO

                // Apply clipping
                // TODO

                // Apply viewport transformation
                // TODO

                // Apply rasterization
                // TODO

            }

            // Solid
            if (mesh->type == 1){
                // Apply perspective division
                // TODO

                // Apply viewport transformation
                // TODO

                // Apply rasterization
                // TODO

            }
        }
    }
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
            auto translationMatrix = Matrix4(
                    (double[4][4]){{1,0,0,x},
                                   {0,1,0,y},
                                   {0,0,1,z},
                                   {0,0,0,1}});
            modelingTransformation = multiplyMatrixWithMatrix(translationMatrix, modelingTransformation);
        }

        // Scaling
        if (mesh->transformationTypes[i] == 's')
        {
            auto scaling = scalings[mesh->transformationIds[i] - 1];
            auto x = scaling->sx;
            auto y = scaling->sy;
            auto z = scaling->sz;
            auto scalingMatrix = Matrix4(
                    (double[4][4]){{x,0,0,0},
                                   {0,y,0,0},
                                   {0,0,z,0},
                                   {0,0,0,1}});
            modelingTransformation = multiplyMatrixWithMatrix(scalingMatrix, modelingTransformation);
        }

        // Rotation
        if (mesh->transformationTypes[i] == 's')
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
            auto rotationMatrix = Matrix4(
                    (double[4][4]){{1,      0,          0,          0},
                                   {0,      cos(ra),    -sin(ra),   0},
                                   {0,      sin(ra),    cos(ra),    0},
                                   {0,      0,          0,          1}});

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

    auto translate = Matrix4((double[4][4]){
            {1, 0, 0, -x},
            {0, 1, 0, -y},
            {0, 0, 1, -z},
            {0, 0, 0, 1}
    });

    auto rotate = Matrix4((double[4][4]){
            {u.x,   u.y,    u.z,    0},
            {v.x,   v.y,    v.z,    0},
            {w.x,   w.y,    w.z,    0},
            {0,     0,      0,      1}
    });

    return multiplyMatrixWithMatrix(rotate, translate);
}

Matrix4 Scene::calculateOrthographicTransformationMatrix(Camera* camera){
    auto l = camera->left;
    auto r = camera->right;
    auto t = camera->top;
    auto b = camera->bottom;
    auto n = camera->near;
    auto f = camera->far;

    auto perspective = Matrix4((double[4][4]){
        {2/(r - l),             0,                  0,                  -((r + l) / (r - l))},
        {0,                     2/(t - b),          0,                  -((t + b) / (t - b))},
        {0,                     0,                  -(2/(f - n)),       -((f + n) / (f - n))},
        {0,                     0,                  0,                  1                   }
    });

    return perspective;
}

Matrix4 Scene::calculatePerspectiveTransformationMatrix(Camera* camera){
    auto l = camera->left;
    auto r = camera->right;
    auto t = camera->top;
    auto b = camera->bottom;
    auto n = camera->near;
    auto f = camera->far;

    auto perspective = Matrix4((double[4][4]){
            {(2*n) / (r - l),   0,                  (r + l) / (r - l),      0                       },
            {0,                 (2*n) / (t - b),    (t + b) / (t - b),      0                       },
            {0,                 0,                  -((f + n) / (f - n)),   -((2 * f * n) / (f - n))},
            {0,                 0,                  -1,                     0                       }
    });

    return perspective;
}

Matrix4 Scene::calculateViewportTransformationMatrix(Camera* camera){
    auto h = camera->horRes;
    auto v = camera->verRes;

    auto viewport = Matrix4((double[4][4]){
        {h/2.0,     0,          0,          (h-1)/2.0},
        {0,         v/2.0,      0,          (v-1)/2.0},
        {0,         0,          0.5,        0.5      },
        {0,         0,          0,          1        }
    });

    return viewport;
}
#pragma endregion

#pragma region Rasterization
void Scene::rasterizeWireframe(Vec4 points[3]){

}

void Scene::rasterizeSolid(Vec4 points[3]){

}
#pragma endregion