class Triangle {
  Triangle(float[] V1, float[] V2, float[] V3) {  // does DEEP copy!!
    v1 = Arrays.copyOf(V1, V1.length); 
    v2 = Arrays.copyOf(V2, V2.length);
    v3 = Arrays.copyOf(V3, V3.length);
  }

  // position data. in 3D space
  float[] v1; // 3 triangle vertices
  float[] v2;
  float[] v3;

  // projected data. On the screen raster
  float[] pv1; // (p)rojected vertices
  float[] pv2;
  float[] pv3;

  // edge vectors
  float[] ev1; //from 1 to 2
  float[] ev2; //from 2 to 3
  float[] ev3; //from 3 to 1

  // vertex normals
  float[] norm;
  
  //vertex colors (phong result per vertex)
  float[] col1;
  float[] col2;
  float[] col3;
  
  float[] colMid; //phong from the middle of the face.


  // add other things as needed, like normals (face, vectors), edge vectors, colors, etc.
}

Triangle[] sphereList;
Triangle[] rotatedList;

void setup() {
  ortho(-320, 320, 320, -320); // hard coded, 640x640 canvas, RHS
  resetMatrix();
  colorMode(RGB, 1.0f);

  sphereList = makeSphere(SPHERE_SIZE, 10);
  rotatedList = new Triangle[sphereList.length];
  announceSettings();

  //float[] v = subtract(new float[]{2,3,4}, new float[]{1,1,1});
  //print("(" + v[X] + "," + v[Y] + "," + v[Z] + ") ");
}

void settings() {
  size(640, 640, P3D); // hard coded 640x640 canvas
}

float theta = 0.0;
float delta = 0.01;
void draw() {
  clear();

  if (rotate)
  {
    theta += delta;
    while (theta > PI*2) theta -= PI*2;
  } 

  if (lineTest)
    lineTest();
  else
  {
    rotateSphere(sphereList, rotatedList, theta);
    drawSphere(rotatedList, lighting, shading);
  }

  //DEBUG DRAWING SPHERE POINTS
}

//////////////////////  MAIN PROGRAM
// creates a sphere made of triangles, centered on 0,0,0, with given radius
//
// also - 
// calculates the 3 edge vectors for each triangle
// calculates the face normal (unit length)
//
// HINT: first setup a loop to calculate all the points around the sphere,
//       store in an array then loop over those points and setup your triangles.
Triangle[] makeSphere(int radius, int divisions)
{

  float[][][] points = new float[divisions][divisions + 1][3];

  float x = 0;
  float y = 0;
  float z = 0;

  //Using the parametric equation for a sphere
  for (int u = 0; u < divisions; u++)//u is the horizontal angle
  {
    for (int v = 0; v < divisions; v++) //v is the vertical angle
    {
      float theta = PI*((float)v/(divisions));
      float phi = TWO_PI*((float)u/(divisions));

      x = radius*cos(phi)*sin(theta);
      y = radius*sin(phi)*sin(theta);
      z = radius*cos(theta);

      //print("(" + x + "," + y + "," + z + ") ");

      //packing into a single array
      points[u][v][X] = x;
      points[u][v][Y] = z;
      points[u][v][Z] = y;
    }
  }   
  points[divisions-1][divisions][X] = 0;
  points[divisions-1][divisions][Y] = -radius;
  points[divisions-1][divisions][Z] = 0;

  Triangle[] tris = new Triangle[divisions*divisions*2];

  for (int i = 0; i < divisions; i++) //vertical
  {
    for (int j = 0; j < divisions; j++) //horizontal
    {
      float[] p1 = new float[3];
      float[] p2 = new float[3];
      float[] p3 = new float[3];
      float[] pn;

      if (j==0)
      {
        p2 = points[i][j];
        p1 = points[i][(j+1)%divisions];
        p3 = points[(i+1)%(divisions)][(j+1)%(divisions)];
      } else if ((j+1)%divisions == 0)
      {
        p2 = points[i][j];
        p1 = points[divisions-1][divisions];
        p3 = points[(i+1)%(divisions)][j];
      } else
      {
        p2 = points[i][j];
        p1 = points[i][(j+1)%divisions];
        p3 = points[(i+1)%(divisions)][j];
      }


      Triangle t = new Triangle(p1, p2, p3);

      //if(j < 20 && i < 1)tris[(j + i*divisions)] = setupTriangle(t);
      //else tris[(j + i*divisions)] = new Triangle(new float[]{0,0,0}, new float[]{0,0,0}, new float[]{0,0,0});

      tris[(j + i*divisions)*2] = setupTriangle(t);
      ///////

      if (j != 0 && (j+1)%divisions != 0)
      {
        p1 = points[(i+1)%(divisions)][j];
        p3 = points[i][(j+1)%(divisions)];
        p2 = points[(i+1)%(divisions)][(j+1)%(divisions)];
        t = new Triangle(p2, p3, p1);
      } else
      {
        t = new Triangle(new float[]{0, 0, 0}, new float[]{0, 0, 0}, new float[]{0, 0, 0});
      }




      tris[(j + i*divisions)*2 + 1] = setupTriangle(t);
    }
  }



  //return new Triangle[0];
  return tris;
}

// takes a new triangle, and calculates it's normals and edge vectors
Triangle setupTriangle(Triangle t)
{
  t.ev1 = subtract3(t.v2, t.v1);
  t.ev2 = subtract3(t.v3, t.v2);
  t.ev3 = subtract3(t.v1, t.v3);

  t.norm = cross3(t.ev1, t.ev2);
  normalize(t.norm);

  //print("(" +  t.ev1[X] + "," +  t.ev1[Y] + "," +  t.ev1[Z] + ") ");

  return t;
}

// This function draws the 2D, already projected triangle, on the raster
// - it culls degenerate or back-facing triangles
//
// - it calls fillTriangle to do the actual filling, and bresLine to
// make the triangle outline. 
//
// - implements the specified lighting model (using the global enum type)
// to calculate the vertex colors before calling fill triangle. Doesn't do shading
//
// - if needed, it draws the outline and normals (check global variables)
//
// HINT: write it first using the gl LINES/TRIANGLES calls, then replace
// those with your versions once it works.
void draw2DTriangle(Triangle t, Lighting lighting, Shading shading)
{
  //The center of the face of the triangle
  float[] faceMid = scale3(add(t.v3, add(t.v1, t.v2)), 1/3f);

  //This sees if the face is even pointing at the camera, if not then just cull
  if (dot(t.norm, subtract3(EYE, faceMid)) > 0)
  {
    //calculate phong
    t.col1 = phong(t.v1, t.v1, EYE, LIGHT, MATERIAL, FILL_COLOR, PHONG_SPECULAR);
    t.col2 = phong(t.v2, t.v2, EYE, LIGHT, MATERIAL, FILL_COLOR, PHONG_SPECULAR);
    t.col3 = phong(t.v3, t.v3, EYE, LIGHT, MATERIAL, FILL_COLOR, PHONG_SPECULAR);
    t.colMid = phong(faceMid, t.norm, EYE, LIGHT, MATERIAL, FILL_COLOR, PHONG_SPECULAR);
    
    //fill the triangle here
    fillTriangle(t, shading, lighting);

    //Do the outline/overlay stuff overtop if needed here
    if (doOutline)
    {    
      stroke(OUTLINE_COLOR[R], OUTLINE_COLOR[G], OUTLINE_COLOR[B]);
      bresLine((int)t.pv1[X], (int)t.pv1[Y], (int)t.pv2[X], (int)t.pv2[Y]);
      bresLine((int)t.pv2[X], (int)t.pv2[Y], (int)t.pv3[X], (int)t.pv3[Y]);
      bresLine((int)t.pv3[X], (int)t.pv3[Y], (int)t.pv1[X], (int)t.pv1[Y]);
    }
    if (normals) //finally draw the normals on top if we need to
    {
      //calculate the average of the face vertex vectors for the middle
      float[] normalPos = add(faceMid, scale3(t.norm, NORMAL_SCALE));

      //project the normal itself and then then face's midpoint
      float[] proj1 = project(normalPos);
      float[] proj2 = project(faceMid);

      //draw the small line with the correct color :)
      stroke(NORMAL_COLOR[R], NORMAL_COLOR[G], NORMAL_COLOR[B]);
      if (proj1 != null) bresLine((int)proj1[X], (int)proj1[Y], (int)proj2[X], (int)proj2[Y]);
    }
  }
}

// uses a scanline algorithm to fill the 2D on-raster triangle
// - implements the specified shading algorithm to set color as specified
// in the global variable shading. Note that for NONE, this function simply
// returns without doing anything
// - uses POINTS to draw on the raster
void fillTriangle(Triangle t, Shading shading, Lighting lighting)
{
  if (shading != Shading.NONE)
  {
    //we implement scanlines here
    float[] projEV1 = subtract2(t.pv2, t.pv1);
    float[] projEV2 = subtract2(t.pv3, t.pv2);
    float[] projEV3 = subtract2(t.pv1, t.pv3);

    //the bounds of the single projected triangle
    float maxX = ceil(max(new float[]{t.pv1[X], t.pv2[X], t.pv3[X]}));
    float minX = floor(min(new float[]{t.pv1[X], t.pv2[X], t.pv3[X]}));
    float maxY = ceil(max(new float[]{t.pv1[Y], t.pv2[Y], t.pv3[Y]}));
    float minY = floor(min(new float[]{t.pv1[Y], t.pv2[Y], t.pv3[Y]}));


    int xDelta = (int)(maxX - minX);
    int yDelta = (int)(maxY - minY);            

    //Keeping in mind to watch when we use floats to avoid floating point error.
    for (int y = 0; y < yDelta; y++)
    {
      for (int x = 0; x < xDelta; x++)
      {
        
        //Get our point and look at cross products to determine if valid 
        float[] p = new float[]{x + minX, y + minY};
        float c1 = cross2(projEV1, subtract2(p, t.pv1));
        float c2 = cross2(projEV2, subtract2(p, t.pv2));
        float c3 = cross2(projEV3, subtract2(p, t.pv3));

        //if the point is inside, draw the fragment
        if ((c1 < 0 && c2 < 0 && c3 < 0) || (c1 > 0 && c2 > 0 && c3 > 0))
        {
          //Barycentric shadin takes over all as it doesn't need lighting
          if (shading == Shading.BARYCENTRIC)
          {
            float invArea = 2/cross2(projEV1, projEV2);
            float r = c1 * invArea / 2;
            float g = c2 * invArea / 2;
            float b = c3 * invArea / 2;

            stroke(r, g, b);
            beginShape(POINTS);
            vertex((int)p[X], (int)p[Y]);
            endShape();
          }
          else if (shading == Shading.FLAT) //Then flat shading takes over if we want it
          {
            float[] outColor = FILL_COLOR;
            //We look at the lighting data as specified by the user
            if(lighting == Lighting.PHONG_FACE) outColor = t.colMid;
            if(lighting == Lighting.PHONG_VERTEX) 
              outColor = scale3(add(t.col1, add(t.col2, t.col3)), 1/3f);
            
            stroke(outColor[R], outColor[G], outColor[B]);
            beginShape(POINTS);
              vertex((int)p[X], (int)p[Y]);
            endShape();
          } 
      
          else if (shading == Shading.GOURAUD)
          {
            //If we use gouraud shading we implement something very similar to
            // our barycentric visualization
            float invArea = 2/cross2(projEV1, projEV2);
            float u = c1 * invArea / 2;
            float v = c2 * invArea / 2;
            float w = c3 * invArea / 2;

            float[] c = add(scale3(t.col1, v), add(scale3(t.col2, w), scale3(t.col3, u)));

            stroke(c[X], c[Y], c[Z]);
            beginShape(POINTS);
              vertex((int)p[X], (int)p[Y]);
            endShape();
          }
          else if (shading == Shading.PHONG)
          {
            //Finally the bonus phong shading requires us to directly calculate
            //our lighting data here within the fragment itself, instead of on
            //a per-triangle basis.
           
            float invArea = 2/cross2(projEV1, projEV2);
            float u = c1 * invArea / 2;
            float v = c2 * invArea / 2;
            float w = c3 * invArea / 2;
            
            float[] v1Normed = {t.v1[X], t.v1[Y], t.v1[Z]};
            float[] v2Normed = {t.v2[X], t.v2[Y], t.v2[Z]};
            float[] v3Normed = {t.v3[X], t.v3[Y], t.v3[Z]};

            //Interpolating normals
            float[] interNorm = 
                  add(scale3(v1Normed, v), add(scale3(v2Normed, w), scale3(v3Normed, u)));
            float[] interPoint = 
                  add(scale3(t.v1, v), add(scale3(t.v2, w), scale3(t.v3, u)));
                  
            //Calling phong!
            float[] phong = phong(interPoint, interNorm, EYE, LIGHT, MATERIAL, FILL_COLOR, PHONG_SPECULAR);

            stroke(phong[X], phong[Y], phong[Z]);
            beginShape(POINTS);
              vertex((int)p[X], (int)p[Y]);
            endShape();
          }
        }
      }
    }
  }
}

// given point p, normal n, eye location, light location, calculates phong
// - material represents ambient, diffuse, specular as defined at the top of the file
// - calculates the diffuse, specular, and multiplies it by the material and
// - fillcolor values
float[] phong(float[] p, float[] n, float[] eye, float[] light, 
  float[] material, float[] fillColor, float s)
{
  
  //Watch out! We were told to create and in-place normalization, so I copy
  //Our data inside in order to normalize and not mess up other data (and possibly cause bugs)
  float[] point = {p[X], p[Y], p[Z]};
  normalize(point);
  float[] norm = {n[X], n[Y], n[Z]};
  normalize(norm);
  float[] eyeCopy = {eye[X], eye[Y], eye[Z]};
  normalize(eyeCopy);
  float[] lightCopy = {light[X], light[Y], light[Z]};
  normalize(lightCopy);
  
  //ambient
  float ambient = material[M_AMBIENT];
  
  //diffuse amount
  float diffuse = dot(lightCopy, norm) * material[M_DIFFUSE];
  
  //specularity
  float[] viewVec = subtract3(eye,p); //order matters !
  normalize(viewVec);
  float[] reflectedLight = subtract3(scale3(norm, 2*dot(norm,lightCopy)), lightCopy);
  normalize(reflectedLight);
  
  float specular = 0;
  //This if removes the "double shine" spot due to how the dot product
  //is symmetric and even powers will make it pos, or odd powers will do weird negitive things.
  if(dot(reflectedLight, viewVec) > 0)
   specular = pow(dot(reflectedLight, viewVec),s) * material[M_SPECULAR];
  
  //mix it all together
  float[] finalColor = {0,0,0};
  finalColor = add(finalColor, scale3(FILL_COLOR, ambient));
  finalColor = add(finalColor, scale3(FILL_COLOR, diffuse));
  finalColor = add(finalColor, scale3(FILL_COLOR, specular));
    
  return finalColor;
}



// implements Bresenham's line algorithm
void bresLine(int fromX, int fromY, int toX, int toY)
{
  //Initalizations and our error value here
  float error = 0.5f;
  int dx = toX - fromX;
  int dy = toY - fromY;
  float m;

  int xPos = fromX;
  int yPos = fromY;

  int xStep = 1;
  int yStep = 1;

  //Pre checks here (analyizing line slopes
  if (dx < 0)
  {
    dx = -dx;
    xStep = -1;
  }
  if (dy < 0)
  {
    dy = -dy;
    yStep = -1;
  }

  //Initalize the slope properly here
  m = dy/((float)dx); //float cast to avoid integer division trunction

  //Running the alg
  if (dx > dy)
  { 
    //The true alg
    while (xPos != toX) //We can approach from either side
    {
      beginShape(POINTS);
      vertex(xPos, yPos);
      endShape();

      xPos += xStep;
      error += m;
      if (error >= .5f)
      {
        yPos += yStep;
        error--;
      }
    }
  } else
  {
    //invert m cause we're stepping differently in this case
    m = 1/m;

    while (yPos != toY) //So we can approach from either side!
    {
      beginShape(POINTS);
      vertex(xPos, yPos);
      endShape();

      yPos += yStep;
      error += m;
      if (error > .5f)
      {
        xPos += xStep;
        error--;
      }
    }
  }
}
