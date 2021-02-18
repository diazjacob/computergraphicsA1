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
  divisions *= 2;
  
  float[][][] spherePoints = new float[divisions][divisions][3];
  
  float x = 0;
  float y = 0;
  float z = 0;
    
  //Using the parametric equation for a sphere
  for(int u = 0; u < divisions; u++)//u is the horizontal angle
  {
     for(int v = 0; v < divisions; v++) //v is the vertical angle
     {
       float vFrac = (v+1f) / divisions;
       float uFrac = (u+1f) / divisions;
       
        x = radius*cos(TWO_PI*uFrac)*sin(TWO_PI*vFrac);
        y = radius*sin(TWO_PI*uFrac)*sin(TWO_PI*vFrac);
        z = radius*cos(TWO_PI*vFrac);
        
        //print("(" + x + "," + y + "," + z + ") ");
                        
        //packing into a single array
        spherePoints[u][v][X] = x;
        spherePoints[u][v][Y] = z;
        spherePoints[u][v][Z] = y;
     }
  }
  
  Triangle[] tris = new Triangle[divisions*divisions];
  
  for(int i = 0; i < divisions; i++) //vertical
  {
    for(int j = 0; j < divisions; j++) //horizontal
    {
       float[] p1 = new float[3];
       float[] p2 = new float[3];
       float[] p3 = new float[3];
       
       if(i%2==0)
       {
         p1 = spherePoints[j][i];
         p2 = spherePoints[(j+1)%divisions][i];
         p3 = spherePoints[j][(i+1)%divisions];
       }
       else
       {
         p1 = spherePoints[j][i];
         p2 = spherePoints[(j+1)%divisions][i];
         p3 = spherePoints[(j+1)%divisions][(i+1)%divisions];
       }


       
       
       Triangle t = new Triangle(p1,p2,p3);
       
       tris[(j + i*divisions)] = setupTriangle(t);
       
       //if((j + i*divisions) % 2 == 0)tris[(j + i*divisions)] = setupTriangle(t);
       //else tris[(j + i*divisions)] = new Triangle(new float[]{0,0,0}, new float[]{0,0,0}, new float[]{0,0,0});

    }
    
  }
  
  //return new Triangle[0];
  return tris;
}

// takes a new triangle, and calculates it's normals and edge vectors
Triangle setupTriangle(Triangle t)
{
  t.ev1 = subtract(t.v1, t.v2);
  t.ev2 = subtract(t.v2, t.v3);
  t.ev3 = subtract(t.v3, t.v1);
  
  t.norm = cross3(t.ev1, t.ev2);

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
  if(doOutline && t.norm[Z] < 0)
  {    
     bresLine((int)t.pv1[X], (int)t.pv1[Y], (int)t.pv2[X], (int)t.pv2[Y]);
     bresLine((int)t.pv2[X], (int)t.pv2[Y], (int)t.pv3[X], (int)t.pv3[Y]);
     bresLine((int)t.pv3[X], (int)t.pv3[Y], (int)t.pv1[X], (int)t.pv1[Y]);
  }
  if(normals)
  {
    float[] avgPos = scale(add(t.v3, add(t.v1, t.v2)),1/3f);
   
    float[] normalPos = add(avgPos, t.norm);
    
    float[] proj1 = project(normalPos);
    float[] proj2 = project(avgPos);
    
    beginShape(POINTS);
      vertex(proj2[X], proj2[Y]);
      //vertex(proj1[X], proj1[Y]);
    endShape();
    
    //bresLine((int)proj1[X], (int)proj1[Y], (int)proj2[X], (int)proj2[Y]);
    
  }
}

// uses a scanline algorithm to fill the 2D on-raster triangle
// - implements the specified shading algorithm to set color as specified
// in the global variable shading. Note that for NONE, this function simply
// returns without doing anything
// - uses POINTS to draw on the raster
void fillTriangle(Triangle t, Shading shading)
{
  
}

// given point p, normal n, eye location, light location, calculates phong
// - material represents ambient, diffuse, specular as defined at the top of the file
// - calculates the diffuse, specular, and multiplies it by the material and
// - fillcolor values
float[] phong(float[] p, float[] n, float[] eye, float[] light, 
  float[] material, float[] fillColor, float s)
{
  return new float[]{0, 0, 0};
}



// implements Bresenham's line algorithm
void bresLine(int fromX, int fromY, int toX, int toY)
{
  float error = 0.5;
  int dx = toX - fromX;
  int dy = toY - fromY;
  float m;

  int xPos = fromX;
  int yPos = fromY;
  
  int xStep = 1;
  int yStep = 1;

  //Pre checks here
  if(dx < 0)
  {
    dx = -dx;
    xStep = -1;
  }
  if(dy < 0)
  {
     dy = -dy;
     yStep = -1;
  }
  
  //Initalize the slope properly here
  m = dy/((float)dx); //float cast to avoid integer division trunction
  
  //Running the alg
  if(dx > dy)
  { 
    //The true alg
    while (xPos != toX) //We can approach from either side
    {
       beginShape(POINTS);
         vertex(xPos, yPos);
        endShape();
        
        xPos += xStep;
        error += m;
        if(error >= .5f)
        {
          yPos += yStep;
          error--;
        }
    }
  }
  else
  {
    //invert m cause we're stepping differently
    m = 1/m;
    
    while (yPos != toY) //We can approach from either side
    {
       beginShape(POINTS);
         vertex(xPos, yPos);
        endShape();
        
        yPos += yStep;
        error += m;
        if(error > .5f)
        {
          xPos += xStep;
          error--;
        }
    }
  }
}
