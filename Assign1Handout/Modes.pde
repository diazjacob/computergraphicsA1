final char KEY_LIGHTING = 'l';
final char KEY_SHADING = 's';
final char KEY_OUTLINE = 'o';
final char KEY_ROTATE = 'r';
final char KEY_NORMALS = 'n';
final char KEY_ACCELERATED = '!';
final char KEY_LINE_TEST = 't';

enum Lighting {
  FLAT, 
    PHONG_FACE, 
    PHONG_VERTEX
}
Lighting lighting = Lighting.FLAT;

enum Shading {
  NONE, 
    BARYCENTRIC, 
    FLAT, 
    GOURAUD, 
    PHONG
}
Shading shading = Shading.NONE;

boolean doOutline = true; // to draw, or not to draw (outline).. is the question
boolean rotate = false;
boolean normals = false;
boolean accelerated = false;
boolean lineTest = false;

void keyPressed()
{
  if (key == KEY_SHADING)
  {
    int next = (shading.ordinal()+1)%Shading.values().length;
    shading = Shading.values()[next];
  }
  if (key == KEY_LIGHTING)
  {
    int next = (lighting.ordinal()+1)%Lighting.values().length;
    lighting = Lighting.values()[next];
  }
  if (key == KEY_OUTLINE)
    doOutline = !doOutline;

  if (key == KEY_ROTATE)
    rotate = !rotate;

  if (key == KEY_NORMALS)
    normals = !normals;

  if (key == KEY_ACCELERATED)
    accelerated = !accelerated;

  if (key == KEY_LINE_TEST)
    lineTest = !lineTest;

  announceSettings();
}

void announceSettings()
{
  String msg = "";
  if (lineTest)
    msg += "Line Test";
  else
    msg += "Shading: "+shading +" Lighting: "+lighting;

  if (accelerated)
    msg += "(accelerated) ";
  println(msg);
}
