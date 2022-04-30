import FreeCAD
import CfdCaseWriterFoam
import CfdTools
import CfdConsoleProcess
import Draft
import Part
import PartDesignGui
import Keel

def createBulbSplines(filepath):
    body = FreeCAD.ActiveDocument #.addObject("PartDesign::Body", "Bulb")
    FreeCAD.ActiveDocument.recompute()
    slices = Keel.readSlices(filepath)
    splines = []

    for slice in slices[:-1]:
        points = [tuple(point) for point in slice]
        spline = Draft.makeBSpline(points, closed=True, face=False, support=None)
        splines.append(spline)
        
    #body.addObjects(splines)

    return splines

def createLoft(splines, name, solid = True, ruled = True):
    FreeCAD.ActiveDocument.addObject('Part::Loft',name)
    FreeCAD.ActiveDocument.ActiveObject.Sections= splines
    FreeCAD.ActiveDocument.ActiveObject.Solid=solid
    FreeCAD.ActiveDocument.ActiveObject.Ruled=ruled
    FreeCAD.ActiveDocument.ActiveObject.Closed=False

def createBulb(filepath):
    splines = createBulbSplines(filepath)
    FreeCAD.ActiveDocument.recompute()
    createLoft(splines=splines,name='BulbLoft')
    FreeCAD.ActiveDocument.recompute()

def createKeel():
    #body = FreeCAD.ActiveDocument.addObject("PartDesign::Body", "Keel")
    FreeCAD.ActiveDocument.recompute()

    keelslices = Keel.keelSlices()
    splines = []
    for slice in keelslices:
        points = [(x,y,z) for x,y,z in zip(slice[0],slice[1],slice[2])]
        spline = Draft.makeBSpline(points, closed=True, face=False, support=None)
        splines.append(spline)

    #body.addObjects(splines)
    FreeCAD.ActiveDocument.recompute()
    createLoft(splines=splines,name='KeelLoft')
    FreeCAD.ActiveDocument.recompute()
    
def createCube(Length, Width, Height, x, y, z):
    FreeCAD.ActiveDocument.addObject("Part::Box","Cube")
    FreeCAD.ActiveDocument.ActiveObject.Length = str(Length) + ' mm'
    FreeCAD.ActiveDocument.ActiveObject.Width = str(Width) + ' mm'
    FreeCAD.ActiveDocument.ActiveObject.Height = str(Height) + ' mm'
    FreeCAD.ActiveDocument.ActiveObject.Placement = FreeCAD.Placement(FreeCAD.Vector(x,y,z),FreeCAD.Rotation(FreeCAD.Vector(0,0,1),0))
    FreeCAD.ActiveDocument.recompute()

def createCut(onlyKeel= False):


    if not onlyKeel:

        FreeCAD.ActiveDocument.addObject("Part::Cut",'CutKeel')
        FreeCAD.ActiveDocument.CutKeel.Base = FreeCAD.ActiveDocument.Cube
        FreeCAD.ActiveDocument.CutKeel.Tool = FreeCAD.ActiveDocument.KeelLoft
        FreeCAD.ActiveDocument.Cube.Visibility=False
        FreeCAD.ActiveDocument.KeelLoft.Visibility=False
        FreeCAD.ActiveDocument.recompute()

        FreeCAD.ActiveDocument.addObject("Part::Cut",'CutKeelBulb')
        FreeCAD.ActiveDocument.CutKeelBulb.Base = FreeCAD.ActiveDocument.CutKeel
        FreeCAD.ActiveDocument.CutKeelBulb.Tool = FreeCAD.ActiveDocument.BulbLoft
        FreeCAD.ActiveDocument.CutKeel.Visibility=False
        FreeCAD.ActiveDocument.BulbLoft.Visibility=False
        FreeCAD.ActiveDocument.recompute()
    else:
        FreeCAD.ActiveDocument.addObject("Part::Cut",'CutKeelBulb')
        FreeCAD.ActiveDocument.CutKeelBulb.Base = FreeCAD.ActiveDocument.Cube
        FreeCAD.ActiveDocument.CutKeelBulb.Tool = FreeCAD.ActiveDocument.KeelLoft
        FreeCAD.ActiveDocument.Cube.Visibility=False
        FreeCAD.ActiveDocument.KeelLoft.Visibility=False
        FreeCAD.ActiveDocument.recompute()


def createKeelBulbCut(filepath):
    L_cube, W_cube, H_cube = 4000, 500, 1200
    x_cube,y_cube,z_cube = -1000,0,-400

    createBulb(filepath)
    createKeel()
    createCube(L_cube,W_cube,H_cube,x_cube,y_cube,z_cube)
    createCut()


