import FreeCAD
import SimulationCFD
import SolidGenerator
import CfdAnalysis
import CfdTools
import CfdPhysicsSelection
import CfdFluidMaterial
import CfdInitialiseFlowField
import CfdSolverFoam
import CfdMesh
import CfdMeshTools
import CfdConsoleProcess
import CfdMeshRefinement
import CfdFluidBoundary
import CfdCaseWriterFoam
import Draft
import Part
import PartDesignGui
import shutil
from msilib.schema import Class
import numpy as np
import os

class Airfoil:

    def y_naca(x,t):
      return 5*t*(0.2969*(x**0.5) - 0.126*x - 0.3516*(x**2) + 0.2843*(x**3) - 0.1015*(x**4) )

    def arange(start,end,N):
      step = (end-start)/(N-1)
      return [start + i*step for i in range(N) ]

    def genNACAxyz(xx,c,z=0,N=200):
      t = xx/100

      xtab = []
      y_pos = []
      y_neg = []

      for x in Airfoil.arange(0,1.0,N):
        xtab.append(x*c)
        y_ = Airfoil.y_naca(x,t)*c
        y_pos.append(y_)
        y_neg.append(-y_)

      xtab = xtab + xtab[::-1]
      ytab = y_pos + y_neg[::-1]
      ztab = [z for i in range(len(xtab))]
      return (xtab,ytab,ztab)

    def find_closest(y_find,y_array):
      y0 = 0
      i = 0
      for y in y_array:
        if y_find == y:
          return (y,i)
        elif (y_find >= y0)and(y_find < y):
          return (y,i)
        i += 1

class Automate:
        
    def genGeometry(bulbpath):
        SolidGenerator.createBulb(bulbpath)
        SolidGenerator.createKeel()
        SolidGenerator.createCube(4000, 700, 1700, -1000, 0, -480)
        SolidGenerator.createCut()


    class Design:

        default_mesh_size = 23 #mm
        default_mesh_refin = 0.5
        default_inlet_v = 5000 #mm/s

        def __init__(self, bulbpath, outputfolder, inlet_v = default_inlet_v):

            self.bulbpath = bulbpath
            self.outputfolder = outputfolder
            self.inlet_v = inlet_v

        def start(self):
            self.genGeometry()
            self.startSimulation()
            self.saveDocument()
            self.clearAll()


        def genGeometry(self):
            SolidGenerator.createBulb(self.bulbpath)
            SolidGenerator.createKeel()
            SolidGenerator.createCube(4000, 700, 1700, -1000, 0, -480)
            SolidGenerator.createCut()
        
        def startSimulation(self):
            mesh_size = Automate.Design.default_mesh_size

            SimulationCFD.createAnalysis()
            SimulationCFD.createMeshCase(size=mesh_size,refinement=Automate.Design.default_mesh_refin)
            SimulationCFD.runMesh()
            SimulationCFD.createBoundaries(inlet_v=self.inlet_v)

            cfd_valid = SimulationCFD.runSolver()

            if not cfd_valid:
                for n in range(5):
                    mesh_size -= 1
                    SimulationCFD.createMeshCase(size=mesh_size,refinement=Automate.Design.default_mesh_refin)
                    SimulationCFD.runMesh()
                    cfd_valid = SimulationCFD.runSolver()
                    if cfd_valid:
                        break
            
            if not cfd_valid:
                raise Exception('Error running simualtion')
            else:
                SimulationCFD.createResults(self.outputfolder)

        def saveDocument(self):
            doc = FreeCAD.ActiveDocument
            path = self.outputfolder + '\\copyDoc.FCStd'
            doc.saveCopy(path)

        def clearAll(self):
            doc = FreeCAD.ActiveDocument
            for obj in doc.Objects:
                doc.removeObject(obj.Name)

class Bezier:
        
    class BezierCurve:
        def __init__(self,c0,z0,c1,z1,z2,c3,z3,I):
            self.c0 = c0
            self.z0 = z0
            self.c1 = c1
            self.z1 = z1
            self.z2 = z2
            self.c3 = c3
            self.z3 = z3
            self.I = I

            self.c2 = self.get_c2()
        
        def Bezier(self,t):

            c = ((1-t)**3)*self.c0 + 3*((1-t)**2)*t*self.c1 + 3*((1-t)*t**2)*self.c2 + self.c3*t**3
            z = ((1-t)**3)*self.z0 + 3*((1-t)**2)*t*self.z1 + 3*((1-t)*t**2)*self.z2 + self.z3*t**3
            return (c, z)

        def getBezierCurve(self, resolution = 0.0001):

            c_tab , z_tab = [self.c0],[self.z0]

            for t in np.arange(0,1,resolution):
                c,z = self.Bezier(t)
                c_tab.append(c)
                z_tab.append(z)

            c_tab.append(self.c3)
            z_tab.append(self.z3)

            self.bezierCurve = (c_tab,z_tab)


            return (c_tab,z_tab)

        def get_c2(self):
            c0 = self.c0
            z0 = self.z0
            c1 = self.c1
            z1 = self.z1
            z2 = self.z2
            c3 = self.c3
            z3 = self.z3
            I = self.I

            sqrt = lambda x: x**0.5
            solution = (-10*c0*z0 + 6*c0*z2 + 4*c0*z3 - 15*c1*z0 - 9*c1*z1 + 9*c1*z2 + 15*c1*z3 - 5*c3*z0 - 15*c3*z1 - 15*c3*z2 + 35*c3*z3 - sqrt(3)*sqrt(-1120*I*z0 - 1680*I*z1 + 2800*I*z3 - 340*c0**2*z0**2 - 280*c0**2*z0*z1 + 40*c0**2*z0*z2 + 920*c0**2*z0*z3 + 420*c0**2*z1**2 + 120*c0**2*z1*z2 - 680*c0**2*z1*z3 + 12*c0**2*z2**2 - 184*c0**2*z2*z3 - 28*c0**2*z3**2 - 180*c0*c1*z0**2 - 240*c0*c1*z0*z1 + 600*c0*c1*z0*z3 + 180*c0*c1*z1**2 + 144*c0*c1*z1*z2 - 264*c0*c1*z1*z3 + 36*c0*c1*z2**2 - 216*c0*c1*z2*z3 - 60*c0*c1*z3**2 + 20*c0*c3*z0**2 + 72*c0*c3*z0*z1 + 88*c0*c3*z0*z2 - 200*c0*c3*z0*z3 - 12*c0*c3*z1**2 - 48*c0*c3*z1*z2 - 60*c0*c3*z2**2 + 80*c0*c3*z2*z3 + 60*c0*c3*z3**2 - 45*c1**2*z0**2 - 90*c1**2*z0*z1 - 18*c1**2*z0*z2 + 198*c1**2*z0*z3 + 27*c1**2*z1**2 + 54*c1**2*z1*z2 - 18*c1**2*z1*z3 + 27*c1**2*z2**2 - 90*c1**2*z2*z3 - 45*c1**2*z3**2 + 18*c1*c3*z0**2 + 84*c1*c3*z0*z1 + 120*c1*c3*z0*z2 - 240*c1*c3*z0*z3 + 18*c1*c3*z1**2 - 120*c1*c3*z1*z3 - 90*c1*c3*z2**2 + 60*c1*c3*z2*z3 + 150*c1*c3*z3**2 - 5*c3**2*z0**2 - 50*c3**2*z0*z1 - 230*c3**2*z0*z2 + 290*c3**2*z0*z3 - 45*c3**2*z1**2 - 270*c3**2*z1*z2 + 410*c3**2*z1*z3 + 75*c3**2*z2**2 + 350*c3**2*z2*z3 - 525*c3**2*z3**2))/(6*(2*z0 + 3*z1 - 5*z3))
            return solution.real

class Keel:

    def genBulbSlices(c1,z1,z2,c3,z3, xx=12, nSlices=10, slope = False, error = 0.07):
        
        '''
        c0,z0 = 810,0
        c1,z1 = 990,20
        z2 = 180
        c3,z3 = 0,226
        I = 1.956e08
        '''
        c0, z0 = 836, 0

        M  = 181.44     #(kg) -> 400 lbs
        density = 11343 #(kg/m3) Lead
        naca_factor = 8221.0/12000.0 # = A/(t*c^2)

        h = z3

        v_intersection = naca_factor * (12/100) * (0.0417951312013333*h**3 + 296.025928*h**2 + 698896.0*h)
        v_total = (M/density) * 1e9 #Volume in mm3
        v_bulb = (v_total + v_intersection)*(1-error)

        I = v_bulb/(naca_factor*(xx/100.0))

        bCurve = Bezier.BezierCurve(
            c0,z0,
            c1,z1,
            z2,
            c3,z3,
            I
        )

        ccurve, zcurve = bCurve.getBezierCurve()

        x_z_func = lambda z : -0.354098*z  

        ztab, ctab = [],[]
        slices = []

        for z in Airfoil.arange(z0,z3,nSlices):
            z_find,i = Airfoil.find_closest(z,zcurve)
            c = ccurve[i]
            xyz = Airfoil.genNACAxyz(xx,c,z)
            if slope:
                dx = x_z_func(z)
                xyz = Keel.movexyz(xyz,dx=dx)
            slices.append(xyz)

        return slices

    def genBulbFile(folderpath,c1,z1,z2,c3,z3,xx=12,nSlices=10,slope=False, error = 0.07):

        slices = Keel.genBulbSlices(c1,z1,z2,c3,z3,xx,nSlices,slope,error)

        csv = ''
        for n in range(len(slices[0][0])):
            row = ''
            for slice in slices:
                slice_csv = ';'.join([ str(slice[0][n]) , str(slice[1][n]) , str(slice[2][n]) ])
                row += slice_csv + '|'
            row = row[:-1]
            csv += row + '\n'

        filename = '\\' +  '_'.join([
            str(xx),
            str(c1),
            str(z1),
            str(z2),
            str(c3),
            str(z3)
        ]) 

        f = open(folderpath + filename + '.csv','w')
        f.write(csv)
        f.close()

    def readSlices(filepath):

        f = open(filepath,'r')
        lines = f.readlines()
        nSlices = len(lines[0].split('|'))
        slices = [[] for n in range(nSlices)]
        for line in lines:
            slices_row = line.replace('\n','').split('|')
            for n in range(nSlices):
                xyz = [float(i) for i in slices_row[n].split(';')]
                slices[n].append(xyz)

        return slices

    def movexyz(xyz_array,dx=0,dy=0,dz=0):
            x_array = [j + dx for j in xyz_array[0]]
            y_array = [j + dy for j in xyz_array[1]]
            z_array = [j + dz for j in xyz_array[2]]

            return (x_array,y_array,z_array)

    def keelSlices():

        keel_sections = [
                    [838 , 0     ],
                    [940 , 305   ],
                    [1054, 305*2 ],
                    [1165, 305*3 ],
                    [1289, 305*4 ]
        ]

        keel_slices = []
        top_chord = keel_sections[-1][0]
        bottom_chord = keel_sections[0][0]
        for c,z in keel_sections:
            keel_naca = Airfoil.genNACAxyz(xx=12,c=c,z=z)
            keel_slices.append(Keel.movexyz(keel_naca,dx=(bottom_chord-c)))
        
        return keel_slices

class SimulationCFD:
        
    def createAnalysis():

        analysis = CfdAnalysis.makeCfdAnalysis('CfdAnalysis')
        CfdTools.setActiveAnalysis(analysis)
        analysis.addObject(CfdPhysicsSelection.makeCfdPhysicsSelection())
        analysis.addObject(CfdFluidMaterial.makeCfdFluidMaterial('FluidProperties'))
        analysis.addObject(CfdInitialiseFlowField.makeCfdInitialFlowField())
        analysis.addObject(CfdSolverFoam.makeCfdSolverFoam())

        obj = FreeCAD.ActiveDocument.PhysicsModel
        obj.Time = 'Steady'
        obj.Phase = 'Single'
        obj.Flow = 'Incompressible'
        obj.Thermal = 'None'
        obj.Turbulence = 'RANS'
        obj.TurbulenceModel = 'kOmegaSST'
        obj.gx = '0 mm/s^2'
        obj.gy = '-9,8e+03 mm/s^2'
        obj.gz = '0 mm/s^2'
        
        FreeCAD.ActiveDocument.FluidProperties.Material = {'CardName': 'WaterIsothermal', 'AuthorAndLicense': 'Water', 'Name': 'Water', 'Type': 'Isothermal', 'Description': 'Standard distilled water properties at 20 Degrees Celsius and 1 atm', 'Density': '998 kg/m^3', 'DynamicViscosity': '1.003e-3 kg/m/s'}

        init = FreeCAD.ActiveDocument.InitialiseFields
        init.PotentialFlow = True
        init.UseInletUValues = False
        init.Ux = '0.0 mm/s'
        init.Uy = '0.0 mm/s'
        init.Uz = '0.0 mm/s'
        init.UseOutletPValue = False
        init.PotentialFlowP = True
        init.Pressure = '0.0 kg/(mm*s^2)'
        init.VolumeFractions = {}
        init.UseInletTemperatureValue = False
        init.Temperature = '290.0 K'
        init.UseInletTurbulenceValues = False
        init.omega = '57.0 deg/s'
        init.k = '10000.0 mm^2/s^2'
        init.BoundaryU = None
        init.BoundaryP = None
        init.BoundaryT = None
        init.BoundaryTurb = None

        FreeCAD.ActiveDocument.recompute()

    def createMeshCase(size=25.0,refinement=0.5):

        char_length = str(size) + ' mm'

        CfdMesh.makeCfdMesh('CutKeelBulb_Mesh')
        FreeCAD.ActiveDocument.ActiveObject.Part = FreeCAD.ActiveDocument.CutKeelBulb
        CfdTools.getActiveAnalysis().addObject(FreeCAD.ActiveDocument.ActiveObject)

        FreeCAD.ActiveDocument.CutKeelBulb_Mesh.CharacteristicLengthMax = char_length
        FreeCAD.ActiveDocument.CutKeelBulb_Mesh.MeshUtility = 'cfMesh'
        FreeCAD.ActiveDocument.CutKeelBulb_Mesh.ElementDimension = '3D'
        FreeCAD.ActiveDocument.CutKeelBulb_Mesh.CellsBetweenLevels = 3
        FreeCAD.ActiveDocument.CutKeelBulb_Mesh.EdgeRefinement = 1.0
        FreeCAD.ActiveDocument.CutKeelBulb_Mesh.PointInMesh = {'x': '0.0 mm', 'y': '0.0 mm', 'z': '0.0 mm'}

        FreeCAD.ActiveDocument.recompute()

        CfdMeshRefinement.makeCfdMeshRefinement(FreeCAD.ActiveDocument.CutKeelBulb_Mesh)

        FreeCAD.ActiveDocument.MeshRefinement.RelativeLength = refinement
        FreeCAD.ActiveDocument.MeshRefinement.RefinementThickness = '0.0 mm'
        FreeCAD.ActiveDocument.MeshRefinement.NumberLayers = 1
        FreeCAD.ActiveDocument.MeshRefinement.ExpansionRatio = 1.2
        FreeCAD.ActiveDocument.MeshRefinement.FirstLayerHeight = '0.0 mm'
        FreeCAD.ActiveDocument.MeshRefinement.RegionEdgeRefinement = 1.0
        FreeCAD.ActiveDocument.MeshRefinement.Internal = True
        FreeCAD.ActiveDocument.MeshRefinement.ShapeRefs = [
        (FreeCAD.ActiveDocument.getObject('KeelLoft'), ('Solid1',)),
        (FreeCAD.ActiveDocument.getObject('BulbLoft'), ('Solid1',))]

        FreeCAD.ActiveDocument.recompute()

    def runMesh():
        cart_mesh = CfdMeshTools.CfdMeshTools(FreeCAD.ActiveDocument.CutKeelBulb_Mesh)
        FreeCAD.ActiveDocument.CutKeelBulb_Mesh.Proxy.cart_mesh = cart_mesh
        cart_mesh.writeMesh()
        cart_mesh = CfdMeshTools.CfdMeshTools(FreeCAD.ActiveDocument.CutKeelBulb_Mesh)
        proxy = FreeCAD.ActiveDocument.CutKeelBulb_Mesh.Proxy
        proxy.cart_mesh = cart_mesh
        cart_mesh.error = False
        cmd = CfdTools.makeRunCommand('./Allmesh', cart_mesh.meshCaseDir, source_env=False)
        FreeCAD.Console.PrintMessage('Executing: ' + ' '.join(cmd) + '\n')
        env_vars = CfdTools.getRunEnvironment()
        '''
        proxy.running_from_macro = True
        if proxy.running_from_macro:
            mesh_process = CfdConsoleProcess.CfdConsoleProcess()
            mesh_process.start(cmd, env_vars=env_vars)
            mesh_process.waitForFinished()
        else:
            proxy.mesh_process.start(cmd, env_vars=env_vars)
        '''
        mesh_process = CfdConsoleProcess.CfdConsoleProcess()
        mesh_process.start(cmd, env_vars=env_vars)
        mesh_process.waitForFinished()

        return mesh_process.exitCode() == 0

    def createBoundaries(inlet_v = 5000.0, surface_slip = True):

        # INLET

        CfdTools.getActiveAnalysis().addObject(CfdFluidBoundary.makeCfdFluidBoundary())

        bc = FreeCAD.ActiveDocument.CfdFluidBoundary
        bc.BoundaryType = 'inlet'
        bc.BoundarySubType = 'uniformVelocityInlet'
        bc.ThermalBoundaryType = 'fixedValue'
        bc.VelocityIsCartesian = True
        bc.Ux = str(inlet_v) + ' mm/s'
        bc.Uy = '0.0 mm/s'
        bc.Uz = '0.0 mm/s'
        bc.VelocityMag = '0.0 mm/s'
        bc.DirectionFace = ''
        bc.ReverseNormal = True
        bc.MassFlowRate = '0.0 kg/s'
        bc.VolFlowRate = '0.0 mm^3/s'
        bc.Pressure = '0.0 kg/(mm*s^2)'
        bc.SlipRatio = '0.0'
        bc.Temperature = '290.0 K'
        bc.HeatFlux = '0.0 kg/s^3'
        bc.HeatTransferCoeff = '0.0 kg/(s^3*K)'
        bc.TurbulenceInletSpecification = 'intensityAndLengthScale'
        bc.TurbulentKineticEnergy = '10000.0 mm^2/s^2'
        bc.SpecificDissipationRate = '57.0 deg/s'
        bc.TurbulenceIntensity = '0.1'
        bc.TurbulenceLengthScale = '100.0 mm'
        bc.VolumeFractions = {}
        bc.PorousBaffleMethod = 'porousCoeff'
        bc.PressureDropCoeff = '0.0'
        bc.ScreenWireDiameter = '0.2 mm'
        bc.ScreenSpacing = '2.0 mm'
        FreeCAD.ActiveDocument.CfdFluidBoundary.Label = 'inlet'
        FreeCAD.ActiveDocument.CfdFluidBoundary.ShapeRefs = [
        (FreeCAD.ActiveDocument.getObject('CutKeelBulb'), ('Face1',))]
        bc.DefaultBoundary = False

        FreeCAD.ActiveDocument.recompute()

        # OUTLET

        outlet_face = 'Face' + str(len(FreeCAD.ActiveDocument.CutKeelBulb.Shape.Faces)-1)

        CfdTools.getActiveAnalysis().addObject(CfdFluidBoundary.makeCfdFluidBoundary())
        bc = FreeCAD.ActiveDocument.CfdFluidBoundary001
        bc.BoundaryType = 'outlet'
        bc.BoundarySubType = 'staticPressureOutlet'
        bc.ThermalBoundaryType = 'fixedValue'
        bc.VelocityIsCartesian = True
        bc.Ux = '0.0 mm/s'
        bc.Uy = '0.0 mm/s'
        bc.Uz = '0.0 mm/s'
        bc.VelocityMag = '0.0 mm/s'
        bc.DirectionFace = ''
        bc.ReverseNormal = False
        bc.MassFlowRate = '0.0 kg/s'
        bc.VolFlowRate = '0.0 mm^3/s'
        bc.Pressure = '0.0 kg/(mm*s^2)'
        bc.SlipRatio = '0.0'
        bc.Temperature = '290.0 K'
        bc.HeatFlux = '0.0 kg/s^3'
        bc.HeatTransferCoeff = '0.0 kg/(s^3*K)'
        bc.TurbulenceInletSpecification = 'intensityAndLengthScale'
        bc.TurbulentKineticEnergy = '10000.0 mm^2/s^2'
        bc.SpecificDissipationRate = '57.0 deg/s'
        bc.TurbulenceIntensity = '0.1'
        bc.TurbulenceLengthScale = '100.0 mm'
        bc.VolumeFractions = {}
        bc.PorousBaffleMethod = 'porousCoeff'
        bc.PressureDropCoeff = '0.0'
        bc.ScreenWireDiameter = '0.2 mm'
        bc.ScreenSpacing = '2.0 mm'
        FreeCAD.ActiveDocument.CfdFluidBoundary001.Label = 'outlet'
        FreeCAD.ActiveDocument.CfdFluidBoundary001.ShapeRefs = [
        (FreeCAD.ActiveDocument.getObject('CutKeelBulb'), (outlet_face,))]
        bc.DefaultBoundary = False
        FreeCAD.ActiveDocument.recompute()

        # WALLS

        CfdTools.getActiveAnalysis().addObject(CfdFluidBoundary.makeCfdFluidBoundary())

        bc = FreeCAD.ActiveDocument.CfdFluidBoundary002
        bc.BoundaryType = 'wall'
        bc.BoundarySubType = 'slipWall'
        bc.ThermalBoundaryType = 'fixedValue'
        bc.VelocityIsCartesian = True
        bc.Ux = '0.0 mm/s'
        bc.Uy = '0.0 mm/s'
        bc.Uz = '0.0 mm/s'
        bc.VelocityMag = '0.0 mm/s'
        bc.DirectionFace = ''
        bc.ReverseNormal = False
        bc.MassFlowRate = '0.0 kg/s'
        bc.VolFlowRate = '0.0 mm^3/s'
        bc.Pressure = '0.0 kg/(mm*s^2)'
        bc.SlipRatio = '0.0'
        bc.Temperature = '290.0 K'
        bc.HeatFlux = '0.0 kg/s^3'
        bc.HeatTransferCoeff = '0.0 kg/(s^3*K)'
        bc.TurbulenceInletSpecification = 'intensityAndLengthScale'
        bc.TurbulentKineticEnergy = '10000.0 mm^2/s^2'
        bc.SpecificDissipationRate = '57.0 deg/s'
        bc.TurbulenceIntensity = '0.1'
        bc.TurbulenceLengthScale = '100.0 mm'
        bc.VolumeFractions = {}
        bc.PorousBaffleMethod = 'porousCoeff'
        bc.PressureDropCoeff = '0.0'
        bc.ScreenWireDiameter = '0.2 mm'
        bc.ScreenSpacing = '2.0 mm'
        FreeCAD.ActiveDocument.CfdFluidBoundary002.Label = 'wall'
        FreeCAD.ActiveDocument.CfdFluidBoundary002.ShapeRefs = [
        (FreeCAD.ActiveDocument.getObject('CutKeelBulb'), ('Face3', 'Face4', 'Face5'))]
        bc.DefaultBoundary = False

        FreeCAD.ActiveDocument.recompute()

        # SYMMETRY

        CfdTools.getActiveAnalysis().addObject(CfdFluidBoundary.makeCfdFluidBoundary())

        bc = FreeCAD.ActiveDocument.CfdFluidBoundary003
        bc.BoundaryType = 'constraint'
        bc.BoundarySubType = 'symmetry'
        bc.ThermalBoundaryType = 'fixedValue'
        bc.VelocityIsCartesian = True
        bc.Ux = '0.0 mm/s'
        bc.Uy = '0.0 mm/s'
        bc.Uz = '0.0 mm/s'
        bc.VelocityMag = '0.0 mm/s'
        bc.DirectionFace = ''
        bc.ReverseNormal = False
        bc.MassFlowRate = '0.0 kg/s'
        bc.VolFlowRate = '0.0 mm^3/s'
        bc.Pressure = '0.0 kg/(mm*s^2)'
        bc.SlipRatio = '0.0'
        bc.Temperature = '290.0 K'
        bc.HeatFlux = '0.0 kg/s^3'
        bc.HeatTransferCoeff = '0.0 kg/(s^3*K)'
        bc.TurbulenceInletSpecification = 'intensityAndLengthScale'
        bc.TurbulentKineticEnergy = '10000.0 mm^2/s^2'
        bc.SpecificDissipationRate = '57.0 deg/s'
        bc.TurbulenceIntensity = '0.1'
        bc.TurbulenceLengthScale = '100.0 mm'
        bc.VolumeFractions = {}
        bc.PorousBaffleMethod = 'porousCoeff'
        bc.PressureDropCoeff = '0.0'
        bc.ScreenWireDiameter = '0.2 mm'
        bc.ScreenSpacing = '2.0 mm'
        FreeCAD.ActiveDocument.CfdFluidBoundary003.Label = 'constraint'
        FreeCAD.ActiveDocument.CfdFluidBoundary003.ShapeRefs = [
        (FreeCAD.ActiveDocument.getObject('CutKeelBulb'), ('Face2',))]
        bc.DefaultBoundary = False
        FreeCAD.ActiveDocument.recompute()

    def addPostProcessForces():
        forces = r'C:\Users\zeemarquez\Documents\FreeCAD\OpenFoam\forces.txt'
        f = open(forces,'r')
        txt = f.read()
        f.close()
        
        output_path = r'C:\Users\ZEEMAR~1\AppData\Local\Temp\case\system\forces'
        o = open(output_path,'w')
        o.write(txt)
        o.close()
        
        
        controlDict = open(r'C:\Users\zeemarquez\Documents\FreeCAD\OpenFoam\controlDict.txt','r')
        txt = controlDict.read()
        controlDict.close()
        
        output_path = r'C:\Users\ZEEMAR~1\AppData\Local\Temp\case\system\controlDict'
        o = open(output_path,'w')
        o.write(txt)
        o.close()

    def runSolver():
        
        FreeCAD.ActiveDocument.getObject('CfdSolver').MaxIterations = 1

        FreeCAD.ActiveDocument.getObject('CfdSolver').MaxIterations = 10

        FreeCAD.ActiveDocument.getObject('CfdSolver').MaxIterations = 100

        FreeCAD.ActiveDocument.getObject('CfdSolver').MaxIterations = 1000

        FreeCAD.ActiveDocument.CfdSolver.Proxy.case_writer = CfdCaseWriterFoam.CfdCaseWriterFoam(FreeCAD.ActiveDocument.CfdAnalysis)
        writer = FreeCAD.ActiveDocument.CfdSolver.Proxy.case_writer
        writer.writeCase()
        
        SimulationCFD.addPostProcessForces()

        analysis_object = FreeCAD.ActiveDocument.CfdAnalysis
        solver_object = FreeCAD.ActiveDocument.CfdSolver
        working_dir = CfdTools.getOutputPath(analysis_object)
        case_name = solver_object.InputCaseName
        solver_directory = os.path.abspath(os.path.join(working_dir, case_name))
        proxy = FreeCAD.ActiveDocument.CfdSolver.Proxy
        proxy.running_from_macro = True

        import CfdRunnableFoam
        solver_runner = CfdRunnableFoam.CfdRunnableFoam(analysis_object, solver_object)

        cmd = solver_runner.get_solver_cmd(solver_directory)
        FreeCAD.Console.PrintMessage(' '.join(cmd)+'\n')
        env_vars = solver_runner.getRunEnvironment()

        solver_process = CfdConsoleProcess.CfdConsoleProcess(stdoutHook=solver_runner.process_output)
        solver_process.start(cmd, env_vars=env_vars)
        solver_process.waitForFinished()
        
        return solver_process.exitCode() == 0
        
    def createResults(folderpath = r"C:\Users\zeemarquez\KeelBulb"):
        source_dir = r"C:\Users\ZEEMAR~1\AppData\Local\Temp\case"
        destination_dir = folderpath + "\\Results\\case"
        shutil.copytree(source_dir, destination_dir)

        forces_path = destination_dir + r'\postProcessing\forces\0\force.dat'
        f = open(forces_path,'r')
        lines = f.readlines()

        total_xyz       = [float(x) for x in lines[-1].split('\t')[1].replace('(','').replace(')','').replace('\n','').split(' ')]
        pressure_xyz    = [float(x) for x in lines[-1].split('\t')[2].replace('(','').replace(')','').replace('\n','').split(' ')]
        visc_xyz        = [float(x) for x in lines[-1].split('\t')[3].replace('(','').replace(')','').replace('\n','').split(' ')]

        #f.close()

        o = open(folderpath + "\\Results\\forces.csv",'w')

        o.write(';x;y;z\n')
        o.write('Total (-);' + ';'.join([str(k) for k in total_xyz]) + '\n')
        o.write('Total (N);' + ';'.join([str(k*1000) for k in total_xyz]) + '\n')
        o.write('Pressure (-);' + ';'.join([str(k) for k in pressure_xyz]) + '\n')
        o.write('Pressure (N);' + ';'.join([str(k*1000) for k in pressure_xyz]) + '\n')
        o.write('Viscous (-);' + ';'.join([str(k) for k in visc_xyz]) + '\n')
        o.write('Viscous (N);' + ';'.join([str(k*1000) for k in visc_xyz]) + '\n')
        o.write('Drag (N);' + str(total_xyz[0]*1000) + '\n')
        o.write('Lift (N);' + str(total_xyz[2]*1000) + '\n')

        o.close()

class SolidGenerator:

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
        splines = SolidGenerator.createBulbSplines(filepath)
        FreeCAD.ActiveDocument.recompute()
        SolidGenerator.createLoft(splines=splines,name='BulbLoft')
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
        SolidGenerator.createLoft(splines=splines,name='KeelLoft')
        FreeCAD.ActiveDocument.recompute()
        
    def createCube(Length, Width, Height, x, y, z):
        FreeCAD.ActiveDocument.addObject("Part::Box","Cube")
        FreeCAD.ActiveDocument.ActiveObject.Length = str(Length) + ' mm'
        FreeCAD.ActiveDocument.ActiveObject.Width = str(Width) + ' mm'
        FreeCAD.ActiveDocument.ActiveObject.Height = str(Height) + ' mm'
        FreeCAD.ActiveDocument.ActiveObject.Placement = FreeCAD.Placement(FreeCAD.Vector(x,y,z),FreeCAD.Rotation(FreeCAD.Vector(0,0,1),0))
        FreeCAD.ActiveDocument.recompute()

    def createCut():
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

    def createKeelBulbCut(filepath):
        L_cube, W_cube, H_cube = 4000, 500, 1200
        x_cube,y_cube,z_cube = -1000,0,-400

        SolidGenerator.createBulb(filepath)
        SolidGenerator.createKeel()
        SolidGenerator.createCube(L_cube,W_cube,H_cube,x_cube,y_cube,z_cube)
        SolidGenerator.createCut()




def runAllDesigns(folderpath = r'C:\Users\zeemarquez\Documents\FreeCAD\BulbDataFiles'):
    
    files = [x for x in os.listdir(folderpath) if ('.csv' == x[-4:])]
    for f in files:
        bulb_path = os.path.join(folderpath,f)
        bulb_name = f[:-4]
        out_folder = r'C:\Users\zeemarquez\KeelBulb' + '\\' + bulb_name 

        try:
            os.mkdir(out_folder)
        except:
            pass

        design = Automate.Design(bulb_path,out_folder)
        design.start()





    