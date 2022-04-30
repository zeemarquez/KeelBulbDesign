import FreeCAD
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
import shutil
import os


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

def createMeshCase(size=25.0,refinement=0.5, onlyKeel = False):

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
    if onlyKeel:
        FreeCAD.ActiveDocument.MeshRefinement.ShapeRefs = [(FreeCAD.ActiveDocument.getObject('KeelLoft'), ('Solid1',))]
    else:
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

    fvSolution = open(r'C:\Users\zeemarquez\Documents\FreeCAD\OpenFoam\fvSolution.txt','r')
    txt = fvSolution.read()
    controlDict.close()
    
    output_path = r'C:\Users\ZEEMAR~1\AppData\Local\Temp\case\system\fvSolution'
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
    
    addPostProcessForces()

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
    
def getForces():

    forces_path =  r'C:\Users\ZEEMAR~1\AppData\Local\Temp\case\postProcessing\forces\0\force.dat'
    f = open(forces_path,'r')
    lines = f.readlines()

    total_xyz       = [float(x) for x in lines[-1].split('\t')[1].replace('(','').replace(')','').replace('\n','').split(' ')]
    pressure_xyz    = [float(x) for x in lines[-1].split('\t')[2].replace('(','').replace(')','').replace('\n','').split(' ')]
    visc_xyz        = [float(x) for x in lines[-1].split('\t')[3].replace('(','').replace(')','').replace('\n','').split(' ')]

    forces = {
        'TotalX':total_xyz[0]*1000,
        'TotalY':total_xyz[1]*1000,
        'TotalZ':total_xyz[2]*1000,
        'PressureX':pressure_xyz[0]*1000,
        'PressureY':pressure_xyz[1]*1000,
        'PressureZ':pressure_xyz[2]*1000,
        'ViscousX':visc_xyz[0]*1000,
        'ViscousY':visc_xyz[1]*1000,
        'ViscousZ':visc_xyz[2]*1000
    }

    return forces


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

    FreeCAD.Console.PrintMessage('\n\n--------------------------------------')
    FreeCAD.Console.PrintMessage('Drag (N):' + str(total_xyz[0]*1000))
    FreeCAD.Console.PrintMessage('--------------------------------------\n\n\n')

    return total_xyz[0]*1000

def changeInletVelocity(vx, vy = 0.0 ,vz = 0.0):

    bc = FreeCAD.ActiveDocument.CfdFluidBoundary
    bc.BoundaryType = 'inlet'
    bc.BoundarySubType = 'uniformVelocityInlet'
    bc.ThermalBoundaryType = 'fixedValue'
    bc.VelocityIsCartesian = True
    bc.Ux = str(vx) + ' mm/s'
    bc.Uy = str(vy) + ' mm/s'
    bc.Uz = str(vz) + ' mm/s'
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
