import FreeCAD
import SimulationCFD
import SolidGenerator
import shutil
import os
import time
import numpy as np
import FreeCADGui as Gui

class Design:

    default_mesh_size = 23 #mm
    default_mesh_refin = 0.5
    default_inlet_v = 5000 #mm/s

    def __init__(self, bulbpath, outputfolder, inlet_v = default_inlet_v, bulbid = None):

        self.bulbpath = bulbpath
        self.inlet_v = inlet_v
        self.final_mesh_size = None
        self.drag_result = None
        self.elapsed_time = None
        

        if bulbid == None:
            self.bulbId = bulbpath.split('\\')[-1][:-4]
        else:
            self.bulbId = bulbid

        self.outputfolder = outputfolder + '\\' + self.bulbId
        self.docpath = self.outputfolder + '\\' + 'FreeCADdoc.FCStd'
        self.doc = None
        
    def start(self):
        self.createFolder(self.outputfolder)
        FreeCAD.Console.PrintMessage('\n\n-----------------------------')
        print(self.bulbId)
        FreeCAD.Console.PrintMessage('-----------------------------\n')
        FreeCAD.Console.PrintMessage('Opening new file...')
        self.openNewDocument()
        FreeCAD.Console.PrintMessage('Genearting geometry...')
        self.genGeometry()
        self.doc.save()
        FreeCAD.Console.PrintMessage('Starting simulation...')
        self.start_time = time.time()
        self.startSimulation()
        FreeCAD.Console.PrintMessage('Simulation ended!\n')
        FreeCAD.Console.PrintMessage('Saving Document...')
        self.doc.save()
        FreeCAD.closeDocument(self.doc.Name)
        FreeCAD.Console.PrintMessage('Finished !')
        FreeCAD.Console.PrintMessage('-----------------------------')
        FreeCAD.Console.PrintMessage('-----------------------------\n\n')

    def createFolder(self, folder):
        try:
            os.mkdir(folder)
        except:
            pass
    def openNewDocument(self):
        self.doc = FreeCAD.newDocument()
        self.doc.saveAs(self.docpath)
    def genGeometry(self, onlyKeel = False):

        if not onlyKeel:
            SolidGenerator.createBulb(self.bulbpath)

        SolidGenerator.createKeel()
        SolidGenerator.createCube(4000, 700, 1700, -1000, 0, -480)
        SolidGenerator.createCut(onlyKeel)

    def startSimulation(self):

        mesh_size = Design.default_mesh_size

        SimulationCFD.createAnalysis()
        SimulationCFD.createMeshCase(size=mesh_size,refinement=Design.default_mesh_refin)
        SimulationCFD.runMesh()
        SimulationCFD.createBoundaries(inlet_v=self.inlet_v)

        cfd_valid = SimulationCFD.runSolver()

        if not cfd_valid:
            for n in range(5):
                mesh_size -= 1
                SimulationCFD.createMeshCase(size=mesh_size,refinement=Design.default_mesh_refin)
                SimulationCFD.runMesh()
                cfd_valid = SimulationCFD.runSolver()
                if cfd_valid:
                    FreeCAD.Console.PrintMessage("Simulation ended successfully | mesh size:" + str(mesh_size))
                    self.final_mesh_size = mesh_size
                    break
        
        if not cfd_valid:
            raise Exception('Error running simualtion')
        else:
            drag = SimulationCFD.createResults(self.outputfolder)
            
            self.end_time = time.time()
            self.elapsed_time = self.end_time - self.start_time

            o = open(self.outputfolder + "\\Results\\simulation_time.txt",'w')
            o.write("Elapsed time (s):" + str(self.elapsed_time))
            o.close()
            self.drag_result = drag

    def runSimulation(self, mesh_size):

        cfd_valid = SimulationCFD.runSolver()

        result = None

        if not cfd_valid:
            for n in range(5):
                mesh_size -= 1
                SimulationCFD.createMeshCase(size=mesh_size,refinement=Design.default_mesh_refin)
                SimulationCFD.runMesh()
                cfd_valid = SimulationCFD.runSolver()
                if cfd_valid:
                    FreeCAD.Console.PrintMessage("Simulation ended successfully | mesh size:" + str(mesh_size))
                    self.final_mesh_size = mesh_size
                    break
        
        if not cfd_valid:
            raise Exception('Error running simualtion')
        else:
            forces = SimulationCFD.getForces()
            drag = forces['TotalX']
            self.drag_result = drag
            result = forces

        return result
        


def runAllDesigns(bulbsfolder = r'C:\Users\zeemarquez\Documents\FreeCAD\BulbDataFiles', designsfolder = r'C:\Users\zeemarquez\Documents\FreeCAD\Simulations Results\Simulations01', delete_case = True):

    bulbfiles = [ file for file in os.listdir(bulbsfolder) if file[-4:] == '.csv']

    resultspath = designsfolder + "\\results.csv"
    resultsfile = open(resultspath,'w')
    resultsfile.write("n;ID;Target mesh;Actual mesh;Simulation time (s);Drag Force (N)")
    resultsfile.close()
    n = 0

    for file in bulbfiles:

        bulbfile = bulbsfolder + '\\' + file
        design = Design(bulbfile, designsfolder)
        design.start()

        resultsfile = open(resultspath,'a')
        resultsfile.write(";".join([str(n), design.bulbId ,str(design.default_mesh_size),str(design.final_mesh_size),str(design.elapsed_time),str(design.drag_result)]))
        resultsfile.write('\n')
        resultsfile.close()

        if delete_case:
            case_folder = design.outputfolder + "\\Results\\case"
            shutil.rmtree(case_folder)

        n += 1

    FreeCAD.Console.PrintMessage('\n\n\n*******************************')
    FreeCAD.Console.PrintMessage('FINISHED CORRECTLY !!!')
    FreeCAD.Console.PrintMessage('*******************************')

def runAllDesigns(bulbsfolder, designsfolder, delete_case = True):

    bulbfiles = [ file for file in os.listdir(bulbsfolder) if file[-4:] == '.csv']

    resultspath = designsfolder + "\\results.csv"
    resultsfile = open(resultspath,'w')
    resultsfile.write("n;ID;Target mesh;Actual mesh;Simulation time (s);Drag Force (N)")
    resultsfile.close()
    n = 0

    for file in bulbfiles:

        bulbfile = bulbsfolder + '\\' + file
        design = Design(bulbfile, designsfolder)
        design.start()

        resultsfile = open(resultspath,'a')

        resultsfile.write(";".join(
            [str(n), design.bulbId,
            str(design.default_mesh_size),
            str(design.final_mesh_size),
            str(design.elapsed_time),
            str(design.drag_result)]))

        resultsfile.write('\n')
        resultsfile.close()

        if delete_case:
            case_folder = design.outputfolder + "\\Results\\case"
            shutil.rmtree(case_folder)

        n += 1

    FreeCAD.Console.PrintMessage('\n\n\n*******************************')
    FreeCAD.Console.PrintMessage('FINISHED CORRECTLY !!!')
    FreeCAD.Console.PrintMessage('*******************************')

def runMeshConvergence(bulbfile, outputfolder, meshSizes):

    resultspath = outputfolder + "\\mesh_convergence.csv"
    resultsfile = open(resultspath,'w')
    resultsfile.write("n;ID;Target mesh;Actual mesh;Simulation time (s);Drag Force (N)")
    resultsfile.close()

    n = 0
    for meshsize in meshSizes:

        Design.default_mesh_size = meshsize
        design = Design(bulbfile, outputfolder, bulbid= "mesh_" + str(meshsize))
        design.start()
        final_mesh = design.final_mesh_size
        drag = design.drag_result

        resultsfile = open(resultspath,'a')
        resultsfile.write(";".join([str(n), design.bulbId ,str(meshsize),str(final_mesh),str(design.elapsed_time),str(drag)]))
        resultsfile.write('\n')
        resultsfile.close()

        case_folder = design.outputfolder + "\\Results\\case"
        shutil.rmtree(case_folder)

        n+=1

    FreeCAD.Console.PrintMessage('\n\n\n*******************************')
    FreeCAD.Console.PrintMessage('FINISHED CORRECTLY !!!')
    FreeCAD.Console.PrintMessage('*******************************')

def runVelocityAnalysis(bulbfile, outputfolder, onlyKeel = False, startV = 500, endV=5000, step=500):

    velocities = list(range(startV,endV,step))

    resultspath = outputfolder + "\\velocity_analysis.csv"
    resultsfile = open(resultspath,'w')
    resultsfile.write("n;v(mm/s);Simulation time (s);Drag Force (N)\n")
    resultsfile.close()

    design = Design(bulbfile, outputfolder, bulbid= "velocity_analysis")
    design.createFolder(design.outputfolder)

    design.openNewDocument()
    design.genGeometry(onlyKeel=onlyKeel)
    design.doc.save()
    design.start_time = time.time()
    Gui.SendMsgToActiveView("ViewFit")


    mesh_size = 20

    SimulationCFD.createAnalysis()
    SimulationCFD.createMeshCase(size=mesh_size,refinement=Design.default_mesh_refin,onlyKeel=onlyKeel)
    SimulationCFD.runMesh()
    SimulationCFD.createBoundaries()

    
     
    n = 0

    for v in velocities:
        
        start_time = time.time()
        
        SimulationCFD.changeInletVelocity(v)

        forces = design.runSimulation(mesh_size)
        design.doc.save()

        elapsed_time = time.time() - start_time
        drag = forces['TotalX']

        resultsfile = open(resultspath,'a')
        resultsfile.write(";".join( [str(n),str(v),str(elapsed_time),str(drag)] ))
        resultsfile.write('\n')
        resultsfile.close()

        '''
        FreeCAD.Console.PrintMessage('\n\n')
        FreeCAD.Console.PrintMessage('-----------------------------')
        FreeCAD.Console.PrintMessage('Velocity (mm/s):'+str(v))
        FreeCAD.Console.PrintMessage('-----------------------------\n\n')


        for f in forces:
            FreeCAD.Console.PrintMessage(f+':\t'+ forces[f])

        FreeCAD.Console.PrintMessage('-----------------------------')
        FreeCAD.Console.PrintMessage('-----------------------------\n\n')

        '''
        

        n+=1






def genGeometry(bulbpath):
    SolidGenerator.createBulb(bulbpath)
    SolidGenerator.createKeel()
    SolidGenerator.createCube(4000, 700, 1700, -1000, 0, -480)
    SolidGenerator.createCut()

