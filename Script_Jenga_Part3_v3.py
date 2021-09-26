
# Import all required modules in Abaqus CAE
from abaqus import *
from abaqusConstants import *
from caeModules import *
import numpy as np

def jenga_model (jobname='Jenga',L=60.,H=10.,W=20.,elem_size=10.,N=6,density=2e-09,shape=1,
                  modules=100,poisson=.3,mu_bricks=0.2,mu_bricks_ground=.8,
                  solve_model=True):
      
    
    # New model
    Mdb()
    mymodel = mdb.models['Model-1']

    p=brick(mymodel,L,H,W,shape,elem_size)

     # Material
    mymaterial=mymodel.Material(name='Material-1')
    mymaterial.Density(table=((density, ), ))
    mymaterial.Elastic(table=((modules, poisson), ))

     # Section
    mysection = mymodel.HomogeneousSolidSection(name='Section-Brick', material=mymaterial.name, thickness=None)

     # Section assignment
    p.SectionAssignment(region=p.sets['ALL'], sectionName=mysection.name, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
        
    # Assembly
    a = mymodel.rootAssembly
    
    i = 0
    for level in range(N):

        # Instance 1
        b1 = a.Instance(name='Brick-'+str(i+1), part=p, dependent=ON)
        # Instance 2
        b2 = a.Instance(name='Brick-'+str(i+2), part=p, dependent=ON)
        
        if level % 2 == 0:  # Even level (par)   
            
            a.translate(instanceList=(b2.name, ), vector=(L-W, 0.0, 0.0))
        else: # odd level
             # Instance 3
            b1=a.Instance(name='Brick-'+str(i+1), part=p, dependent=ON)
            a.rotate(instanceList=(b1.name, ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(
                0.0, 1.0, 0.0), angle=90.0)
            a.translate(instanceList=(b1.name,), vector=(0.0, 0, W))
            
             # Instance 4
            b2=a.Instance(name='Brick-'+str(i+2), part=p, dependent=ON) 
            a.rotate(instanceList=(b2.name,), axisPoint=(0.0, 0.0, 0.0), axisDirection=(
                0.0, 1.0, 0.0), angle=90.0)
            a.translate(instanceList=(b2.name, ), vector=(0.0, 0, L))

        a.translate(instanceList=(b1.name,b2.name), vector=(0.0, level*H, 0.0))
        # i = i + 2
        i += 2
      
      
    # Ground  
    W_ground = 2.*N*L 
    s1 = mymodel.ConstrainedSketch(name='ground', sheetSize=200.0)
    s1.rectangle(point1=(0.0, 0.0), point2=(W_ground, W_ground))
    pg = mymodel.Part(name='Ground', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    pg.BaseShell(sketch=s1)

    # Set
    pg.Set(faces=pg.faces, name='ALL')
    pg.Surface(side1Faces=pg.faces, name='SURF')

    # Mesh
    pg.seedPart(size=W_ground, deviationFactor=0.1, minSizeFactor=0.1)
    pg.generateMesh()

    # Section shell
    mymodel.HomogeneousShellSection(name='Section-Shell', 
        preIntegrate=OFF, material='Material-1', thicknessType=UNIFORM, 
        thickness=1.0, thicknessField='', nodalThicknessField='', 
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
        integrationRule=SIMPSON, numIntPts=5)
        
    pg.SectionAssignment(region=p.sets['ALL'], sectionName='Section-Shell', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

    # Assembly: Ground
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    ag = a.Instance(name='Ground', part=pg, dependent=ON)

    a.rotate(instanceList=(ag.name, ), axisPoint=(0.0, 0.0, 0.0), 
        axisDirection=(1.0, 0.0, 0.0), angle=90.0)
    a.translate(instanceList=(ag.name, ), vector=(-W_ground/2, 0.0, -W_ground/2))

    # STEP

    mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', 
        timePeriod=2.0, improvedDtMethod=ON)

    # Interaction 


    int_prop_bricks=mymodel.ContactProperty('IntProp-Bricks')
    int_prop_bricks.NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)
    int_prop_bricks.TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((mu_bricks, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
     
    int_prop_bricks_ground=mymodel.ContactProperty('IntProp-Bricks_ground')
    int_prop_bricks_ground.NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)
    int_prop_bricks_ground.TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((mu_bricks_ground, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)


     # Interaction definition    
    Contact_exp=mymodel.ContactExp(name='Int-1', createStepName='Step-1')
    Contact_exp.includedPairs.setValuesInStep(
        stepName='Step-1', useAllstar=ON)

    Contact_exp.contactPropertyAssignments.appendInStep(
        stepName='Step-1', assignments=((GLOBAL, SELF, int_prop_bricks.name), ))

    pairs=[]
    for i in range(2*N):
        aux=(a.instances['Brick-'+str(i+1)].surfaces['Surf'],a.instances['Ground'].surfaces['SURF'],int_prop_bricks_ground.name)
        pairs.append(aux)

    Contact_exp.contactPropertyAssignments.appendInStep(stepName='Step-1', assignments=pairs)

     #  Boundary Conditions

    mymodel.Gravity(name='Gravity', createStepName='Step-1', 
        comp2=-9800.0, distributionType=UNIFORM, field='')

    region = a.instances['Ground'].sets['ALL']
    mymodel.EncastreBC(name='Encastre-brick', 
        createStepName='Step-1', region=region, localCsys=None)
    myamplitude=mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((
        0.0, 0.0), (0.05, 1.0), (1.05, 1.0), (1.5, 0.0)))
    a = mymodel.rootAssembly
    region = a.instances['Brick-3'].sets['ALL']
    mymodel.VelocityBC(name='Velocity', createStepName='Step-1', 
        region=region, v1=L, v2=UNSET, v3=UNSET, vr1=UNSET, vr2=UNSET, 
        vr3=UNSET, amplitude=myamplitude.name, localCsys=None, distributionType=UNIFORM, 
        fieldName='')

     # Results
    mymodel.fieldOutputRequests['F-Output-1'].setValues(variables=(
        'S', 'E', 'U', 'V', 'A'), numIntervals=200)
    mymodel.historyOutputRequests['H-Output-1'].setValues(variables=(
        'ALLFD', 'ALLKE', 'ALLSE'))
        
      # Visualization
    session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
    session.viewports['Viewport: 1'].enableMultipleColors()
    session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
    cmap=session.viewports['Viewport: 1'].colorMappings['Part instance']
    session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
    session.viewports['Viewport: 1'].disableMultipleColors()
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        adaptiveMeshConstraints=OFF)
        
       # Job  
        
    myjob=mdb.Job(name=jobname, model=mymodel.name, description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB)
       # Save model
    mdb.saveAs(pathName='C:/temp/2.Practice/Jenga')


    if solve_model:
        myjob.submit()
        myjob.waitForCompletion()
    return myjob
    
def brick(mymodel,L,H,W,shape,elem_size):

    if (shape==0):
    # rectangular
        s = mymodel.ConstrainedSketch(name='brick', sheetSize=200.0)
        s.rectangle(point1=(0.0, 0.0), point2=(W, H))

    elif (shape==1):
    # L shape
        s = mymodel.ConstrainedSketch(name='brick', sheetSize=200.0)
        s.Line(point1=(0.0, 0.0), point2=(W, 0.0))
        s.Line(point1=(W, 0.0), point2=(W, H/2))
        s.Line(point1=(W, H/2), point2=(W/2, H/2))
        s.Line(point1=(W/2, H/2), point2=(W/2, H))
        s.Line(point1=(W/2, H), point2=(0.0, H))
        s.Line(point1=(0.0, H), point2=(0.0, 0.0))
        
    p = mymodel.Part(name='Brick', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=L)

     # Set and surface 
    p.Set(cells=p.cells, name='ALL')
    p.Surface(side1Faces=p.faces, name='Surf')

     # Mesh
    p.seedPart(size=elem_size, deviationFactor=0.1, minSizeFactor=0.1)
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, 
        kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
        hourglassControl=DEFAULT, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    p.setElementType(regions=p.sets['ALL'], elemTypes=(elemType1, elemType2, 
        elemType3))
    p.generateMesh()

    return p
  

def export_displacement(odbName,instancename,node=1):

    session.mdbData.summary()
    odb = session.openOdb(name=odbName+'.odb')
    myview = session.viewports[session.currentViewportName]
    odb = session.openOdb(name='Jenga_L.odb')
    session.viewports['Viewport: 1'].setValues(displayedObject=odb)
    session.viewports['Viewport: 1'].makeCurrent()
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)

    datalist = session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
        NODAL, ((COMPONENT, 'U1'), (COMPONENT, 'U2'), (COMPONENT, 'U3'), )), ), 
        nodeLabels=((instancename, (node, )), ))

    data_ux,data_uy,data_uz=datalist

    t,ux=zip(*data_ux.data)
    t,uy=zip(*data_uy.data)
    t,uz=zip(*data_uz.data)

    np.savetxt(odbName+'Displacement_'+instancename+'_1.txt',np.transpose([ux,uy,uz]),fmt='%12.3f')
    return odb

def Screenshots(odbName):

    odb = session.openOdb(name=odbName+'.odb')
    n_frames=len(odb.steps['Step-1'].frames)    
    
    myview = session.viewports[session.currentViewportName]
    a = mdb.models['Model-1'].rootAssembly
    myview.setValues(displayedObject=a)
    session.mdbData.summary()

    myview.setValues(displayedObject=odb)
    myview.makeCurrent()
    myview.odbDisplay.display.setValues(plotState=(DEFORMED, ))
    myview.enableMultipleColors()
    myview.setColor(initialColor='#BDBDBD')
    cmap=session.viewports['Viewport: 1'].colorMappings['Part instance']
    myview.setColor(colorMapping=cmap)
    myview.disableMultipleColors()

    for i in range(n_frames):    

        myview.odbDisplay.setFrame(
            step=0, frame=i)
        session.printToFile(fileName=odbName+'Frame_%03d' % (i), format=PNG, canvasObjects=(
            session.viewports[session.currentViewportName], ))
    return odb
  
mus=[0.1, 0.3, 0.5, 0.7]

for i in range(len(mus)):
    name='Jenga_L_'+str(i+1)
    jenga_model (jobname='Jenga_L',L=40.,H=10.,W=10.,elem_size=5.,N=6,density=2e-09,shape=1,
                      modules=100,poisson=.3,mu_bricks=mus[i],mu_bricks_ground=.8,
                      solve_model=True)
     
     
     # Post-process: Read Displacement
     export_displacement(name,instancename='Brick-10',node=1)
     # Post-process: Print Frames
     odb=Screenshots(odbName)
     
     odb.close()