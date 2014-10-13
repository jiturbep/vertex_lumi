theApp.EvtMax = 10

# Necessary and sufficient to read Athena pool files!
import AthenaPoolCnvSvc.ReadAthenaPool
import os

svcMgr.EventSelector.InputCollections = ['vertex_aod_pu20.root']
svcMgr.EventSelector.SkipEvents = 0

from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from GaudiSvc.GaudiSvcConf import THistSvc
svcMgr += THistSvc()

svcMgr.THistSvc.Output += ["outfile DATAFILE='VtxTmTest.root' OPT='RECREATE'"]
svcMgr.THistSvc.PrintAll = True

from InDetTruthVertexValidation.InDetTruthVertexValidationConf import VertexTruthMatchAlgorithm as VertexTruthMatchAlgorithm
from InDetTruthVertexValidation.InDetTruthVertexValidationConf import VertexTruthMatchAthAlgTool as VertexTruthMatchAthAlgTool

# Set up the algorithm
alg = VertexTruthMatchAlgorithm( 'VtxTmTest' )
alg.OutputLevel = DEBUG

# Setup for the Vtx truth match tool
ToolSvc = Service( 'ToolSvc' )
ToolSvc += VertexTruthMatchAthAlgTool( 'MyVtxTmTool' )
ToolSvc.MyVtxTmTool.OutputLevel = DEBUG
ToolSvc.MyVtxTmTool.SelectTruthTracks=False
ToolSvc.MyVtxTmTool.PtMin = 350.
ToolSvc.MyVtxTmTool.EtaMax = 2.6
ToolSvc.MyVtxTmTool.RemoveEmptyEvents = True
ToolSvc.MyVtxTmTool.RemoveDummyEvents = True
ToolSvc.MyVtxTmTool.AddOnlyFirstVertex = True


# Assign the tool to the algorithm
alg.VtxTm = ToolSvc.MyVtxTmTool

# Add the algorithm to the sequence
job += alg

THistSvc.Dump = True
