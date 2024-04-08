from ROOT import TFile

inFile = TFile('/home/daniel/an/DstarK/2_luuksel/GenCFDebug_nopc_kStarBW50MeV_fromq_bs10000syst.root')
oFile = TFile('/home/daniel/an/DstarK/2_luuksel/GenCFDebugPatched_nopc_kStarBW50MeV_fromq_bs10000syst.root', 'recreate')

for comb in ['sc', 'oc']:
    print(comb)
    print()
    
    oFile.mkdir(comb)
    oFile.cd(comb)

    for key in inFile.Get(comb).GetListOfKeys():
        name = key.GetName()
        obj = inFile.Get(f'{comb}/{name}')


        if name == 'gCFGenSyst':
            for iPoint in range(obj.GetN() - 1):
                if obj.GetErrorY(iPoint) < 0.001:
                    x = obj.GetPointX(iPoint)
                    
                    print(f'patch error: x={x:.0f}')

                    relunc = obj.GetErrorY(iPoint + 1) / obj.GetPointY(iPoint + 1) 
                    obj.SetPointError(iPoint, obj.GetErrorX(iPoint), relunc * obj.GetPointY(iPoint))

        obj.Write()
oFile.Close()
