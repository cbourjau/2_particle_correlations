from ROOT import TProof, TFileCollection, TChain, gROOT, TFile


gROOT.LoadMacro('DebugClassesMultESA2013.C+')
fn = '/home/christian/msc/analysis/pyproof/ppb.dat'
with open(fn, 'read') as f:
    chain = TChain("tree")
    lines = f.readlines()
    print len(lines)
    i = 0
    for l in lines[]:
        chain.AddFromFile(l)
        i += 1
        if i > 20:
            f = 
            chain.Merge(f)




















