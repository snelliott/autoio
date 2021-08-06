import automol

ich = automol.smiles.inchi('C(=O)=O')
geo = automol.inchi.geometry(ich)
zma = automol.geom.zmatrix(geo)
zma_str = automol.zmat.string(zma)
print(zma_str)
