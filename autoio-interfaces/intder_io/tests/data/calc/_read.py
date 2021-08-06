import automol
import autofile
import elstruct
import intder_io


with open('output.dat') as fobj:
    OUT_STR = fobj.read()

zma = elstruct.reader.opt_zmatrix('gaussian09', OUT_STR)
geo = elstruct.reader.opt_geometry('gaussian09', OUT_STR)
hess = elstruct.reader.hessian('gaussian09', OUT_STR)
geox = automol.zmat.geometry(zma, dummy=True)

zma_str = autofile.data_types.swrite.zmatrix(zma)
geo_str = autofile.data_types.swrite.geometry(geo)
hess_str = autofile.data_types.swrite.hessian(hess)
geox_str = autofile.data_types.swrite.geometry(geox)

ithess_str = intder_io.writer.cart_hess_file(hess)


with open('a.zmat', 'w') as fobj:
    fobj.write(zma_str)
with open('a.geo', 'w') as fobj:
    fobj.write(geo_str)
with open('a.hess', 'w') as fobj:
    fobj.write(hess_str)
with open('a.geox', 'w') as fobj:
    fobj.write(geox_str)
with open('a.ithess', 'w') as fobj:
    fobj.write(ithess_str)
