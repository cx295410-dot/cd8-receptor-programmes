"""Use rdata's XDR object parser with streamed, disk-backed numeric vectors."""
from pathlib import Path
import sys,gzip,json,time
import numpy as np

ROOT=Path('D:/CD8_Immunotherapy_20261023/analysis_extensions')
sys.path.insert(0,str(ROOT/'packages'))
import rdata
from rdata.parser._xdr import ParserXDR
from rdata.parser._parser import Parser,RObjectType

class StreamParser(ParserXDR):
    def __init__(self,stream,out,threshold=1000000):
        Parser.__init__(self)
        self.file=stream; self.out=Path(out); self.out.mkdir(parents=True,exist_ok=True)
        self.threshold=threshold; self.large=[]
    def _parse_array_values(self,dtype,length):
        dtype=np.dtype(dtype)
        if length<0: raise ValueError('Negative array length unsupported')
        if length<self.threshold:
            payload=self.file.read(length*dtype.itemsize)
            if len(payload)!=length*dtype.itemsize: raise EOFError('Incomplete XDR array')
            return np.frombuffer(payload,dtype=dtype.newbyteorder('>')).astype(dtype,copy=True)
        target=self.out/f'large_array_{len(self.large):02d}.npy'
        arr=np.lib.format.open_memmap(target,mode='w+',dtype=dtype,shape=(length,))
        self.large.append(dict(path=str(target),dtype=str(dtype),length=length))
        for start in range(0,length,262144):
            stop=min(length,start+262144)
            payload=self.file.read((stop-start)*dtype.itemsize)
            if len(payload)!=(stop-start)*dtype.itemsize: raise EOFError('Incomplete large array')
            arr[start:stop]=np.frombuffer(payload,dtype=dtype.newbyteorder('>'))
        arr.flush()
        print('PARSED_NUMERIC',json.dumps(self.large[-1]),flush=True)
        return arr

def parse(path,out,threshold=1000000):
    with open(path,'rb') as probe: compressed=probe.read(2)==b'\x1f\x8b'
    opener=gzip.open if compressed else open
    with opener(path,'rb') as f:
        if f.read(2)!=b'X\n': raise ValueError('Only XDR RDS supported')
        parser=StreamParser(f,out,threshold)
        result=parser.parse_all(); parser.check_complete()
    return result,parser.large

def resolve(obj):
    while obj.info.type==RObjectType.REF: obj=obj.referenced_object
    return obj

def attrs(obj):
    node=obj.attributes; result={}
    while node is not None and node.info.type not in {RObjectType.NILVALUE,RObjectType.NIL}:
        if node.info.type!=RObjectType.LIST: raise TypeError(node.info.type)
        tag=resolve(node.tag)
        name=resolve(tag.value).value.decode()
        result[name]=resolve(node.value[0])
        node=resolve(node.value[1])
    return result

def string_vector(obj):
    return [resolve(x).value.decode('utf-8') for x in obj.value]

def selftest():
    from numpy.testing import assert_equal
    out=ROOT/'independent_GSE179994'/'reader_test'
    inputs=['test_vector__xdr__version_3.rds','test_full_named_matrix__xdr__version_3.rds','test_s4__xdr__version_3.rds','test_s4__xdr__version_2.rds']
    generated=Path(rdata.TESTDATA_PATH)/'generated'
    for index,name in enumerate(inputs):
        path=generated/name
        standard=rdata.parser.parse_file(path)
        streamed,_=parse(path,out/str(index),threshold=2)
        assert_equal(str(standard),str(streamed))
    # An independent roundtrip exercises large disk-backed numeric values and NAs.
    path=out/'roundtrip.rds'
    original={'ints':np.array([0,3,2147483647,-2147483648],dtype=np.int32),'floats':np.array([0.,1.5,np.nan,1e20]),'labels':np.array(['g1','g2'])}
    rdata.write_rds(path,original)
    streamed,large=parse(path,out/'roundtrip',threshold=2)
    conventional=rdata.read_rds(path)
    converted=rdata.conversion.convert(streamed)
    for key in original: assert_equal(converted[key],conventional[key])
    report=dict(reference_XDR_files_passed=len(inputs),numeric_NA_roundtrip_passed=True,disk_backed_test_vectors=len(large))
    (out/'verification.json').write_text(json.dumps(report,indent=2))
    print('READER_TESTS_PASS',json.dumps(report),flush=True)

def main():
    selftest()
    out=ROOT/'independent_GSE179994'/'matrix'
    source=ROOT/'download_repair'/'fresh_download'/'GSE179994_all.Tcell.rawCounts.rds.gz'
    out.mkdir(exist_ok=True)
    start=time.monotonic()
    result,large=parse(source,out)
    if result.object.info.type!=RObjectType.S4: raise TypeError('Expected S4 dgCMatrix')
    slots=attrs(result.object)
    cls=string_vector(slots['class'])
    if cls!=['dgCMatrix']: raise TypeError(cls)
    dims=slots['Dim'].value.astype(int).tolist()
    labels=slots['Dimnames'].value
    genes=string_vector(labels[0]); cells=string_vector(labels[1])
    assert dims==[len(genes),len(cells)]
    assert len(set(genes))==len(genes) and len(set(cells))==len(cells)
    for key in ['i','p','x']:
        arr=slots[key].value
        if isinstance(arr,np.ma.MaskedArray): raise ValueError('Missing sparse matrix values')
        if isinstance(arr,np.memmap): path=str(arr.filename)
        else:
            dest=out/(key+'.npy'); np.save(dest,arr); path=str(dest)
        large_slot=dict(path=path,length=len(arr),dtype=str(arr.dtype))
        if key=='p': assert len(arr)==dims[1]+1 and arr[0]==0 and np.all(np.diff(arr)>=0)
        if key=='i':
            for begin in range(0,len(arr),1000000):
                part=arr[begin:begin+1000000]; assert np.all((part>=0)&(part<dims[0]))
        if key=='x':
            for begin in range(0,len(arr),1000000):
                part=arr[begin:begin+1000000]; assert np.all(np.isfinite(part)) and np.all(part>=0) and np.all(part==np.floor(part))
        slots[key]=large_slot
    assert slots['i']['length']==slots['x']['length'] and int(np.load(slots['p']['path'])[-1])==slots['x']['length']
    info=dict(source=str(source),R_class=cls,shape=dims,genes=genes,cells=cells,sparse_arrays={k:slots[k] for k in ['i','p','x']},nonzero_entries=slots['x']['length'],parse_seconds=time.monotonic()-start,finite_nonnegative_integer_counts=True,unique_gene_and_cell_names=True,parser='rdata 1.1.0 XDR parser; streamed numeric vectors to memmap; reference tests passed',E_access_performed=False)
    (out/'matrix.json').write_text(json.dumps(info))
    print('MATRIX_READY',json.dumps({k:v for k,v in info.items() if k not in ['genes','cells','sparse_arrays']}),flush=True)

if __name__=='__main__': main()
