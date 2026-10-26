"""D-only resumed download and complete gzip validation. No E access or automatic archive resume."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.request,json,hashlib,gzip,os,time,shutil,threading

ROOT=Path('D:/CD8_Immunotherapy_20261023/analysis_extensions')
PAUSE_FLAG=ROOT/'E_DRIVE_REPAIR_PAUSE.json'
if not PAUSE_FLAG.exists() or not json.loads(PAUSE_FLAG.read_text()).get('D_only_work_authorized'):
    raise SystemExit('D-only resume authorization required')
STAGE=ROOT/'download_repair'/'fresh_download'
STAGE.mkdir(parents=True,exist_ok=True)
NAME='GSE179994_all.Tcell.rawCounts.rds.gz'
URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/'+NAME
SIZE=441470188; CHUNK=4*1024**2
DEST=STAGE/NAME
STATUS=ROOT/'independent_download_status.json'
LOCK=threading.Lock(); DONE={}; START=time.monotonic()

def status(state,**extra):
    with LOCK:
        n=sum(DONE.values())
        row=dict(state=state,file=NAME,pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),downloaded_bytes=n,expected_bytes=SIZE,percent=round(n/SIZE*100,2),temporary_directory=str(STAGE),destination=str(DEST),verified_chunks=len(DONE),session_MiB_per_second=round(n/max(1,time.monotonic()-START)/1024**2,3),**extra)
        tmp=STATUS.with_suffix('.json.tmp'); tmp.write_text(json.dumps(row,indent=2)); os.replace(tmp,STATUS)

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024**2),b''): h.update(block)
    return h.hexdigest()

def fetch(index):
    start=index*CHUNK; end=min(SIZE,start+CHUNK)-1
    output=STAGE/f'chunk_{index:04d}.bin'; receipt=output.with_suffix('.json')
    if output.exists() and receipt.exists():
        old=json.loads(receipt.read_text())
        if output.stat().st_size==end-start+1 and sha(output)==old['sha256']:
            with LOCK: DONE[index]=end-start+1
            status('redownloading_to_D_verified_chunks'); return
    for attempt in range(1,9):
        try:
            request=urllib.request.Request(URL+f'?fresh_start={start}&fresh_end={end}&attempt={attempt}',headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity','Cache-Control':'no-cache'})
            h=hashlib.sha256(); count=0; tmp=output.with_suffix('.tmp')
            with urllib.request.urlopen(request,timeout=60) as response:
                if response.status!=206 or response.headers.get('Content-Range')!=f'bytes {start}-{end}/{SIZE}': raise ValueError('Official server returned a different range')
                with tmp.open('wb') as f:
                    while block:=response.read(256*1024):
                        count+=len(block)
                        if count>end-start+1: raise ValueError('Too many bytes')
                        h.update(block); f.write(block)
                    f.flush(); os.fsync(f.fileno())
            if count!=end-start+1: raise ValueError('Incomplete range')
            if sha(tmp)!=h.hexdigest(): raise IOError('D readback mismatch')
            os.replace(tmp,output)
            receipt.write_text(json.dumps(dict(start=start,end=end,sha256=h.hexdigest(),source=URL)))
            with LOCK: DONE[index]=count
            status('redownloading_to_D_verified_chunks')
            print(f'CHUNK {index} verified',flush=True)
            return
        except Exception as error:
            print(f'CHUNK {index} attempt {attempt}: {error!r}',flush=True)
            if attempt==8: raise
            time.sleep(min(15,2*attempt))

try:
    if shutil.disk_usage(STAGE).free < SIZE*2+512*1024**2:
        raise OSError('Insufficient D free space for safe download and assembly')
    status('resuming_D_only_download_E_remains_paused',E_io_allowed=False)
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(fetch,i) for i in range((SIZE+CHUNK-1)//CHUNK)]
        for future in as_completed(futures): future.result()
    assembled=DEST
    status('assembling_verified_chunks_on_D',E_io_allowed=False)
    with assembled.open('wb') as out:
        for i in range((SIZE+CHUNK-1)//CHUNK):
            with (STAGE/f'chunk_{i:04d}.bin').open('rb') as src:
                shutil.copyfileobj(src,out,4*1024**2)
        out.flush(); os.fsync(out.fileno())
    if assembled.stat().st_size!=SIZE: raise ValueError('Whole file size mismatch')
    status('verifying_full_gzip_crc_on_D',E_io_allowed=False)
    with gzip.open(assembled,'rb') as f:
        for block in iter(lambda:f.read(4*1024**2),b''): pass
    actual=sha(assembled)
    row=dict(file=NAME,path=str(assembled),url=URL,bytes=SIZE,sha256=actual,gzip_CRC_passed=True,storage_drive='D',E_io_performed=False)
    (STAGE/'D_only_final_verification.json').write_text(json.dumps(row,indent=2))
    (ROOT/'independent_expression_D_manifest.json').write_text(json.dumps(row,indent=2))
    status('validation_expression_complete_on_D_E_remains_paused',gzip_CRC_passed=True,sha256=actual,E_io_allowed=False,auto_resume_archive=False)
    print('COMPLETE: D-only download and full gzip validation; E untouched; no automatic archive resume',flush=True)
except Exception as error:
    status('D_only_download_failed_preserved_all_files',error=repr(error),E_io_allowed=False)
    print(repr(error),flush=True)
    raise
