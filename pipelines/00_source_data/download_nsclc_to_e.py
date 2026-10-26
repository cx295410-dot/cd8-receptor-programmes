"""Download one verified public research archive directly to E:, with resume."""
from pathlib import Path
import hashlib
import json
import msvcrt
import os
import re
import shutil
import time
import urllib.request
from datetime import datetime, timezone


# Coordinated E drive repair pause: no E access until explicit resume.
if Path('D:/CD8_Immunotherapy_20261023/analysis_extensions/E_DRIVE_REPAIR_PAUSE.json').exists():
    raise SystemExit('E drive repair pause active; do not resume without repair-complete notification.')

ROOT=Path('E:/CD8_Data/NSCLC_2026')
DEST=ROOT/'data.zip'
PART=ROOT/'data.zip.part'
STATUS=ROOT/'download_status.json'
URL='https://zenodo.org/api/records/19112009/files/data.zip/content'
SIZE=15761732140
MD5='3d27af5617f0a140a859233ed89c7a05'
ROOT.mkdir(parents=True,exist_ok=True)
START=time.monotonic()
BASE=PART.stat().st_size if PART.exists() else 0

def update(state,**extra):
    received=PART.stat().st_size if PART.exists() else (DEST.stat().st_size if DEST.exists() else 0)
    info=dict(state=state,pid=os.getpid(),updated_utc=datetime.now(timezone.utc).isoformat(),
              destination=str(DEST),partial_file=str(PART),url=URL,
              expected_bytes=SIZE,downloaded_bytes=received,percent=round(100*received/SIZE,3),
              session_MiB_per_second=round(max(0,received-BASE)/max(.01,time.monotonic()-START)/2**20,3),
              expected_md5=MD5)
    info.update(extra)
    temp=STATUS.with_suffix('.json.tmp')
    temp.write_text(json.dumps(info,indent=2,ensure_ascii=False),encoding='utf-8')
    os.replace(temp,STATUS)

def verify(path):
    if path.stat().st_size!=SIZE: raise ValueError('Unexpected archive size')
    update('verifying_md5')
    digest=hashlib.md5()
    with path.open('rb') as source:
        for block in iter(lambda:source.read(8*1024*1024),b''): digest.update(block)
    value=digest.hexdigest()
    if value!=MD5: raise ValueError('MD5 mismatch; keep file for investigation: '+value)
    return value

def run():
    if DEST.exists():
        value=verify(DEST)
        update('complete',actual_md5=value,verified=True)
        return
    offset=PART.stat().st_size if PART.exists() else 0
    if offset>SIZE: raise ValueError('Existing partial archive exceeds expected size')
    if shutil.disk_usage(ROOT).free<SIZE-offset+1024**3:
        raise OSError('Insufficient E: free space for archive and 1GiB reserve')
    update('checking_source')
    with urllib.request.urlopen('https://zenodo.org/api/records/19112009',timeout=60) as response:
        record=json.load(response)
    item=next(f for f in record['files'] if f['key']=='data.zip')
    if item['size']!=SIZE or item['checksum']!='md5:'+MD5:
        raise ValueError('Official file metadata changed; review before downloading')
    (ROOT/'source_manifest.json').write_text(json.dumps(dict(record='https://zenodo.org/records/19112009',
        file=item,retrieved_utc=datetime.now(timezone.utc).isoformat()),indent=2),encoding='utf-8')
    retries=0
    while offset<SIZE:
        end=min(SIZE-1,offset+128*1024**2-1)
        try:
            update('downloading',retry_count=retries)
            # Distinct range URL prevents cached responses for a different Range header.
            req=urllib.request.Request(URL+f'?download=1&range_start={offset}&range_end={end}',
                headers={'User-Agent':'CD8-public-data-research/1.0','Range':f'bytes={offset}-{end}',
                         'Accept-Encoding':'identity','Cache-Control':'no-cache'})
            with urllib.request.urlopen(req,timeout=90) as response:
                match=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',response.headers.get('Content-Range',''))
                if response.status!=206 or not match or tuple(map(int,match.groups()))!=(offset,end,SIZE):
                    raise IOError('Server returned an unexpected byte range; partial file preserved')
                last=time.monotonic()
                with PART.open('ab') as output:
                    while offset<=end:
                        block=response.read(min(1024**2,end-offset+1))
                        if not block: raise IOError('Connection closed before requested range completed')
                        output.write(block)
                        offset+=len(block)
                        if time.monotonic()-last>=30:
                            output.flush()
                            update('downloading',retry_count=retries)
                            last=time.monotonic()
                    output.flush()
            retries=0
            print(f'{offset}/{SIZE} bytes ({100*offset/SIZE:.2f}%)',flush=True)
            update('downloading',retry_count=retries)
        except Exception as error:
            offset=PART.stat().st_size if PART.exists() else 0
            retries+=1
            update('retrying',retry_count=retries,last_error=f'{type(error).__name__}: {error}')
            print(f'Retry {retries}: {type(error).__name__}: {error}',flush=True)
            if retries>=10: raise
            time.sleep(min(60,3*2**min(retries,4)))
    value=verify(PART)
    os.replace(PART,DEST)
    update('complete',actual_md5=value,verified=True)
    print('COMPLETE: size and official MD5 verified. '+str(DEST),flush=True)

if __name__=='__main__':
    with (ROOT/'download.lock').open('a+b') as lock:
        if lock.tell()==0:
            lock.write(b'0'); lock.flush()
        lock.seek(0)
        try:
            msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
        except OSError:
            raise SystemExit('Another downloader already holds this project lock.')
        try:
            run()
        except Exception as error:
            update('failed',last_error=f'{type(error).__name__}: {error}')
            raise
