from pathlib import Path
import concurrent.futures,json,urllib.request
ROOT=Path(r'D:\CD8_Immunotherapy_20261023'); OUT=ROOT/'data'/'raw'/'publication_supplements'
BASE='https://ars.els-cdn.com/content/image/1-s2.0-S0092867418313941-'
def get(name):
    url=BASE+name; dest=OUT/('SadeFeldman_'+name)
    try:
        b=urllib.request.urlopen(url,timeout=45).read()
        if b[:2]!=b'PK': raise ValueError('Not an xlsx file')
        dest.write_bytes(b); print('CELL_SUPPLEMENT',name,len(b),flush=True)
        return {'url':url,'path':str(dest),'bytes':len(b),'status':'downloaded'}
    except Exception as e:
        print('CELL_SUPPLEMENT_ERROR',name,str(e)[:160],flush=True)
        return {'url':url,'status':'error','error':str(e)}
rows=list(concurrent.futures.ThreadPoolExecutor(3).map(get,['mmc1.xlsx','mmc2.xlsx','mmc4.xlsx']))
(ROOT/'data'/'metadata'/'SadeFeldman_supplement_sources.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
