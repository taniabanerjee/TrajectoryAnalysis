import numpy as np
import time
import multiprocessing as mp
import pandas as pd

def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

def main():
    total = 0
    i=0
    spatdf = pd.read_csv('spat.csv')
    df = pd.DataFrame()
    df['timestamp'] = pd.to_datetime(spatdf['timestamp'])
    df['phase'] = spatdf['hexphase']
    df = df.drop_duplicates('phase')
    ts = np.datetime64('2019-07-15 10:46:20.2')
    mask = df['timestamp'] >= ts
    pos = np.flatnonzero(mask)
    row = df.iloc[pos[0]-1]
    while (i < total):
        results = []
        pool = mp.Pool(10)

        # Prepare data
        np.random.RandomState(100)
        arrsize = np.random.randint(100, 5000)
        arr = np.random.randint(0, 10, size=[arrsize, 5])
        data = arr.tolist()
    
        start = time.time()
        for row in data:
            results.append(pool.apply(howmany_within_range, args=(row, 4, 8)))
        end = time.time()
        pool.close()
        pool.join()
        pooltime = time.time()

        print(end-start, pooltime-end)
        i = i + 1

if __name__=='__main__':
    main()
