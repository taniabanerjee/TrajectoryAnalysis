def convertToPixelsPerSec(spkmph):
    return 50 * spkmph / (3 * 2.23694)

pps = convertToPixelsPerSec(60)
print (pps * pps)
