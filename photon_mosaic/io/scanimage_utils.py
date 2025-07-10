# Utility functions for working with ScanImage tiff files
# Adapted from https://github.com/hjmh/fly2p/blob/master/fly2p/preproc/scanImageUtils.py

# TODO discard flyback: https://github.com/T-E-Lab/GULP/blob/6847021222ea2ace0dff11fb847d205c568283a3/gulp2p/preproc/tiff.py#L19
import json
import numpy as np
from skimage import io
import tifffile

def extractSIbasicMetadata(metadat):

    #initilize dict
    metadict = {}

    if type(metadat) == dict:
        nCh = metadat['SI.hChannels.channelSave']
        fpsscan = metadat['SI.hRoiManager.scanFrameRate']
        discardFBFrames = metadat['SI.hFastZ.discardFlybackFrames']
        nDiscardFBFrames = metadat['SI.hFastZ.numDiscardFlybackFrames']
        fpv = metadat['SI.hFastZ.numFramesPerVolume']
        nVols = metadat['SI.hFastZ.numVolumes']
        stackZStepSize = metadat['SI.hStackManager.stackZStepSize']
        scanVolumeRate = metadat['SI.hRoiManager.scanVolumeRate']
        [p00, p10, p01, p11] = metadat['SI.hRoiManager.imagingFovUm']

    else:
        for i, line in enumerate(metadat.split('\n')):

            if not 'SI.' in line: continue
            # extract version
            if 'VERSION_' in line: print(line)

            # get channel info
            if 'channelSave' in line:
                #print(line)
                if not '[' in line:
                    nCh = 1
                else:
                    strchanlist = line.split('=')[-1].strip()
                    try: chanlist = [int(i) for i in strchanlist.strip('][').split(' ')]    
                    except ValueError:  chanlist = [int(i) for i in strchanlist.strip('][').split(';')] 
                    nCh = len(chanlist)

            if 'scanFrameRate' in line:
                fpsscan = float(line.split('=')[-1].strip())


            #if 'hFastZ' in line:
            if 'discardFlybackFrames' in line:
                discardFBFrames = line.split('=')[-1].strip()

            if 'numDiscardFlybackFrames' in line:
                nDiscardFBFrames = int(line.split('=')[-1].strip())

            if 'numFramesPerVolume' in line:
                fpv = int(line.split('=')[-1].strip())


            if 'numVolumes' in line:
                nVols = int(line.split('=')[-1].strip())

            if 'hStackManager.stackZStepSize' in line:
                stackZStepSize = float(line.split('=')[-1].strip())

            if 'hRoiManager.scanVolumeRate' in line:
                scanVolumeRate = float(line.split('=')[-1].strip())

            if 'SI.hRoiManager.imagingFovUm' in line:
                imagingFovUm = line.split('=')[-1].strip()
                p00 = np.fromstring(imagingFovUm[1:-1].split(';')[0], dtype=float, count=2, sep=' ')
                p10 = np.fromstring(imagingFovUm[1:-1].split(';')[1], dtype=float, count=2, sep=' ')
                p01 = np.fromstring(imagingFovUm[1:-1].split(';')[2], dtype=float, count=2, sep=' ')
                p11 = np.fromstring(imagingFovUm[1:-1].split(';')[3], dtype=float, count=2, sep=' ')

    metadict["nCh"] = nCh
    metadict["fpsscan"] = fpsscan
    metadict["discardFBFrames"] = discardFBFrames
    metadict["nDiscardFBFrames"] = nDiscardFBFrames
    metadict["fpv"] = fpv
    metadict["nVols"] = nVols
    metadict["stackZStepSize"] = stackZStepSize
    metadict["scanVolumeRate"] = scanVolumeRate
    metadict["fovCoords"] = {'p00':list(p00),'p10':list(p01),
            'p01':list(p10),'p11':list(p11)}
    metadict["xrange_um"] = p01[0]-p00[0]
    metadict["yrange_um"] = p11[1]-p00[1]

    return metadict

def get_si_metadata(tiff_path):
    """
    Reads metadata from a Scanimage TIFF and return a dictionary with
    specified key values.

    Currently can only extract numerical data.

    Args:
        tiff_path: path to TIFF or directory containing tiffs

    Returns:
        dict: dictionary of SI parameters

    """
    if tiff_path.suffix != ".tif":
        tiffs = [tiff_path / tiff for tiff in sorted(tiff_path.glob("*.tif"))]
    else:
        tiffs = [
            tiff_path,
        ]
    if tiffs:
        return extractSIbasicMetadata(tifffile.TiffFile(tiffs[0]).scanimage_metadata["FrameData"])
    else:
        return None
    
def reshape_si_stack(data:np.array, metadata:dict):
    size_y = size_x = 512
    size_t = (np.size(stack) // (metadata["fpv"] * metadata["nCh"] * size_y * size_x))
    stack = data.reshape(size_t, metadata["fpv"], metadata["nCh"], size_y, size_x)
    metadata['dimension_order'] = 'TZCYX'
    return stack, metadata