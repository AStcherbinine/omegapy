
## Saving OMEGA data

## Loading OMEGA data

## Finding wavelength indexes

## Generating masks

???+ abstract "Note from the *SOFT10_README.txt*"
    The last 4 scans (16 pixels), 2 scans (32 pixels) or 1 scan (64, 128
    pixels) of `idat` and `jdat` have only IR data (spectels 0 to 255).

    There is calibration data at the beginning of each cube for the visible
    channel, at the beginning of the *ORBNNNN_0.QUB* cubes for the IR channels
    The number of calibration scans depend on the pixel length (16 to 128)
    and the summation (1, 2 or 4 for 128 pixel modes):

    * for every cube, the first scan (128 pixels x 4), 3 scans (128 x 2) 
    7 scans (128 x 1), 14 scans (64 pixels), 
    28 scans (32 pixels), 56 scans (16 pixels) of the visible channel
    (spectels 256 and above) correspond to an internal calibration

    * for cubes with names *NNNN_0* (first cube in a sequence), the first
    6 scans (128 x 4), 12 scans (128 x 2), 24 scans (128 x 1)
    48 scans (64 pixels), 96 scans (32 pixels) or 192 scans (16 pixels)
    of the IR channel (spectels 0 to 255) correspond to an internal 
    calibration (closed shutter, lamp on at 6 different levels, 
    in order 0,4,3,2,1,0

### Handling corrupted 128-px wide cubes

## Observation search

 * JMARS
 * `find_cube`
