#!/bin/usr/env python
### this function is used to update the ptuse header in the search file

from astropy.io import fits

def update_ptuse_header(search_file, ra, dec):
    """
    Update the PTUSE header in the search file with new coordinates.

    Parameters:
    search_file (str): Path to the search FITS file.
    ra (str): New right ascension to update in the header.
    dec (str): New declination to update in the header.

    Returns:
    None
    """
    # Open the FITS file in update mode
    assert ":" in ra and ":" in dec, "RA and DEC must be in sexagesimal format (e.g., '12:34:56.78')"
    with fits.open(search_file, mode='update') as hdul:
        # Update the PTUSE header with new coordinates
        hdul[0].header["RA"] = ra
        hdul[0].header["DEC"] = dec
        hdul[0].header["STT_CRD1"] = ra
        hdul[0].header["STT_CRD2"] = dec
        hdul[0].header["STP_CRD1"] = ra
        hdul[0].header["STP_CRD2"] = dec
        
        # Save changes to the FITS file
        # hdul.flush()

if __name__ == "__main__":
    ra = "16:46:22.98"
    dec = "-44:05:40.8"

    import glob
    searchfiles = glob.glob("./*.sf")
    for file in searchfiles:
        print(f"Updating PTUSE header in {file} with RA={ra}, DEC={dec}")
        update_ptuse_header(file, ra, dec)