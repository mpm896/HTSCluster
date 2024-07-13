"""
Helper function to save a Polars DataFrame as an xlsx file with embedded images
@author: Matthew Martinez

Adapted from rdkit.Chem.PandasTools -> SaveXlsxFromFrame
https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/PandasTools.py
"""
# %%
from io import BytesIO
from os import PathLike
import time
from typing import Optional, Tuple

import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw  # type: ignore
import polars as pl


InteractiveRenderer = None
drawOptions = None
if hasattr(rdkit, 'IPythonConsole'):
    try:
      from rdkit.Chem.Draw.IPythonConsole import InteractiveRenderer, drawOptions
    except ImportError:
      pass

def xlsx_from_polarsdf(
    df: pl.DataFrame,
    outFile: str | PathLike,
    molCol: str='Molecule',
    size: Tuple[int,int]=(300,300),
    formats: Optional[str]=None
) -> None:
    """
    Save Polars DataFrame as an xlsx file with embedded images.
    Use the Dataframe write_excel method to first save the file,
    then open the file and embed the images

    :param df: Polars DataFrame
    :param outFile: Name of output file
    :param molCol: Name of the column containing Mol objects
    :param size: Image size
    :param format: Column format
    """
    import xlsxwriter  # type: ignore
    from .utils import get_mols, mols_to_img  # Added here to avoid circular import
    from .parallel import write_parallel_imgs

    assert 'SMILES' in df.columns, "Dataframe must have a column named 'SMILES'"

    print()
    print('Getting images...')
    imgs = mols_to_img(get_mols(df, column_name='SMILES'))
    col_num = len(df.columns)

    with xlsxwriter.Workbook(outFile) as workbook:
        # Create a worksheet
        worksheet = workbook.add_worksheet(name="Sheet1")
        worksheet.set_column(col_num, col_num, width=size[0] / 6.)
        df.write_excel(workbook=workbook, worksheet="Sheet1")
        
        # Insert the images
        for i, img in enumerate(imgs):
            if i > 0:
                worksheet.set_row(i, height=size[1])
            img_data = BytesIO()
            img.save(img_data, format='PNG')
            worksheet.insert_image(i+1, col_num, "f", {'image_data': img_data, 'x_offset': 10, 'y_offset': 10})