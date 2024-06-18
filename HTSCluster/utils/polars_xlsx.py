"""
Helper function to save a Polars DataFrame as an xlsx file with embedded images
@author: Matthew Martinez

Adapted from rdkit.Chem.PandasTools -> SaveXlsxFromFrame
https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/PandasTools.py
"""
# %%
from io import BytesIO
from os import PathLike
from typing import Optional, Tuple

import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
import polars as pl

InteractiveRenderer = None
drawOptions = None
if hasattr(rdkit, 'IPythonConsole'):
    try:
      from rdkit.Chem.Draw.IPythonConsole import InteractiveRenderer, drawOptions
    except ImportError:
      pass

# %%
def xlsx_from_polarsdf(
    df: pl.DataFrame,
    outFile: str | PathLike,
    molCol: str='ROMol',
    size: Tuple[int,int]=(300,300),
    formats: Optional[str]=None
) -> None:
    """
    Saves Polars DataFrame as xlsx from with embedded images.
    For now, only support a single column of Mols
    
    :param df: Polars DataFrame
    :param outFile: Name of output file
    :param molCol: Name of the column containing Mol objects
    :param size: Image size
    :param format: Column format
    """
    import xlsxwriter

    # Setup params from the DataFrame
    cols = list(df.columns)
    df_types = list(df.dtypes)
    if isinstance(molCol, str):
        molCol = [molCol]
    molCol = list(set(molCol))
    dataTypes = dict(zip(cols, df_types))
    molCol_indices = [cols.index(mc) for mc in molCol]

    # New Excel workbook
    workbook = xlsxwriter.Workbook(outFile)
    cell_formats = {}
    formats = formats or {}
    for key in ['write_string', 'write_number', 'write_datetime']:
        format = formats.get(key, None)
        if format is not None:
            format = workbook.add_format(format)
        cell_formats[key] = format
    worksheet = workbook.add_worksheet()  # New worksheet

    # Write first row with column names
    for col_idx, col in enumerate(cols):
        worksheet.write_string(0, col_idx, col)

    for row_idx, row in enumerate(df.iter_rows()):
        row_idx_actual = row_idx + 1

        worksheet.set_row(row_idx_actual, height=size[1])  # looks like height is not in px?

        for col_idx, col in enumerate(cols):
            if col_idx in molCol_indices:
                image_data = BytesIO()
                m = row[col_idx]
                img = Draw.MolToImage(
                    m if isinstance(m, Chem.Mol) else Chem.Mol(), size=size,
                                    options=drawOptions
                )
                img.save(image_data, format='PNG')
                worksheet.insert_image(
                    row_idx_actual,
                    col_idx,
                    "f",
                    {'image_data': image_data}
                )
                worksheet.set_column(
                    col_idx,
                    col_idx,
                    width=size[0] / 6.
                )  # looks like height is not in px?
            elif str(dataTypes[col]) == "String":
              # string length is limited in xlsx
                worksheet.write_string(
                    row_idx_actual,
                    col_idx,
                    str(row[col_idx])[:32000], cell_formats['write_string']
                )
            elif (('float' in str(dataTypes[col]).lower()) 
                    or ('int' in str(dataTypes[col]).lower())):
                if (row[col_idx] != np.nan) or (row[col_idx] != np.inf):
                  worksheet.write_number(
                      row_idx_actual,
                      col_idx,
                      row[col_idx],
                      cell_formats['write_number']
                  )
            elif 'datetime' in str(dataTypes[col]):
                worksheet.write_datetime(
                    row_idx_actual, 
                    col_idx,
                    row[col_idx],
                    cell_formats['write_datetime']
                )

    workbook.close()
    image_data.close()