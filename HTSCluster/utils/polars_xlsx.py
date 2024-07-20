"""
Helper function to save a Polars DataFrame as an xlsx file with embedded images
@author: Matthew Martinez

Adapted from rdkit.Chem.PandasTools -> SaveXlsxFromFrame
https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/PandasTools.py
"""

# %%
from io import BytesIO
import math
from os import PathLike
from typing import Optional, Tuple

import polars as pl


def xlsx_from_polarsdf(
    df: pl.DataFrame,
    outFile: str | PathLike,
    molCol: str = "Molecule",
    size: Tuple[int, int] = (300, 300),
    formats: Optional[str] = None,
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

    assert "SMILES" in df.columns, "Dataframe must have a column named 'SMILES'"

    print()
    print("Getting images...")
    imgs = mols_to_img(get_mols(df, column_name="SMILES"))
    col_num = len(df.columns)

    with xlsxwriter.Workbook(outFile) as workbook:
        # Create a worksheet
        worksheet = workbook.add_worksheet(name="Sheet1")
        worksheet.set_column(col_num, col_num, width=size[0] / 6.0)
        df.write_excel(workbook=workbook, worksheet="Sheet1")

        # Insert the images
        for i, img in enumerate(imgs):
            if i > 0:
                worksheet.set_row(i, height=size[1])
            img_data = BytesIO()
            img.save(img_data, format="PNG")
            worksheet.insert_image(
                i + 1,
                col_num,
                "f",
                {"image_data": img_data, "x_offset": 10, "y_offset": 10},
            )


def xlsx_from_polarsdf_chunks(
    df: pl.DataFrame,
    outFile: str | PathLike,
    molCol: str = "Molecule",
    size: Tuple[int, int] = (150, 150),
    formats: Optional[str] = None,
    chunksize: int = 1000,
) -> None:
    """
    Save Polars DataFrame as an xlsx file with embedded images.
    Convert SMILES to Mols to Images in chunks
    Use the Dataframe write_excel method to first save the file,
    then open the file and embed the images

    :param df: Polars DataFrame
    :param outFile: Name of output file
    :param molCol: Name of the column containing Mol objects
    :param size: Image size
    :param format: Column format
    :param chunksize: Chunk size
    """
    import xlsxwriter  # type: ignore
    from .utils import get_mols, mols_to_img  # Added here to avoid circular import

    assert "SMILES" in df.columns, "Dataframe must have a column named 'SMILES'"

    chunks = math.ceil(len(df) / chunksize)
    col_num = len(df.columns)
    
    with xlsxwriter.Workbook(outFile) as workbook:
        # Create a worksheet
        worksheet = workbook.add_worksheet(name="Sheet1")
        worksheet.set_column(col_num, col_num, width=size[0] / 6.0)
        df.write_excel(workbook=workbook, worksheet="Sheet1")

        # Process data in chunks to avoid having a large amount of images in memory at once
        for i in range(chunks):
            print(f"Chunk {i + 1} of {chunks} -----")
            print("Getting images...")
            mols = get_mols(df[i * chunksize:(i + 1) * chunksize], column_name="SMILES")
            imgs = mols_to_img(mols)
            startrow = i * chunksize + 1
            # Insert the images
            for i, img in enumerate(imgs):
                worksheet.set_row(i + startrow, height=size[1])
                img_data = BytesIO()
                img.save(img_data, format="JPEG")
                worksheet.insert_image(
                    i + startrow,
                    col_num,
                    "f",
                    {"image_data": img_data, "x_offset": 10, "y_offset": 10},
                )
        
