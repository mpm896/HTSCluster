# Benchmark times for different functions
import time
from typing import Optional, Tuple
from io import BytesIO

import numpy as np
import polars as pl
from rdkit import Chem  # type: ignore
from rdkit.Chem import Draw, PandasTools  # type: ignore

from ..prepare import file_to_df
from ..utils.polars_xlsx import xlsx_from_polarsdf
from ..utils.utils import get_mols, insert_mols


# Writing an xlsx file with and without images, using Polars, PandasTools, and xlsx_from_polarsdf
def xlsx_from_polarsdf_deprecated(
    df: pl.DataFrame,
    outFile: str,
    molCol: str = "Molecule",
    size: Tuple[int, int] = (300, 300),
    formats: Optional[str] = None,
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
    import xlsxwriter  # type: ignore

    drawOptions = None

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
    for key in ["write_string", "write_number", "write_datetime"]:
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

        worksheet.set_row(
            row_idx_actual, height=size[1]
        )  # looks like height is not in px?

        for col_idx, col in enumerate(cols):
            if col_idx in molCol_indices:
                image_data = BytesIO()
                m = row[col_idx]
                img = Draw.MolToImage(
                    m if isinstance(m, Chem.Mol) else Chem.Mol(),
                    size=size,
                    options=drawOptions,
                )
                img.save(image_data, format="PNG")
                worksheet.insert_image(
                    row_idx_actual, col_idx, "f", {"image_data": image_data}
                )
                worksheet.set_column(
                    col_idx, col_idx, width=size[0] / 6.0
                )  # looks like height is not in px?
            elif str(dataTypes[col]) == "String":
                # string length is limited in xlsx
                worksheet.write_string(
                    row_idx_actual,
                    col_idx,
                    str(row[col_idx])[:32000],
                    cell_formats["write_string"],
                )
            elif ("float" in str(dataTypes[col]).lower()) or (
                "int" in str(dataTypes[col]).lower()
            ):
                if (row[col_idx] != np.nan) or (row[col_idx] != np.inf):
                    worksheet.write_number(
                        row_idx_actual,
                        col_idx,
                        row[col_idx],
                        cell_formats["write_number"],
                    )
            elif "datetime" in str(dataTypes[col]):
                worksheet.write_datetime(
                    row_idx_actual,
                    col_idx,
                    row[col_idx],
                    cell_formats["write_datetime"],
                )

    workbook.close()
    image_data.close()


def benchmark_write_xlsx() -> None:
    """Benchmark writing xlsx files with and without images"""
    # Read in 2 CSVs: 1 for small dataset, 1 for large dataset
    OUT_DIR = "tests/data/output/"
    small_df = file_to_df("tests/data/HITS.csv")
    large_df = file_to_df("tests/data/Chembrigde_Div.csv")
    small_df_pandas = small_df.to_pandas()
    large_df_pandas = large_df.to_pandas()

    # Get Mols to use for images
    small_mols = get_mols(small_df)
    large_mols = get_mols(large_df)

    # Inset Mols into Pandas DataFrames
    small_df_pandas_with_mols = small_df_pandas
    small_df_pandas_with_mols["Molecule"] = small_mols
    large_df_pandas_with_mols = large_df_pandas
    large_df_pandas_with_mols["Molecule"] = large_mols

    # Inset Mols into Polars DataFrames
    small_df_with_mols = insert_mols(small_mols, small_df)
    large_df_with_mols = insert_mols(large_mols, large_df)

    # Write XLSX files and benchmark times
    # First: CSVs without images
    print("--------------------")
    print("Writing small xlsx without images")
    print("--------------------")
    print()

    start = time.time()
    small_df_pandas.to_excel(f"{OUT_DIR}/small-pandas.xlsx", index=False)
    pd_write_excel_small = time.time() - start

    start = time.time()
    small_df.write_excel(
        f"{OUT_DIR}/small-pl.xlsx",
    )
    pl_write_excel_small = time.time() - start

    print("--------------------")
    print("Writing large xlsx without images")
    print("--------------------")
    print()

    start = time.time()
    large_df_pandas.to_excel(f"{OUT_DIR}/small-pandas.xlsx", index=False)
    pd_write_excel_large = time.time() - start

    start = time.time()
    large_df.write_excel(f"{OUT_DIR}/small-pl.xlsx")
    pl_write_excel_large = time.time() - start

    print()
    print("-------- TIMES -> NO IMAGES --------")
    print()

    print("SMALL DATAFRAME: 315 COMPOUNDS")
    print(f"Pandas: {pd_write_excel_small}")
    print(f"Polars: {pl_write_excel_small}")
    print()

    print("LARGE DATAFRAME: 43,000 COMPOUNDS")
    print(f"Pandas: {pd_write_excel_large}")
    print(f"Polars: {pl_write_excel_large}")
    print()
    print()

    print("--------------------")
    print("Writing small xlsx with images")
    print("--------------------")
    print()

    start = time.time()
    PandasTools.SaveXlsxFromFrame(
        small_df_pandas_with_mols,
        f"{OUT_DIR}/small-pandas-with-images.xlsx",
        molCol="Molecule",
        size=(300, 300),
    )
    pd_write_excel_with_images_small = time.time() - start

    start = time.time()
    xlsx_from_polarsdf_deprecated(
        small_df_with_mols,
        f"{OUT_DIR}/small-pl-with-images.xlsx",
        molCol="Molecule",
        size=(300, 300),
    )
    pl_likepandastools_write_excel_with_images_small = time.time() - start

    start = time.time()
    xlsx_from_polarsdf(
        small_df,
        f"{OUT_DIR}/small-pl-with-images.xlsx",
        molCol="Molecule",
        size=(300, 300),
    )
    pl_write_excel_with_images_small = time.time() - start

    print()
    print("-------- TIMES -> WITH IMAGES --------")
    print()

    print("SMALL DATAFRAME: 315 COMPOUNDS")
    print(f"PandasTools: {pd_write_excel_with_images_small}")
    print(
        f"Polars (like Pandas Tools Algo): {pl_likepandastools_write_excel_with_images_small}"
    )
    print(f"Polars: {pl_write_excel_with_images_small}")
    print()

    print("--------------------")
    print("Writing large xlsx with images")
    print("--------------------")
    print()

    start = time.time()
    PandasTools.SaveXlsxFromFrame(
        large_df_pandas_with_mols,
        f"{OUT_DIR}/small-pandas-with-images.xlsx",
        molCol="Molecule",
        size=(300, 300),
    )
    pd_write_excel_with_images_large = time.time() - start

    start = time.time()
    xlsx_from_polarsdf_deprecated(
        large_df_with_mols,
        f"{OUT_DIR}/small-pl-with-images.xlsx",
        molCol="Molecule",
        size=(300, 300),
    )
    pl_likepandastools_write_excel_with_images_large = time.time() - start

    start = time.time()
    xlsx_from_polarsdf(
        large_df,
        f"{OUT_DIR}/small-pl-with-images.xlsx",
        molCol="Molecule",
        size=(300, 300),
    )
    pl_write_excel_with_images_large = time.time() - start

    print("LARGE DATAFRAME: 43,000 COMPOUNDS")
    print(f"PandasTools: {pd_write_excel_with_images_large}")
    print(
        f"Polars (like Pandas Tools Algo): {pl_likepandastools_write_excel_with_images_large}"
    )
    print(f"Polars: {pl_write_excel_with_images_large}")
    print()


def benchmark_parallel() -> None:
    from ..prepare import file_to_df
    from ..utils.utils import get_mols

    df_large = file_to_df("tests/data/Chembrigde_Div.csv")
    large_mols = get_mols(df_large)

    small_df = file_to_df("tests/data/HITS.csv")

    print()
    print("STARTING TO WRITE XLSX FILE")
    start = time.time()
    xlsx_from_polarsdf(
        small_df,
        "tests/data/output/test-parallel.xlsx",
        molCol="Molecule",
        size=(300, 300),
    )
    print(f"Took {time.time() - start} seconds to write concurrently")


if __name__ == "__main__":
    # benchmark_write_xlsx()
    benchmark_parallel()
