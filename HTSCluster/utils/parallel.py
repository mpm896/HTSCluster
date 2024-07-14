# %%
import functools
from io import BytesIO
import multiprocessing as mp
import time

from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map

from .polars_xlsx import xlsx_from_polarsdf


# %%
def img_writer(worksheet, col_num, imgs, size=(300, 300)) -> None:
    """Write images to an xlsx worksheet"""

    for i, img in tqdm(enumerate(imgs)):
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


# %%
def write_parallel_imgs(imgs, worksheet, col_num, size=(300, 300)) -> None:
    """Write an xlsx file with images using multiprocessing"""

    n_workers = 2 * mp.cpu_count()
    batch = round(len(imgs) / n_workers)

    all_args = [worksheet, col_num]

    print()
    print("----- STARTING CONCURRENT WRITING -----")
    writer = process_map(
        functools.partial(img_writer, *all_args),
        imgs,
        max_workers=n_workers,
        chunksize=batch,
    )


if __name__ == "__main__":
    from ..prepare import file_to_df
    from .utils import get_mols

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
