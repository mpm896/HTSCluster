from pathlib import Path
from typing import TextIO
import pytest

TEST_DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def screen_file_csv() -> Path:
    """CSV file containing screening data"""
    return TEST_DATA_DIR / "HTS_SG.csv"


@pytest.fixture
def screen_file_csv_handle() -> TextIO:
    """File handle for the CSV file containing screening data"""
    return open(screen_file_csv, "rt")


@pytest.fixture
def screen_file_xlsx() -> Path:
    """Excel file containing screening data"""
    pass


@pytest.fixture
def screen_file_xlsx_handle() -> TextIO:
    "File handle for the Excel file containing screening data"
    return open(screen_file_xlsx)


@pytest.fixture
def query_file_csv() -> Path:
    """Text file for query SMILES"""
    return TEST_DATA_DIR / "query.csv"


@pytest.fixture
def query_file_csv_handle() -> TextIO:
    """Handle for text file with query SMILES"""
    return open(query_file_csv)


@pytest.fixture
def query_file_xlsx() -> Path:
    """Excel file for query SMILES"""
    return TEST_DATA_DIR / "query.xlsx"


@pytest.fixture
def query_file_xlsx_handle() -> TextIO:
    """Handle for Excel file with query SMILES"""
    return open(query_file_xlsx)


@pytest.fixture
def query_file_multicol() -> Path:
    """Excel file for query SMILES with multiple columns"""
    return TEST_DATA_DIR / "query_multicol.xlsx"


@pytest.fixture
def query_file_error() -> Path:
    """Excel file, multiple columns, no SMILES heading, should throw error"""
    return TEST_DATA_DIR / "query_error.xlsx"
