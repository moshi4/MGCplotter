import os
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def testdata_dir() -> Path:
    """Testdata directory fixture"""
    return Path(__file__).parent / "testdata"


@pytest.fixture(scope="session")
def small_dataset_dir(testdata_dir: Path) -> Path:
    """Small dataset directory fixture"""
    return testdata_dir / "small_dataset"


@pytest.fixture(scope="session")
def reference_file(small_dataset_dir: Path) -> Path:
    """Query file fixture"""
    return small_dataset_dir / "reference" / "GCF_000286675.1_ASM28667v1_genomic.gbff"


@pytest.fixture(scope="session")
def query_gbff_dir(small_dataset_dir: Path) -> Path:
    """Query gbff directory fixture"""
    return small_dataset_dir / "query_gbff"


@pytest.fixture(scope="session")
def query_faa_dir(small_dataset_dir: Path) -> Path:
    """Query faa directory fixture"""
    return small_dataset_dir / "query_faa"


@pytest.fixture(scope="session")
def cog_color_json_file(testdata_dir: Path) -> Path:
    """COG color json file fixture"""
    return testdata_dir / "cog_color.json"


@pytest.fixture(scope="session")
def invalid_cog_color_json_file(testdata_dir: Path) -> Path:
    """Invalid COG color json file fixture"""
    return testdata_dir / "invalid_cog_color.json"


@pytest.fixture(scope="session")
def use_thread_num() -> int:
    """Use thread number fixture"""
    cpu_num = os.cpu_count()
    return 1 if cpu_num is None or cpu_num == 1 else cpu_num - 1
