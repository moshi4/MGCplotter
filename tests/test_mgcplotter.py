import os
import subprocess as sp
from pathlib import Path


def test_circos_installation():
    """Test Circos installation"""
    res = sp.run("circos -modules", shell=True, capture_output=True, text=True)
    lines = [line for line in res.stdout.split("\n") if line != ""]
    for line in lines:
        assert line.startswith("ok")


def test_minimal_option_run(
    reference_file: Path,
    tmp_path: Path,
):
    """Test minimal option run"""
    cmd = f"MGCplotter -r {reference_file} -o {tmp_path}"
    res = sp.run(cmd, shell=True, capture_output=True)

    assert res.returncode == 0
    assert (tmp_path / "circos.png").exists()


def test_normal_option_run(
    reference_file: Path,
    query_gbff_dir: Path,
    tmp_path: Path,
    use_thread_num: int,
):
    """Test normal option run"""
    cmd = (
        f"MGCplotter -r {reference_file} --query_files {query_gbff_dir / '*.gbff'} "
        + f"--assign_cog_color --thread_num {use_thread_num} -o {tmp_path} "
    )
    res = sp.run(cmd, shell=True, capture_output=True)

    assert res.returncode == 0
    assert (tmp_path / "circos.png").exists()


def test_full_option_run(
    reference_file: Path,
    query_faa_dir: Path,
    tmp_path: Path,
    use_thread_num: int,
    cog_color_json_file: Path,
):
    """Test full option run"""
    cmd = (
        f"MGCplotter -r {reference_file} --query_files {query_faa_dir / '*.faa'} "
        + f"--assign_cog_color --thread_num {use_thread_num} -o {tmp_path} "
        + "--mmseqs_evalue 1e-10 --cog_evalue 1e-5 --ticks_labelsize 40 "
        + "--forward_cds_r 0.06 --reverse_cds_r 0.06 --rrna_r 0.06 --trna_r 0.06 "
        + "--conserved_cds_r 0.03 --gc_content_r 0.20 --gc_skew_r 0.20 "
        + "--forward_cds_color blue --reverse_cds_color red --rrna_color magenta "
        + "--trna_color green --conserved_cds_color crimson --gc_content_p_color grey "
        + "--gc_content_n_color black --gc_skew_p_color purple --gc_skew_n_color olive "
        + f"--cog_color_json {cog_color_json_file}"
    )
    res = sp.run(cmd, shell=True, capture_output=True)

    assert res.returncode == 0
    assert (tmp_path / "circos.png").exists()


def test_invalid_suffix_error(reference_file: Path, tmp_path: Path):
    """Test invalid suffix error"""
    cmd = f"MGCplotter -r {reference_file} -o {tmp_path} --query_files test.dummy"
    res = sp.run(cmd, shell=True, capture_output=True)

    assert res.returncode != 0
    assert b"is invalid file suffix" in res.stderr


def test_invalid_color_like_string_error(reference_file: Path, tmp_path: Path):
    """Test invalid color like string error"""
    cmd = f"MGCplotter -r {reference_file} -o {tmp_path} --rrna_color invalidcolor"
    res = sp.run(cmd, shell=True, capture_output=True)

    assert res.returncode != 0
    assert b"invalid color like string" in res.stderr


def test_invalid_radius_value_range_error(reference_file: Path, tmp_path: Path):
    """Test invalid radius value range error"""
    cmd = f"MGCplotter -r {reference_file} -o {tmp_path} --trna_r 0.5"
    res = sp.run(cmd, shell=True, capture_output=True)

    assert res.returncode != 0
    assert b"is invalid value range" in res.stderr


def test_cog_color_file_not_found_error(reference_file: Path, tmp_path: Path):
    """Test cog color file not found error"""
    cmd = f"MGCplotter -r {reference_file} -o {tmp_path} --cog_color_json noexist.json"
    res = sp.run(cmd, shell=True, capture_output=True)

    assert res.returncode != 0
    assert b"File not found" in res.stderr


def test_invalid_cog_color_file_error(
    reference_file: Path,
    tmp_path: Path,
    invalid_cog_color_json_file: Path,
):
    """Test cog color file not found error"""
    cmd = (
        f"MGCplotter -r {reference_file} -o {tmp_path} "
        + f"--cog_color_json {invalid_cog_color_json_file}"
    )
    res = sp.run(cmd, shell=True, capture_output=True)

    assert res.returncode != 0
    assert b"is not COG letter" in res.stderr
    assert b"is not color like string" in res.stderr


def test_generate_cog_color_template(tmp_path: Path):
    """Test generate_cog_color_template"""
    os.chdir(tmp_path)
    res = sp.run("generate_cog_color_template", capture_output=True)

    cog_color_template_json_file = Path(os.getcwd()) / "cog_color_template.json"
    assert cog_color_template_json_file.exists()
    assert res.returncode == 0
