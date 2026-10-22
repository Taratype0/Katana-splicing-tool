from scripts.verify_annotated_external_product import write_changed


def test_generated_input_rerun_preserves_fingerprint(tmp_path):
    path = tmp_path/'samples.tsv'
    write_changed(path, 'sample\tcondition\r\ns1\tWT\r\n')
    before = path.stat().st_mtime_ns
    write_changed(path, 'sample\tcondition\r\ns1\tWT\r\n')
    assert path.stat().st_mtime_ns == before
    assert path.read_bytes() == b'sample\tcondition\ns1\tWT\n'
    path.write_bytes(b'sample\tcondition\r\r\ns1\tWT\r\r\n')
    before = path.stat().st_mtime_ns
    write_changed(path, 'sample\tcondition\r\ns1\tWT\r\n')
    assert path.stat().st_mtime_ns == before
