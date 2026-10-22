from zipfile import ZipFile

from scripts.audit_rbfox2_supplements import read_sheets


def test_strict_ooxml_sparse_cells_and_shared_strings(tmp_path):
    path = tmp_path / 'strict.xlsx'
    ns = 'http://purl.oclc.org/ooxml/spreadsheetml/main'
    with ZipFile(path, 'w') as z:
        z.writestr('xl/sharedStrings.xml', f'<sst xmlns="{ns}"><si><t>IncLevel1</t></si></sst>')
        z.writestr('xl/worksheets/sheet1.xml', f'<worksheet xmlns="{ns}"><sheetData>'
                   '<row r="1"><c r="A1" t="s"><v>0</v></c><c r="C1" t="inlineStr"><is><t>PSI</t></is></c></row>'
                   '<row r="2"><c r="A2"><v>0</v></c><c r="C2"><v>0.5</v></c></row>'
                   '</sheetData></worksheet>')
    assert list(read_sheets(path)) == [('xl/worksheets/sheet1.xml', [['IncLevel1', '', 'PSI'], ['0', '', '0.5']])]
