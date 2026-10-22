"""English report output policy; archived source results are never modified."""
from pathlib import Path
import re
import shutil


def is_retired_report(path):
    """Identify language-specific report products, not user sample/gene names."""
    name = Path(path).name.lower()
    if name.endswith(('.html', '.htm', '.pdf', '.md', '.txt')):
        return bool(re.search(r'(?:^|[_-])(?:report[_-]zh|zh[_-]report)(?:\.|[_-])', name)) or bool(
            re.search(r'(?:report|interpretation|summary|conclusion).*[_-]zh\.', name))
    return name.endswith('_bilingual_complete_results.zip')


def copy_english_reader(source, destination):
    """Copy an existing English report without its obsolete language switch."""
    from bs4 import BeautifulSoup

    source, destination = Path(source).resolve(), Path(destination)
    if not (source / 'REPORT_EN.html').is_file():
        raise ValueError('Approved English report input missing: REPORT_EN.html')
    if destination.exists():
        raise FileExistsError('Use a new report directory; existing reports are preserved')
    shutil.copytree(source, destination, ignore=lambda directory, names: [
        name for name in names if (Path(directory) / name).is_file() and is_retired_report(name)])
    for path in destination.rglob('*.html'):
        soup = BeautifulSoup(path.read_text(encoding='utf-8'), 'html.parser')
        changed = False
        for anchor in soup.find_all('a', href=True):
            if is_retired_report(anchor['href'].split('#', 1)[0].split('?', 1)[0]):
                anchor.decompose()
                changed = True
        if changed:
            path.write_text(str(soup), encoding='utf-8')
    return destination / 'REPORT_EN.html'
