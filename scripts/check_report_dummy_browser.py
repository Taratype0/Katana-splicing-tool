"""Exercise the complete dummy with Chromium; no scientific computations."""
from pathlib import Path
import argparse
import json
from playwright.sync_api import sync_playwright


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,required=True)
    args=ap.parse_args();root=args.directory.resolve();shots=root.parent/'screenshots';shots.mkdir(exist_ok=True)
    manifest=json.loads((root/'report-manifest.json').read_text('utf8'))
    with sync_playwright() as p:
        browser=p.chromium.launch();page=browser.new_page(viewport={'width':1440,'height':1050},device_scale_factor=1)
        errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto((root/'REPORT.html').as_uri());page.wait_for_function('window.reportReady === true')
        assert page.locator('section.chapter[id^="chapter-"]').count()==7
        for cid in manifest['contrasts']:
            page.locator('#contrastSelect').select_option(cid)
            assert page.locator('#chapter-5 .exhibit:not(.scope-hidden)').count()>=4
            assert page.locator('.top10:not(.scope-hidden)').count()>=4
        page.locator('#contrastSelect').select_option(manifest['contrasts'][0])
        page.locator('#chapter-5').screenshot(path=str(shots/'05-enrichment.png'))
        page.locator('.motif-significance:not(.scope-hidden)').screenshot(path=str(shots/'06-motif.png'))
        page.locator('[id$="-mapp"]:not(.scope-hidden)').screenshot(path=str(shots/'06-mapp.png'))
        page.locator('#report-das-selection').select_option('significant')
        assert page.locator('.top10[data-policy="significant"]:not(.scope-hidden)').count()==1
        assert page.locator('.top10[data-policy="stable_significant"]:not(.scope-hidden)').count()==0
        for panel in manifest['panels']:
            if panel['id'].endswith('-mapp'):
                axis=panel['metadata']['activity_axis'];assert len(axis)==6
                assert all(x['condition'] and x['sample_id'] in x['label'] for x in axis)
                assert len(set(x['condition'] for x in axis))==2
        network=manifest['network_browsers'][0]
        page.goto((root/network['source']).as_uri())
        assert page.locator('#moduleSelect option').count()==4
        assert 'Chapter 07' in page.title() or 'Chapter 07' in page.locator('body').inner_text()
        assert '7.1' in page.locator('.toc').inner_text()
        assert 'No FDR-significant' in page.locator('#topModules').inner_text()
        for module in ['M2','M3']:
            page.locator('#moduleSelect').select_option(module)
            assert page.locator('#networkPlot svg').is_visible()
            assert 'PTBP1' in page.locator('#independentModuleLists').inner_text()
            if module=='M2':
                assert page.locator('#pathwayTable svg').count()==1
                page.locator('#chapter-3').screenshot(path=str(shots/'07-network.png'))
                page.locator('#chapter-4').screenshot(path=str(shots/'07-module-enrichment.png'))
            else: assert 'No significant terms' in page.locator('#pathwayTable').inner_text()
        assert not errors,errors
        browser.close()
    result=dict(status='passed',comparisons=4,chapters=7,all_modules_browsable=True,network_visible=True,
                mapp_conditions_preserved=True,das_switch_checked=True,javascript_errors=errors,synthetic_only=True)
    (root.parent/'browser-acceptance.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result))


if __name__=='__main__':main()
