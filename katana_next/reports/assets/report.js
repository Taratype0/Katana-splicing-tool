(() => {
  const $ = id => document.getElementById(id);
  const figures = [...document.querySelectorAll('.exhibit')];
  const chapters = [...document.querySelectorAll('section.chapter')];
  const settings = JSON.parse($('report-settings')?.textContent || '{}');
  let eventSet = settings.event_set || 'stable_significant';
  let expandAll = false;
  const initialComparison = settings.initial_comparison || $('contrastSelect')?.options[0]?.value || 'all';
  if ($('contrastSelect')) $('contrastSelect').value=initialComparison;

  const policySelectors = [...document.querySelectorAll('.das-select')];
  policySelectors.forEach(select => {
    select.value = eventSet;
    select.addEventListener('change', () => {
      eventSet = select.value;
      policySelectors.forEach(s => { s.value = eventSet; });
      applyScope();
    });
  });

  if (location.protocol === 'file:') {
    document.querySelectorAll('a[href="REPORT.zip"]').forEach(a => { a.hidden = true; });
  }

  function syncPills(contrast) {
    document.querySelectorAll('.comparison-pill').forEach(btn => {
      btn.classList.toggle('active', !expandAll && btn.dataset.contrastValue === contrast);
    });
    const expandBtn = $('expandAllScope');
    if (expandBtn) {
      expandBtn.classList.toggle('active-pill', expandAll);
      expandBtn.textContent = expandAll ? 'Focused comparison view' : 'Show all comparisons';
    }
  }

  function applyScope() {
    const contrast = expandAll ? 'all' : ($('contrastSelect')?.value || 'all');
    const celltype = $('celltypeSelect')?.value || 'all';
    let n = 0;
    for (const figure of figures) {
      const c = figure.dataset.contrast || 'all';
      const t = figure.dataset.celltype || 'all';
      const group = JSON.parse(figure.dataset.comparisons || '[]');
      const policy = figure.dataset.policy || 'common';
      const show = (contrast === 'all' || c === contrast || (c === 'all' && (!group.length || group.includes(contrast))))
        && (celltype === 'all' || t === 'all' || t === celltype)
        && (policy === 'common' || policy === eventSet);
      figure.classList.toggle('scope-hidden', !show);
      if (show) n++;
      for (const frame of figure.querySelectorAll('iframe[data-src]')) {
        if (show && !frame.getAttribute('src')) frame.src = frame.dataset.src;
        if(!show&&frame.getAttribute('src'))frame.removeAttribute('src');
      }
    }
    chapters.forEach(ch => {
      ch.classList.toggle(
        'no-visible-figures',
        !!ch.querySelector('.exhibit') && ![...ch.querySelectorAll('.exhibit')].some(f => !f.classList.contains('scope-hidden'))
      );
    });
    if ($('scopeSummary')) {
      $('scopeSummary').textContent = expandAll
        ? `${n} exhibits · all declared comparisons displayed`
        : `${n} exhibits · active comparison & study-wide QC`;
    }
    syncPills($('contrastSelect')?.value || initialComparison);
    window.visibleReportFigures = n;
    window.reportEventSet = eventSet;
    window.reportReady = true;
    updateNavigation();
  }

  for (const id of ['contrastSelect', 'celltypeSelect']) {
    $(id)?.addEventListener('change', () => {
      expandAll = false;
      applyScope();
    });
  }

  document.querySelectorAll('.comparison-pill').forEach(btn => {
    btn.addEventListener('click', () => {
      expandAll = false;
      if ($('contrastSelect')) $('contrastSelect').value = btn.dataset.contrastValue || initialComparison;
      applyScope();
    });
  });

  $('expandAllScope')?.addEventListener('click', () => {
    expandAll = !expandAll;
    applyScope();
  });

  $('resetScope')?.addEventListener('click', () => {
    expandAll = false;
    if ($('contrastSelect')) $('contrastSelect').value = initialComparison;
    if ($('celltypeSelect')) $('celltypeSelect').value = 'all';
    applyScope();
  });

  $('printReport')?.addEventListener('click', () => print());

  document.querySelectorAll('.plot-button').forEach(button => {
    button.addEventListener('click', () => {
      const im = button.querySelector('img');
      if (!im) return;
      $('modalTitle').textContent = im.alt;
      $('modalImage').src = im.src;
      $('modalImage').alt = im.alt;
      $('downloadFigure').href = im.src;
      $('modalScroll').classList.remove('actual');
      $('figureModal').showModal();
    });
  });

  $('fitFigure')?.addEventListener('click', () => $('modalScroll').classList.remove('actual'));
  $('actualFigure')?.addEventListener('click', () => $('modalScroll').classList.add('actual'));
  $('closeFigure')?.addEventListener('click', () => $('figureModal').close());
  $('figureModal')?.addEventListener('click', e => {
    if (e.target === $('figureModal')) $('figureModal').close();
  });

  function updateNavigation() {
    let selected = chapters[0];
    for (const chapter of chapters) {
      if (chapter.getBoundingClientRect().top <= 165) selected = chapter;
      else break;
    }
    document.querySelectorAll('.toc a').forEach(a => {
      a.classList.toggle('active', a.getAttribute('href') === '#' + selected?.id);
    });
  }

  let scheduled = false;
  addEventListener('scroll', () => {
    if (!scheduled) {
      scheduled = true;
      requestAnimationFrame(() => { scheduled = false; updateNavigation(); });
    }
  }, { passive: true });
  addEventListener('resize', updateNavigation);
  addEventListener('load', updateNavigation);
  applyScope();
})();
