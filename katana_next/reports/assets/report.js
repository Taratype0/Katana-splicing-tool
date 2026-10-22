(()=>{
  const $=id=>document.getElementById(id);
  const figures=[...document.querySelectorAll('.exhibit')];
  const chapters=[...document.querySelectorAll('section.chapter')];
  let eventSet=JSON.parse($('report-settings')?.textContent||'{}').event_set||'stable_significant';
  const policySelectors=[...document.querySelectorAll('.das-select')];
  policySelectors.forEach(select=>{select.value=eventSet;select.addEventListener('change',()=>{eventSet=select.value;policySelectors.forEach(s=>s.value=eventSet);applyScope()})});
  if(location.protocol==='file:')document.querySelectorAll('a[href="REPORT.zip"]').forEach(a=>a.hidden=true);
  function applyScope(){
    const contrast=$('contrastSelect')?.value||'all',celltype=$('celltypeSelect')?.value||'all';
    let n=0;
    for(const figure of figures){
      const c=figure.dataset.contrast||'all',t=figure.dataset.celltype||'all';
      const policy=figure.dataset.policy||'common';
      const show=(contrast==='all'||c==='all'||c===contrast)&&(celltype==='all'||t==='all'||t===celltype)&&(policy==='common'||policy===eventSet);
      figure.classList.toggle('scope-hidden',!show);if(show)n++;
    }
    chapters.forEach(ch=>ch.classList.toggle('no-visible-figures',!!ch.querySelector('.exhibit')&&![...ch.querySelectorAll('.exhibit')].some(f=>!f.classList.contains('scope-hidden'))));
    if($('scopeSummary'))$('scopeSummary').textContent=n+' exhibits in this view';
    window.visibleReportFigures=n;window.reportEventSet=eventSet;window.reportReady=true;updateNavigation();
  }
  for(const id of ['contrastSelect','celltypeSelect'])$(id)?.addEventListener('change',applyScope);
  $('resetScope')?.addEventListener('click',()=>{
    for(const id of ['contrastSelect','celltypeSelect'])if($(id))$(id).value='all';
    applyScope();
  });
  $('printReport')?.addEventListener('click',()=>print());
  document.querySelectorAll('.plot-button').forEach(button=>button.addEventListener('click',()=>{
    const im=button.querySelector('img');
    $('modalTitle').textContent=im.alt;$('modalImage').src=im.src;$('modalImage').alt=im.alt;
    $('downloadFigure').href=im.src;$('modalScroll').classList.remove('actual');$('figureModal').showModal();
  }));
  $('fitFigure')?.addEventListener('click',()=>$('modalScroll').classList.remove('actual'));
  $('actualFigure')?.addEventListener('click',()=>$('modalScroll').classList.add('actual'));
  $('closeFigure')?.addEventListener('click',()=>$('figureModal').close());
  $('figureModal')?.addEventListener('click',e=>{if(e.target===$('figureModal'))$('figureModal').close()});
  function updateNavigation(){
    // Whole chapters span many screens; use the latest heading passed.
    let selected=chapters[0];
    for(const chapter of chapters){if(chapter.getBoundingClientRect().top<=160)selected=chapter;else break;}
    document.querySelectorAll('.toc a').forEach(a=>a.classList.toggle('active',a.getAttribute('href')==='#'+selected?.id));
  }
  let scheduled=false;
  addEventListener('scroll',()=>{if(!scheduled){scheduled=true;requestAnimationFrame(()=>{scheduled=false;updateNavigation()})}},{passive:true});
  addEventListener('resize',updateNavigation);addEventListener('load',updateNavigation);
  applyScope();
})();
