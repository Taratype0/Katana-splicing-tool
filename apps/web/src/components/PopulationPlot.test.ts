import {expect,test} from 'vitest';
import {vegaLite,vega} from 'vega-embed';
import {volcanoSpec,type Population} from './PopulationPlot';

test('application volcano compiles with symmetric zero-centred axes, threshold rules and original rows',()=>{
 const population:Population={mode:'points',total:3,plotted:3,unplottable:0,x_min:-2,x_max:7,y_max:8,points:[[-2,3,'r1','Fos'],[7,8,'r2','Ranbp2'],[.1,.2,'r3','Nonsignificant']]};
 const spec=volcanoSpec(population,'log2fc','deg_fdr',.05,1,true);
 expect(spec.layer[0].data.values).toHaveLength(3);
 expect(spec.layer[0].encoding!.x!.scale.domain).toEqual([-7,7]);
 expect(spec.layer[0].encoding!.x!.axis.values).toContain(0);
 expect(spec.layer[0].encoding!.y!.scale.domain[0]).toBe(0);
 expect(spec.layer[2].data.values).toEqual([{x:-1},{x:1}]);
 expect(spec.layer[3].data.values).toEqual([{y:-Math.log10(.05)}]);
 expect(()=>vega.parse(vegaLite.compile(spec as never).spec)).not.toThrow();
});

test('PS-DAS keeps unsupported rows in the background and preserves exact zero FDR',()=>{
 const population:Population={mode:'points',total:2,plotted:2,unplottable:0,x_min:-.5,x_max:.5,y_max:300,points:[[.2,300,'r1','Fos',true,0],[-.2,3,'r2','Other',false,.001]]};
 const spec=volcanoSpec(population,'dpsi','FDR',.05,.1,true,true,true);
 expect(spec.layer[0].data.values).toEqual(expect.arrayContaining([expect.objectContaining({row_id:'r1',p:0,status:'Higher in numerator'}),expect.objectContaining({row_id:'r2',status:'Other measured rows'})]));
 expect(()=>vega.parse(vegaLite.compile(spec as never).spec)).not.toThrow();
});
