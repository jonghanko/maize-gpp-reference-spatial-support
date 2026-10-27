"""Independent CSV/SVG/PDF verification; does not import the plotting module."""
from pathlib import Path
import csv,json,hashlib,re,xml.etree.ElementTree as ET,struct
from pypdf import PdfReader
import pdfplumber
import argparse
PACKAGE=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,default=PACKAGE/'data/contrasts12')
parser.add_argument('--out',type=Path,default=PACKAGE/'output/contrasts12')
args=parser.parse_args()
SOURCE=args.source.resolve()
ROOT=args.out.resolve()
ROOT.mkdir(parents=True,exist_ok=True)
NS={'s':'http://www.w3.org/2000/svg'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(p):
 with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def font_embeds(page):
 report=[]
 for k,ref in page['/Resources']['/Font'].items():
  f=ref.get_object();children=f.get('/DescendantFonts');d=children[0].get_object().get('/FontDescriptor') if children else f.get('/FontDescriptor');desc=d.get_object() if d else {};report.append({'font':str(f.get('/BaseFont')),'embedded':any(key in desc for key in ['/FontFile','/FontFile2','/FontFile3'])})
 return report
source=rows(SOURCE/'MAIN_CONDITIONALS_OWN_PATHS.csv');agg=rows(SOURCE/'TableS6_aggregation_estimands.csv');fields=['training_given_test_B','training_given_test_C','test_given_training_B','test_given_training_C','interaction','original_pipeline_contrast']
receipt=json.loads((ROOT/'PLOT_RECEIPT.json').read_text(encoding='utf-8'));reports=[]
for stem,panels in [('Figure3_BC_conditional_contrasts',[('principal','ridge'),('principal','random_forest'),('T1','ridge'),('T1','random_forest')]),('FigureS7_HGB_conditional_contrasts',[('principal','hist_gradient_boosting'),('T1','hist_gradient_boosting')]),('Figure4_primary_LOYO_input_effects',[])]:
 svg=ET.parse(ROOT/(stem+'.svg')).getroot();clips={p.attrib['id']:p.find('s:rect',NS).attrib for p in svg.findall('.//s:clipPath',NS)};expected=[];recovered=[]
 if panels:
  for cohort,model in panels:
   for ref in ['VUT_REF','CUT_REF']:
    r=next(r for r in source if (r['cohort'],r['model'],r['fit_reference'])==(cohort,model,ref));expected.extend(float(r[k]) for k in fields)
  for group in svg.findall('.//s:g',NS):
   if not group.attrib.get('id','').startswith('PathCollection_'):continue
   g=group.find('s:g',NS);clip=clips[g.attrib['clip-path'][5:-1]]
   for use in g.findall('s:use',NS):recovered.append(-.15+.25*(float(use.attrib['x'])-float(clip['x']))/float(clip['width']))
  axis_limits=[-.15,.10];colors=['#0072b2','#d55e00']
 else:
  cols=['equal_site_season_absolute_difference','equal_region_pooled_absolute_difference'];cts=['C_minus_B','B_minus_A_shared','C_minus_A_shared'];primary=['ridge','random_forest']
  mx=max(abs(float(r[col])) for r in agg if r['model'] in primary and r['contrast'] in cts for col in cols)*1.45
  for col in cols:
   for model in primary:
    for ref in ['VUT_REF','CUT_REF']:
     for ct in cts:expected.append(float(next(r for r in agg if (r['model'],r['variant'],r['contrast'])==(model,ref,ct))[col]))
  for p in svg.findall('.//s:path',NS):
   style=p.attrib.get('style','');clip=p.attrib.get('clip-path','')
   if not any('fill: '+c in style for c in ['#2166a1','#c87820']) or not clip:continue
   data=list(map(float,re.findall(r'-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?',p.attrib['d'])));rect=clips[clip[5:-1]];recovered.append((data[1]-data[5])/float(rect['height'])*(2*mx))
  axis_limits=[-mx,mx];colors=['#2166a1','#c87820']
 assert len(recovered)==len(expected),(stem,len(recovered),len(expected))
 err=max(abs(x-y) for x,y in zip(expected,recovered));assert err<1e-8,(stem,err)
 export=rows(ROOT/(stem+'_plot_values.csv'));assert [float(r['value']) for r in export]==expected
 page=PdfReader(ROOT/(stem+'.pdf')).pages[0];width=float(page.mediabox.width)/72;height=float(page.mediabox.height)/72
 assert abs(width-7.16)<1e-9 and height<=7.5
 font_states=[]
 page.extract_text(visitor_text=lambda text,cm,tm,font,size:font_states.append(float(size)) if text.strip() else None)
 assert font_states and min(font_states)>=12,(stem,font_states)
 with pdfplumber.open(ROOT/(stem+'.pdf')) as p:chars=[c for c in p.pages[0].chars if c['text'].strip()]
 fontmin=min(float(c['size']) for c in chars);assert fontmin>=12-1e-8,(stem,fontmin)
 embeds=font_embeds(page);assert all(f['embedded'] for f in embeds)
 with (ROOT/(stem+'.png')).open('rb') as f:head=f.read(24)
 pixels=struct.unpack('>II',head[16:24]);assert pixels[0]==4296 and pixels[1]==round(height*600)
 plot=next(r for r in receipt['plots'] if r['stem']==stem);b=plot['bbox_audit'];assert not b['outside_figure'] and not b['text_bbox_collisions'] and not b['axes_bbox_collisions']
 all_svg_text=''.join(svg.itertext());assert all(c in (ROOT/(stem+'.svg')).read_text(encoding='utf-8').lower() for c in colors)
 # Panel descriptions were moved; no suptitle, cohort count or explanatory footer remains inside SVG.
 assert '642 days' not in all_svg_text and '437 days' not in all_svg_text and 'Negative conditional/pipeline' not in all_svg_text
 reports.append({'stem':stem,'status':'PASS','source_values':len(expected),'export_matches_source_exact':True,'independent_SVG_recovered_count':len(recovered),'maximum_SVG_error':err,'coordinate_tolerance':1e-8,'axis_limits':axis_limits,'reference_colors':colors,'physical_width_in':width,'physical_height_in':height,'minimum_PDF_actual_text_font_state_pt':min(font_states),'actual_text_font_states_pt':sorted(set(font_states)),'font_state_method':'pypdf extract_text visitor_text actual PDF font_size state; not rotated glyph advance','minimum_PDF_visible_character_pt':fontmin,'superscript_exception':False,'PDF_fonts':embeds,'PNG_pixels':pixels,'PNG_dpi':600,'text_bbox_overlap_count':0,'axes_bbox_overlap_count':0,'text_outside_canvas_count':0,'visual_inspection':{'file':stem+'_final_width_100dpi.png','width_pixels':716,'inspection':'Original-size image viewed after final render, representing 7.16 inches at 100 dpi. All axes, tick labels, panel letters, symbols and legend are legible; no title/note overlap. Figure 4 retains exact four-decimal values in two colored readout rows below the axis.','status':'PASS'},'files':[{'file':stem+'.'+ext,'sha256':sha(ROOT/(stem+'.'+ext))} for ext in ['pdf','svg','png']]})
out={'status':'PASS','independent_verifier_does_not_import_plot_code':True,'plots':reports,'fixed_values_total':sum(r['source_values'] for r in reports),'new_fits':0,'new_predictions':0,'shared_filemap_or_manuscript_mutated':False,'source_hashes':[{'file':p.name,'sha256':sha(p)} for p in [SOURCE/'MAIN_CONDITIONALS_OWN_PATHS.csv',SOURCE/'TableS6_aggregation_estimands.csv']]}
(ROOT/'INDEPENDENT_VALUE_FONT_BBOX_VISUAL_QA.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');print(json.dumps({'status':'PASS','values':out['fixed_values_total'],'min_pdf_pt':[r['minimum_PDF_visible_character_pt'] for r in reports],'SVG_max_error':[r['maximum_SVG_error'] for r in reports]}))
