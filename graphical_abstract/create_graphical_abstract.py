"""Draw a deterministic, one-page scientific graphical abstract.

Requires reportlab, Pillow, pypdf and Poppler's pdftoppm. Fonts are embedded
in the PDF but the font software is not distributed. No models/data are run.
"""
import argparse
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

from PIL import Image
from pypdf import PdfReader
from pypdf.generic import ContentStream
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas

TITLE = 'Training-Test Representation and Reference-Path Sensitivity in a Satellite-Based Maize Gross Primary Production Benchmark'
CAPTION = (
    'Graphical Abstract. A matched-case maize gross primary production (GPP) benchmark compares '
    'uniform-cell MODIS inputs (B) with normalized monthly retained-core weights (C). B and C share '
    'cells, accepted scenes, meteorology, observations and validation roles. Stored B- and C-trained '
    'models are tested on both representations without refitting. In the matrix, the first letter '
    'denotes training and the second test input; fitting/scoring references and aggregation are '
    'held fixed within each four-cell comparison. The principal random-forest (RF) callout shows '
    '2.20-2.51% lower equal-site-season root mean square error (RMSE) for CC versus BB across the '
    'dependent own-fit/own-score VUT/CUT reference paths. However, C versus B training lowers RMSE '
    'at fixed B test input and raises it at fixed C test input. The 642-case principal cohort and '
    'separate 437-case Mead hindcast are not pooled; the 120 saved models comprise 100 principal '
    'and 20 hindcast objects. These are conditional historical comparisons, not physical daily '
    'footprint correction or independent forecast/transfer validation.'
)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--font', type=Path, required=True)
    p.add_argument('--pdftoppm', required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    assert args.font.is_file(), 'Supply a licensed Arial-compatible sans-serif TTF via --font.'
    pdfmetrics.registerFont(TTFont('GA-Sans', str(args.font)))
    W, H = 660, 295
    colors = {'ink': '#17293D', 'muted': '#4B5563', 'B': '#245891',
              'C': '#096B63', 'adverse': '#963D25', 'line': '#919AA4'}
    pdf = args.out / 'gagraphic.pdf'
    c = canvas.Canvas(str(pdf), pagesize=(W, H), invariant=1, pageCompression=1)
    c.setTitle(TITLE + ' - Graphical Abstract')
    c.setAuthor('Jonghan Ko')
    c.setSubject('Matched inputs, saved-model cross-evaluation, conditional RF contrasts')
    c.setCreator('Deterministic vector artwork; create_graphical_abstract.py')
    c.setFillColorRGB(1, 1, 1)
    c.rect(0, 0, W, H, fill=1, stroke=0)
    text = []

    def label(value, x, y, size=16, color='ink', align='left', region=None):
        """Coordinates y are measured downward from page top."""
        width = pdfmetrics.stringWidth(value, 'GA-Sans', size)
        x0 = x if align == 'left' else x-width/2
        ascent = pdfmetrics.getAscent('GA-Sans') * size / 1000
        descent = abs(pdfmetrics.getDescent('GA-Sans')) * size / 1000
        bbox = [x0, y-ascent, x0+width, y+descent]
        assert min(bbox[0], bbox[1]) >= 8 and bbox[2] <= W-8 and bbox[3] <= H-8, (value, bbox)
        if region:
            assert bbox[0] >= region[0] and bbox[2] <= region[1], (value, bbox, region)
        c.setFont('GA-Sans', size)
        c.setFillColor(colors[color])
        c.drawString(x0, H-y, value)
        text.append({'text': value, 'font_pt': size, 'color': colors[color], 'bbox': bbox})

    def line(x1, y1, x2, y2, color='line', width=1):
        c.setStrokeColor(colors[color])
        c.setLineWidth(width)
        c.line(x1, H-y1, x2, H-y2)

    def arrow(x1, x2, y):
        line(x1, y, x2, y, 'muted', 1.5)
        line(x2-5, y-4, x2, y, 'muted', 1.5)
        line(x2-5, y+4, x2, y, 'muted', 1.5)

    label('Maize GPP: training × test inputs', 18, 32, 22)
    label('Matched cases', 18, 77, region=(18, 175))
    label('MODIS + AgERA5', 18, 108, region=(18, 175))
    label('Same cells/scenes', 18, 133, region=(18, 175), color='muted')
    label('B: Uniform cells', 18, 174, color='B', region=(18, 175))
    label('C: Core weights', 18, 207, color='C', region=(18, 175))
    arrow(181, 201, 157)

    label('2 × 2 evaluation', 323, 77, align='center', region=(208, 438))
    label('Test B', 321, 106, align='center')
    label('Test C', 397, 106, align='center')
    label('Train B', 208, 145, color='B', region=(208, 278))
    label('Train C', 208, 194, color='C', region=(208, 278))
    cells = [('BB', 285, 118, 'B'), ('BC', 361, 118, 'line'),
             ('CB', 285, 167, 'line'), ('CC', 361, 167, 'C')]
    for name, x, y, color in cells:
        c.setStrokeColor(colors[color])
        c.setLineWidth(1.6 if name in ['BB', 'CC'] else 0.8)
        c.rect(x, H-y-45, 72, 45, fill=0, stroke=1)
        label(name, x+36, y+29, 18, color=color if name in ['BB', 'CC'] else 'ink', align='center')
    label('Saved models', 323, 231, align='center')
    label('No refitting', 323, 251, align='center', color='muted')
    arrow(442, 462, 157)

    label('Principal RF', 478, 77, region=(478, 644))
    label('2.20-2.51%', 478, 108, 23, 'C', region=(478, 644))
    label('lower pipeline RMSE', 478, 133, region=(478, 644))
    label('C vs B training', 478, 166, region=(478, 644))
    label('B test: lower RMSE', 478, 194, color='C', region=(478, 644))
    label('C test: higher RMSE', 478, 220, color='adverse', region=(478, 644))
    label('Pipeline gains can mask opposite training contrasts', W/2, 282, 18, align='center')
    c.showPage()
    c.save()

    # Render the vector master, then write pixel-exact IEEE and larger review variants.
    temp = args.out / '_render'
    temp.mkdir(exist_ok=True)
    for scale, name in [(1, 'gagraphic.png'), (4, 'gagraphic_preview.png')]:
        prefix = temp / f'ga_{scale}'
        subprocess.run([args.pdftoppm, '-png', '-singlefile', '-scale-to-x', str(W*scale),
                        '-scale-to-y', str(H*scale), str(pdf), str(prefix)], check=True)
        with Image.open(str(prefix)+'.png') as im:
            assert im.size == (W*scale, H*scale)
            rendered = im.convert('RGB')
            if scale == 1:
                # An indexed palette meets the IEEE size recommendation while
                # retaining anti-aliased labels and the fixed vector geometry.
                rendered = rendered.quantize(colors=256, method=Image.Quantize.MEDIANCUT,
                                             dither=Image.Dither.NONE)
            rendered.save(args.out/name, optimize=True, compress_level=9, dpi=(300, 300))
    shutil.rmtree(temp)
    (args.out/'GRAPHICAL_ABSTRACT_CAPTION_EN.txt').write_text(CAPTION+'\n', encoding='utf-8', newline='\n')
    (args.out/'GRAPHICAL_ABSTRACT_CAPTION_EN.md').write_text('# Graphical Abstract caption\n\n'+CAPTION+'\n', encoding='utf-8', newline='\n')

    # Machine preflight complements visual inspection; it is not a substitute for it.
    reader = PdfReader(pdf)
    sizes = [float(values[1]) for values, op in ContentStream(reader.pages[0].get_contents(), reader).operations if op == b'Tf']
    expected_labels = [r['text'] for r in text]
    extracted = reader.pages[0].extract_text()
    for label_text in expected_labels:
        assert label_text in extracted, label_text
    overlaps = []
    for i, a in enumerate(text):
        for b in text[i+1:]:
            ax1, ay1, ax2, ay2 = a['bbox']; bx1, by1, bx2, by2 = b['bbox']
            if min(ax2,bx2)-max(ax1,bx1) > 0.01 and min(ay2,by2)-max(ay1,by1) > 0.01:
                overlaps.append([a['text'], b['text']])
    assert not overlaps, overlaps

    def luminance(code):
        values = [int(code[i:i+2],16)/255 for i in [1,3,5]]
        lin = [v/12.92 if v <= 0.04045 else ((v+0.055)/1.055)**2.4 for v in values]
        return sum(a*b for a,b in zip(lin,[0.2126,0.7152,0.0722]))

    ratios = {key: 1.05/(luminance(colors[key])+0.05) for key in ['ink','muted','B','C','adverse']}
    assert min(ratios.values()) >= 4.5
    png = args.out/'gagraphic.png'
    assert png.stat().st_size < 45000, png.stat().st_size
    with Image.open(png) as im:
        dpi = im.info['dpi']
        assert min(dpi) >= 299.9
    qa = {'status':'PASS', 'pages':len(reader.pages), 'pdf_points':[W,H],
          'IEEE_png_pixels':[W,H], 'preview_pixels':[4*W,4*H],
          'IEEE_png_bytes':png.stat().st_size, 'IEEE_png_dpi':list(dpi),
          'minimum_PDF_text_font_pt':min(sizes), 'visible_text_minimum_pixels_at_660px':16,
          'font_embedded':all('/FontDescriptor' in f.get_object() and
                              '/FontFile2' in f.get_object()['/FontDescriptor']
                              for f in reader.pages[0]['/Resources']['/Font'].values()
                              if f.get_object().get('/Subtype') == '/TrueType'),
          'visible_words':sum(len(r['text'].split()) for r in text),
          'text_overlap_pairs':overlaps, 'text_contrast_ratios_on_white':ratios,
          'text_runs':text, 'caption':CAPTION, 'models_executed':0,
          'new_results':0, 'public_uploads':0, 'journal_submissions':0,
          'visual_QA':'PENDING'}
    (args.out/'PREFLIGHT.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in qa.items() if k not in ['text_runs','caption']},ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
