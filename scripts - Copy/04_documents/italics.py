import re, copy
GENERA=['Escherichia coli','E. coli','Escherichia','Moraxella','Faucicola','Aeromonas','Shewanella',
 'Buttiauxella','Cupriavidus','Citrobacter','C. freundii','Klebsiella pneumoniae','Klebsiella','Neisseria meningitidis',
 'N. meningitidis','Neisseria','Campylobacter jejuni','C. jejuni','Campylobacter','Enterobacterales','Enterobacteriaceae',
 'Moraxellaceae','Pseudomonas','Edwardsiella tarda','Edwardsiella','Hafnia alvei','Hafnia','Salmonella enterica','Salmonella',
 'Acinetobacter','Burkholderia','Vibrio','Kluyvera','Ewingella','Rahnella','Plesiomonas shigelloides','Plesiomonas',
 'B. cochleicola','F. osloensis','A. media','S. baltica','Lwoffella lincolnii','Lwoffella','Serratia','Yersinia','Proteus',
 'Providencia','Pantoea','Pectobacterium','Photorhabdus','Xenorhabdus','Erwinia','Cronobacter','Cedecea','Leclercia',
 'Lelliottia','Morganella','Phytobacter','Pluralibacter','Tatumella','Yokenella','Franconibacter','Kosakonia','Dickeya']
GENES=[r'mcr-\d+(?:\.\d+)?',r'\beptA\b',r'\beptB\b',r'\beptC\b',r'\bcptA\b',r'\bopgE\b',r'\bopgB\b',
 r'\bnikAB\b',r'\bnikA\b',r'\bnikB\b',r'\bglnD\b',r'\bhemA\b',r'\blolB\b',r'\btyrS\b',r'\btnpA\b']
SPECIES=[r'\b(?:osloensis|mancuniensis|boevrei|atlantae|lincolnii|media|baltica|hydrophila|jandaei|sobria|rivipollensis|cochleicola|freundii|coli|enterica|pneumoniae|meningitidis|jejuni|tarda|alvei|shigelloides|frigidimarina|porci|catarrhalis|algae|tetraodonis|mytilicola|bicestrii|khirikhana|salmonicida|perminowiae|kookii|universalis|cryocrescens|parvum|gergoviae|fluvialis|cowanii|ursingii|pulveris|morganii|variigena|dadantii|dispersa|ptyseos|psychrotolerans|respiraculi|sedlakii|pasteurii|roggenkampii|kobei|hangzhouensis|pseudotuberculosis)\b']
PAT=re.compile('|'.join([r'\b'+re.escape(g)+r'\b' for g in GENERA]+GENES+SPECIES))
def apply(doc):
    for p in doc.paragraphs:
        _fix(p)
    for tb in doc.tables:
        for row in tb.rows:
            for c in row.cells:
                for p in c.paragraphs: _fix(p)
def _fix(p):
    runs=list(p.runs)
    if not runs: return
    pieces=[]
    for r in runs:
        txt=r.text
        if not txt: continue
        if r.italic:
            pieces.append((txt,r,True)); continue
        last=0
        for m in PAT.finditer(txt):
            if m.start()>last: pieces.append((txt[last:m.start()],r,False))
            pieces.append((m.group(0),r,True))
            last=m.end()
        if last<len(txt): pieces.append((txt[last:],r,False))
    if len(pieces)==len(runs) and all(pieces[i][2]==bool(runs[i].italic) for i in range(len(runs))): return
    for r in runs:
        r._element.getparent().remove(r._element)
    for txt,src,ital in pieces:
        nr=p.add_run(txt)
        nr.bold=src.bold; nr.underline=src.underline; nr.italic=ital
        if src.font.size: nr.font.size=src.font.size
