"""Create a separate QC DOCX with logged mechanical edits; preserve supplied original."""
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED
import json
import re
import html
import csv
import platform
from copy import deepcopy
from xml.etree import ElementTree as E
from figure1_bootstrap import METHODS

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'GDF15_evidence_curve_DRAFT_v0.1_2-1.docx'
TARGET=ROOT/'GDF15_evidence_curve_QC.docx'
W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'
E.register_namespace('w',W)
NS={'w':W}

def plain(xml):
    return ''.join(html.unescape(m.group(1)) for m in re.finditer(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>',xml,re.S))

def replace_text(raw,old,new):
    aggregate=plain(raw);start=aggregate.index(old);end=start+len(old);offset=0;done=False
    def change(match):
        nonlocal offset,done
        text=html.unescape(match.group(2));a=offset;b=a+len(text);offset=b
        if b<=start or a>=end:return match.group()
        left=text[:max(0,start-a)];right=text[max(0,end-a):] if end<b else ''
        middle='' if done else new;done=True
        return match.group(1)+html.escape(left+middle+right,quote=False)+match.group(3)
    out=re.sub(r'(<w:t(?:\s[^>]*)?>)(.*?)(</w:t>)',change,raw,flags=re.S)
    assert plain(out)==aggregate[:start]+new+aggregate[end:]
    return out


def main():
    result=json.loads((ROOT/'results/results.json').read_text());m=result['models']['IPAH+HPAH']
    verified=json.loads((ROOT/'results/centerpiece_verification.json').read_text());assert verified['status']=='PASS'
    with ZipFile(SOURCE) as z:files={n:z.read(n) for n in z.namelist()}
    xml=files['word/document.xml'].decode();body=re.search(r'<w:body>(.*?)</w:body>',xml,re.S).group(1)
    blocks=list(re.finditer(r'<w:tbl(?:\s[^>]*)?>.*?</w:tbl>|<w:p(?:\s[^>]*)?>.*?</w:p>',body,re.S))
    assert plain(blocks[14].group()).startswith('A recent research letter')
    changes=[];replacements={}
    def edit(i,old,new,reason):
        original=replacements.get(i,blocks[i].group());updated=replace_text(original,old,new)
        replacements[i]=updated;changes.append({'body_block':i,'before':old,'after':new,'reason':reason})
    edit(2,plain(blocks[2].group()),'QC revision for author review. Verified mechanical corrections are logged in qc/manuscript_changes.json. Figure 1 uses the approved finite-fit bootstrap handling. Scientific framing is preserved; remaining author items are recorded in qc/published_data_freeze.md.','QC status')
    for i in (14,40):edit(i,'0.990',f"{m['auc_ci'][1]:.3f}",'Final 4000-replicate AUC CI independently reproduced')
    edit(41,'0.965',f"{m['internal_validation']['corrected_calibration_slope']:.3f}",'Final bootstrap calibration slope')
    before='0.000 against 0.000, 0.002 against 0.000, 0.012 against 0.030, 0.088 against 0.121 and 0.665 against 0.618'
    rows=m['calibration_quintiles'];parts=[f"{r['predicted_mean']:.3f} against {r['observed_proportion']:.3f}" for r in rows]
    edit(42,before,', '.join(parts[:-1])+' and '+parts[-1],'Independent calibration quintiles')
    t=m['threshold_comparison'];round10=lambda x:round(x/10)*10
    edit(14,'440 pg/mL',f"{round10(t['lr10_ci_width']):,} pg/mL",'Width rounded from unrounded interval')
    edit(56,'1,910 to 2,350',f"{round10(t['lr10_ci_low']):,} to {round10(t['lr10_ci_high']):,}",'Final stratified bootstrap interval')
    edit(56,'440',str(round10(t['lr10_ci_width'])),'Width rounded from unrounded interval')
    edit(58,'Interval endpoints are rounded to the nearest 10 pg/mL because bootstrap variation between seeds exceeds the final digit.',
         'Interval endpoints and widths are rounded separately to the nearest 10 pg/mL because bootstrap variation between seeds exceeds the final digit. Widths and the width ratio are calculated from unrounded endpoints.','Explain rounded endpoints versus rounded widths')
    # Table 1 cells are populated directly from independently verified rows.
    table=blocks[47].group();wrapper='<root xmlns:w="'+W+'">'+table+'</root>'
    tbl=E.fromstring(wrapper)[0];trs=tbl.findall('w:tr',NS);assert len(trs)==3
    def cell_text(cell,value):
        ts=cell.findall('.//w:t',NS);assert ts
        ts[0].text=value
        for el in ts[1:]:el.text=''
    cell_text(trs[0].findall('w:tc',NS)[0],'Row')
    for i,r in enumerate(verified['adjacent_observations'],1):
        values=[str(i),r['diagnosis'],f"{float(r['raw_concentration']):,.2f} pg/mL",f"{r['continuous_lr']:.2f}",
                f"{r['continuous_post_test_probability']:.0%}",f"{r['rule_lr']:.2f}" if r['rule_lr']<1 else f"{r['rule_lr']:.1f}",f"{r['rule_post_test_probability']:.0%}"]
        cells=trs[i].findall('w:tc',NS);assert len(cells)==len(values)
        for c,v in zip(cells,values):cell_text(c,v)
    replacements[47]=E.tostring(tbl,encoding='unicode')
    changes.append({'body_block':47,'before':plain(table),'after':plain(replacements[47]),'reason':'Programmatic verified Table 1; row numbers are not clinical IDs'})
    edit(49,'Post-test probabilities assume a pre-test probability of 10%.','Post-test probabilities assume a pre-test probability of 10%. Table row numbers are not reconstructed Patient_IDs. The observed values are from Fig_1a!C36 and Fig_1a!B27; their full-precision gap is '+verified['absolute_gap']+' pg/mL.','Source provenance')
    ethics='The original study reported approval of the UK National Cohort Study protocol (13/EE/0203) by the East of England Research Ethics Committee and of the UCSF PH Biobank and Lung Transplant Registry protocols (22-36931 and 13-10738, respectively) by the UCSF Institutional Review Board. The NCS gene-identification protocol (08/H0802/32) was approved by the UK Research Ethics Committee; the supplement reports written informed consent for whole-blood collection for those gene-identification studies [1].'
    edit(87,'[A3: AKIKO OR AMIT. The study approval and consent paragraph from the JCI supplemental material, with approving bodies named.]',ethics,'Verified supplement pp 5–6; consent scoped to source')
    funding='NIH (T32HL007731 to JTD; R01AG058659 to MAS; K24HL174231 to JPS; HL151552, VA CX002011, and CFF HAYS19AB3 to JRG; and R01HL164581 and R01HL153915 to AH); Gates Cambridge Trust (OPP1144 to EMDDDB); National Institute for Health and Care Research (NIHR) (NIHR203321 to DGK); United Kingdom National Cohort Study of Idiopathic and Heritable Pulmonary Arterial Hypertension Consortium grant to LSH; Cardiorespiratory NIHR Cambridge Bioresource Centre to MT; British Heart Foundation (SP/12/12/29836 and SP/18/10/33975 to the UK Pulmonary Arterial Hypertension Cohort Study Consortium); and Nina Ireland Program for Lung Health to PJW [1]'
    edit(88,'[funding list to be lifted verbatim from the source letter]',funding,'Published original-cohort funding, separately attributed')
    edit(88,'This secondary analysis received no specific funding [A10: GIORGIO to confirm student support].',
         'Jordan’s secondary-analysis work was supported through the 2026 Responsible AI Summer Research Program at Santa Clara University under the award “Transparent AI for Rare Disease Diagnosis”. The programme had no role in the analysis.','Giorgio funding confirmation')
    note='[A8: this sentence rests on a citation we have not yet been able to verify against PubMed. Either the source is confirmed and cited, or the sentence is deleted. Do not submit with it unverified.]'
    edit(77,note,'[6]','A8 verified, PMID 42275164; exact Discussion limitation')
    refs={115:'1. Prabhakar A, De Bie EMDD, DesJardin JT, et al. GDF15 is a putative biomarker for distinguishing pulmonary veno-occlusive disease and pulmonary arterial hypertension. J Clin Invest. 2026;136(2):e199159. doi:10.1172/JCI199159.',
          116:'2. Aguado B, Grynblat J, Budhram B, et al. Pulmonary veno-occlusive disease: a paradigm of diagnosis and therapeutic challenges in pulmonary hypertension. Curr Opin Pulm Med. 2025;31(5):443–455. doi:10.1097/MCP.0000000000001184.',
          120:'6. Andruska AM, Cruz-Utrilla A, Nossent EJ, et al. Identification of pulmonary arterial hypertension patients with venous or capillary involvement. Am J Respir Crit Care Med. 2026;212(9):2199–2210. doi:10.1093/ajrccm/aamag289.',
          121:'7. Bossuyt PM, Reitsma JB, Bruns DE, et al.; STARD Group. STARD 2015: an updated list of essential items for reporting diagnostic accuracy studies. BMJ. 2015;351:h5527. doi:10.1136/bmj.h5527.',
          122:'8. Steyerberg EW, Harrell FE Jr, Borsboom GJ, Eijkemans MJC, Vergouwe Y, Habbema JDF. Internal validation of predictive models: efficiency of some procedures for logistic regression analysis. J Clin Epidemiol. 2001;54(8):774–781. doi:10.1016/S0895-4356(01)00341-9.'}
    for i,text in refs.items():edit(i,plain(blocks[i].group()),text,'Verified reference identity')
    edit(124,plain(blocks[124].group()),'Bibliographic identities have been checked against published records. Sentence-level support and limits of full-text access are recorded in qc/reference_audit.csv. Author-owned and scientific interpretation issues remain open.','Avoid blanket unsupported reference-verification claim')
    edit(37,'[A9: JORDAN, version and package versions]',platform.python_version()+'; package versions are specified in requirements-lock.txt','Recorded environment; URL/DOI remain pending')
    for i,old,new in [(16,'genotyping.','genotyping [2,4].'),(17,'specificity.','specificity [1].'),(18,'two decades.','two decades [3].'),(22,'senior author,','senior author [1],'),(33,'standard bootstrap correction.','standard bootstrap correction [8].'),(36,'diagnostic accuracy components;','diagnostic accuracy components [5,7];')]:
        edit(i,old,new,'Explicit reference call with support recorded in audit')
    edit(33, 'Uncertainty was estimated by stratified case bootstrap with 4,000 replicates, resampling PVOD and comparator patients separately at their observed group sizes and re-deriving every quantity inside each replicate, so that quantities derived under different procedures remain comparable.',
         'The Figure 1 band was estimated by stratified case bootstrap, resampling PVOD and comparator patients separately at their observed group sizes. ' + METHODS + ' The existing AUC and paired threshold/landmark comparisons use stratified case resampling; comparator-specific bands and the retained supplementary intervals use the original parametric bootstrap.',
         'Approved Figure 1 finite-fit accounting and explicit scope; existing estimators unchanged')
    edit(98, 'the band is a 95% stratified case bootstrap interval over 4,000 replicates.',
         'the band is a pointwise 95% stratified case bootstrap interval from 3,991 finite fits out of 4,000 requested resamples; 9 resamples with complete separation were excluded. With these exclusions the band narrows in the tail, with its upper edge effectively unchanged at and below the published cutoff.',
         'Approved finite-fit Figure 1 caption')
    edit(98, 'it is drawn broken there and should not be read as evidence.',
         'it is drawn broken there and should not be read as evidence. The upper comparator support limit is also marked at 1,809.5 pg/mL; beyond it the curve is broken and the band is lighter.',
         'Caption matches both observed support markers in integrated Figure 1')
    edit(36,'the completed checklists are provided as supplementary material.','reporting mappings accompany this QC revision, with outstanding author items identified.','Checklists are not submission-complete')
    # Close only source-resolved rows in the author task table.
    for old,new in [('The study approval and consent paragraph from the JCI supplemental material, with approving bodies named','A3 source passage verified; see Declarations and qc/sources.md'),
                    ('We could not retrieve this citation from PubMed. Until it is verified the sentence that rests on it is unsupported and has been marked.','VERIFIED: Andruska et al. 2026; PMID 42275164; Discussion explicitly identifies dichotomization information loss.'),
                    ('Expected to be "no specific funding". Confirm student support arrangements.','Confirmed: 2026 Responsible AI Summer Research Program at Santa Clara University; “Transparent AI for Rare Disease Diagnosis”; programme had no role in analysis.')]:
        edit(6,old,new,'Resolve only confirmed A3/A8/A10 item')
    # Resolved source insertions inherit red placeholder formatting unless normalized.
    resolved_starts = ('The original study reported approval', 'NIH (T32HL007731',
                       '6. Andruska', '7. Bossuyt', '8. Steyerberg', platform.python_version(), '[6]')
    for i, raw in list(replacements.items()):
        def normalize(match):
            run = match.group()
            if any(plain(run).startswith(prefix) for prefix in resolved_starts):
                run = re.sub(r'<w:color\b[^>]*/>', '', run)
                run = re.sub(r'<w:b(?:\s[^>]*)?/>', '', run)
                run = re.sub(r'<w:bCs(?:\s[^>]*)?/>', '', run)
                run = run.replace('</w:rPr>', '<w:b w:val="0"/><w:bCs w:val="0"/><w:color w:val="000000"/></w:rPr>')
            return run
        replacements[i] = re.sub(r'<w:r(?:\s[^>]*)?>.*?</w:r>', normalize, raw, flags=re.S)
    # Existing source typo; preserve prose while restoring the missing word space.
    if 'pg/mL.They' in plain(blocks[45].group()):
        edit(45, 'pg/mL.They', 'pg/mL. They', 'Missing sentence space')
    for i in sorted(replacements,reverse=True):
        match=blocks[i];body=body[:match.start()]+replacements[i]+body[match.end():]
    xml=re.sub(r'<w:body>.*?</w:body>',lambda _: '<w:body>'+body+'</w:body>',xml,flags=re.S)
    xml = xml.replace('<w:t>', '<w:t xml:space="preserve">')
    E.fromstring(xml)
    files['word/document.xml']=xml.encode()
    images = [n for n in files if n.startswith('word/media/') and n.endswith('.png')]
    assert len(images) == 1, 'Expected exactly one supplied Figure 1 image'
    files[images[0]] = (ROOT / 'figures/figure1_lr_curve.png').read_bytes()
    with ZipFile(TARGET,'w',compression=ZIP_DEFLATED) as z:
        for name,content in files.items():
            info=ZipInfo(name,(2026,9,7,0,0,0));info.compress_type=ZIP_DEFLATED;z.writestr(info,content)
    (ROOT/'qc/manuscript_changes.json').write_text(json.dumps(changes,indent=2,ensure_ascii=False)+'\n')
    print(TARGET.name,len(changes),'logged edits')

if __name__=='__main__':main()
