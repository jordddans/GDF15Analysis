Hi Giorgio,

I would like to submit an abstract to BayLearn. I attached a two-page draft so you
have something concrete to react to.

The framing presents the published analysis and our model as two stages. Here is a
quick description of the full process I explained in our last meeting. I believe it is
worth submitting, because it contains a partial answer to the question you asked me to
work out — whether the other clinical variables can separate PVOD patients on their
own, and whether GDF15 adds anything on top of them.

**Stage one — the published result.** The paper found that GDF15 distinguishes PVOD
from similar forms of pulmonary hypertension and identified 1,658 pg/mL as a cutoff:
above that value patients were more likely to have PVOD, below it less likely. One
thing I want to flag about the data itself. I have been calling it our "practice
dataset," but that undersells it. It is the individual-level source data published
with the paper — 322 real patients, including 26 confirmed PVOD cases from three
registries. It is real published patient data that anyone can download and check,
which I think is the strongest thing about this submission.

**Stage two — our reanalysis.** We built a second model to examine how well the
underlying data support that cutoff. Rather than sorting patients into two groups by a
single dividing line, our model treats GDF15 as continuous and asks how much a
patient's exact level should raise or lower confidence in a PVOD diagnosis.

**On the parametric analysis.** We have partly answered this already. Each clinical
variable was tested individually against PVOD, and none discriminates well alone:
six-minute walk distance is the best at 0.67, with mean pulmonary artery pressure,
pulmonary vascular resistance, cardiac index, NT-proBNP and age all between 0.52 and
0.62. GDF15 alone reaches 0.96, and it is close to uncorrelated with the others, so it
appears to be adding information rather than restating the right-heart
catheterization. What we have not done is fit those variables jointly, or formally
test the improvement from adding GDF15 to a combined model. That is limited by missing
data rather than effort — requiring all five clinical variables at once leaves only 25
patients, and the ones who survive that filter are disproportionately PVOD cases,
which would bias the result rather than just weaken it. It may be more tractable in
the institutional cohort, and I think it is the natural next analysis there.

**What the reanalysis shows.** Our preliminary results support the paper's broader
conclusion that GDF15 is a strong and useful biomarker. They also suggest its meaning
changes gradually rather than suddenly at exactly 1,658. Values well below the cutoff
give strong evidence against PVOD, which is the most reliable single finding in the
data. Very high values give strong evidence in favor. But measurements near the cutoff
are much less decisive than a positive-or-negative result implies — at 1,658 itself the
uncertainty interval reaches down to "no information at all."

**A limitation I did not expect.** There are plenty of PVOD patients at high GDF15
levels: 18 of the 26 are above 2,500 pg/mL. What is missing is comparison patients.
Above roughly 2,000 pg/mL there are none, because the common subtypes never reach
those values in this data. Without both groups present you cannot measure the ratio
between them, so our model detects when it has left the range where both were observed
and reports a defensible lower bound instead of a falsely precise estimate.

**Why a machine learning audience would care.** The standard statistic used to check
overfitting barely moves here — corrected and uncorrected values differ in the fourth
decimal place, which would tell most readers that 26 cases are plenty. But removing a
single patient and refitting shifts the model's answer at the high end by 15 percent.
The aggregate statistic is stable exactly where the number a clinician would act on is
not. I think that generalizes well beyond this biomarker.

I would frame the whole thing as a reanalysis and extension of the published work, not
as a claim that its central finding is wrong. I would also state clearly that the
institutional patient-data analysis is ongoing and that nothing from it appears in the
abstract.

Two things I would like your guidance on: who should be included as authors, and
whether the original authors should be notified as a courtesy. This is secondary use
of data they published openly, which I understand is normal, but you will know the
convention better than I do.

Happy to jump on a call today if that is easier than going through the draft.

Best,
Jordan
