Hi Giorgio,

I would like to submit an abstract to BayLearn. A two-page draft is attached so you
have something concrete to react to — I know I am past the window you asked for, so
please do not feel obliged to give it a full read. Even a quick yes/no on the framing
and the author list would be enough for me to submit before tomorrow's deadline.

One terminology point first. I have been calling the GDF15 data our "practice
dataset," but that undersells what it is. It is the individual-level source data
published with Prabhakar et al. in *J Clin Invest* this year — 322 real patients,
including 26 confirmed PVOD cases from three registries. It is real published patient
data that anyone can download and check, which I think is the single strongest thing
about this submission. I would rather describe it that way in the abstract.

The framing I settled on presents the published analysis and our model as two stages.

The published paper found that GDF15 distinguishes PVOD from similar forms of
pulmonary hypertension, and identified 1,658 pg/mL as a cutoff: above that value
patients were more likely to have PVOD, below it less likely.

We then built a second model to examine how well the underlying data support that
cutoff. Rather than sorting patients into two groups by a single dividing line, our
model treats GDF15 as continuous and asks how much a patient's exact level should
raise or lower confidence in a PVOD diagnosis.

Our preliminary results support the paper's broader conclusion that GDF15 is a strong
and useful biomarker. It separates PVOD from idiopathic and heritable PAH with an AUC
of 0.96, against 0.67 for six-minute walk distance, 0.62 for mean pulmonary artery
pressure, 0.60 for pulmonary vascular resistance, and 0.56 for NT-proBNP. It is also
close to uncorrelated with those markers, so it is supplying genuinely new information
rather than restating the right-heart catheterization.

At the same time, the results suggest the meaning of a GDF15 result changes gradually
rather than suddenly at exactly 1,658. Values well below the cutoff give strong
evidence against PVOD — this is the most reliable single finding in the data. Very
high values give strong evidence in favor. But measurements near the cutoff are much
less decisive than a positive-or-negative result implies: at 1,658 itself the evidence
is weak enough that its uncertainty interval reaches down to "no information at all."

The model also exposes a specific limitation of a 26-case sample, and it is not the
one I expected. There are plenty of PVOD patients at high GDF15 levels — 18 of the 26
are above 2,500 pg/mL. What is missing is comparison patients: above roughly 2,000
pg/mL there are none, because the common subtypes simply never reach those values in
this data. Without both groups present you cannot measure the ratio between them, so
our model detects when it has left the range where both were observed and reports a
defensible lower bound instead of a falsely precise estimate.

One further result is the reason I think this is interesting to a machine learning
audience and not only a clinical one. The standard summary statistic used to check
whether a model is overfitting barely moves here — corrected and uncorrected values
differ in the fourth decimal place, which would tell most readers that 26 cases are
plenty. But when we remove a single patient and refit, the model's answer at the high
end shifts by 15 percent. The aggregate statistic is stable exactly where the number a
clinician would actually act on is not. I think that generalizes well beyond this
biomarker.

I would frame the whole thing as a reanalysis and extension of the published work, not
as a claim that its central finding is wrong. I would also state clearly that the
institutional patient-data analysis is ongoing and that nothing from it appears in the
abstract.

Two things I would like your guidance on. First, who should be included as authors.
Second, whether the original authors should be notified as a courtesy — this is
secondary use of data they published openly, which I understand is normal, but you
will know the convention here better than I do.

Happy to jump on a call today if that is easier than going through the draft.

Best,
Jordan
