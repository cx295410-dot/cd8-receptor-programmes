import unittest
import numpy as np
import pandas as pd
from evidence_core import donor_summary, paired_inference, bh_adjust, scoped_clone, post_clone_history, pareto_front

class EvidenceCoreTests(unittest.TestCase):
    def test_donor_equal_summary_rejects_cell_weighting(self):
        frame=pd.DataFrame({'patient':['a']*100+['b']*20,'time':['post']*120,'target':[1]*100+[0]*20})
        out=donor_summary(frame,['target'],min_cells=20)
        self.assertAlmostEqual(out['target'].mean(),.5)
        self.assertNotAlmostEqual(frame.target.mean(),.5)

    def test_eligibility_is_per_patient_time(self):
        frame=pd.DataFrame({'patient':['a']*19+['b']*20,'time':['post']*39,'target':[1]*39})
        out=donor_summary(frame,['target'],min_cells=20)
        self.assertEqual(out.patient.tolist(),['b'])

    def test_exact_paired_test_uses_patients(self):
        f=pd.DataFrame({'patient':list('abcd')*2,'time':['pre']*4+['post']*4,'target':[0]*4+[1]*4})
        result=paired_inference(f,'target',seed=1,n_boot=500)
        self.assertEqual(result['n_pairs'],4)
        self.assertAlmostEqual(result['p_signflip'],.125)
        self.assertAlmostEqual(result['mean_delta'],1.)

    def test_missing_pairs_are_not_zero_imputed(self):
        f=pd.DataFrame({'patient':['a','a','b'],'time':['pre','post','post'],'target':[0,1,1]})
        result=paired_inference(f,'target',seed=1,n_boot=500)
        self.assertEqual(result['n_pairs'],1)

    def test_bh_matches_known_reference_and_preserves_missingness(self):
        q=bh_adjust([.01,.04,.03,np.nan])
        np.testing.assert_allclose(q[:3],[.03,.04,.04])
        self.assertTrue(np.isnan(q[3]))

    def test_fixed_family_keeps_missing_hypotheses_in_denominator(self):
        q=bh_adjust([.01,.04,.03,np.nan],family_size=4)
        np.testing.assert_allclose(q[:3],[.04,.0533333333333,.0533333333333])
        self.assertTrue(np.isnan(q[3]))

    def test_identical_sequences_in_different_donors_are_different_clones(self):
        self.assertNotEqual(scoped_clone('a','TRB:ACG'),scoped_clone('b','TRB:ACG'))
        self.assertIsNone(scoped_clone('a','TRB:ACG;TRB:AAA'))

    def test_clone_history_requires_both_times_in_same_patient(self):
        f=pd.DataFrame({'patient':['a','a','a','b'],'time':['pre','post','post','post'],
                        'clone':['a|TRB:1','a|TRB:1','a|TRB:2','b|TRB:1']})
        h=post_clone_history(f)
        self.assertEqual(h.tolist(),['pre','previously_detected','newly_detected','no_paired_baseline'])

    def test_pareto_marks_only_nondominated_candidates(self):
        x=np.array([[1.,1.],[.5,.5],[1.,0.],[0.,1.],[np.nan,.5]])
        self.assertEqual(pareto_front(x).tolist(),[True,False,False,False,False])

if __name__=='__main__': unittest.main(verbosity=2)
