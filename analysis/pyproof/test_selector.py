import unittest
from mock import patch, Mock
from selector import MyPySelector
from math import pi
import numpy as np
import logging


class Testselector(unittest.TestCase):

    def load_selector(self, **kwargs):
        args = {"assoc_inter": [0.5, 1.0], "trigger_inter": [1.0, 2.0], "mc_and_recon_valid": False, "only_charge": None, "data_type": "MC", "track_cut": "golden", "pt_threshold": 0.0, "workers": None, "allow_0_trig": False, "events_have_recon_triggers": False, "append_to_dir": "", "events_have_one_mc_trigger": False, "use_chain": False, "exclude_tpc_border": False, "directory": "assocs0.5-1.0_trigs1.0-2.0_thresh0.0_cut_golden_type_MC", "data": "MC.dat"}
        for (k, v) in kwargs.items():
            args[k] = v
        # overwrite whatever old pickle was found
        for (k, v) in args.items():
            setattr(self, k, v)

        self.sel = MyPySelector()
        self.sel.GetOutputList = Mock()
        self.sel.GetOutputList.return_value = Mock()
        self.sel.GetOutputList().Add.return_value = True
        self.sel.SlaveBegin('tree')


    def load_hard_event(self, option='valid'):
        """fChain needs extra care, calling function needs to be decorated with
        @patch.object(MyPySelector, 'fChain', Mock())"""
        e = self.sel.fChain.event
        e.trig = 1; e.vtxstatus = 1; e.zvtxMC=0; e.ptmaxMC = 1 + self.pt_threshold
        e.cent = 50; e.zvtx=0; e.ptmax = 1 + self.pt_threshold;
        t = Mock()
        t.qMC=3; t.ptMC = 1.5; t.phiMC = 0.12 - pi/18; t.etaMC = 0.0
        t.q=3; t.pt = 1.5; t.phi = 0.12 - pi/18; t.eta = 0; t.filter=1; t.primary=1;
        t.dcaxy=0.0; t.dcaz = 0.0
        self.sel.get_tracks_generator = Mock(return_value=[t,])
        if option == 'valid':
            return None
        if option == 'bad threshold':
            e.ptmaxMC = e.ptmax = self.pt_threshold * .5
            return None
        
    @patch.object(MyPySelector, 'fChain', Mock())
    def test_process_event_without_trigger(self):
        """using asso 0.5 1 trigs 1 2"""
        self.load_selector()
        self.load_hard_event('valid')
        self.sel.get_triggers_assocs = Mock(return_value=([], [[0,0, .75, 1]]))
        self.sel.added_to_pool = Mock(return_value=False)
        self.sel.replace_pool = Mock(return_value=True)
        self.sel.get_background_tracks = Mock(return_value=np.asarray([[0,0, .75, 1]]))
        self.assertEqual(1, self.sel.Process(0))


    @patch.object(MyPySelector, 'fChain', Mock())
    def test_validate_mc_and_recon_mc_run(self):
        """Check that both types are check if mc_and_recon_valid is True"""
        self.load_selector()
        self.load_hard_event('valid')
        self.sel.mc_and_recon_valid = True
        self.sel.data_type = 'MC'
        self.sel.fChain.event.ptmaxMC = self.sel.fChain.event.ptmax = 5
        self.assertTrue(self.sel.validate_event_mc() and self.sel.validate_event_rec() and self.sel.validate_event_wrt_n_triggers())
        # event does not meet threshold:
        self.sel.fChain.event.ptmaxMC = self.sel.fChain.event.ptmax = 1
        self.sel.pt_threshold = 4
        self.assertFalse(self.sel.validate_event_mc() and self.sel.validate_event_rec() and self.sel.validate_event_wrt_n_triggers())


    @patch.object(MyPySelector, 'fChain', Mock())
    def test_validate_mc_and_recon_recon_run(self):
        """Check that both types are check if mc_and_recon_valid is True
        and its running on reconstructed data"""
        self.load_selector(datatype='reconstructed')
        self.load_hard_event('valid')
        self.sel.mc_and_recon_valid = True
        self.sel.fChain.event.ptmaxMC = self.sel.fChain.event.ptmax = 5
        self.assertTrue(self.sel.validate_event_mc() and self.sel.validate_event_rec() and self.sel.validate_event_wrt_n_triggers())
        # event does not meet threshold:
        self.sel.fChain.event.ptmaxMC = self.sel.fChain.event.ptmax = 1
        self.sel.pt_threshold = 4
        self.assertFalse(self.sel.validate_event_mc() and self.sel.validate_event_rec() and self.sel.validate_event_wrt_n_triggers())
        # Reconstructed does not meet requirement
        self.sel.fChain.event.ptmaxMC = 5
        self.sel.fChain.event.ptmax = 1
        self.sel.pt_threshold = 4
        self.assertFalse(self.sel.validate_event_mc() and self.sel.validate_event_rec() and self.sel.validate_event_wrt_n_triggers())
        # MC does not meet requirement
        self.sel.fChain.event.ptmaxMC = 1
        self.sel.fChain.event.ptmax = 5
        self.sel.pt_threshold = 4
        self.assertFalse(self.sel.validate_event_mc() and self.sel.validate_event_rec() and self.sel.validate_event_wrt_n_triggers())
        

    def test__validated_mc_track(self):
        """Test the mc track validator for various scenarios"""
        self.load_selector()
        self.sel.data_type = 'MC'
        track = Mock()
        track.etaMC = 0  # dont care about eta here

        # only pos charge, no TPC border treatment!
        self.sel.only_charge = 3; self.sel.exclude_tpc_border = False
        track.qMC=3; track.ptMC = 2; track.phiMC = 0.12 - pi/18        
        self.assertTrue(self.sel._validate_mc_track(track))
        track.qMC=-3
        self.assertFalse(self.sel._validate_mc_track(track))

        # only pos charge, with TPC border exclusion
        self.sel.only_charge = 3; self.sel.exclude_tpc_border = True
        track.qMC=3; track.ptMC = 2; track.phiMC = 0.12 - pi/18        
        self.assertFalse(self.sel._validate_mc_track(track))
        track.qMC = 3; track.phiMC = 0.3  - pi/18
        self.assertTrue(self.sel._validate_mc_track(track))

        # exclude TPC borders, include both charges
        self.sel.only_charge = None; self.sel.exclude_tpc_border = True
        track.qMC = 3; track.ptMC = 2; track.phiMC = 0.12 - pi/18        
        self.assertFalse(self.sel._validate_mc_track(track))
        track.qMC = -3; track.phiMC = 0.05  - pi/18
        self.assertTrue(self.sel._validate_mc_track(track))


    def test__validated_recon_track(self):
        """Test the recon track validator for various scenarios"""
        self.load_selector()
        self.sel.data_type = 'reconstructed'
        track = Mock()
        track.eta = 0  # dont care about eta here
        track.filter = 1; self.sel.track_filter = 1
        track.primary = 1; track.dcaxy = 0; track.dcaz = 0

        # only pos charge, no TPC border treatment!
        self.sel.only_charge = 3; self.sel.exclude_tpc_border = False
        track.q=3; track.pt = 2; track.phi = 0.12 - pi/18        
        self.assertTrue(self.sel._validate_recon_track(track))
        track.q=-3
        self.assertFalse(self.sel._validate_recon_track(track))

        # only pos charge, with TPC border exclusion
        self.sel.only_charge = 3; self.sel.exclude_tpc_border = True
        track.q=3; track.pt = 2; track.phi = 0.12 - pi/18        
        self.assertFalse(self.sel._validate_recon_track(track))
        track.q = 3; track.phi = 0.3  - pi/18
        self.assertTrue(self.sel._validate_recon_track(track))

        # exclude TPC borders, include both charges
        self.sel.only_charge = None; self.sel.exclude_tpc_border = True
        track.q= 3; track.pt = 2; track.phi = 0.12 - pi/18
        self.assertFalse(self.sel._validate_recon_track(track))
        track.q = -3; track.phi = 0.05  - pi/18
        self.assertTrue(self.sel._validate_recon_track(track))


    def test_is_track_at_tpc_border_for_MC(self):
        self.load_selector()
        self.sel.data_type = 'MC'
        track = Mock()
        track.qMC=3; track.ptMC = 2; track.phiMC = 0.12 - pi/18
        self.assertTrue(self.sel.is_track_at_tpc_border(track))
        track.qMC = -3; track.phiMC=.23 - pi/18
        self.assertTrue(self.sel.is_track_at_tpc_border(track))
        track.qMC = 3; track.phiMC = 0.3  - pi/18
        self.assertFalse(self.sel.is_track_at_tpc_border(track))

    def test_is_track_at_tpc_border_for_recon(self):
        self.load_selector()
        self.sel.data_type = 'reconstructed'
        track = Mock()
        track.q=3; track.pt = 2; track.phi = 0.12 - pi/18
        self.assertTrue(self.sel.is_track_at_tpc_border(track))
        track.q = -3; track.phi=.23 - pi/18
        self.assertTrue(self.sel.is_track_at_tpc_border(track))
        track.q = 3; track.phi = 0.3  - pi/18
        self.assertFalse(self.sel.is_track_at_tpc_border(track))


    # def test_get_zvtx(self):
    #     sel = MyPySelector()
    #     with patch('aapje.MyPySelector.fChain') as mock_fChain:
    #         mock_fChain.event.zvtx = -9.9
    #         self.assertEqual(sel.get_zvtx(), 0)
    #         mock_fChain.event.zvtx = -10
    #         self.assertEqual(sel.get_zvtx(), 0)
    #         mock_fChain.event.zvtx = 9.9
    #         self.assertEqual(sel.get_zvtx(), 9)
    #         mock_fChain.event.zvtx = 10
    #         self.assertEqual(sel.get_zvtx(), 9)

    # def test_get_bg(self):
    #     sel = MyPySelector()
    #     sel.pool = Mock()
    #     sel.pool = [[(i, i)] for i in range(0,10)]  # i equals z index
    #     with patch('aapje.MyPySelector.fChain') as mock_fChain:
    #         mock_fChain.event.zvtx = -9.9
    #         self.assertEqual(sel.get_bg(), (0, 0))
    #         mock_fChain.event.zvtx = 9.9
    #         self.assertEqual(sel.get_bg(), (9, 9))

    # class Fake_track(object):
    #     def __init__(self, val):
    #         self.pt = val
    #         self.phi = val
    #         self.eta = val    

    # def test_it_assocs(self):
    #     """mocking up a few tracks with float pt values"""
    #     sel = MyPySelector()
    #     sel.assoc_inter = Mock(spec=list)
    #     with patch('aapje.MyPySelector.fChain') as mock_fChain:
    #         mock_fChain.track = [self.Fake_track(i/2.0) for i in range(-1, 5)]
    #         sel.assoc_inter = [0.5, 1.5]
    #         self.assertEqual(len(list(sel.it_assocs())), 2)  # exclude upper limit

    # def test_it_triggerss(self):
    #     """mocking up a few tracks with float pt values"""
    #     sel = MyPySelector()
    #     sel.assoc_inter = Mock(spec=list)
    #     with patch('aapje.MyPySelector.fChain') as mock_fChain:
    #         mock_fChain.track = [self.Fake_track(i/2.0) for i in range(-1, 10)]
    #         sel.assoc_inter = [2.0, 4.0]
    #         self.assertEqual(len(list(sel.it_assocs())), 4)  # exclude upper limit

    # def test_validate_event(self):
    #     sel = MyPySelector()
    #     with patch('aapje.MyPySelector.fChain') as mock_fChain:
    #         mock_fChain.event.trig = 1
    #         mock_fChain.event.zvtx = 999.0
    #         self.assertFalse(sel.validate_event())

if __name__ == "__main__":
    unittest.main()
