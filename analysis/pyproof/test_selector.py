import unittest
from mock import patch, Mock
from selector import MyPySelector
from ROOT import TH1I


class Testselector(unittest.TestCase):
    def test_get_zvtx(self):
        sel = MyPySelector()
        with patch('selector.MyPySelector.fChain') as mock_fChain:
            mock_fChain.event.zvtx = -9.9
            self.assertEqual(sel.get_zvtx(), 0)
            mock_fChain.event.zvtx = -10
            self.assertEqual(sel.get_zvtx(), 0)
            mock_fChain.event.zvtx = 9.9
            self.assertEqual(sel.get_zvtx(), 9)
            mock_fChain.event.zvtx = 10
            self.assertEqual(sel.get_zvtx(), 9)

    def test_get_bg(self):
        sel = MyPySelector()
        sel.pool = Mock()
        sel._bin_zvtx = TH1I('bin_zvtx', 'used to find z vtx',
                             10, -10, 10)

        sel.pool = [[(i, i)] for i in range(0, 10)]  # i equals z index
        with patch('selector.MyPySelector.fChain') as mock_fChain:
            iterator = sel.it_background_tracks(1)
            mock_fChain.event.zvtx = -9.9
            self.assertEqual(iterator.next(), (0, 0))
            iterator = sel.it_background_tracks(1)
            mock_fChain.event.zvtx = 9.9
            self.assertEqual(iterator.next(), (9, 9))

    class Fake_track(object):
        def __init__(self, val):
            self.pt = val
            self.phi = val
            self.eta = val

    def test_it_assocs(self):
        """mocking up a few tracks with float pt values"""
        sel = MyPySelector()
        sel.assoc_inter = Mock(spec=list)
        with patch('selector.MyPySelector.fChain') as mock_fChain:
            mock_fChain.track = [self.Fake_track(i/2.0) for i in range(-1, 5)]
            sel.assoc_inter = [0.5, 1.5]
            self.assertEqual(len(list(sel.it_assocs())), 2)  # exclude upper limit

    def test_it_triggerss(self):
        """mocking up a few tracks with float pt values"""
        sel = MyPySelector()
        sel.assoc_inter = Mock(spec=list)
        with patch('selector.MyPySelector.fChain') as mock_fChain:
            mock_fChain.track = [self.Fake_track(i/2.0) for i in range(-1, 10)]
            sel.assoc_inter = [2.0, 4.0]
            self.assertEqual(len(list(sel.it_assocs())), 4)  # exclude upper limit

    def test_validate_event(self):
        sel = MyPySelector()
        with patch('selector.MyPySelector.fChain') as mock_fChain:
            mock_fChain.event.trig = 1
            mock_fChain.event.zvtx = 999.0
            self.assertFalse(sel.validate_event())




















