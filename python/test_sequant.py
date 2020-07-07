import unittest
import _sequant as sequant
from _sequant import Tensor, Sum, Product, Expr

sequant.IndexSpace.occupied = "i"

class TestSequant(unittest.TestCase):

  def test_expressions(self):
    t = Tensor("T", ["i_1","i_2"], ["i_1"])

    s = t + t - t
    self.assertEqual(type(s), Sum)
    self.assertEqual([ type(e) for e in s.summands], [ Tensor, Tensor, Product ])
    self.assertFalse(s.factors)

    p = s*t
    self.assertEqual(type(p), Product)
    self.assertEqual([type(e) for e in p.factors], [ Sum, Tensor ])
    self.assertFalse(p.summands)

    self.assertEqual(t.latex, "{T^{{i_1}}_{{i_1}{i_2}}}")
    self.assertTrue((p+s).latex)

if __name__ == '__main__':
  unittest.main()
