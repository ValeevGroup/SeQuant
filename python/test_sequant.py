import unittest
import _sequant as sequant
from _sequant import Tensor, Sum, Product, Expr

#sequant.IndexSpace.occupied = "i"

def visit(visitor, expr):
  if isinstance(expr,Sum): return visitor.Sum(*expr.summands)
  if isinstance(expr,Product): return visitor.Product(*expr.factors)
  if isinstance(expr,Tensor): return visitor.Tensor(expr)
  assert False, expr

class TestSequant(unittest.TestCase):

  def _test_expressions(self):
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

  def test_ccsd(self):
    from _sequant.mbpt import A,H,T,VacuumAverage
    ccd = VacuumAverage( A(2) * H() * T(2) * T(2), [(1, 2), (1, 3)] );
    print (ccd.latex)

    class String:
      def Sum(self,*args):
        return "(%s)" % " + ".join(visit(self,a) for a in args)
      def Product(self,*args):
        return " * ".join(visit(self,a) for a in args)
      def Tensor(self,arg):
        return '%s[%s]' % (arg.label,",".join('"%s"' % a for a in arg.braket))

    s = visit(String(), ccd)
    print (s)


if __name__ == '__main__':
  unittest.main()
