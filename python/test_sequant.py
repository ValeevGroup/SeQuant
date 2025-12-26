import unittest
import _sequant as sq
from _sequant import Tensor, Sum, Product, Constant, Expr, zRational

#sq.IndexSpace.occupied = "i"

def visit(visitor, expr):
  if isinstance(expr,Sum): return visitor.Sum(*expr.summands)
  if isinstance(expr,Product): return visitor.Product(expr.scalar,*expr.factors)
  if isinstance(expr,Tensor): return visitor.Tensor(expr)
  if isinstance(expr,Constant): return visitor.Constant(expr)
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
    from _sequant.mbpt import A,H,T,T_,VacuumAverage
    ccd = VacuumAverage( A(-2) * H() * T_(2) * T_(2), [("h", "t")] );  # explicit list of operator connections ..
    ccsd = VacuumAverage( A(-2) * H() * T(2) * T(2));                  # .. is not needed since H and T are connected by default
    print (ccsd.latex)

    class String:
      def Sum(self,*args):
        return "(%s)" % " + ".join(visit(self,a) for a in args)
      def Product(self,pfac,*factors):
        pfac_str = (str(pfac) + " * ") if pfac != 1 else ""
        return pfac_str + " * ".join(visit(self,a) for a in factors)
      def Tensor(self,arg):
        return '%s[%s]' % (arg.label,",".join('"%s"' % a for a in arg.braket))
      def Constant(self,arg):
        return '%s' % arg
      def zRational(self,arg):
        return '%s' % arg

    s = visit(String(), ccd)
    print (s)


if __name__ == '__main__':
  unittest.main()
