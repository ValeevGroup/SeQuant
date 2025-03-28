Coupled-Cluster Class
==========================

SeQuant's coupled-cluster (CC) class (see :class:`sequant::mbpt::CC`) supports derivations of ground and excited state methods,
using traditional, unitary and orbital-optimized ansatze (see :enum:`sequant::mbpt::CC::Ansatz`).
Expressions are generated in spin-orbital basis, and they can post-processed using SeQuant's spin-tracing capabilities.

Examples
--------
Here are some examples of using the CC class for deriving expressions.

CC amplitude equations
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

    using namespace sequant::mbpt;
    // Traditional CCSD
    auto t_eqs = CC{2}.t();
    std::wcout << "R[1]: " << to_latex(t_eqs[1]) << "\n"
               << "R[2]: " << to_latex(t_eqs[2]) << "\n";

    auto l_eqs = CC{2}.λ();
    std::wcout << "L[1]: " << to_latex(l_eqs[1]) << "\n"
               << "L[2]: " << to_latex(l_eqs[2]) << "\n";

    // Unitary CCSD
    auto Ut_eqs = CC{2, CC::Ansatz::U}.t();
    std::wcout << "R[1]: " << to_latex(Ut_eqs[1]) << "\n"
               << "R[2]: " << to_latex(Ut_eqs[2]) << "\n";


EOM-CC equations
^^^^^^^^^^^^^^^^^^

.. code-block:: c++

    using namespace sequant::mbpt;
    // EE-EOM-CCSD
    auto r_eqs = CC{2}.eom_r(nₚ(2), nₕ(2));
    std::wcout << "R[1]: " << to_latex(r_eqs[1]) << "\n"
               << "R[2]: " << to_latex(r_eqs[2]) << "\n";
