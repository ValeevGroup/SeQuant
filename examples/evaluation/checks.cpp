  // CHECKS
  // auto t1 = std::make_shared<sequant::Tensor>(
  //     sequant::Tensor(L"g", {L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
  // auto t2 = std::make_shared<sequant::Tensor>(
  //     sequant::Tensor(L"g", {L"a_1", L"a_2"}, {L"i_1", L"i_2"}));

  // auto sum1 = std::make_shared<Sum>(Sum{});
  // sum1->append(t1);
  // sum1->append(t2);

  // auto t3 = std::make_shared<sequant::Tensor>(
  //     sequant::Tensor(L"g", {L"a_4", L"a_3"}, {L"i_1", L"i_2"}));
  // auto t4 = std::make_shared<sequant::Tensor>(
  //     sequant::Tensor(L"g", {L"i_1", L"i_2"}, {L"a_4", L"a_3"}));

  // auto sum2 = std::make_shared<Sum>(Sum{});
  // sum2->append(t3);
  // sum2->append(t4);

  // auto evt_ptr1 = std::make_shared<EvalTensor>(EvalTensor(t1));
  // auto evt_ptr2 = std::make_shared<EvalTensor>(EvalTensor(t2));
  // print_evtensor(evt_ptr1);
  // print_evtensor(evt_ptr2);

  // print_evtensor(std::make_shared<EvalTensor>(EvalTensor(sum1)));
  // wcout << endl << endl;
  // print_evtensor(std::make_shared<EvalTensor>(EvalTensor(sum2)));
  //
  // auto fac1 = std::make_shared<Tensor>(
  //     Tensor(L"g", {L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
  // auto fac2 = std::make_shared<Tensor>(
  //     Tensor(L"g", {L"a_1", L"i_5"}, {L"i_3", L"i_4"}));
  // auto prod1 = std::make_shared<Product>(Product{});
  // prod1->append(1.0, fac1);
  // prod1->append(2.5, fac2);
  // print_evtensor(std::make_shared<EvalTensor>(EvalTensor(prod1)));
  // std::wcout << "\n\n";

  // auto fac3 = std::make_shared<Tensor>(
  //     Tensor(L"g", {L"a_1", L"a_2"}, {L"i_3", L"i_4"}));
  // auto fac4 = std::make_shared<Tensor>(
  //     Tensor(L"g", {L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
  // auto prod2 = std::make_shared<Product>(Product{});
  // prod2->append(3.5, fac1);
  // prod2->append(1.0, fac2);

  // print_evtensor(std::make_shared<EvalTensor>(EvalTensor(prod2)));
  // std::wcout << "\n\n";

  // auto fac5 = std::make_shared<Tensor>(Tensor(L"g", {L"a_3"}, {L"a_4"}));
  // auto prod3 = std::make_shared<Product>(Product{});
  // prod3->append(1, fac1);
  // prod3->append(1, fac2);
  // prod3->append(1, fac5);

  // print_evtensor(std::make_shared<EvalTensor>(EvalTensor(prod3)));
  // std::wcout << "\n\n";
  // auto prod4 = std::make_shared<Product>(Product{});
  // prod4->append(7, fac1);
  // prod4->append(7, fac5);
  // prod4->append(1, fac2);
  // auto prod4_evt_ptr = std::make_shared<EvalTensor>(EvalTensor(prod4));
  // prod4_evt_ptr->fill_btas_indices();
  // print_evtensor(prod4_evt_ptr);
  // std::wcout << "\n\n";

  // counting how many times different evaluations occur
  // auto factorized_evptr =
  //     std::make_shared<EvalTensor>(EvalTensor(factorized_expr));

  // hash_count_map counts;
  // count_hashes(factorized_evptr, counts);
  // for (const auto& items : counts) {
  //   if (items.second > 1)
  //     std::wcout
  //       << "hash: " << items.first
  //       << "    counts: " << items.second
  //       << "\n";
  // }
  // // //
