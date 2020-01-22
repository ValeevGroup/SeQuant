   using str_to_tensor_pair = std::pair<std::wstring, btas::Tensor<double> const *>;
   using str_to_tensor_map = std::map<str_to_tensor_pair::first_type,
                                       str_to_tensor_pair::second_type>;
  btas::Tensor<double> Fock_oo_raw(*Fock_oo),
                       Fock_ov_raw(*Fock_ov),
                       Fock_vv_raw(*Fock_vv),
                       G_oooo_raw(*G_oooo),
                       G_vvvv_raw(*G_vvvv),
                       G_ovvv_raw(*G_ovvv),
                       G_ooov_raw(*G_ooov),
                       G_oovv_raw(*G_oovv),
                       G_ovov_raw(*G_ovov),
                       T_ov_raw(*T_ov),
                       T_oovv_raw(*T_oovv);

   str_to_tensor_map btensor_map;
   btensor_map.insert(str_to_tensor_pair(L"f_oo",   &Fock_oo_raw));
   btensor_map.insert(str_to_tensor_pair(L"f_ov",   &Fock_ov_raw));
   btensor_map.insert(str_to_tensor_pair(L"f_vv",   &Fock_vv_raw));
   btensor_map.insert(str_to_tensor_pair(L"g_oooo", &G_oooo_raw));
   btensor_map.insert(str_to_tensor_pair(L"g_vvvv", &G_vvvv_raw));
   btensor_map.insert(str_to_tensor_pair(L"g_ovvv", &G_ovvv_raw));
   btensor_map.insert(str_to_tensor_pair(L"g_ooov", &G_ooov_raw));
   btensor_map.insert(str_to_tensor_pair(L"g_oovv", &G_oovv_raw));
   btensor_map.insert(str_to_tensor_pair(L"g_ovov", &G_ovov_raw));
   btensor_map.insert(str_to_tensor_pair(L"t_ov",   &T_ov_raw));
   btensor_map.insert(str_to_tensor_pair(L"t_oovv", &T_oovv_raw));



   auto result2 = interpret::eval_equation(factorized_expr, btensor_map);

   std::wcout  << "norm(result2) = " << std::sqrt(btas::dot(result2.tensor(), result2.tensor())) << "\n";
