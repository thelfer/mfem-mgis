#include<MFEMMGIS/Parameters.hxx>
#include<cassert>


int main()
{
	using namespace mfem_mgis;

  /** some values */
  real real_value = real(666.666);
  size_type size_type_value = 42;
  std::string string_value = "Merci à toi jeune aventurier de bien vouloir regarder des tests unitaires";

	/** test parameter types */ 
	Parameter param_real(real_value);
	Parameter param_size_type(size_type_value);
	Parameter param_string(string_value);

	std::vector<Parameter> param_vector = {param_real, param_size_type, param_string};

	/** test others constructors */
	Parameter param_real2(param_real);
  Parameter param_real3({real_value});

	Parameter param_size_type2(param_size_type);
	Parameter param_string2(param_string);
	Parameter param_vector2(param_vector); 
//  std::string msg_constructor_1 = "Parameter& Parameter::operator=(std::string_view src)";
//  Parameter test_constructor_parameter_1(msg_constructor_1);
  std::string_view msg_constructor_2 = "Parameter& Parameter::operator=(std::string_view src)";
  Parameter test_constructor_parameter_2(msg_constructor_2);
  Parameter& test_equal_parameter_string = test_constructor_parameter_2;
  test_equal_parameter_string = msg_constructor_2;

	/** test parameter functions */
	Parameter equal = param_vector2;

	/** test parameter constructors */
	Parameters params = Parameters();            // test default constructor
	Parameters params2 = Parameters();

  /** insert test */ 
	params.insert("real", param_real);           // fill a parameters 
	params.insert("size_type", param_size_type);
	params2.insert("string", param_string);
	params2.insert("vector", param_vector);

  /** Parameters contructors */
	Parameters params_copy(params);              // init a parameters with a parameters
  Parameters move(params);
	Parameters params_copy2(std::move(params));              // init a parameters with a parameters
  Parameter param_params(params_copy);         // init a parameter with a paramerters

  /** map */

  Parameters params3;
  std::map<std::string, Parameter> my_map;
  my_map["real"] = param_real;
  params3.insert(my_map); // test insert std::map<std::string, Parameter>&
  params3.insert({{"real", param_real}}); // test insert std::map<std::string, Parameter>&
  auto equal_test_ = Parameters{{"real", param_real}};

  /** trigger warning msg (doublon) */
  Parameters params4;
  params4.insert("real", param_real);
  params4.insert("real", param_real);

	/** test functions */
  auto test_get_parameter = get(params_copy,"real");
	auto test_get_parameters = params_copy.get("real"); // get
	if(std::get<real>(test_get_parameters) != real_value) 
	{
		std::cout << "Error: test_get" << std::endl;
		std::abort();
  }

  Parameter status = "failed";
  auto test_get_if_parameter_ok = get_if(params_copy,"real", status);
  auto test_get_if_parameter_no_ok = get_if(params_copy,"should no exist", status);

  /** check function */ 
  checkParameters(params, {"real", "size_type"});

	/** should not contain string parameter */
	bool contain_string = params_copy.contains("string");
	if( contain_string == true ) 
	{
		std::cout << "Error: contains part 1" << std::endl;
		std::abort();
  }
	params_copy.insert(params2);
	/** should contain string parameter */
	contain_string = params_copy.contains("string");
	if( contain_string == false ) 
	{
		std::cout << "Error: contains part 2" << std::endl;
		std::abort();
  }

	Parameters extractor = extract(params_copy, {"string","vector"});
	/** should not contain parameters in params (real and size_type) */
	auto test_extractor = extractor.contains("string") && extractor.contains("vector")
		&& !extractor.contains("real") && !extractor.contains("size_type");
	if( !test_extractor ) 
	{
		std::cout << "Error: test_extractor" << std::endl;
		std::abort();
  }

  /** test iterator */
  /** Juste test the number of item */
  int ccount(0), count(0);
  for(auto it = params2.cbegin() ; it != params2.cend() ; it++) { ccount++; }
  for(auto it = params2.begin() ; it != params2.end() ; it++) { count++; }
  if( ccount != count || count != 2) /* two items */
  {
		std::cout << "Error: test_iterator" << std::endl;
		std::abort();
  }

	std::cout << "SUCCESS" << std::endl;
  return 0;
}
