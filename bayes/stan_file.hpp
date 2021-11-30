
// Code generated by stanc v2.28.1
#include <stan/model/model_header.hpp>
namespace stan_file_model_namespace {

using stan::io::dump;
using stan::model::assign;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using namespace stan::math;


stan::math::profile_map profiles__;
static constexpr std::array<const char*, 15> locations_array__ = 
{" (found before start of program)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 8, column 3 to column 10)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 9, column 3 to column 10)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 18, column 3 to column 57)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 13, column 4 to column 35)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 14, column 3 to column 21)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 15, column 3 to column 21)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 2, column 1 to column 16)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 3, column 10 to column 11)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 3, column 3 to column 15)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 4, column 10 to column 11)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 4, column 3 to column 15)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 5, column 8 to column 9)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 5, column 1 to column 14)",
 " (in '/Users/jakob/Documents/Studium/master_thesis/bayes/stan_file.stan', line 18, column 9 to column 10)"};



class stan_file_model final : public model_base_crtp<stan_file_model> {

 private:
  int n;
  Eigen::Matrix<double, -1, 1> x__;
  Eigen::Matrix<double, -1, 1> y__;
  Eigen::Matrix<double, -1, 1> dy__; 
  Eigen::Map<Eigen::Matrix<double, -1, 1>> x{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double, -1, 1>> y{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double, -1, 1>> dy{nullptr, 0};
 
 public:
  ~stan_file_model() { }
  
  inline std::string model_name() const final { return "stan_file_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.28.1", "stancflags = "};
  }
  
  
  stan_file_model(stan::io::var_context& context__,
                  unsigned int random_seed__ = 0,
                  std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static constexpr const char* function__ = "stan_file_model_namespace::stan_file_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 7;
      context__.validate_dims("data initialization","n","int",
           std::vector<size_t>{});
      n = std::numeric_limits<int>::min();
      
      current_statement__ = 7;
      n = context__.vals_i("n")[(1 - 1)];
      current_statement__ = 7;
      check_greater_or_equal(function__, "n", n, 0);
      current_statement__ = 8;
      validate_non_negative_index("x", "n", n);
      current_statement__ = 9;
      context__.validate_dims("data initialization","x","double",
           std::vector<size_t>{static_cast<size_t>(n)});
      x__ = Eigen::Matrix<double, -1, 1>(n);
      new (&x) Eigen::Map<Eigen::Matrix<double, -1, 1>>(x__.data(), n);
      
      {
        std::vector<local_scalar_t__> x_flat__;
        current_statement__ = 9;
        x_flat__ = context__.vals_r("x");
        current_statement__ = 9;
        pos__ = 1;
        current_statement__ = 9;
        for (int sym1__ = 1; sym1__ <= n; ++sym1__) {
          current_statement__ = 9;
          assign(x, x_flat__[(pos__ - 1)],
            "assigning variable x", index_uni(sym1__));
          current_statement__ = 9;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 10;
      validate_non_negative_index("y", "n", n);
      current_statement__ = 11;
      context__.validate_dims("data initialization","y","double",
           std::vector<size_t>{static_cast<size_t>(n)});
      y__ = Eigen::Matrix<double, -1, 1>(n);
      new (&y) Eigen::Map<Eigen::Matrix<double, -1, 1>>(y__.data(), n);
      
      {
        std::vector<local_scalar_t__> y_flat__;
        current_statement__ = 11;
        y_flat__ = context__.vals_r("y");
        current_statement__ = 11;
        pos__ = 1;
        current_statement__ = 11;
        for (int sym1__ = 1; sym1__ <= n; ++sym1__) {
          current_statement__ = 11;
          assign(y, y_flat__[(pos__ - 1)],
            "assigning variable y", index_uni(sym1__));
          current_statement__ = 11;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 12;
      validate_non_negative_index("dy", "n", n);
      current_statement__ = 13;
      context__.validate_dims("data initialization","dy","double",
           std::vector<size_t>{static_cast<size_t>(n)});
      dy__ = Eigen::Matrix<double, -1, 1>(n);
      new (&dy) Eigen::Map<Eigen::Matrix<double, -1, 1>>(dy__.data(), n);
      
      {
        std::vector<local_scalar_t__> dy_flat__;
        current_statement__ = 13;
        dy_flat__ = context__.vals_r("dy");
        current_statement__ = 13;
        pos__ = 1;
        current_statement__ = 13;
        for (int sym1__ = 1; sym1__ <= n; ++sym1__) {
          current_statement__ = 13;
          assign(dy, dy_flat__[(pos__ - 1)],
            "assigning variable dy", index_uni(sym1__));
          current_statement__ = 13;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 14;
      validate_non_negative_index("y_tilde", "n", n);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1 + 1;
    
  }
  
  template <bool propto__, bool jacobian__ , typename VecR, typename VecI, 
  stan::require_vector_like_t<VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "stan_file_model_namespace::log_prob";
    (void) function__;  // suppress unused var warning
    
    try {
      local_scalar_t__ a;
      a = DUMMY_VAR__;
      
      current_statement__ = 1;
      a = in__.template read<local_scalar_t__>();
      local_scalar_t__ b;
      b = DUMMY_VAR__;
      
      current_statement__ = 2;
      b = in__.template read<local_scalar_t__>();
      {
        current_statement__ = 4;
        lp_accum__.add(
          normal_lpdf<propto__>(y, add(b, multiply(a, stan::math::cos(x))),
            dy));
        current_statement__ = 5;
        lp_accum__.add(normal_lpdf<propto__>(a, 0, 100));
        current_statement__ = 6;
        lp_accum__.add(normal_lpdf<propto__>(b, 0, 100));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, 
  stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, 
  stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr> 
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    (void) propto__;
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    int current_statement__ = 0; 
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    constexpr bool jacobian__ = false;
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "stan_file_model_namespace::write_array";
    (void) function__;  // suppress unused var warning
    
    try {
      double a;
      a = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      a = in__.template read<local_scalar_t__>();
      double b;
      b = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      b = in__.template read<local_scalar_t__>();
      out__.write(a);
      out__.write(b);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
      std::vector<double> y_tilde;
      y_tilde = std::vector<double>(n, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 3;
      assign(y_tilde,
        normal_rng(add(b, multiply(a, stan::math::cos(x))), dy, base_rng__),
        "assigning variable y_tilde");
      out__.write(y_tilde);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, 
  stan::require_std_vector_t<VecVar>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline void transform_inits_impl(VecVar& params_r__, VecI& params_i__,
                                   VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      local_scalar_t__ a;
      a = DUMMY_VAR__;
      
      a = in__.read<local_scalar_t__>();
      out__.write(a);
      local_scalar_t__ b;
      b = DUMMY_VAR__;
      
      b = in__.read<local_scalar_t__>();
      out__.write(b);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"a", "b", "y_tilde"};
    
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
      std::vector<size_t>{}, std::vector<size_t>{static_cast<size_t>(n)}};
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "a");
    param_names__.emplace_back(std::string() + "b");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= n; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "y_tilde" + '.' + std::to_string(sym1__));
        } 
      }
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "a");
    param_names__.emplace_back(std::string() + "b");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= n; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "y_tilde" + '.' + std::to_string(sym1__));
        } 
      }
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"a\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"b\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"y_tilde\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(n) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"}]");
    
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"a\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"b\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"y_tilde\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(n) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"}]");
    
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ = 
  (1 + 1);
      const size_t num_transformed = 0;
      const size_t num_gen_quantities = n;
      std::vector<double> vars_vec(num_params__
       + (emit_transformed_parameters * num_transformed)
       + (emit_generated_quantities * num_gen_quantities));
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        vars_vec.data(), vars_vec.size());
    }

    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ = 
  (1 + 1);
      const size_t num_transformed = 0;
      const size_t num_gen_quantities = n;
      vars.resize(num_params__
        + (emit_transformed_parameters * num_transformed)
        + (emit_generated_quantities * num_gen_quantities));
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }

    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }


    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits(context, params_i, params_r_vec, pstream);
      params_r = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        params_r_vec.data(), params_r_vec.size());
    }

  inline void transform_inits(const stan::io::var_context& context,
                              std::vector<int>& params_i,
                              std::vector<double>& vars,
                              std::ostream* pstream__ = nullptr) const {
     constexpr std::array<const char*, 2> names__{"a", "b"};
      const std::array<Eigen::Index, 2> constrain_param_sizes__{1, 1};
      const auto num_constrained_params__ = std::accumulate(
        constrain_param_sizes__.begin(), constrain_param_sizes__.end(), 0);
    
     std::vector<double> params_r_flat__(num_constrained_params__);
     Eigen::Index size_iter__ = 0;
     Eigen::Index flat_iter__ = 0;
     for (auto&& param_name__ : names__) {
       const auto param_vec__ = context.vals_r(param_name__);
       for (Eigen::Index i = 0; i < constrain_param_sizes__[size_iter__]; ++i) {
         params_r_flat__[flat_iter__] = param_vec__[i];
         ++flat_iter__;
       }
       ++size_iter__;
     }
     vars.resize(num_params_r__);
     transform_inits_impl(params_r_flat__, params_i, vars, pstream__);
    } // transform_inits() 
    
};
}
using stan_model = stan_file_model_namespace::stan_file_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

stan::math::profile_map& get_stan_profile_data() {
  return stan_file_model_namespace::profiles__;
}

#endif


