
// Code generated by stanc v2.28.1
#include <stan/model/model_header.hpp>
namespace toyMC_stan_model_namespace {

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
static constexpr std::array<const char*, 56> locations_array__ = 
{" (found before start of program)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 50, column 1 to column 31)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 51, column 1 to column 13)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 52, column 4 to column 16)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 54, column 4 to column 38)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 55, column 1 to column 17)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 56, column 4 to column 20)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 63, column 4 to column 122)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 62, column 13 to line 64, column 1)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 62, column 0 to line 64, column 1)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 67, column 4 to column 66)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 66, column 13 to line 68, column 1)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 66, column 0 to line 68, column 1)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 72, column 4 to column 23)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 73, column 4 to column 23)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 74, column 4 to column 27)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 75, column 4 to column 27)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 71, column 13 to line 76, column 1)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 71, column 0 to line 76, column 1)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 77, column 25 to column 26)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 77, column 22 to column 24)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 77, column 0 to column 28)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 78, column 29 to column 30)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 78, column 26 to column 28)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 78, column 0 to column 32)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 36, column 1 to column 16)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 37, column 8 to column 9)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 37, column 1 to column 21)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 38, column 8 to column 9)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 38, column 1 to column 21)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 40, column 4 to column 19)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 41, column 11 to column 12)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 41, column 4 to column 23)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 42, column 11 to column 12)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 42, column 4 to column 23)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 44, column 4 to column 28)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 4, column 4 to column 18)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 5, column 4 to column 18)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 7, column 4 to column 13)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 8, column 4 to column 13)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 11, column 8 to column 22)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 12, column 8 to column 22)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 10, column 17 to line 13, column 5)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 10, column 4 to line 13, column 5)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 15, column 4 to column 41)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 16, column 4 to column 56)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 18, column 4 to column 29)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 19, column 4 to column 20)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 20, column 4 to column 22)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 22, column 8 to column 32)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 23, column 8 to column 34)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 21, column 17 to line 24, column 5)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 21, column 4 to line 24, column 5)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 25, column 4 to column 72)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 31, column 4 to column 30)",
 " (in '/home/jakob/Documents/Studium/master_thesis/bayes/event_based_fit/toyMC_stan.stan', line 3, column 62 to line 32, column 1)"};


template <typename T0__, typename T1__, typename T2__, typename T3__,
typename T4__>
stan::promote_args_t<T0__, T1__, T2__, stan::value_type_t<T3__>,
stan::value_type_t<T4__>>
mylpdf(const T0__& phi, const T1__& pol, const T2__& sigma,
       const T3__& a_arg__, const T4__& b_arg__, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__,
          T1__,
          T2__,
          stan::value_type_t<T3__>,
          stan::value_type_t<T4__>>;
  int current_statement__ = 0; 
  const auto& a = to_ref(a_arg__);
  const auto& b = to_ref(b_arg__);
  static constexpr bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    Eigen::Matrix<local_scalar_t__, -1, 1> m_a;
    m_a = Eigen::Matrix<local_scalar_t__, -1, 1>(5);
    stan::math::fill(m_a, DUMMY_VAR__);
    
    Eigen::Matrix<local_scalar_t__, -1, 1> m_b;
    m_b = Eigen::Matrix<local_scalar_t__, -1, 1>(5);
    stan::math::fill(m_b, DUMMY_VAR__);
    
    current_statement__ = 38;
    assign(m_a, 0, "assigning variable m_a", index_uni(1));
    current_statement__ = 39;
    assign(m_b, 1, "assigning variable m_b", index_uni(1));
    current_statement__ = 43;
    for (int k = 1; k <= 4; ++k) {
      current_statement__ = 40;
      assign(m_a, rvalue(a, "a", index_uni(k)),
        "assigning variable m_a", index_uni((k + 1)));
      current_statement__ = 41;
      assign(m_b, rvalue(b, "b", index_uni(k)),
        "assigning variable m_b", index_uni((k + 1)));
    }
    local_scalar_t__ enmrtr;
    enmrtr = DUMMY_VAR__;
    
    current_statement__ = 44;
    enmrtr = (1 - (((0.5 * rvalue(m_a, "m_a", index_uni(3))) * pol) * sigma));
    local_scalar_t__ nmrtr;
    nmrtr = DUMMY_VAR__;
    
    current_statement__ = 45;
    nmrtr = (1 +
              ((pol * sigma) *
                stan::math::cos(
                  (((2 * (-45 - phi)) * stan::math::pi()) / 180.))));
    local_scalar_t__ angle;
    angle = DUMMY_VAR__;
    
    current_statement__ = 46;
    angle = ((phi * stan::math::pi()) / 180.);
    Eigen::Matrix<local_scalar_t__, -1, 1> sines;
    sines = Eigen::Matrix<local_scalar_t__, -1, 1>(5);
    stan::math::fill(sines, DUMMY_VAR__);
    
    Eigen::Matrix<local_scalar_t__, -1, 1> cosines;
    cosines = Eigen::Matrix<local_scalar_t__, -1, 1>(5);
    stan::math::fill(cosines, DUMMY_VAR__);
    
    current_statement__ = 52;
    for (int k = 0; k <= 4; ++k) {
      current_statement__ = 49;
      assign(sines, stan::math::sin((k * angle)),
        "assigning variable sines", index_uni((k + 1)));
      current_statement__ = 50;
      assign(cosines, stan::math::cos((k * angle)),
        "assigning variable cosines", index_uni((k + 1)));
    }
    local_scalar_t__ nnmrtr;
    nnmrtr = DUMMY_VAR__;
    
    current_statement__ = 53;
    nnmrtr = (nmrtr * (dot_product(sines, m_a) + dot_product(cosines, m_b)));
    current_statement__ = 54;
    return stan::math::log((nnmrtr / enmrtr));
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
  
}

struct mylpdf_functor__ {
template <typename T0__, typename T1__, typename T2__, typename T3__,
typename T4__>
stan::promote_args_t<T0__, T1__, T2__, stan::value_type_t<T3__>,
stan::value_type_t<T4__>>
operator()(const T0__& phi, const T1__& pol, const T2__& sigma,
           const T3__& a, const T4__& b, std::ostream* pstream__)  const 
{
return mylpdf(phi, pol, sigma, a, b, pstream__);
}
};

class toyMC_stan_model final : public model_base_crtp<toyMC_stan_model> {

 private:
  int N;
  Eigen::Matrix<double, -1, 1> phi_prmpt__;
  Eigen::Matrix<double, -1, 1> pol_prmpt__;
  int M;
  Eigen::Matrix<double, -1, 1> phi_side__;
  Eigen::Matrix<double, -1, 1> pol_side__;
  double f; 
  Eigen::Map<Eigen::Matrix<double, -1, 1>> phi_prmpt{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double, -1, 1>> pol_prmpt{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double, -1, 1>> phi_side{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double, -1, 1>> pol_side{nullptr, 0};
 
 public:
  ~toyMC_stan_model() { }
  
  inline std::string model_name() const final { return "toyMC_stan_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.28.1", "stancflags = "};
  }
  
  
  toyMC_stan_model(stan::io::var_context& context__,
                   unsigned int random_seed__ = 0,
                   std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static constexpr const char* function__ = "toyMC_stan_model_namespace::toyMC_stan_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 25;
      context__.validate_dims("data initialization","N","int",
           std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      
      current_statement__ = 25;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 25;
      check_greater_or_equal(function__, "N", N, 0);
      current_statement__ = 26;
      validate_non_negative_index("phi_prmpt", "N", N);
      current_statement__ = 27;
      context__.validate_dims("data initialization","phi_prmpt","double",
           std::vector<size_t>{static_cast<size_t>(N)});
      phi_prmpt__ = Eigen::Matrix<double, -1, 1>(N);
      new (&phi_prmpt) Eigen::Map<Eigen::Matrix<double, -1, 1>>(phi_prmpt__.data(), N);
      
      
      {
        std::vector<local_scalar_t__> phi_prmpt_flat__;
        current_statement__ = 27;
        phi_prmpt_flat__ = context__.vals_r("phi_prmpt");
        current_statement__ = 27;
        pos__ = 1;
        current_statement__ = 27;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 27;
          assign(phi_prmpt, phi_prmpt_flat__[(pos__ - 1)],
            "assigning variable phi_prmpt", index_uni(sym1__));
          current_statement__ = 27;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 28;
      validate_non_negative_index("pol_prmpt", "N", N);
      current_statement__ = 29;
      context__.validate_dims("data initialization","pol_prmpt","double",
           std::vector<size_t>{static_cast<size_t>(N)});
      pol_prmpt__ = Eigen::Matrix<double, -1, 1>(N);
      new (&pol_prmpt) Eigen::Map<Eigen::Matrix<double, -1, 1>>(pol_prmpt__.data(), N);
      
      
      {
        std::vector<local_scalar_t__> pol_prmpt_flat__;
        current_statement__ = 29;
        pol_prmpt_flat__ = context__.vals_r("pol_prmpt");
        current_statement__ = 29;
        pos__ = 1;
        current_statement__ = 29;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 29;
          assign(pol_prmpt, pol_prmpt_flat__[(pos__ - 1)],
            "assigning variable pol_prmpt", index_uni(sym1__));
          current_statement__ = 29;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 30;
      context__.validate_dims("data initialization","M","int",
           std::vector<size_t>{});
      M = std::numeric_limits<int>::min();
      
      current_statement__ = 30;
      M = context__.vals_i("M")[(1 - 1)];
      current_statement__ = 30;
      check_greater_or_equal(function__, "M", M, 0);
      current_statement__ = 31;
      validate_non_negative_index("phi_side", "M", M);
      current_statement__ = 32;
      context__.validate_dims("data initialization","phi_side","double",
           std::vector<size_t>{static_cast<size_t>(M)});
      phi_side__ = Eigen::Matrix<double, -1, 1>(M);
      new (&phi_side) Eigen::Map<Eigen::Matrix<double, -1, 1>>(phi_side__.data(), M);
      
      
      {
        std::vector<local_scalar_t__> phi_side_flat__;
        current_statement__ = 32;
        phi_side_flat__ = context__.vals_r("phi_side");
        current_statement__ = 32;
        pos__ = 1;
        current_statement__ = 32;
        for (int sym1__ = 1; sym1__ <= M; ++sym1__) {
          current_statement__ = 32;
          assign(phi_side, phi_side_flat__[(pos__ - 1)],
            "assigning variable phi_side", index_uni(sym1__));
          current_statement__ = 32;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 33;
      validate_non_negative_index("pol_side", "M", M);
      current_statement__ = 34;
      context__.validate_dims("data initialization","pol_side","double",
           std::vector<size_t>{static_cast<size_t>(M)});
      pol_side__ = Eigen::Matrix<double, -1, 1>(M);
      new (&pol_side) Eigen::Map<Eigen::Matrix<double, -1, 1>>(pol_side__.data(), M);
      
      
      {
        std::vector<local_scalar_t__> pol_side_flat__;
        current_statement__ = 34;
        pol_side_flat__ = context__.vals_r("pol_side");
        current_statement__ = 34;
        pos__ = 1;
        current_statement__ = 34;
        for (int sym1__ = 1; sym1__ <= M; ++sym1__) {
          current_statement__ = 34;
          assign(pol_side, pol_side_flat__[(pos__ - 1)],
            "assigning variable pol_side", index_uni(sym1__));
          current_statement__ = 34;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 35;
      context__.validate_dims("data initialization","f","double",
           std::vector<size_t>{});
      f = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 35;
      f = context__.vals_r("f")[(1 - 1)];
      current_statement__ = 35;
      check_greater_or_equal(function__, "f", f, 0);
      current_statement__ = 35;
      check_less_or_equal(function__, "f", f, 1);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1 + 4 + 4 + 1 + 4 + 4;
    
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
    static constexpr const char* function__ = "toyMC_stan_model_namespace::log_prob";
    (void) function__;  // suppress unused var warning
    
    try {
      local_scalar_t__ sigma;
      sigma = DUMMY_VAR__;
      
      current_statement__ = 1;
      sigma = in__.template read_constrain_lub<local_scalar_t__, jacobian__>(
                -1, 1, lp__);
      Eigen::Matrix<local_scalar_t__, -1, 1> a;
      a = Eigen::Matrix<local_scalar_t__, -1, 1>(4);
      stan::math::fill(a, DUMMY_VAR__);
      
      current_statement__ = 2;
      a = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(4);
      Eigen::Matrix<local_scalar_t__, -1, 1> b;
      b = Eigen::Matrix<local_scalar_t__, -1, 1>(4);
      stan::math::fill(b, DUMMY_VAR__);
      
      current_statement__ = 3;
      b = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(4);
      local_scalar_t__ sigma_bkg;
      sigma_bkg = DUMMY_VAR__;
      
      current_statement__ = 4;
      sigma_bkg = in__.template read_constrain_lub<local_scalar_t__, jacobian__>(
                    -1, 1, lp__);
      Eigen::Matrix<local_scalar_t__, -1, 1> a_bkg;
      a_bkg = Eigen::Matrix<local_scalar_t__, -1, 1>(4);
      stan::math::fill(a_bkg, DUMMY_VAR__);
      
      current_statement__ = 5;
      a_bkg = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(4);
      Eigen::Matrix<local_scalar_t__, -1, 1> b_bkg;
      b_bkg = Eigen::Matrix<local_scalar_t__, -1, 1>(4);
      stan::math::fill(b_bkg, DUMMY_VAR__);
      
      current_statement__ = 6;
      b_bkg = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(4);
      {
        current_statement__ = 9;
        for (int k = 1; k <= N; ++k) {
          current_statement__ = 7;
          lp_accum__.add(
            ((f *
               mylpdf(rvalue(phi_prmpt, "phi_prmpt", index_uni(k)),
                 rvalue(pol_prmpt, "pol_prmpt", index_uni(k)), sigma, a,
                 b, pstream__)) +
              ((1 - f) *
                mylpdf(rvalue(phi_prmpt, "phi_prmpt", index_uni(k)),
                  rvalue(pol_prmpt, "pol_prmpt", index_uni(k)), sigma_bkg,
                  a_bkg, b_bkg, pstream__))));
        }
        current_statement__ = 12;
        for (int k = 1; k <= M; ++k) {
          current_statement__ = 10;
          lp_accum__.add(
            mylpdf(rvalue(phi_side, "phi_side", index_uni(k)),
              rvalue(pol_side, "pol_side", index_uni(k)), sigma_bkg, a_bkg,
              b_bkg, pstream__));
        }
        current_statement__ = 18;
        for (int k = 1; k <= 4; ++k) {
          current_statement__ = 13;
          lp_accum__.add(
            normal_lpdf<propto__>(rvalue(a, "a", index_uni(k)), 0, 1));
          current_statement__ = 14;
          lp_accum__.add(
            normal_lpdf<propto__>(rvalue(b, "b", index_uni(k)), 0, 1));
          current_statement__ = 15;
          lp_accum__.add(
            normal_lpdf<propto__>(rvalue(a_bkg, "a_bkg", index_uni(k)), 0, 1));
          current_statement__ = 16;
          lp_accum__.add(
            normal_lpdf<propto__>(rvalue(b_bkg, "b_bkg", index_uni(k)), 0, 1));
        }
        current_statement__ = 20;
        if (logical_lt(sigma, -1)) {
          current_statement__ = 20;
          lp_accum__.add(stan::math::negative_infinity());
        } else {
          current_statement__ = 19;
          if (logical_gt(sigma, 1)) {
            current_statement__ = 19;
            lp_accum__.add(stan::math::negative_infinity());
          } else {
            current_statement__ = 19;
            lp_accum__.add(
              -log_diff_exp(normal_cdf_log(1, 0, 1),
                 normal_cdf_log(-1, 0, 1)));
          }
        }
        current_statement__ = 21;
        lp_accum__.add(normal_lpdf<propto__>(sigma, 0, 1));
        current_statement__ = 23;
        if (logical_lt(sigma_bkg, -1)) {
          current_statement__ = 23;
          lp_accum__.add(stan::math::negative_infinity());
        } else {
          current_statement__ = 22;
          if (logical_gt(sigma_bkg, 1)) {
            current_statement__ = 22;
            lp_accum__.add(stan::math::negative_infinity());
          } else {
            current_statement__ = 22;
            lp_accum__.add(
              -log_diff_exp(normal_cdf_log(1, 0, 1),
                 normal_cdf_log(-1, 0, 1)));
          }
        }
        current_statement__ = 24;
        lp_accum__.add(normal_lpdf<propto__>(sigma_bkg, 0, 1));
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
    static constexpr const char* function__ = "toyMC_stan_model_namespace::write_array";
    (void) function__;  // suppress unused var warning
    
    try {
      double sigma;
      sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      sigma = in__.template read_constrain_lub<local_scalar_t__, jacobian__>(
                -1, 1, lp__);
      Eigen::Matrix<double, -1, 1> a;
      a = Eigen::Matrix<double, -1, 1>(4);
      stan::math::fill(a, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 2;
      a = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(4);
      Eigen::Matrix<double, -1, 1> b;
      b = Eigen::Matrix<double, -1, 1>(4);
      stan::math::fill(b, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 3;
      b = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(4);
      double sigma_bkg;
      sigma_bkg = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      sigma_bkg = in__.template read_constrain_lub<local_scalar_t__, jacobian__>(
                    -1, 1, lp__);
      Eigen::Matrix<double, -1, 1> a_bkg;
      a_bkg = Eigen::Matrix<double, -1, 1>(4);
      stan::math::fill(a_bkg, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 5;
      a_bkg = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(4);
      Eigen::Matrix<double, -1, 1> b_bkg;
      b_bkg = Eigen::Matrix<double, -1, 1>(4);
      stan::math::fill(b_bkg, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 6;
      b_bkg = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(4);
      out__.write(sigma);
      out__.write(a);
      out__.write(b);
      out__.write(sigma_bkg);
      out__.write(a_bkg);
      out__.write(b_bkg);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
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
      local_scalar_t__ sigma;
      sigma = DUMMY_VAR__;
      
      sigma = in__.read<local_scalar_t__>();
      out__.write_free_lub(-1, 1, sigma);
      Eigen::Matrix<local_scalar_t__, -1, 1> a;
      a = Eigen::Matrix<local_scalar_t__, -1, 1>(4);
      stan::math::fill(a, DUMMY_VAR__);
      
      for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
        assign(a, in__.read<local_scalar_t__>(),
          "assigning variable a", index_uni(sym1__));
      }
      out__.write(a);
      Eigen::Matrix<local_scalar_t__, -1, 1> b;
      b = Eigen::Matrix<local_scalar_t__, -1, 1>(4);
      stan::math::fill(b, DUMMY_VAR__);
      
      for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
        assign(b, in__.read<local_scalar_t__>(),
          "assigning variable b", index_uni(sym1__));
      }
      out__.write(b);
      local_scalar_t__ sigma_bkg;
      sigma_bkg = DUMMY_VAR__;
      
      sigma_bkg = in__.read<local_scalar_t__>();
      out__.write_free_lub(-1, 1, sigma_bkg);
      Eigen::Matrix<local_scalar_t__, -1, 1> a_bkg;
      a_bkg = Eigen::Matrix<local_scalar_t__, -1, 1>(4);
      stan::math::fill(a_bkg, DUMMY_VAR__);
      
      for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
        assign(a_bkg, in__.read<local_scalar_t__>(),
          "assigning variable a_bkg", index_uni(sym1__));
      }
      out__.write(a_bkg);
      Eigen::Matrix<local_scalar_t__, -1, 1> b_bkg;
      b_bkg = Eigen::Matrix<local_scalar_t__, -1, 1>(4);
      stan::math::fill(b_bkg, DUMMY_VAR__);
      
      for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
        assign(b_bkg, in__.read<local_scalar_t__>(),
          "assigning variable b_bkg", index_uni(sym1__));
      }
      out__.write(b_bkg);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"sigma", "a", "b", "sigma_bkg",
      "a_bkg", "b_bkg"};
    
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
      std::vector<size_t>{static_cast<size_t>(4)},
      std::vector<size_t>{static_cast<size_t>(4)}, std::vector<size_t>{
      }, std::vector<size_t>{static_cast<size_t>(4)},
      std::vector<size_t>{static_cast<size_t>(4)}};
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "sigma");
    for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "a" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "b" + '.' + std::to_string(sym1__));
      } 
    }
    param_names__.emplace_back(std::string() + "sigma_bkg");
    for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "a_bkg" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "b_bkg" + '.' + std::to_string(sym1__));
      } 
    }
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "sigma");
    for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "a" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "b" + '.' + std::to_string(sym1__));
      } 
    }
    param_names__.emplace_back(std::string() + "sigma_bkg");
    for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "a_bkg" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= 4; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "b_bkg" + '.' + std::to_string(sym1__));
      } 
    }
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"a\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(4) + "},\"block\":\"parameters\"},{\"name\":\"b\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(4) + "},\"block\":\"parameters\"},{\"name\":\"sigma_bkg\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"a_bkg\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(4) + "},\"block\":\"parameters\"},{\"name\":\"b_bkg\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(4) + "},\"block\":\"parameters\"}]");
    
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"a\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(4) + "},\"block\":\"parameters\"},{\"name\":\"b\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(4) + "},\"block\":\"parameters\"},{\"name\":\"sigma_bkg\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"a_bkg\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(4) + "},\"block\":\"parameters\"},{\"name\":\"b_bkg\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(4) + "},\"block\":\"parameters\"}]");
    
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
  (((((1 + 4) + 4) + 1) + 4) + 4);
      const size_t num_transformed = 0;
      const size_t num_gen_quantities = 0;
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
  (((((1 + 4) + 4) + 1) + 4) + 4);
      const size_t num_transformed = 0;
      const size_t num_gen_quantities = 0;
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
     constexpr std::array<const char*, 6> names__{"sigma", "a", "b",
      "sigma_bkg", "a_bkg", "b_bkg"};
      const std::array<Eigen::Index, 6> constrain_param_sizes__{1, 4, 
       4, 1, 4, 4};
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
using stan_model = toyMC_stan_model_namespace::toyMC_stan_model;

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
  return toyMC_stan_model_namespace::profiles__;
}

#endif


