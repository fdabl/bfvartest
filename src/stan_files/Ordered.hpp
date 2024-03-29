/*
    bfvartest is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bfvartest is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bfvartest.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0

#include <stan/model/model_header.hpp>

namespace model_Ordered_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_Ordered");
    reader.add_event(56, 54, "end", "model_Ordered");
    return reader;
}

class model_Ordered
  : public stan::model::model_base_crtp<model_Ordered> {
private:
        int k;
        double alpha;
        int nr_equal;
        int nr_ordered;
        std::vector<int> index_vector;
        vector_d s2;
        vector_d N;
        int priors_only;
        vector_d n;
        vector_d b;
        double nplus;
public:
    model_Ordered(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_Ordered(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;

        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "model_Ordered_namespace::model_Ordered";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {
            // initialize data block variables from context__
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "k", "int", context__.to_vec());
            k = int(0);
            vals_i__ = context__.vals_i("k");
            pos__ = 0;
            k = vals_i__[pos__++];

            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "alpha", "double", context__.to_vec());
            alpha = double(0);
            vals_r__ = context__.vals_r("alpha");
            pos__ = 0;
            alpha = vals_r__[pos__++];

            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "nr_equal", "int", context__.to_vec());
            nr_equal = int(0);
            vals_i__ = context__.vals_i("nr_equal");
            pos__ = 0;
            nr_equal = vals_i__[pos__++];

            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "nr_ordered", "int", context__.to_vec());
            nr_ordered = int(0);
            vals_i__ = context__.vals_i("nr_ordered");
            pos__ = 0;
            nr_ordered = vals_i__[pos__++];

            current_statement_begin__ = 7;
            validate_non_negative_index("index_vector", "k", k);
            context__.validate_dims("data initialization", "index_vector", "int", context__.to_vec(k));
            index_vector = std::vector<int>(k, int(0));
            vals_i__ = context__.vals_i("index_vector");
            pos__ = 0;
            size_t index_vector_k_0_max__ = k;
            for (size_t k_0__ = 0; k_0__ < index_vector_k_0_max__; ++k_0__) {
                index_vector[k_0__] = vals_i__[pos__++];
            }

            current_statement_begin__ = 8;
            validate_non_negative_index("s2", "k", k);
            context__.validate_dims("data initialization", "s2", "vector_d", context__.to_vec(k));
            s2 = Eigen::Matrix<double, Eigen::Dynamic, 1>(k);
            vals_r__ = context__.vals_r("s2");
            pos__ = 0;
            size_t s2_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < s2_j_1_max__; ++j_1__) {
                s2(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "s2", s2, 0);

            current_statement_begin__ = 9;
            validate_non_negative_index("N", "k", k);
            context__.validate_dims("data initialization", "N", "vector_d", context__.to_vec(k));
            N = Eigen::Matrix<double, Eigen::Dynamic, 1>(k);
            vals_r__ = context__.vals_r("N");
            pos__ = 0;
            size_t N_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < N_j_1_max__; ++j_1__) {
                N(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "N", N, 0);

            current_statement_begin__ = 10;
            context__.validate_dims("data initialization", "priors_only", "int", context__.to_vec());
            priors_only = int(0);
            vals_i__ = context__.vals_i("priors_only");
            pos__ = 0;
            priors_only = vals_i__[pos__++];


            // initialize transformed data variables
            current_statement_begin__ = 14;
            validate_non_negative_index("n", "k", k);
            n = Eigen::Matrix<double, Eigen::Dynamic, 1>(k);
            stan::math::fill(n, DUMMY_VAR__);

            current_statement_begin__ = 15;
            validate_non_negative_index("b", "k", k);
            b = Eigen::Matrix<double, Eigen::Dynamic, 1>(k);
            stan::math::fill(b, DUMMY_VAR__);

            current_statement_begin__ = 16;
            nplus = double(0);
            stan::math::fill(nplus, DUMMY_VAR__);

            // execute transformed data statements
            current_statement_begin__ = 19;
            stan::math::assign(n, divide(subtract(N, 1.0), 2.0));
            current_statement_begin__ = 20;
            stan::math::assign(b, elt_multiply(s2, N));
            current_statement_begin__ = 21;
            stan::math::assign(nplus, sum(n));

            // validate transformed data
            current_statement_begin__ = 14;
            check_greater_or_equal(function__, "n", n, 0);

            current_statement_begin__ = 15;
            check_greater_or_equal(function__, "b", b, 0);


            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 26;
            num_params_r__ += 1;
            current_statement_begin__ = 27;
            validate_non_negative_index("lambda_unconstrained", "(k - nr_equal)", (k - nr_equal));
            num_params_r__ += (k - nr_equal);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_Ordered() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        current_statement_begin__ = 26;
        if (!(context__.contains_r("tau")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable tau missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("tau");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "tau", "double", context__.to_vec());
        double tau(0);
        tau = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, tau);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable tau: ") + e.what()), current_statement_begin__, prog_reader__());
        }

        current_statement_begin__ = 27;
        if (!(context__.contains_r("lambda_unconstrained")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable lambda_unconstrained missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("lambda_unconstrained");
        pos__ = 0U;
        validate_non_negative_index("lambda_unconstrained", "(k - nr_equal)", (k - nr_equal));
        context__.validate_dims("parameter initialization", "lambda_unconstrained", "vector_d", context__.to_vec((k - nr_equal)));
        Eigen::Matrix<double, Eigen::Dynamic, 1> lambda_unconstrained((k - nr_equal));
        size_t lambda_unconstrained_j_1_max__ = (k - nr_equal);
        for (size_t j_1__ = 0; j_1__ < lambda_unconstrained_j_1_max__; ++j_1__) {
            lambda_unconstrained(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.positive_ordered_unconstrain(lambda_unconstrained);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable lambda_unconstrained: ") + e.what()), current_statement_begin__, prog_reader__());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        typedef T__ local_scalar_t__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);

            // model parameters
            current_statement_begin__ = 26;
            local_scalar_t__ tau;
            (void) tau;  // dummy to suppress unused var warning
            if (jacobian__)
                tau = in__.scalar_lb_constrain(0, lp__);
            else
                tau = in__.scalar_lb_constrain(0);

            current_statement_begin__ = 27;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lambda_unconstrained;
            (void) lambda_unconstrained;  // dummy to suppress unused var warning
            if (jacobian__)
                lambda_unconstrained = in__.positive_ordered_constrain((k - nr_equal), lp__);
            else
                lambda_unconstrained = in__.positive_ordered_constrain((k - nr_equal));

            // transformed parameters
            current_statement_begin__ = 31;
            validate_non_negative_index("sds", "k", k);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> sds(k);
            stan::math::initialize(sds, DUMMY_VAR__);
            stan::math::fill(sds, DUMMY_VAR__);

            current_statement_begin__ = 32;
            validate_non_negative_index("rho", "k", k);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> rho(k);
            stan::math::initialize(rho, DUMMY_VAR__);
            stan::math::fill(rho, DUMMY_VAR__);

            current_statement_begin__ = 33;
            validate_non_negative_index("lambda", "k", k);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lambda(k);
            stan::math::initialize(lambda, DUMMY_VAR__);
            stan::math::fill(lambda, DUMMY_VAR__);

            current_statement_begin__ = 34;
            validate_non_negative_index("prec", "k", k);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> prec(k);
            stan::math::initialize(prec, DUMMY_VAR__);
            stan::math::fill(prec, DUMMY_VAR__);

            // transformed parameters block statements
            current_statement_begin__ = 35;
            for (int i = 1; i <= k; ++i) {
                current_statement_begin__ = 35;
                stan::model::assign(lambda, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            get_base1(lambda_unconstrained, get_base1(index_vector, i, "index_vector", 1), "lambda_unconstrained", 1), 
                            "assigning variable lambda");
            }
            current_statement_begin__ = 37;
            stan::math::assign(rho, divide(lambda, sum(lambda)));
            current_statement_begin__ = 38;
            stan::math::assign(prec, multiply(multiply(rho, tau), k));
            current_statement_begin__ = 39;
            stan::math::assign(sds, elt_divide(1.0, stan::math::sqrt(multiply(multiply(rho, tau), k))));

            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            current_statement_begin__ = 31;
            size_t sds_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < sds_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(sds(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: sds" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable sds: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 32;
            size_t rho_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < rho_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(rho(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: rho" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable rho: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            stan::math::check_simplex(function__, "rho", rho);

            current_statement_begin__ = 33;
            size_t lambda_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < lambda_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(lambda(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: lambda" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable lambda: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 34;
            size_t prec_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < prec_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(prec(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: prec" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable prec: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }

            // model body

            current_statement_begin__ = 43;
            lp_accum__.add(-(stan::math::log(tau)));
            current_statement_begin__ = 44;
            lp_accum__.add(gamma_lpdf<propto__>(lambda_unconstrained, alpha, 1));
            current_statement_begin__ = 47;
            lp_accum__.add(stan::math::lgamma((nr_ordered + 1)));
            current_statement_begin__ = 49;
            if (as_bool(logical_negation(logical_eq(priors_only, 1)))) {

                current_statement_begin__ = 51;
                lp_accum__.add(dot_product(n, stan::math::log(prec)));
                current_statement_begin__ = 52;
                lp_accum__.add((-(0.5) * dot_product(prec, b)));
            }

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("tau");
        names__.push_back("lambda_unconstrained");
        names__.push_back("sds");
        names__.push_back("rho");
        names__.push_back("lambda");
        names__.push_back("prec");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((k - nr_equal));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(k);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(k);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(k);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(k);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;

        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_Ordered_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning

        // read-transform, write parameters
        double tau = in__.scalar_lb_constrain(0);
        vars__.push_back(tau);

        Eigen::Matrix<double, Eigen::Dynamic, 1> lambda_unconstrained = in__.positive_ordered_constrain((k - nr_equal));
        size_t lambda_unconstrained_j_1_max__ = (k - nr_equal);
        for (size_t j_1__ = 0; j_1__ < lambda_unconstrained_j_1_max__; ++j_1__) {
            vars__.push_back(lambda_unconstrained(j_1__));
        }

        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        if (!include_tparams__ && !include_gqs__) return;

        try {
            // declare and define transformed parameters
            current_statement_begin__ = 31;
            validate_non_negative_index("sds", "k", k);
            Eigen::Matrix<double, Eigen::Dynamic, 1> sds(k);
            stan::math::initialize(sds, DUMMY_VAR__);
            stan::math::fill(sds, DUMMY_VAR__);

            current_statement_begin__ = 32;
            validate_non_negative_index("rho", "k", k);
            Eigen::Matrix<double, Eigen::Dynamic, 1> rho(k);
            stan::math::initialize(rho, DUMMY_VAR__);
            stan::math::fill(rho, DUMMY_VAR__);

            current_statement_begin__ = 33;
            validate_non_negative_index("lambda", "k", k);
            Eigen::Matrix<double, Eigen::Dynamic, 1> lambda(k);
            stan::math::initialize(lambda, DUMMY_VAR__);
            stan::math::fill(lambda, DUMMY_VAR__);

            current_statement_begin__ = 34;
            validate_non_negative_index("prec", "k", k);
            Eigen::Matrix<double, Eigen::Dynamic, 1> prec(k);
            stan::math::initialize(prec, DUMMY_VAR__);
            stan::math::fill(prec, DUMMY_VAR__);

            // do transformed parameters statements
            current_statement_begin__ = 35;
            for (int i = 1; i <= k; ++i) {
                current_statement_begin__ = 35;
                stan::model::assign(lambda, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            get_base1(lambda_unconstrained, get_base1(index_vector, i, "index_vector", 1), "lambda_unconstrained", 1), 
                            "assigning variable lambda");
            }
            current_statement_begin__ = 37;
            stan::math::assign(rho, divide(lambda, sum(lambda)));
            current_statement_begin__ = 38;
            stan::math::assign(prec, multiply(multiply(rho, tau), k));
            current_statement_begin__ = 39;
            stan::math::assign(sds, elt_divide(1.0, stan::math::sqrt(multiply(multiply(rho, tau), k))));

            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            current_statement_begin__ = 32;
            stan::math::check_simplex(function__, "rho", rho);

            // write transformed parameters
            if (include_tparams__) {
                size_t sds_j_1_max__ = k;
                for (size_t j_1__ = 0; j_1__ < sds_j_1_max__; ++j_1__) {
                    vars__.push_back(sds(j_1__));
                }
                size_t rho_j_1_max__ = k;
                for (size_t j_1__ = 0; j_1__ < rho_j_1_max__; ++j_1__) {
                    vars__.push_back(rho(j_1__));
                }
                size_t lambda_j_1_max__ = k;
                for (size_t j_1__ = 0; j_1__ < lambda_j_1_max__; ++j_1__) {
                    vars__.push_back(lambda(j_1__));
                }
                size_t prec_j_1_max__ = k;
                for (size_t j_1__ = 0; j_1__ < prec_j_1_max__; ++j_1__) {
                    vars__.push_back(prec(j_1__));
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    std::string model_name() const {
        return "model_Ordered";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        size_t lambda_unconstrained_j_1_max__ = (k - nr_equal);
        for (size_t j_1__ = 0; j_1__ < lambda_unconstrained_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda_unconstrained" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            size_t sds_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < sds_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "sds" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t rho_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < rho_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "rho" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t lambda_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < lambda_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "lambda" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t prec_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < prec_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "prec" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        size_t lambda_unconstrained_j_1_max__ = (k - nr_equal);
        for (size_t j_1__ = 0; j_1__ < lambda_unconstrained_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda_unconstrained" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            size_t sds_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < sds_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "sds" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t rho_j_1_max__ = (k - 1);
            for (size_t j_1__ = 0; j_1__ < rho_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "rho" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t lambda_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < lambda_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "lambda" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t prec_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < prec_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "prec" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__) return;
    }

}; // model

}  // namespace

typedef model_Ordered_namespace::model_Ordered stan_model;

#ifndef USING_R

stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

#endif


#endif
