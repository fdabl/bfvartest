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
// Code generated by Stan version 2.18.0

#include <stan/model/model_header.hpp>

namespace model_Mixed_namespace {

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
    reader.add_event(0, 0, "start", "model_Mixed");
    reader.add_event(56, 54, "end", "model_Mixed");
    return reader;
}

#include <meta_header.hpp>
 class model_Mixed : public prob_grad {
private:
    int k;
    double alpha;
    int nr_equal;
    vector<int> index_vector;
    vector_d s2;
    vector_d N;
    int priors_only;
    vector_d n;
    vector_d b;
    double nplus;
public:
    model_Mixed(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_Mixed(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
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

        static const char* function__ = "model_Mixed_namespace::model_Mixed";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "k", "int", context__.to_vec());
            k = int(0);
            vals_i__ = context__.vals_i("k");
            pos__ = 0;
            k = vals_i__[pos__++];
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "alpha", "double", context__.to_vec());
            alpha = double(0);
            vals_r__ = context__.vals_r("alpha");
            pos__ = 0;
            alpha = vals_r__[pos__++];
            current_statement_begin__ = 7;
            context__.validate_dims("data initialization", "nr_equal", "int", context__.to_vec());
            nr_equal = int(0);
            vals_i__ = context__.vals_i("nr_equal");
            pos__ = 0;
            nr_equal = vals_i__[pos__++];
            current_statement_begin__ = 8;
            validate_non_negative_index("index_vector", "k", k);
            context__.validate_dims("data initialization", "index_vector", "int", context__.to_vec(k));
            validate_non_negative_index("index_vector", "k", k);
            index_vector = std::vector<int>(k,int(0));
            vals_i__ = context__.vals_i("index_vector");
            pos__ = 0;
            size_t index_vector_limit_0__ = k;
            for (size_t i_0__ = 0; i_0__ < index_vector_limit_0__; ++i_0__) {
                index_vector[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("s2", "k", k);
            context__.validate_dims("data initialization", "s2", "vector_d", context__.to_vec(k));
            validate_non_negative_index("s2", "k", k);
            s2 = vector_d(static_cast<Eigen::VectorXd::Index>(k));
            vals_r__ = context__.vals_r("s2");
            pos__ = 0;
            size_t s2_i_vec_lim__ = k;
            for (size_t i_vec__ = 0; i_vec__ < s2_i_vec_lim__; ++i_vec__) {
                s2[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 10;
            validate_non_negative_index("N", "k", k);
            context__.validate_dims("data initialization", "N", "vector_d", context__.to_vec(k));
            validate_non_negative_index("N", "k", k);
            N = vector_d(static_cast<Eigen::VectorXd::Index>(k));
            vals_r__ = context__.vals_r("N");
            pos__ = 0;
            size_t N_i_vec_lim__ = k;
            for (size_t i_vec__ = 0; i_vec__ < N_i_vec_lim__; ++i_vec__) {
                N[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 11;
            context__.validate_dims("data initialization", "priors_only", "int", context__.to_vec());
            priors_only = int(0);
            vals_i__ = context__.vals_i("priors_only");
            pos__ = 0;
            priors_only = vals_i__[pos__++];

            // validate, data variables
            current_statement_begin__ = 5;
            current_statement_begin__ = 6;
            current_statement_begin__ = 7;
            current_statement_begin__ = 8;
            current_statement_begin__ = 9;
            check_greater_or_equal(function__,"s2",s2,0);
            current_statement_begin__ = 10;
            check_greater_or_equal(function__,"N",N,0);
            current_statement_begin__ = 11;
            // initialize data variables
            current_statement_begin__ = 15;
            validate_non_negative_index("n", "k", k);
            n = vector_d(static_cast<Eigen::VectorXd::Index>(k));
            stan::math::fill(n,DUMMY_VAR__);
            current_statement_begin__ = 16;
            validate_non_negative_index("b", "k", k);
            b = vector_d(static_cast<Eigen::VectorXd::Index>(k));
            stan::math::fill(b,DUMMY_VAR__);
            current_statement_begin__ = 17;
            nplus = double(0);
            stan::math::fill(nplus,DUMMY_VAR__);

            current_statement_begin__ = 20;
            stan::math::assign(n, divide(subtract(N,1.0),2.0));
            current_statement_begin__ = 21;
            stan::math::assign(b, elt_multiply(s2,N));
            current_statement_begin__ = 22;
            stan::math::assign(nplus, sum(n));

            // validate transformed data
            current_statement_begin__ = 15;
            check_greater_or_equal(function__,"n",n,0);
            current_statement_begin__ = 16;
            check_greater_or_equal(function__,"b",b,0);
            current_statement_begin__ = 17;

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 28;
            ++num_params_r__;
            current_statement_begin__ = 29;
            validate_non_negative_index("lambda_unconstrained", "(k - nr_equal)", (k - nr_equal));
            num_params_r__ += (k - nr_equal);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_Mixed() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("tau")))
            throw std::runtime_error("variable tau missing");
        vals_r__ = context__.vals_r("tau");
        pos__ = 0U;
        context__.validate_dims("initialization", "tau", "double", context__.to_vec());
        double tau(0);
        tau = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,tau);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable tau: ") + e.what());
        }

        if (!(context__.contains_r("lambda_unconstrained")))
            throw std::runtime_error("variable lambda_unconstrained missing");
        vals_r__ = context__.vals_r("lambda_unconstrained");
        pos__ = 0U;
        validate_non_negative_index("lambda_unconstrained", "(k - nr_equal)", (k - nr_equal));
        context__.validate_dims("initialization", "lambda_unconstrained", "double", context__.to_vec((k - nr_equal)));
        std::vector<double> lambda_unconstrained((k - nr_equal),double(0));
        for (int i0__ = 0U; i0__ < (k - nr_equal); ++i0__)
            lambda_unconstrained[i0__] = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < (k - nr_equal); ++i0__)
            try {
            writer__.scalar_lb_unconstrain(0,lambda_unconstrained[i0__]);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable lambda_unconstrained: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        typedef T__ local_scalar_t__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);

            local_scalar_t__ tau;
            (void) tau;  // dummy to suppress unused var warning
            if (jacobian__)
                tau = in__.scalar_lb_constrain(0,lp__);
            else
                tau = in__.scalar_lb_constrain(0);

            vector<local_scalar_t__> lambda_unconstrained;
            size_t dim_lambda_unconstrained_0__ = (k - nr_equal);
            lambda_unconstrained.reserve(dim_lambda_unconstrained_0__);
            for (size_t k_0__ = 0; k_0__ < dim_lambda_unconstrained_0__; ++k_0__) {
                if (jacobian__)
                    lambda_unconstrained.push_back(in__.scalar_lb_constrain(0,lp__));
                else
                    lambda_unconstrained.push_back(in__.scalar_lb_constrain(0));
            }


            // transformed parameters
            current_statement_begin__ = 33;
            validate_non_negative_index("sds", "k", k);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  sds(static_cast<Eigen::VectorXd::Index>(k));
            (void) sds;  // dummy to suppress unused var warning

            stan::math::initialize(sds, DUMMY_VAR__);
            stan::math::fill(sds,DUMMY_VAR__);
            current_statement_begin__ = 34;
            validate_non_negative_index("rho", "k", k);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  rho(static_cast<Eigen::VectorXd::Index>(k));
            (void) rho;  // dummy to suppress unused var warning

            stan::math::initialize(rho, DUMMY_VAR__);
            stan::math::fill(rho,DUMMY_VAR__);
            current_statement_begin__ = 35;
            validate_non_negative_index("lambda", "k", k);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lambda(static_cast<Eigen::VectorXd::Index>(k));
            (void) lambda;  // dummy to suppress unused var warning

            stan::math::initialize(lambda, DUMMY_VAR__);
            stan::math::fill(lambda,DUMMY_VAR__);
            current_statement_begin__ = 36;
            validate_non_negative_index("prec", "k", k);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  prec(static_cast<Eigen::VectorXd::Index>(k));
            (void) prec;  // dummy to suppress unused var warning

            stan::math::initialize(prec, DUMMY_VAR__);
            stan::math::fill(prec,DUMMY_VAR__);


            current_statement_begin__ = 37;
            for (int i = 1; i <= k; ++i) {
                current_statement_begin__ = 37;
                stan::model::assign(lambda, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            get_base1(lambda_unconstrained,get_base1(index_vector,i,"index_vector",1),"lambda_unconstrained",1), 
                            "assigning variable lambda");
            }
            current_statement_begin__ = 39;
            stan::math::assign(rho, divide(lambda,sum(lambda)));
            current_statement_begin__ = 40;
            stan::math::assign(prec, multiply(multiply(rho,tau),k));
            current_statement_begin__ = 41;
            stan::math::assign(sds, elt_divide(1.0,stan::math::sqrt(multiply(multiply(rho,tau),k))));

            // validate transformed parameters
            for (int i0__ = 0; i0__ < k; ++i0__) {
                if (stan::math::is_uninitialized(sds(i0__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: sds" << '[' << i0__ << ']';
                    throw std::runtime_error(msg__.str());
                }
            }
            for (int i0__ = 0; i0__ < k; ++i0__) {
                if (stan::math::is_uninitialized(rho(i0__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: rho" << '[' << i0__ << ']';
                    throw std::runtime_error(msg__.str());
                }
            }
            for (int i0__ = 0; i0__ < k; ++i0__) {
                if (stan::math::is_uninitialized(lambda(i0__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: lambda" << '[' << i0__ << ']';
                    throw std::runtime_error(msg__.str());
                }
            }
            for (int i0__ = 0; i0__ < k; ++i0__) {
                if (stan::math::is_uninitialized(prec(i0__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: prec" << '[' << i0__ << ']';
                    throw std::runtime_error(msg__.str());
                }
            }

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 33;
            current_statement_begin__ = 34;
            stan::math::check_simplex(function__,"rho",rho);
            current_statement_begin__ = 35;
            current_statement_begin__ = 36;

            // model body

            current_statement_begin__ = 46;
            lp_accum__.add(-(stan::math::log(tau)));
            current_statement_begin__ = 47;
            lp_accum__.add(gamma_lpdf<propto__>(lambda_unconstrained, alpha, 1));
            current_statement_begin__ = 49;
            if (as_bool(logical_negation(logical_eq(priors_only,1)))) {

                current_statement_begin__ = 51;
                lp_accum__.add(dot_product(n,stan::math::log(prec)));
                current_statement_begin__ = 52;
                lp_accum__.add((-(0.5) * dot_product(prec,b)));
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
        stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
        static const char* function__ = "model_Mixed_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double tau = in__.scalar_lb_constrain(0);
        vector<double> lambda_unconstrained;
        size_t dim_lambda_unconstrained_0__ = (k - nr_equal);
        for (size_t k_0__ = 0; k_0__ < dim_lambda_unconstrained_0__; ++k_0__) {
            lambda_unconstrained.push_back(in__.scalar_lb_constrain(0));
        }
        vars__.push_back(tau);
            for (int k_0__ = 0; k_0__ < (k - nr_equal); ++k_0__) {
            vars__.push_back(lambda_unconstrained[k_0__]);
            }

        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {
            current_statement_begin__ = 33;
            validate_non_negative_index("sds", "k", k);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  sds(static_cast<Eigen::VectorXd::Index>(k));
            (void) sds;  // dummy to suppress unused var warning

            stan::math::initialize(sds, DUMMY_VAR__);
            stan::math::fill(sds,DUMMY_VAR__);
            current_statement_begin__ = 34;
            validate_non_negative_index("rho", "k", k);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  rho(static_cast<Eigen::VectorXd::Index>(k));
            (void) rho;  // dummy to suppress unused var warning

            stan::math::initialize(rho, DUMMY_VAR__);
            stan::math::fill(rho,DUMMY_VAR__);
            current_statement_begin__ = 35;
            validate_non_negative_index("lambda", "k", k);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lambda(static_cast<Eigen::VectorXd::Index>(k));
            (void) lambda;  // dummy to suppress unused var warning

            stan::math::initialize(lambda, DUMMY_VAR__);
            stan::math::fill(lambda,DUMMY_VAR__);
            current_statement_begin__ = 36;
            validate_non_negative_index("prec", "k", k);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  prec(static_cast<Eigen::VectorXd::Index>(k));
            (void) prec;  // dummy to suppress unused var warning

            stan::math::initialize(prec, DUMMY_VAR__);
            stan::math::fill(prec,DUMMY_VAR__);


            current_statement_begin__ = 37;
            for (int i = 1; i <= k; ++i) {
                current_statement_begin__ = 37;
                stan::model::assign(lambda, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            get_base1(lambda_unconstrained,get_base1(index_vector,i,"index_vector",1),"lambda_unconstrained",1), 
                            "assigning variable lambda");
            }
            current_statement_begin__ = 39;
            stan::math::assign(rho, divide(lambda,sum(lambda)));
            current_statement_begin__ = 40;
            stan::math::assign(prec, multiply(multiply(rho,tau),k));
            current_statement_begin__ = 41;
            stan::math::assign(sds, elt_divide(1.0,stan::math::sqrt(multiply(multiply(rho,tau),k))));

            // validate transformed parameters
            current_statement_begin__ = 33;
            current_statement_begin__ = 34;
            stan::math::check_simplex(function__,"rho",rho);
            current_statement_begin__ = 35;
            current_statement_begin__ = 36;

            // write transformed parameters
            if (include_tparams__) {
            for (int k_0__ = 0; k_0__ < k; ++k_0__) {
            vars__.push_back(sds[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < k; ++k_0__) {
            vars__.push_back(rho[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < k; ++k_0__) {
            vars__.push_back(lambda[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < k; ++k_0__) {
            vars__.push_back(prec[k_0__]);
            }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities



            // validate generated quantities

            // write generated quantities
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
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "model_Mixed";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= (k - nr_equal); ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda_unconstrained" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            for (int k_0__ = 1; k_0__ <= k; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "sds" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_0__ = 1; k_0__ <= k; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "rho" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_0__ = 1; k_0__ <= k; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "lambda" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_0__ = 1; k_0__ <= k; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "prec" << '.' << k_0__;
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
        for (int k_0__ = 1; k_0__ <= (k - nr_equal); ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda_unconstrained" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            for (int k_0__ = 1; k_0__ <= k; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "sds" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_0__ = 1; k_0__ <= (k - 1); ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "rho" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_0__ = 1; k_0__ <= k; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "lambda" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_0__ = 1; k_0__ <= k; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "prec" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
        }


        if (!include_gqs__) return;
    }

}; // model

}

typedef model_Mixed_namespace::model_Mixed stan_model;


#endif
