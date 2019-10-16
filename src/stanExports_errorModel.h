// Generated by rstantools.  Do not edit by hand.

/*
    deconR is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    deconR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with deconR.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.19.1
#include <stan/model/model_header.hpp>
namespace model_errorModel_namespace {
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
    reader.add_event(0, 0, "start", "model_errorModel");
    reader.add_event(104, 102, "end", "model_errorModel");
    return reader;
}
#include <stan_meta_header.hpp>
class model_errorModel : public prob_grad {
private:
        int numCellTypes;
        int numSamples;
        vector_d genotype;
        vector_d y;
        vector_d measProp;
        vector_d sd2;
public:
    model_errorModel(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_errorModel(stan::io::var_context& context__,
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
        static const char* function__ = "model_errorModel_namespace::model_errorModel";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 15;
            context__.validate_dims("data initialization", "numCellTypes", "int", context__.to_vec());
            numCellTypes = int(0);
            vals_i__ = context__.vals_i("numCellTypes");
            pos__ = 0;
            numCellTypes = vals_i__[pos__++];
            check_greater_or_equal(function__, "numCellTypes", numCellTypes, 0);
            current_statement_begin__ = 17;
            context__.validate_dims("data initialization", "numSamples", "int", context__.to_vec());
            numSamples = int(0);
            vals_i__ = context__.vals_i("numSamples");
            pos__ = 0;
            numSamples = vals_i__[pos__++];
            check_greater_or_equal(function__, "numSamples", numSamples, 1);
            current_statement_begin__ = 19;
            validate_non_negative_index("genotype", "numSamples", numSamples);
            context__.validate_dims("data initialization", "genotype", "vector_d", context__.to_vec(numSamples));
            genotype = Eigen::Matrix<double, Eigen::Dynamic, 1>(numSamples);
            vals_r__ = context__.vals_r("genotype");
            pos__ = 0;
            size_t genotype_j_1_max__ = numSamples;
            for (size_t j_1__ = 0; j_1__ < genotype_j_1_max__; ++j_1__) {
                genotype(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 21;
            validate_non_negative_index("y", "numSamples", numSamples);
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec(numSamples));
            y = Eigen::Matrix<double, Eigen::Dynamic, 1>(numSamples);
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_1_max__ = numSamples;
            for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 24;
            validate_non_negative_index("measProp", "numSamples", numSamples);
            context__.validate_dims("data initialization", "measProp", "vector_d", context__.to_vec(numSamples));
            measProp = Eigen::Matrix<double, Eigen::Dynamic, 1>(numSamples);
            vals_r__ = context__.vals_r("measProp");
            pos__ = 0;
            size_t measProp_j_1_max__ = numSamples;
            for (size_t j_1__ = 0; j_1__ < measProp_j_1_max__; ++j_1__) {
                measProp(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 27;
            validate_non_negative_index("sd2", "numSamples", numSamples);
            context__.validate_dims("data initialization", "sd2", "vector_d", context__.to_vec(numSamples));
            sd2 = Eigen::Matrix<double, Eigen::Dynamic, 1>(numSamples);
            vals_r__ = context__.vals_r("sd2");
            pos__ = 0;
            size_t sd2_j_1_max__ = numSamples;
            for (size_t j_1__ = 0; j_1__ < sd2_j_1_max__; ++j_1__) {
                sd2(j_1__) = vals_r__[pos__++];
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 34;
            num_params_r__ += 1;
            current_statement_begin__ = 37;
            validate_non_negative_index("beta", "numCellTypes", numCellTypes);
            num_params_r__ += numCellTypes;
            current_statement_begin__ = 43;
            validate_non_negative_index("Sigma", "numCellTypes", numCellTypes);
            num_params_r__ += numCellTypes;
            current_statement_begin__ = 46;
            validate_non_negative_index("trueProp", "numSamples", numSamples);
            num_params_r__ += numSamples;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_errorModel() { }
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
        current_statement_begin__ = 34;
        if (!(context__.contains_r("beta0")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta0 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta0");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "beta0", "double", context__.to_vec());
        double beta0(0);
        beta0 = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(beta0);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta0: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 37;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "numCellTypes", numCellTypes);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(numCellTypes));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta(numCellTypes);
        size_t beta_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            beta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 43;
        if (!(context__.contains_r("Sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable Sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("Sigma");
        pos__ = 0U;
        validate_non_negative_index("Sigma", "numCellTypes", numCellTypes);
        context__.validate_dims("parameter initialization", "Sigma", "vector_d", context__.to_vec(numCellTypes));
        Eigen::Matrix<double, Eigen::Dynamic, 1> Sigma(numCellTypes);
        size_t Sigma_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < Sigma_j_1_max__; ++j_1__) {
            Sigma(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_lb_unconstrain(0, Sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable Sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 46;
        if (!(context__.contains_r("trueProp")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable trueProp missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("trueProp");
        pos__ = 0U;
        validate_non_negative_index("trueProp", "numSamples", numSamples);
        context__.validate_dims("parameter initialization", "trueProp", "vector_d", context__.to_vec(numSamples));
        Eigen::Matrix<double, Eigen::Dynamic, 1> trueProp(numSamples);
        size_t trueProp_j_1_max__ = numSamples;
        for (size_t j_1__ = 0; j_1__ < trueProp_j_1_max__; ++j_1__) {
            trueProp(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(trueProp);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable trueProp: ") + e.what()), current_statement_begin__, prog_reader__());
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
            current_statement_begin__ = 34;
            local_scalar_t__ beta0;
            (void) beta0;  // dummy to suppress unused var warning
            if (jacobian__)
                beta0 = in__.scalar_constrain(lp__);
            else
                beta0 = in__.scalar_constrain();
            current_statement_begin__ = 37;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_constrain(numCellTypes, lp__);
            else
                beta = in__.vector_constrain(numCellTypes);
            current_statement_begin__ = 43;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> Sigma;
            (void) Sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                Sigma = in__.vector_lb_constrain(0, numCellTypes, lp__);
            else
                Sigma = in__.vector_lb_constrain(0, numCellTypes);
            current_statement_begin__ = 46;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> trueProp;
            (void) trueProp;  // dummy to suppress unused var warning
            if (jacobian__)
                trueProp = in__.vector_constrain(numSamples, lp__);
            else
                trueProp = in__.vector_constrain(numSamples);
            // transformed parameters
            current_statement_begin__ = 59;
            local_scalar_t__ betaNormal;
            (void) betaNormal;  // dummy to suppress unused var warning
            stan::math::initialize(betaNormal, DUMMY_VAR__);
            stan::math::fill(betaNormal, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 60;
            stan::math::assign(betaNormal, (get_base1(beta, 1, "beta", 1) + get_base1(beta, 3, "beta", 1)));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 59;
            if (stan::math::is_uninitialized(betaNormal)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: betaNormal";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable betaNormal: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            // model body
            current_statement_begin__ = 74;
            lp_accum__.add(normal_log<propto__>(beta0, 0.0, 10));
            current_statement_begin__ = 75;
            lp_accum__.add(normal_log<propto__>(get_base1(beta, 1, "beta", 1), 0.0, 10));
            current_statement_begin__ = 76;
            lp_accum__.add(normal_log<propto__>(get_base1(beta, 2, "beta", 1), 0.0, 10));
            current_statement_begin__ = 77;
            lp_accum__.add(normal_log<propto__>(get_base1(beta, 3, "beta", 1), 0.0, 10));
            current_statement_begin__ = 78;
            lp_accum__.add(gamma_log<propto__>(get_base1(Sigma, 1, "Sigma", 1), (1.0 / 10), (1.0 / 10)));
            current_statement_begin__ = 79;
            lp_accum__.add(gamma_log<propto__>(get_base1(Sigma, 2, "Sigma", 1), (1.0 / 10), (1.0 / 10)));
            current_statement_begin__ = 85;
            lp_accum__.add(normal_log<propto__>(trueProp, .5, 5));
            current_statement_begin__ = 90;
            lp_accum__.add(normal_log<propto__>(measProp, trueProp, sd2));
            current_statement_begin__ = 96;
            for (int i = 1; i <= numSamples; ++i) {
                current_statement_begin__ = 98;
                lp_accum__.add(normal_log<propto__>(get_base1(y, i, "y", 1), (((beta0 + (get_base1(beta, 1, "beta", 1) * get_base1(genotype, i, "genotype", 1))) + (get_base1(beta, 2, "beta", 1) * get_base1(trueProp, i, "trueProp", 1))) + (get_base1(beta, 3, "beta", 1) * (get_base1(trueProp, i, "trueProp", 1) * get_base1(genotype, i, "genotype", 1)))), stan::math::sqrt(((((1 - (2 * get_base1(trueProp, i, "trueProp", 1))) + pow(get_base1(trueProp, i, "trueProp", 1), 2)) / get_base1(Sigma, 1, "Sigma", 1)) + (pow(get_base1(trueProp, i, "trueProp", 1), 2) / get_base1(Sigma, 2, "Sigma", 1))))));
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
        names__.push_back("beta0");
        names__.push_back("beta");
        names__.push_back("Sigma");
        names__.push_back("trueProp");
        names__.push_back("betaNormal");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(numCellTypes);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(numCellTypes);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(numSamples);
        dimss__.push_back(dims__);
        dims__.resize(0);
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
        static const char* function__ = "model_errorModel_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double beta0 = in__.scalar_constrain();
        vars__.push_back(beta0);
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta = in__.vector_constrain(numCellTypes);
        size_t beta_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            vars__.push_back(beta(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> Sigma = in__.vector_lb_constrain(0, numCellTypes);
        size_t Sigma_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < Sigma_j_1_max__; ++j_1__) {
            vars__.push_back(Sigma(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> trueProp = in__.vector_constrain(numSamples);
        size_t trueProp_j_1_max__ = numSamples;
        for (size_t j_1__ = 0; j_1__ < trueProp_j_1_max__; ++j_1__) {
            vars__.push_back(trueProp(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 59;
            double betaNormal;
            (void) betaNormal;  // dummy to suppress unused var warning
            stan::math::initialize(betaNormal, DUMMY_VAR__);
            stan::math::fill(betaNormal, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 60;
            stan::math::assign(betaNormal, (get_base1(beta, 1, "beta", 1) + get_base1(beta, 3, "beta", 1)));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                vars__.push_back(betaNormal);
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
    static std::string model_name() {
        return "model_errorModel";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "beta0";
        param_names__.push_back(param_name_stream__.str());
        size_t beta_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t Sigma_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < Sigma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "Sigma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t trueProp_j_1_max__ = numSamples;
        for (size_t j_1__ = 0; j_1__ < trueProp_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "trueProp" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "betaNormal";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "beta0";
        param_names__.push_back(param_name_stream__.str());
        size_t beta_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t Sigma_j_1_max__ = numCellTypes;
        for (size_t j_1__ = 0; j_1__ < Sigma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "Sigma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t trueProp_j_1_max__ = numSamples;
        for (size_t j_1__ = 0; j_1__ < trueProp_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "trueProp" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "betaNormal";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_errorModel_namespace::model_errorModel stan_model;
#endif