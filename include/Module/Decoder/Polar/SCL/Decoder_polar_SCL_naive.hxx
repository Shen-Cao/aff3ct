#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <algorithm>
#include <sstream>
#include <limits>
#include <cmath>
#include <tuple>

#include "Tools/Exception/exception.hpp"
#include "Tools/Math/utils.h"
#include "Tools/Code/Polar/fb_assert.h"
#include "Module/Decoder/Polar/SCL/Decoder_polar_SCL_naive.hpp"

namespace aff3ct
{
namespace module
{

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
Task& Decoder_polar_SCL_naive<B,R,F,G>
::operator[](const dec::tsk t)
{
	return Module::operator[]((size_t)t);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
Socket& Decoder_polar_SCL_naive<B,R,F,G>
::operator[](const dec::sck::decode_siho_cw_flexible_frozen s)
{
	return Module::operator[]((size_t)dec::tsk::decode_siho_cw_flexible_frozen)[(size_t)s];
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
Socket& Decoder_polar_SCL_naive<B,R,F,G>
::operator[](const dec::sck::set_known_bits s)
{
	return Module::operator[]((size_t)dec::tsk::set_known_bits)[(size_t)s];
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
Decoder_polar_SCL_naive<B,R,F,G>
::Decoder_polar_SCL_naive(const int& K, const int& N, const int& L, const std::vector<bool>& frozen_bits)
: Decoder_SIHO<B,R>(K, N),
  m((int)std::log2(N)),
  metric_init(std::numeric_limits<R>::min()),
  frozen_bits(frozen_bits),
  L(L),
  polar_trees(L, tools::Binary_tree_metric<Contents_SCL<B,R>,R>(this->m + 1, metric_init)), 
  // add a bool flag to indicate whether the decoding was successful or not
  decode_success((bool)true)
{
	const std::string name = "Decoder_polar_SCL_naive";
	this->set_name(name);

	// SCL decoder will not fail to decode, so we set the fail flag to true by default
	this->decode_success = true;

	if (!tools::is_power_of_2(this->N))
	{
		std::stringstream message;
		message << "'N' has to be a power of 2 ('N' = " << N << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	if (this->N != (int)frozen_bits.size())
	{
		std::stringstream message;
		message << "'frozen_bits.size()' has to be equal to 'N' ('frozen_bits.size()' = " << frozen_bits.size()
		        << ", 'N' = " << N << ").";
		throw tools::length_error(__FILE__, __LINE__, __func__, message.str());
	}

	if (this->L <= 0 || !tools::is_power_of_2(this->L))
	{
		std::stringstream message;
		message << "'L' has to be a positive power of 2 ('L' = " << L << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	auto k = 0; for (auto i = 0; i < this->N; i++) if (frozen_bits[i] == 0) k++;
	if (this->K != k)
	{
		std::stringstream message;
		message << "The number of information bits in the frozen_bits is invalid ('K' = " << K << ", 'k' = "
		        << k << ").";
		throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
	}

	this->active_paths.insert(0);
	for (auto i = 0; i < L; i++)
	{
		this->recursive_allocate_nodes_contents(this->polar_trees[i].get_root(), this->N);
		this->recursive_initialize_frozen_bits(this->polar_trees[i].get_root(), frozen_bits);
	}
	for (auto i = 0; i < L; i++)
		leaves_array.push_back(this->polar_trees[i].get_leaves());

	auto &pf = this->create_task("decode_siho_cw_flexible_frozen");
	auto pfs_Y_N = this->template create_socket_in <R     >(pf, "Y_N", this->N); // R for float input
	auto pfs_F_N = this->template create_socket_in <B     >(pf, "F_N", this->N); // B for int input
	auto pfs_CWD = this->template create_socket_out<int8_t>(pf, "CWD",       1);
	auto pfs_V_N = this->template create_socket_out<B     >(pf, "V_N", this->N);
	this->create_codelet(pf, [pfs_Y_N, pfs_F_N, pfs_CWD, pfs_V_N](Module &m, Task &t, const size_t frame_id) -> int
	{
		auto &dec = static_cast<Decoder_polar_SCL_naive<B,R,F,G>&>(m);

		auto ret = dec._decode_siho_cw_flexible_frozen(static_cast<R     *>(t[pfs_Y_N].get_dataptr()),
		                            				   static_cast<B     *>(t[pfs_F_N].get_dataptr()),
		                            				   static_cast<int8_t*>(t[pfs_CWD].get_dataptr()),
		                            				   static_cast<B     *>(t[pfs_V_N].get_dataptr()),
		                            				   frame_id);

		if (dec.is_auto_reset())
			dec._reset(frame_id);

		return ret;
	});

	auto &pk = this->create_task("set_known_bits");
	auto pks_KNOWN_BITS = this->template create_socket_in<B>(pk, "KNOWN_BITS", this->N);
	this->create_codelet(pk, [pks_KNOWN_BITS](Module &m, Task &t, const size_t frame_id) -> int
	{
		auto &dec = static_cast<Decoder_polar_SCL_naive<B,R,F,G>&>(m);

		dec._set_known_bits(static_cast<B*>(t[pks_KNOWN_BITS].get_dataptr()), frame_id);

		return status_t::SUCCESS;
	});
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
Decoder_polar_SCL_naive<B,R,F,G>
::~Decoder_polar_SCL_naive()
{
	for (auto i = 0; i < L; i++)
		this->recursive_deallocate_nodes_contents(this->polar_trees[i].get_root());
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
Decoder_polar_SCL_naive<B,R,F,G>* Decoder_polar_SCL_naive<B,R,F,G>
::clone() const
{
	auto m = new Decoder_polar_SCL_naive(*this);
	m->deep_copy(*this);
	return m;
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::deep_copy(const Decoder_polar_SCL_naive<B,R,F,G> &m)
{
	Module::deep_copy(m);
	this->leaves_array.clear();
	for (auto i = 0; i < L; i++)
	{
		this->recursive_deep_copy(m.polar_trees[i].get_root(), this->polar_trees[i].get_root());
		this->leaves_array.push_back(this->polar_trees[i].get_leaves());
	}
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::recursive_deep_copy(const tools::Binary_node<Contents_SCL<B,R>> *nref, tools::Binary_node<Contents_SCL<B,R>> *nclone)
{
	auto cref = nref->get_contents();
	auto cclone = new Contents_SCL<B,R>(*cref);
	nclone->set_contents(cclone);

	if (!nref->is_leaf() && !nclone->is_leaf())
	{
		if (nref->get_left() != nullptr)
			this->recursive_deep_copy(nref->get_left(), nclone->get_left());

		if (nref->get_right() != nullptr)
			this->recursive_deep_copy(nref->get_right(), nclone->get_right());
	}
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::set_frozen_bits(const std::vector<bool>& fb)
{
	aff3ct::tools::fb_assert(frozen_bits, this->K, this->N);
	std::copy(fb.begin(), fb.end(), this->frozen_bits.begin());
	for (auto i = 0; i < L; i++)
		this->recursive_initialize_frozen_bits(this->polar_trees[i].get_root(), frozen_bits);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
const std::vector<bool>& Decoder_polar_SCL_naive<B,R,F,G>
::get_frozen_bits() const
{
	return this->frozen_bits;
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::_set_known_bits(const B *KNOWN_BITS, const size_t frame_id)
{
	this->known_bits.assign(KNOWN_BITS, KNOWN_BITS + this->N);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::_load(const R *Y_N)
{

	for (auto path = 0; path < this->L; path ++)
	{
		auto *contents = this->polar_trees[path].get_root()->get_contents();
		for (auto i = 0; i < this->N; i++)
			contents->lambda[i] = Y_N[i];

		polar_trees[path].set_path_metric(metric_init);
	}

	// initialization
	active_paths.clear();
	active_paths.insert(0);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::_decode(const size_t frame_id)
{
	std::set<int> last_active_paths;
	int cur_path;

	// tuples to be sorted. <Path,estimated bit,metric>
	std::vector<std::tuple<int,B,R>> metrics_vec;

	// run through each leaf
	for (auto leaf_index = 0 ; leaf_index < this->N; leaf_index++)
	{
		// compute LLR for current leaf
		for (auto path : active_paths)
			this->recursive_compute_llr(leaves_array[path][leaf_index], tools::compute_depth(leaf_index, this->m));

		// if current leaf is a frozen bit
		if (leaves_array[0][leaf_index]->get_c()->is_frozen_bit)
		{
			auto min_phi = std::numeric_limits<R>::max();
			for (auto path : active_paths)
			{
				auto cur_leaf = leaves_array[path][leaf_index];
				cur_leaf->get_c()->s[0] = 0;
				auto phi_cur = tools::phi<R>(polar_trees[path].get_path_metric(), cur_leaf->get_c()->lambda[0], 0);
				this->polar_trees[path].set_path_metric(phi_cur);
				min_phi = std::min<R>(min_phi, phi_cur);
			}

			// normalization
			for (auto path : active_paths)
				this->polar_trees[path].set_path_metric(this->polar_trees[path].get_path_metric() - min_phi);
		}
		else
		{
			
			// metrics vec used to store values of hypothetic path metrics
			metrics_vec.clear();
			auto min_phi = std::numeric_limits<R>::max();
			for (auto path : active_paths)
			{
				auto cur_leaf = leaves_array[path][leaf_index];
				R phi0 = tools::phi<B,R>(polar_trees[path].get_path_metric(), cur_leaf->get_c()->lambda[0],                 (B)0);
				R phi1 = tools::phi<B,R>(polar_trees[path].get_path_metric(), cur_leaf->get_c()->lambda[0], tools::bit_init<B>());
				metrics_vec.push_back(std::make_tuple(path,                 (B)0, phi0));
				metrics_vec.push_back(std::make_tuple(path, tools::bit_init<B>(), phi1));

				min_phi = std::min<R>(min_phi, phi0);
				min_phi = std::min<R>(min_phi, phi1);
			}

			// normalization
			for (auto vec : metrics_vec)
				std::get<2>(vec) -= min_phi;

			if (active_paths.size() <= (unsigned)(L / 2))
			{
				last_active_paths = active_paths;
				for (auto path : last_active_paths)
					this->duplicate_path(path, leaf_index);
			}
			else
			{
				// sort hypothetic metrics
				std::sort(metrics_vec.begin(), metrics_vec.end(),
					[](std::tuple<int,B,R> x, std::tuple<int,B,R> y){
						return std::get<2>(x) < std::get<2>(y);
					});

				// search in worst metrics. If a path is found twice, erase it
				for (auto it = metrics_vec.begin() + metrics_vec.size() / 2; it != metrics_vec.end(); ++it)
				{
					cur_path = std::get<0>(*it);

					auto it_double = std::find_if(it +1, metrics_vec.end(),
						[cur_path](std::tuple<int,B,R> x){
							return std::get<0>(x) == cur_path;
						});

					if (it_double != metrics_vec.end())
						active_paths.erase(std::get<0>(*it));
				}

				// remove worst metrics from list
				metrics_vec.resize(metrics_vec.size() / 2);

				for (auto it = metrics_vec.begin(); it != metrics_vec.end(); ++it)
				{
					cur_path = std::get<0>(*it);

					auto it_double = std::find_if(it +1, metrics_vec.end(),
						[cur_path](std::tuple<int,B,R> x){
							return std::get<0>(x) == cur_path;
						});

					if (it_double != metrics_vec.end())
					{
						// duplicate
						metrics_vec.erase(it_double);
						duplicate_path(std::get<0>(*it), leaf_index);
					}
					else
					{
						// choose
						leaves_array[std::get<0>(*it)][leaf_index]->get_c()->s[0] = std::get<1>(*it);
						polar_trees[std::get<0>(*it)].set_path_metric(std::get<2>(*it));
					}
				}
			}
		}

		// propagate sums
		for (auto path : active_paths)
			this->propagate_sums(leaves_array[path][leaf_index]);
	}

	this->select_best_path(frame_id);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::_decode_flexible_frozen(const B *F_N, const size_t frame_id)
{
	std::set<int> last_active_paths;
	int cur_path;

	// tuples to be sorted. <Path,estimated bit,metric>
	std::vector<std::tuple<int,B,R>> metrics_vec;

	// run through each leaf
	for (auto leaf_index = 0 ; leaf_index < this->N; leaf_index++)
	{
		// compute LLR for current leaf
		for (auto path : active_paths)
			this->recursive_compute_llr(leaves_array[path][leaf_index], tools::compute_depth(leaf_index, this->m));

		// if current leaf is a frozen bit
		if (leaves_array[0][leaf_index]->get_c()->is_frozen_bit)
		{
			auto min_phi = std::numeric_limits<R>::max();
			for (auto path : active_paths)
			{
				auto cur_leaf = leaves_array[path][leaf_index];
                // force u_i to the frozen value f (NOT always 0)
                // cur_leaf->get_c()->s[0] = F_N[leaf_index]; // This implementation is wrong...
				if (F_N[leaf_index] == 0)
				{
					cur_leaf->get_c()->s[0] = 0;
				}
				else
				{
					cur_leaf->get_c()->s[0] = tools::bit_init<B>();
				}
				
				// For frozen bits, we use the provided value directly as it represents the bit value (0 or 1)
				auto phi_cur = tools::phi<R>(polar_trees[path].get_path_metric(), cur_leaf->get_c()->lambda[0], F_N[leaf_index]);
				this->polar_trees[path].set_path_metric(phi_cur);
				min_phi = std::min<R>(min_phi, phi_cur);
			}

			// normalization
			for (auto path : active_paths)
				this->polar_trees[path].set_path_metric(this->polar_trees[path].get_path_metric() - min_phi);
		}
		else
		{
			
			// metrics vec used to store values of hypothetic path metrics
			metrics_vec.clear();
			auto min_phi = std::numeric_limits<R>::max();
			for (auto path : active_paths)
			{
				auto cur_leaf = leaves_array[path][leaf_index];
				R phi0 = tools::phi<B,R>(polar_trees[path].get_path_metric(), cur_leaf->get_c()->lambda[0],                 (B)0);
				R phi1 = tools::phi<B,R>(polar_trees[path].get_path_metric(), cur_leaf->get_c()->lambda[0], tools::bit_init<B>());
				metrics_vec.push_back(std::make_tuple(path,                 (B)0, phi0));
				metrics_vec.push_back(std::make_tuple(path, tools::bit_init<B>(), phi1));

				min_phi = std::min<R>(min_phi, phi0);
				min_phi = std::min<R>(min_phi, phi1);
			}

			// normalization
			for (auto vec : metrics_vec)
				std::get<2>(vec) -= min_phi;

			if (active_paths.size() <= (unsigned)(L / 2))
			{
				last_active_paths = active_paths;
				for (auto path : last_active_paths)
					this->duplicate_path(path, leaf_index);
			}
			else
			{
				// sort hypothetic metrics
				std::sort(metrics_vec.begin(), metrics_vec.end(),
					[](std::tuple<int,B,R> x, std::tuple<int,B,R> y){
						return std::get<2>(x) < std::get<2>(y);
					});

				// search in worst metrics. If a path is found twice, erase it
				for (auto it = metrics_vec.begin() + metrics_vec.size() / 2; it != metrics_vec.end(); ++it)
				{
					cur_path = std::get<0>(*it);

					auto it_double = std::find_if(it +1, metrics_vec.end(),
						[cur_path](std::tuple<int,B,R> x){
							return std::get<0>(x) == cur_path;
						});

					if (it_double != metrics_vec.end())
						active_paths.erase(std::get<0>(*it));
				}

				// remove worst metrics from list
				metrics_vec.resize(metrics_vec.size() / 2);

				for (auto it = metrics_vec.begin(); it != metrics_vec.end(); ++it)
				{
					cur_path = std::get<0>(*it);

					auto it_double = std::find_if(it +1, metrics_vec.end(),
						[cur_path](std::tuple<int,B,R> x){
							return std::get<0>(x) == cur_path;
						});

					if (it_double != metrics_vec.end())
					{
						// duplicate
						metrics_vec.erase(it_double);
						duplicate_path(std::get<0>(*it), leaf_index);
					}
					else
					{
						// choose
						leaves_array[std::get<0>(*it)][leaf_index]->get_c()->s[0] = std::get<1>(*it);
						polar_trees[std::get<0>(*it)].set_path_metric(std::get<2>(*it));
					}
				}
			}
		}

		// propagate sums
		for (auto path : active_paths)
			this->propagate_sums(leaves_array[path][leaf_index]);
	}

	this->select_best_path(frame_id);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
int Decoder_polar_SCL_naive<B,R,F,G>
::_decode_siho(const R *Y_N, B *V_K, const size_t frame_id)
{
//	auto t_load = std::chrono::steady_clock::now(); // ----------------------------------------------------------- LOAD
	this->_load(Y_N);
//	auto d_load = std::chrono::steady_clock::now() - t_load;

//	auto t_decod = std::chrono::steady_clock::now(); // -------------------------------------------------------- DECODE
	this->_decode(frame_id);
//	auto d_decod = std::chrono::steady_clock::now() - t_decod;

//	auto t_store = std::chrono::steady_clock::now(); // --------------------------------------------------------- STORE
	this->_store(V_K);
//	auto d_store = std::chrono::steady_clock::now() - t_store;

//	(*this)[dec::tsk::decode_siho].update_timer(dec::tm::decode_siho::load,   d_load);
//	(*this)[dec::tsk::decode_siho].update_timer(dec::tm::decode_siho::decode, d_decod);
//	(*this)[dec::tsk::decode_siho].update_timer(dec::tm::decode_siho::store,  d_store);

	return this->decode_success ? 0 : -1;
	// return 0;
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
int Decoder_polar_SCL_naive<B,R,F,G>
::_decode_siho_cw(const R *Y_N, B *V_N, const size_t frame_id)
{
//	auto t_load = std::chrono::steady_clock::now(); // ----------------------------------------------------------- LOAD
	this->_load(Y_N);
//	auto d_load = std::chrono::steady_clock::now() - t_load;

//	auto t_decod = std::chrono::steady_clock::now(); // -------------------------------------------------------- DECODE
	this->_decode(frame_id);
//	auto d_decod = std::chrono::steady_clock::now() - t_decod;

//	auto t_store = std::chrono::steady_clock::now(); // --------------------------------------------------------- STORE
	this->_store(V_N, true);
//	auto d_store = std::chrono::steady_clock::now() - t_store;

//	(*this)[dec::tsk::decode_siho_cw].update_timer(dec::tm::decode_siho_cw::load,   d_load);
//	(*this)[dec::tsk::decode_siho_cw].update_timer(dec::tm::decode_siho_cw::decode, d_decod);
//	(*this)[dec::tsk::decode_siho_cw].update_timer(dec::tm::decode_siho_cw::store,  d_store);

	return this->decode_success ? 0 : -1;
	// return 0;
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
int Decoder_polar_SCL_naive<B,R,F,G>
::_decode_siho_cw_flexible_frozen(const R *Y_N, const B *F_N, int8_t *CWD, B *V_N, const size_t frame_id)
{
	const auto status = this->_decode_siho_cw_flexible_frozen(Y_N, F_N, V_N, frame_id);
	std::fill(CWD, CWD + this->get_n_frames_per_wave(), 0);
	return status;
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
int Decoder_polar_SCL_naive<B,R,F,G>
::_decode_siho_cw_flexible_frozen(const R *Y_N, const B *F_N, B *V_N, const size_t frame_id)
{
//	auto t_load = std::chrono::steady_clock::now(); // ----------------------------------------------------------- LOAD
	// Amplify the channel LLRs at the known input-side positions (the BGL
	// shortened and punctured padding bits, fixed to the public constant 0).
	// Uses exactly the same float arithmetic as the former Python-side
	// shortener (add_LLR = 16): y_out = y + ((m * y) * 16) with m in {0,1},
	// so y_out = 17 * y at the masked positions, bit-for-bit identical.
	const R *Y_N_eff = Y_N;
	if ((int)this->known_bits.size() == this->N)
	{
		this->Y_N_known.resize(this->N);
		for (auto i = 0; i < this->N; i++)
		{
			const R m = (R)this->known_bits[i];
			this->Y_N_known[i] = Y_N[i] + ((m * Y_N[i]) * (R)16.0f);
		}
		Y_N_eff = this->Y_N_known.data();
	}
	this->_load(Y_N_eff);
//	auto d_load = std::chrono::steady_clock::now() - t_load;

//	auto t_decod = std::chrono::steady_clock::now(); // -------------------------------------------------------- DECODE
	this->_decode_flexible_frozen(F_N, frame_id);
//	auto d_decod = std::chrono::steady_clock::now() - t_decod;

//	auto t_store = std::chrono::steady_clock::now(); // --------------------------------------------------------- STORE
	this->_store(V_N, true);
//	auto d_store = std::chrono::steady_clock::now() - t_store;

//	(*this)[dec::tsk::decode_siho].update_timer(dec::tm::decode_siho::load,   d_load);
//	(*this)[dec::tsk::decode_siho].update_timer(dec::tm::decode_siho::decode, d_decod);
//	(*this)[dec::tsk::decode_siho].update_timer(dec::tm::decode_siho::store,  d_store);

	return this->decode_success ? 0 : -1;
	// return 0;
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::_store(B *V, bool coded) const
{
	auto *root = (tools::Binary_node<Contents_SCL<B,R>>*)this->polar_trees[*active_paths.begin()].get_root();
	if (!coded)
	{
		auto k = 0;
		this->recursive_store(root, V, k, false);
	}
	else
	{
		auto k = 0;
		this->recursive_store(root, V, k, true);

		// Store the full codeword by collecting bits from leaf nodes
		// Use the leaves_array to get the correct bit values from leaf nodes
		// auto path = *active_paths.begin();
		// for (int i = 0; i < this->N; i++) {
		// 	// Get the bit value from the corresponding leaf node
		// 	auto *leaf_node = leaves_array[path][i];
		// 	auto *contents = leaf_node->get_c();
		// 	// Convert bit values to 0 or 1
		// 	V[i] = contents->s[0] ? 1 : 0;
		// 	std::cout << "Storing Position " << i << " Value = " << contents->s[0] << std::endl;
		// }

		// std::copy(contents_root->s.begin(), contents_root->s.begin() + this->N, V);
	}
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::recursive_compute_llr(tools::Binary_node<Contents_SCL<B, R> >* node_cur, int depth)
{
	if (depth != 0)
		recursive_compute_llr(node_cur->get_father(), --depth);

	if(node_cur->is_left())
		this->apply_f((tools::Binary_node<Contents_SCL<B,R>>*)node_cur->get_father());

	if(node_cur->is_right())
		this->apply_g((tools::Binary_node<Contents_SCL<B,R>>*)node_cur->get_father());
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::propagate_sums(const tools::Binary_node<Contents_SCL<B, R> >* node_cur)
{
	if(!node_cur->is_leaf())
		this->compute_sums((tools::Binary_node<Contents_SCL<B,R>>*)node_cur);

	if (node_cur->is_right() && !node_cur->is_root())
		this->propagate_sums(node_cur->get_father());

}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::recursive_duplicate_tree_llr(tools::Binary_node<Contents_SCL<B,R>>* node_a,
                               tools::Binary_node<Contents_SCL<B,R> >* node_b)
{
	node_b->get_c()->lambda = node_a->get_c()->lambda;

	if(!node_a->get_father()->is_root())
		this->recursive_duplicate_tree_llr(node_a->get_father(), node_b->get_father());
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::recursive_duplicate_tree_sums(tools::Binary_node<Contents_SCL<B,R>>* node_a,
                                tools::Binary_node<Contents_SCL<B,R>>* node_b,
                                tools::Binary_node<Contents_SCL<B,R>>* node_caller)
{
	if (node_a->get_left() != node_caller && !node_a->is_leaf())
		node_b->get_left()->get_c()->s = node_a->get_left()->get_c()->s;

	if(!node_a->is_root())
		this->recursive_duplicate_tree_sums(node_a->get_father(), node_b->get_father(), node_a);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::duplicate_path(int path, int leaf_index)
{
	std::vector<tools::Binary_node<Contents_SCL<B,R>>*> path_leaves, newpath_leaves;
	int newpath = 0;
	while (active_paths.find(newpath++) != active_paths.end()){};
	newpath--;

	active_paths.insert(newpath);

	path_leaves = leaves_array[path];

	newpath_leaves = leaves_array[newpath];

	for (auto i = 0; i < leaf_index; i++)
		newpath_leaves[i]->get_c()->s = path_leaves[i]->get_c()->s;

	recursive_duplicate_tree_sums(leaves_array[path][leaf_index], leaves_array[newpath][leaf_index], nullptr);

	if (leaf_index < this->N - 1)
		recursive_duplicate_tree_llr(leaves_array[path][leaf_index + 1], leaves_array[newpath][leaf_index + 1]);

	leaves_array[newpath][leaf_index]->get_c()->s[0] = tools::bit_init<B>();
	polar_trees[newpath].set_path_metric(tools::phi<B,R>(polar_trees[path].get_path_metric(),
	                                                     leaves_array[path][leaf_index]->get_c()->lambda[0],
	                                                     tools::bit_init<B>()));

	leaves_array[path][leaf_index]->get_c()->s[0] = 0;
	polar_trees[path].set_path_metric(tools::phi<B,R>(polar_trees[path].get_path_metric(),
	                                                  leaves_array[path][leaf_index]->get_c()->lambda[0],
	                                                  0));
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::select_best_path(const size_t frame_id)
{
	int best_path = 0;
	if (active_paths.size() >= 1)
		best_path = *active_paths.begin();

	for (int path : active_paths)
		if(polar_trees[path].get_path_metric() < polar_trees[best_path].get_path_metric())
			best_path = path;

	active_paths.clear();
	active_paths.insert(best_path);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::recursive_allocate_nodes_contents(tools::Binary_node<Contents_SCL<B,R>>* node_curr, const int vector_size)
{
	if (node_curr != nullptr)
	{
		node_curr->set_contents(new Contents_SCL<B,R>(vector_size));

		this->recursive_allocate_nodes_contents(node_curr->get_left() , vector_size / 2);
		this->recursive_allocate_nodes_contents(node_curr->get_right(), vector_size / 2);
	}
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::recursive_initialize_frozen_bits(const tools::Binary_node<Contents_SCL<B,R>>* node_curr,
                                   const std::vector<bool>& frozen_bits)
{
	auto *contents = node_curr->get_contents();

	if (!node_curr->is_leaf()) // stop condition
	{
		this->recursive_initialize_frozen_bits(node_curr->get_left(),  frozen_bits); // recursive call
		this->recursive_initialize_frozen_bits(node_curr->get_right(), frozen_bits); // recursive call
	}
	else
		contents->is_frozen_bit = frozen_bits[node_curr->get_lane_id()];
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::recursive_store(const tools::Binary_node<Contents_SCL<B,R>>* node_curr, B *V_K, int &k, bool contain_frozen) const
{
	auto *contents = node_curr->get_contents();

	if (!node_curr->is_leaf()) // stop condition
	{
		this->recursive_store(node_curr->get_left(),  V_K, k, contain_frozen); // recursive call
		this->recursive_store(node_curr->get_right(), V_K, k, contain_frozen); // recursive call
	}
	else
	{
		if (contain_frozen)
		{
			V_K[k++] = contents->s[0] ? 1 : 0;
			/*std::cout << "Storing position - " << k << " lane_id - " << node_curr->get_lane_id() << " content value: " << contents->s[0] << " V_N value: " << V_K[k-1] << " Is frozen: " << frozen_bits[node_curr->get_lane_id()] << std::endl;*/
		}
		else
		{
			if (!frozen_bits[node_curr->get_lane_id()])
				{
					V_K[k++] = contents->s[0] ? 1 : 0;
					/*std::cout << "Storing position - " << k << " lane_id - " << node_curr->get_lane_id() << " content value: " << contents->s[0] << " V_N value: " << V_K[k-1] << " Is frozen: " << frozen_bits[node_curr->get_lane_id()] << std::endl;*/
				}
		}
	}
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::recursive_deallocate_nodes_contents(tools::Binary_node<Contents_SCL<B,R>>* node_curr)
{
	if (node_curr != nullptr)
	{
		this->recursive_deallocate_nodes_contents(node_curr->get_left()); // recursive call
		this->recursive_deallocate_nodes_contents(node_curr->get_right()); // recursive call

		auto *contents = node_curr->get_contents();
		delete contents;
		node_curr->set_contents(nullptr);
	}
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::apply_f(const tools::Binary_node<Contents_SCL<B,R>>* node_curr)
{
	const auto size = (int)node_curr->get_c()->lambda.size();
	const auto size_2 = size / 2;

	const auto *node_left  = node_curr->get_left(); // get left node

	for (auto i = 0; i < size_2; i++)
		node_left->get_c()->lambda[i] = F(node_curr->get_c()->lambda[        i],  // apply f()
		                                  node_curr->get_c()->lambda[size_2 +i]);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::apply_g(const tools::Binary_node<Contents_SCL<B,R>>* node_curr)
{
	const auto size = (int)node_curr->get_c()->lambda.size();
	const auto size_2 = size / 2;

	const auto *node_left  = node_curr->get_left(); // get left node
	const auto *node_right = node_curr->get_right(); // get right node

	for (auto i = 0; i < size_2; i++)
		node_right->get_c()->lambda[i] = G(node_curr->get_c()->lambda[        i], // apply g()
		                                   node_curr->get_c()->lambda[size_2 +i],
		                                   node_left->get_c()->s     [        i]);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive<B,R,F,G>
::compute_sums(const tools::Binary_node<Contents_SCL<B,R>>* node_curr)
{
	const auto size = (int)node_curr->get_c()->lambda.size();
	const auto size_2 = size / 2;

	const auto *node_left  = node_curr->get_left(); // get left node
	const auto *node_right = node_curr->get_right(); // get right node

	for (auto i = 0; i < size_2; i++)
		node_curr->get_c()->s[i] = node_left->get_c()->s[i] ^ node_right->get_c()->s[i]; // bit xor

	for (auto i = 0; i < size_2; i++)
		node_curr->get_c()->s[size_2 +i] = node_right->get_c()->s[i]; // bit eq
}
}
}
