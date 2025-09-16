#include <string>
#include <sstream>

#include "Tools/Exception/exception.hpp"
#include "Module/Module.hpp"
#include "Module/Decoder/Polar/SCL/CRC/Decoder_polar_SCL_naive_CA.hpp"

namespace aff3ct
{
namespace module
{

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
Task& Decoder_polar_SCL_naive_CA<B,R,F,G>
::operator[](const dec::tsk t)
{
	return Module::operator[]((size_t)t);
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
Socket& Decoder_polar_SCL_naive_CA<B,R,F,G>
::operator[](const dec::sck::set_crc_const s)
{
	return Module::operator[]((size_t)dec::tsk::set_crc_const)[(size_t)s];
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
Decoder_polar_SCL_naive_CA<B,R,F,G>
::Decoder_polar_SCL_naive_CA(const int& K, const int& N, const int& L, const std::vector<bool>& frozen_bits,
                             const CRC<B>& crc)
: Decoder_polar_SCL_naive<B,R,F,G>(K, N, L, frozen_bits), crc(crc.clone())
{
	const std::string name = "Decoder_polar_SCL_naive_CA";
	this->set_name(name);

	if (this->crc->get_size() > K)
	{
		std::stringstream message;
		message << "'crc->get_size()' has to be equal or smaller than 'K' ('crc->get_size()' = " << this->crc->get_size()
		        << ", 'K' = " << K << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	auto &p5 = this->create_task("set_crc_const");
	auto p5s_V_crc = this->template create_socket_in<B>(p5, "V_crc", this->crc->get_size());
	this->create_codelet(p5, [p5s_V_crc](Module &m, Task &t, const size_t frame_id) -> int
	{
		auto &dec = static_cast<Decoder_polar_SCL_naive_CA<B,R,F,G>&>(m);

		dec._set_crc_const(static_cast<B*>(t[p5s_V_crc].get_dataptr()), frame_id);

		return status_t::SUCCESS;
	});
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
Decoder_polar_SCL_naive_CA<B,R,F,G>* Decoder_polar_SCL_naive_CA<B,R,F,G>
::clone() const
{
	auto m = new Decoder_polar_SCL_naive_CA(*this);
	m->deep_copy(*this);
	return m;
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive_CA<B,R,F,G>
::deep_copy(const Decoder_polar_SCL_naive_CA<B,R,F,G> &m)
{
	Decoder_polar_SCL_naive<B,R,F,G>::deep_copy(m);
	if (m.crc != nullptr) this->crc.reset(m.crc->clone());
}

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive_CA<B,R,F,G>
::select_best_path(const size_t frame_id)
{
	std::vector<B> U_test;
	// const int crc_size = crc->get_size();
	// std::vector<B> U_test_crc(this->K + crc_size);
	std::set<int> active_paths_before_crc = this->active_paths;
	for (auto path : active_paths_before_crc)
	{
		U_test.clear();
		// U_test_crc.clear();

		for (auto leaf = 0 ; leaf < this->N ; leaf++)
			// if (!this->frozen_bits[leaf])
			// 	U_test.push_back(this->leaves_array[path][leaf]->get_c()->s[0]);
			// The CRC is applied to the whole string
			U_test.push_back(this->leaves_array[path][leaf]->get_c()->s[0]);


		// Use the customized method check_crc_const
		// Note: You MUST set the crc_const by set_crc_const before decoding!
		// Note: Currently check_crc_const only support one frame per wave!
		bool decode_result = this->crc->check_crc_const(U_test, frame_id);
		// Also, we cannot check U_test's crc according to U_test's crc.
		// crc->build(U_test, U_test_crc, frame_id);
		// bool decode_result = crc->check(U_test_crc, frame_id);
		// The original version is obviously wrong, since crc->check takes a K + crc_size input, but U_test is only K bits long.
		// bool decode_result = crc->check(U_test, frame_id);
		if (!decode_result)
			this->active_paths.erase(path);
	}

	// Unfortunately, we HAVE to setup a new socket to return the decode_success flag. 
	// Hence we straightly disable it currently. 
	// If the decode is unsuccessful, then set the decode flag to false 
	if (this->active_paths.size() < 1)
	{
		this->decode_success=false;
	}
	else
	{
		// Don't forget to set the decode flag to true if the decode is successful.
		this->decode_success=true;
	}

	this->Decoder_polar_SCL_naive<B,R,F,G>::select_best_path(frame_id);
}

// template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
// void Decoder_polar_SCL_naive_CA<B,R,F,G>
// ::set_crc_const(const B *V_crc, const size_t frame_id)
// {
	
// }

template <typename B, typename R, tools::proto_f<R> F, tools::proto_g<B,R> G>
void Decoder_polar_SCL_naive_CA<B,R,F,G>
::_set_crc_const(const B *V_crc, const size_t frame_id)
{
	this->crc->set_crc_const(V_crc, frame_id);
}
}
}
