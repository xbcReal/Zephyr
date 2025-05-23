#pragma once

#include <seal/seal.h>
#include <seal/util/uintarith.h>

#include <vector>

using namespace std;
using namespace seal;
using namespace seal::util;

class CKKSEvaluator {
 private:
  double x_l = 1e-4, x_r = 1e3;
  double k2 = 1.6885913944245319, k1 = 0.08385114624023438;
  double x2 = 339.8933847174716, x1 = 0.12670372040311723;
  double pivot = 0.37277489383454;
  double m1 = -0.5 * k2 * pow(x1, -1.5), c1 = 1.5 * k2 * pow(x1, -0.5);
  double m2 = -0.5 * k2 * pow(x2, -1.5), c2 = 1.5 * k2 * pow(x2, -0.5);
  double sgn_factor = 0.5;
  vector<double> F4_COEFFS = {0, 315, 0, -420, 0, 378, 0, -180, 0, 35};
  // should be divided by (1 << 7)
  int f4_scale = (1 << 7);
  vector<double> G4_COEFFS = {0, 5850, 0, -34974, 0, 97015, 0, -113492, 0, 46623};
  // should be divided by (1 << 10)
  int g4_scale = (1 << 10);
  Plaintext a, b, half, div_b, neg_p, err, M1, M2, C1, C2;

  Ciphertext poly_eval(Ciphertext x, vector<Plaintext> coeff);
  Ciphertext newtonIter(Ciphertext x, Ciphertext res, int iter = 4);
  pair<Ciphertext, Ciphertext> goldSchmidtIter(Ciphertext v, Ciphertext y, int d = 1);
  Ciphertext initGuess(Ciphertext x);
  Ciphertext evalLine(Ciphertext x, Plaintext m, Plaintext c);

 public:
  SEALContext *context = nullptr;
  Encryptor *encryptor = nullptr;
  Decryptor *decryptor = nullptr;
  CKKSEncoder *encoder = nullptr;
  Evaluator *evaluator = nullptr;
  RelinKeys *relin_keys = nullptr;
  GaloisKeys *galois_keys = nullptr;

  double scale;
  size_t N;
  size_t slot_count;
  size_t degree;
  size_t comm = 0;
  size_t round = 0;
  std::vector<std::uint32_t> rots;

  CKKSEvaluator(
      SEALContext &context,
      Encryptor &encryptor,
      Decryptor &decryptor,
      CKKSEncoder &encoder,
      Evaluator &evaluator,
      double scale,
      RelinKeys &relin_keys,
      GaloisKeys &galois_keys) {
    this->scale = scale;
    this->encryptor = &encryptor;
    this->decryptor = &decryptor;
    this->encoder = &encoder;
    this->evaluator = &evaluator;
    this->relin_keys = &relin_keys;
    this->galois_keys = &galois_keys;
    this->context = &context;

    N = encoder.slot_count() * 2;
    degree = N;
    slot_count = encoder.slot_count();
    this->encoder->encode(x_l, scale, a);
    this->encoder->encode(x_r, scale, b);
    this->encoder->encode(1 / (x_r - x_l), scale, div_b);
    this->encoder->encode(-pivot, scale, neg_p);
    this->encoder->encode(0.5, scale, half);
    this->encoder->encode(m1, scale, M1);
    this->encoder->encode(m2, scale, M2);
    this->encoder->encode(c1, scale, C1);
    this->encoder->encode(c2, scale, C2);

    for (int i = 0; i < uint(std::ceil(log2(degree)/2)); i++) {
      rots.push_back((degree / exponentiate_uint(2, 2*i))+1);
    rots.push_back((degree / exponentiate_uint(2, 2*i+1))+1);
    rots.push_back((3*degree / exponentiate_uint(2, 2*i+1))+1);
    }

  }

  void re_encrypt(Ciphertext &ct);
  void print_decrypted_ct(Ciphertext &ct, int nums);
  void print_decoded_pt(Plaintext &pt, int num);

  vector<double> init_vec_with_value(int N, double init_value);
  vector<double> init_mask(int N, int m);
  uint64_t get_modulus(Ciphertext &x, int k);
  void eval_odd_deg9_poly(vector<double> &a, Ciphertext &x, Ciphertext &dest);

  Ciphertext invert_sqrt(Ciphertext x, int d_newt = 5, int d_gold = 2);
  Ciphertext sgn_eval(Ciphertext x, int d_g, int d_f, double sgn_factor = 0.5);
  Ciphertext exp(Ciphertext x);
  Ciphertext inverse(Ciphertext x, int iter = 4);

  double calculateMAE(vector<double> &y_true, Ciphertext &ct, int N);
};
