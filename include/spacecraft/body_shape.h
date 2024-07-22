//
// Created by alex_ on 7/13/2024.
//

#ifndef BODY_SHAPE_H
#define BODY_SHAPE_H

#define MAX_VERTS 100     /* maximum number of polyhedral vertices */
#define MAX_FACES 100     /* maximum number of polyhedral faces */
#define MAX_POLYGON_SZ 10 /* maximum number of verts per polygonal face */

#define X 0
#define Y 1
#define Z 2

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

#include <armadillo>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include <fmt/core.h>

#include "math/volume_integrals.h"

namespace naomi::geometry
{

struct face
{
  std::vector<int> verts;
  double w;
  arma::vec3 norm;

  face(std::vector<int> v, const double w, const arma::vec3& normal):
    verts(std::move(v)), w(w), norm(normal){}
};

class body_shape
{
  std::vector<arma::vec3> m_verts;
  std::vector<face> m_faces;
  double m_mass;
  double m_density;
  int A;   /* alpha */
  int B;   /* beta */
  int C;   /* gamma */

  /* projection integrals */
  double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

  /* face integrals */
  double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;

  /* volume integrals */
  double T0 = 0;
  arma::vec3 T1, T2, TP;
public:
  body_shape(const std::vector<arma::vec3>& verts, const std::vector<std::vector<int>>& faces, const double mass, const double density = 1.0)
  {
    m_verts = verts;
    m_mass = mass;
    m_density = density;
    make_polyhedron(faces);
    compute_volume_integrals();
  }

  [[nodiscard]] arma::vec3 get_center_of_mass() const {
    return T1/T0;
  }

  arma::mat33 get_inertia_tensor()
  {
    auto r = get_center_of_mass();
    arma::mat33 J;
    double mass = m_density * T0;

    /* compute inertia tensor */
    J(X, X) = m_density * (T2[Y] + T2[Z]);
    J(Y, Y) = m_density * (T2[Z] + T2[X]);
    J(Z, Z) = m_density * (T2[X] + T2[Y]);
    J(X, Y) = J(Y, X) = - m_density * TP[X];
    J(Y, Z) = J(Z, Y) = - m_density * TP[Y];
    J(Z, X) = J(X, Z) = - m_density * TP[Z];

    /* translate inertia tensor to center of mass */
    J(X, X) -= mass * (r[Y]*r[Y] + r[Z]*r[Z]);
    J(Y, Y) -= mass * (r[Z]*r[Z] + r[X]*r[X]);
    J(Z, Z) -= mass * (r[X]*r[X] + r[Y]*r[Y]);
    J(X, Y) = J(Y, X) += mass * r[X] * r[Y];
    J(Y, Z) = J(Z, Y) += mass * r[Y] * r[Z];
    J(Z, X) = J(X, Z) += mass * r[Z] * r[X];

    return J;
  }

  static body_shape make_rectangle(double l, double w, double h, const double mass, const arma::vec3& origin = {0, 0, 0})
  {
    double a = l/2;
    double b = w/2;
    double c = h/2;

    std::vector<arma::vec3> verts = {
      {-a, -b, -c},
      {a, -b, -c},
      {a, b, -c},
      {-a, b, -c},
      {-a, -b, c},
      {a, -b, c},
      {a, b, c},
      {-a, b, c}
    };

    std::for_each(verts.begin(), verts.end(), [origin](arma::vec3 v){v += origin;});

    std::vector<std::vector<int>> faces = {
      {0, 3, 2, 1},
      {4, 5, 6, 7},
      {0, 1, 5, 4},
      {6, 2, 3, 7},
      {1, 2, 6, 5},
      {0, 4, 7, 3},
    };

    return {verts, faces, mass};
  }

private:
  void make_polyhedron(const std::vector<std::vector<int>>& faces)
  {
    for (auto f : faces) {
      /* compute face normal and offset w from first 3 vertices */
      const arma::vec3 dv1 = m_verts[f[1]] - m_verts[f[0]];
      const arma::vec3 dv2 = m_verts[f[2]] - m_verts[f[1]];
      const arma::vec3 normal_vec = cross(dv1, dv2);
      const arma::vec3 normal = normalise(normal_vec);
      const double w = -dot(normal, m_verts[f[0]]);
      m_faces.emplace_back(f, w, normal);
    }
  }

  /*
     ============================================================================
     compute mass properties
     ============================================================================
  */

  /* compute various integrations over projection of face */
  void compute_projection_integrals(const face& f)
  {
    double a0, a1, da;
    double b0, b1, db;
    double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
    double a1_2, a1_3, b1_2, b1_3;
    double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
    double Cab, Kab, Caab, Kaab, Cabb, Kabb;
    const std::size_t n_verts = f.verts.size();

    P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;

    for (std::size_t i = 0; i < n_verts; i++) {
      a0 = m_verts[f.verts[i]][A];
      b0 = m_verts[f.verts[i]][B];
      a1 = m_verts[f.verts[(i+1) % n_verts]][A];
      b1 = m_verts[f.verts[(i+1) % n_verts]][B];
      da = a1 - a0;
      db = b1 - b0;
      a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
      b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
      a1_2 = a1 * a1; a1_3 = a1_2 * a1;
      b1_2 = b1 * b1; b1_3 = b1_2 * b1;

      C1 = a1 + a0;
      Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
      Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
      Cab = 3*a1_2 + 2*a1*a0 + a0_2; Kab = a1_2 + 2*a1*a0 + 3*a0_2;
      Caab = a0*Cab + 4*a1_3; Kaab = a1*Kab + 4*a0_3;
      Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
      Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

      P1 += db*C1;
      Pa += db*Ca;
      Paa += db*Caa;
      Paaa += db*Caaa;
      Pb += da*Cb;
      Pbb += da*Cbb;
      Pbbb += da*Cbbb;
      Pab += db*(b1*Cab + b0*Kab);
      Paab += db*(b1*Caab + b0*Kaab);
      Pabb += da*(a1*Cabb + a0*Kabb);
    }

    P1 /= 2.0;
    Pa /= 6.0;
    Paa /= 12.0;
    Paaa /= 20.0;
    Pb /= -6.0;
    Pbb /= -12.0;
    Pbbb /= -20.0;
    Pab /= 24.0;
    Paab /= 60.0;
    Pabb /= -60.0;
  }

  void compute_face_integrals(face& f)
  {
    double k1, k2, k3, k4;

    compute_projection_integrals(f);

    const auto w = f.w;
    const auto n = f.norm;
    k1 = 1 / n[C]; k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

    Fa = k1 * Pa;
    Fb = k1 * Pb;
    Fc = -k2 * (n[A]*Pa + n[B]*Pb + w*P1);

    Faa = k1 * Paa;
    Fbb = k1 * Pbb;
    Fcc = k3 * (SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb
           + w*(2*(n[A]*Pa + n[B]*Pb) + w*P1));

    Faaa = k1 * Paaa;
    Fbbb = k1 * Pbbb;
    Fccc = -k4 * (CUBE(n[A])*Paaa + 3*SQR(n[A])*n[B]*Paab
             + 3*n[A]*SQR(n[B])*Pabb + CUBE(n[B])*Pbbb
             + 3*w*(SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb)
             + w*w*(3*(n[A]*Pa + n[B]*Pb) + w*P1));

    Faab = k1 * Paab;
    Fbbc = -k2 * (n[A]*Pabb + n[B]*Pbbb + w*Pbb);
    Fcca = k3 * (SQR(n[A])*Paaa + 2*n[A]*n[B]*Paab + SQR(n[B])*Pabb
           + w*(2*(n[A]*Paa + n[B]*Pab) + w*Pa));
  }

  void compute_volume_integrals()
  {

    for (auto f: m_faces) {
      const double nx = fabs(f.norm[X]);
      const double ny = fabs(f.norm[Y]);
      const double nz = fabs(f.norm[Z]);
      if (nx > ny && nx > nz) C = X;
      else C = ny > nz ? Y : Z;
      A = (C + 1) % 3;
      B = (A + 1) % 3;

      compute_face_integrals(f);

      T0 += f.norm[X] * ((A == X) ? Fa : ((B == X) ? Fb : Fc));

      T1[A] += f.norm[A] * Faa;
      T1[B] += f.norm[B] * Fbb;
      T1[C] += f.norm[C] * Fcc;
      T2[A] += f.norm[A] * Faaa;
      T2[B] += f.norm[B] * Fbbb;
      T2[C] += f.norm[C] * Fccc;
      TP[A] += f.norm[A] * Faab;
      TP[B] += f.norm[B] * Fbbc;
      TP[C] += f.norm[C] * Fcca;
    }
    T1 /= 2.0;
    T2 /= 3.0;
    TP /= 2.0;
  }
};
}

#endif //BODY_SHAPE_H
