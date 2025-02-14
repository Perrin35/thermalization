OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7409878) q[0];
sx q[0];
rz(2.7317943) q[0];
sx q[0];
rz(8.3538342) q[0];
rz(0.61869705) q[1];
sx q[1];
rz(3.7096042) q[1];
sx q[1];
rz(8.5811442) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7535665) q[0];
sx q[0];
rz(-2.8754398) q[0];
sx q[0];
rz(1.4583336) q[0];
rz(-pi) q[1];
rz(-1.7143519) q[2];
sx q[2];
rz(-1.2782405) q[2];
sx q[2];
rz(2.1168328) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31402097) q[1];
sx q[1];
rz(-2.2306475) q[1];
sx q[1];
rz(0.41484264) q[1];
rz(-pi) q[2];
rz(-2.0655965) q[3];
sx q[3];
rz(-1.1462948) q[3];
sx q[3];
rz(2.692846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7817276) q[2];
sx q[2];
rz(-0.83911506) q[2];
sx q[2];
rz(-0.016999379) q[2];
rz(-1.4012236) q[3];
sx q[3];
rz(-0.56178105) q[3];
sx q[3];
rz(1.9248272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.568999) q[0];
sx q[0];
rz(-1.8091135) q[0];
sx q[0];
rz(-0.78729415) q[0];
rz(0.17838082) q[1];
sx q[1];
rz(-1.2063824) q[1];
sx q[1];
rz(1.23752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84920657) q[0];
sx q[0];
rz(-2.153075) q[0];
sx q[0];
rz(1.6447862) q[0];
x q[1];
rz(2.8493153) q[2];
sx q[2];
rz(-0.9773796) q[2];
sx q[2];
rz(3.127272) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2678821) q[1];
sx q[1];
rz(-2.0036266) q[1];
sx q[1];
rz(-0.3496561) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54333289) q[3];
sx q[3];
rz(-2.6381734) q[3];
sx q[3];
rz(1.0375298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7201207) q[2];
sx q[2];
rz(-1.4482435) q[2];
sx q[2];
rz(3.0778911) q[2];
rz(-1.2740159) q[3];
sx q[3];
rz(-0.84309045) q[3];
sx q[3];
rz(1.564285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.098123) q[0];
sx q[0];
rz(-1.1935357) q[0];
sx q[0];
rz(-2.3980339) q[0];
rz(0.45383635) q[1];
sx q[1];
rz(-1.3498787) q[1];
sx q[1];
rz(-0.41890621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6565257) q[0];
sx q[0];
rz(-0.48975268) q[0];
sx q[0];
rz(-1.9974677) q[0];
rz(-1.9100283) q[2];
sx q[2];
rz(-0.76398173) q[2];
sx q[2];
rz(-2.6126249) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.031098628) q[1];
sx q[1];
rz(-0.95606177) q[1];
sx q[1];
rz(-0.35009274) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7396853) q[3];
sx q[3];
rz(-2.4554376) q[3];
sx q[3];
rz(2.9952733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.087674) q[2];
sx q[2];
rz(-2.3232465) q[2];
sx q[2];
rz(1.2308925) q[2];
rz(-2.6026717) q[3];
sx q[3];
rz(-1.475622) q[3];
sx q[3];
rz(-1.4628791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26381275) q[0];
sx q[0];
rz(-0.75436622) q[0];
sx q[0];
rz(0.60212773) q[0];
rz(1.6944132) q[1];
sx q[1];
rz(-2.4592631) q[1];
sx q[1];
rz(-2.2785861) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38463889) q[0];
sx q[0];
rz(-0.85345972) q[0];
sx q[0];
rz(-1.7216645) q[0];
rz(-2.8406937) q[2];
sx q[2];
rz(-1.931012) q[2];
sx q[2];
rz(0.27733251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26559184) q[1];
sx q[1];
rz(-0.99960281) q[1];
sx q[1];
rz(-0.31293641) q[1];
rz(-pi) q[2];
rz(2.8252271) q[3];
sx q[3];
rz(-0.85815198) q[3];
sx q[3];
rz(2.0848839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0756695) q[2];
sx q[2];
rz(-1.1752335) q[2];
sx q[2];
rz(-2.2389331) q[2];
rz(2.3265808) q[3];
sx q[3];
rz(-0.60324001) q[3];
sx q[3];
rz(0.41316113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0478504) q[0];
sx q[0];
rz(-0.047690064) q[0];
sx q[0];
rz(-0.93511859) q[0];
rz(1.5330261) q[1];
sx q[1];
rz(-2.8159499) q[1];
sx q[1];
rz(-2.8757222) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2926798) q[0];
sx q[0];
rz(-0.32912185) q[0];
sx q[0];
rz(-1.9051244) q[0];
x q[1];
rz(2.3912848) q[2];
sx q[2];
rz(-1.7774095) q[2];
sx q[2];
rz(2.7624102) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7794687) q[1];
sx q[1];
rz(-1.2219795) q[1];
sx q[1];
rz(1.3999983) q[1];
rz(-pi) q[2];
rz(0.10281296) q[3];
sx q[3];
rz(-1.8342557) q[3];
sx q[3];
rz(2.3578532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.461146) q[2];
sx q[2];
rz(-1.6411883) q[2];
sx q[2];
rz(2.6908596) q[2];
rz(-1.9510673) q[3];
sx q[3];
rz(-2.4400986) q[3];
sx q[3];
rz(2.7019971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.180069) q[0];
sx q[0];
rz(-2.838205) q[0];
sx q[0];
rz(0.044483749) q[0];
rz(2.8728409) q[1];
sx q[1];
rz(-2.2046397) q[1];
sx q[1];
rz(0.43089795) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4558894) q[0];
sx q[0];
rz(-2.1834032) q[0];
sx q[0];
rz(1.2489178) q[0];
rz(0.459983) q[2];
sx q[2];
rz(-2.1873133) q[2];
sx q[2];
rz(-3.1107855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2984276) q[1];
sx q[1];
rz(-2.3379605) q[1];
sx q[1];
rz(-0.9982361) q[1];
rz(-0.66029064) q[3];
sx q[3];
rz(-1.5753645) q[3];
sx q[3];
rz(-2.3167603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1162794) q[2];
sx q[2];
rz(-0.37386027) q[2];
sx q[2];
rz(0.6244134) q[2];
rz(-1.1011018) q[3];
sx q[3];
rz(-2.3291984) q[3];
sx q[3];
rz(-2.1541434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2295912) q[0];
sx q[0];
rz(-1.0593375) q[0];
sx q[0];
rz(0.68369317) q[0];
rz(1.4423485) q[1];
sx q[1];
rz(-2.3577299) q[1];
sx q[1];
rz(-0.20194617) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037887427) q[0];
sx q[0];
rz(-1.572084) q[0];
sx q[0];
rz(-1.8989858) q[0];
rz(-0.82345404) q[2];
sx q[2];
rz(-1.5353893) q[2];
sx q[2];
rz(3.0454783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8117042) q[1];
sx q[1];
rz(-2.0913908) q[1];
sx q[1];
rz(1.4938526) q[1];
rz(-pi) q[2];
rz(0.45777623) q[3];
sx q[3];
rz(-1.6094064) q[3];
sx q[3];
rz(0.82070551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.086229) q[2];
sx q[2];
rz(-1.4198885) q[2];
sx q[2];
rz(-1.281338) q[2];
rz(-2.4751439) q[3];
sx q[3];
rz(-2.0446348) q[3];
sx q[3];
rz(2.5950529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9796889) q[0];
sx q[0];
rz(-0.55210102) q[0];
sx q[0];
rz(2.6276278) q[0];
rz(-2.1886096) q[1];
sx q[1];
rz(-2.516808) q[1];
sx q[1];
rz(-2.0228588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361439) q[0];
sx q[0];
rz(-1.4892254) q[0];
sx q[0];
rz(-0.42382999) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0055112) q[2];
sx q[2];
rz(-0.82821199) q[2];
sx q[2];
rz(1.0880053) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9671849) q[1];
sx q[1];
rz(-0.46665114) q[1];
sx q[1];
rz(-1.0789167) q[1];
x q[2];
rz(0.548377) q[3];
sx q[3];
rz(-1.2634687) q[3];
sx q[3];
rz(-0.69068324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.755456) q[2];
sx q[2];
rz(-2.7894207) q[2];
sx q[2];
rz(2.238671) q[2];
rz(1.4593982) q[3];
sx q[3];
rz(-1.7995588) q[3];
sx q[3];
rz(2.3193147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67895401) q[0];
sx q[0];
rz(-2.8215388) q[0];
sx q[0];
rz(-1.1902887) q[0];
rz(-0.32084385) q[1];
sx q[1];
rz(-1.6879993) q[1];
sx q[1];
rz(1.2215337) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71305481) q[0];
sx q[0];
rz(-0.24745169) q[0];
sx q[0];
rz(2.5525981) q[0];
rz(2.3053929) q[2];
sx q[2];
rz(-1.3635907) q[2];
sx q[2];
rz(1.5839603) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62432837) q[1];
sx q[1];
rz(-0.78189865) q[1];
sx q[1];
rz(2.8502591) q[1];
rz(-0.19467312) q[3];
sx q[3];
rz(-1.750573) q[3];
sx q[3];
rz(-2.468022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.90423501) q[2];
sx q[2];
rz(-0.15494896) q[2];
sx q[2];
rz(0.3332738) q[2];
rz(1.0171558) q[3];
sx q[3];
rz(-2.1867496) q[3];
sx q[3];
rz(0.3199544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3453813) q[0];
sx q[0];
rz(-2.6295202) q[0];
sx q[0];
rz(-2.2228125) q[0];
rz(2.9504919) q[1];
sx q[1];
rz(-2.7156576) q[1];
sx q[1];
rz(0.79413116) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0427754) q[0];
sx q[0];
rz(-2.2198027) q[0];
sx q[0];
rz(-0.88943435) q[0];
rz(-pi) q[1];
rz(2.7808519) q[2];
sx q[2];
rz(-1.7684002) q[2];
sx q[2];
rz(0.2195356) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0550784) q[1];
sx q[1];
rz(-0.76755133) q[1];
sx q[1];
rz(-1.961019) q[1];
rz(-pi) q[2];
rz(0.32148949) q[3];
sx q[3];
rz(-1.3703385) q[3];
sx q[3];
rz(-0.050762477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11067757) q[2];
sx q[2];
rz(-2.5258749) q[2];
sx q[2];
rz(-0.75882971) q[2];
rz(-2.2714254) q[3];
sx q[3];
rz(-2.3535959) q[3];
sx q[3];
rz(1.0677392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1398685) q[0];
sx q[0];
rz(-1.682946) q[0];
sx q[0];
rz(-1.2276822) q[0];
rz(2.7036746) q[1];
sx q[1];
rz(-1.7057849) q[1];
sx q[1];
rz(1.5954856) q[1];
rz(0.48702892) q[2];
sx q[2];
rz(-1.9430046) q[2];
sx q[2];
rz(-0.76443048) q[2];
rz(-1.036676) q[3];
sx q[3];
rz(-2.277959) q[3];
sx q[3];
rz(-2.1716519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
