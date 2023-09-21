OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(0.24917319) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(-0.5501774) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43139797) q[0];
sx q[0];
rz(-0.22127998) q[0];
sx q[0];
rz(1.6943323) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1044188) q[2];
sx q[2];
rz(-1.1871157) q[2];
sx q[2];
rz(-1.8431611) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8093811) q[1];
sx q[1];
rz(-0.84377938) q[1];
sx q[1];
rz(-2.0903281) q[1];
rz(-pi) q[2];
rz(0.95991858) q[3];
sx q[3];
rz(-1.5023408) q[3];
sx q[3];
rz(-1.6144891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41574079) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(0.63981167) q[2];
rz(0.85302991) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(-0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2663651) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(0.27045989) q[0];
rz(2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.6289904) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0727901) q[0];
sx q[0];
rz(-1.2501984) q[0];
sx q[0];
rz(-0.94706236) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.011629148) q[2];
sx q[2];
rz(-0.28380576) q[2];
sx q[2];
rz(-1.1139882) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0257033) q[1];
sx q[1];
rz(-1.4374226) q[1];
sx q[1];
rz(2.2504836) q[1];
rz(-0.87725957) q[3];
sx q[3];
rz(-0.3650529) q[3];
sx q[3];
rz(3.0960494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9397395) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(-2.0837636) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(3.1159744) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(-0.29552466) q[0];
rz(-0.23513901) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(-2.3957516) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5999818) q[0];
sx q[0];
rz(-2.9020502) q[0];
sx q[0];
rz(-1.177686) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8638641) q[2];
sx q[2];
rz(-0.54364294) q[2];
sx q[2];
rz(3.1090528) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5513788) q[1];
sx q[1];
rz(-0.56400245) q[1];
sx q[1];
rz(-0.53932921) q[1];
rz(-1.3219222) q[3];
sx q[3];
rz(-2.0119152) q[3];
sx q[3];
rz(0.53019023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0535584) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(-0.68874613) q[2];
rz(3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21928366) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(-3.02137) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(2.8682958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1647987) q[0];
sx q[0];
rz(-1.8935888) q[0];
sx q[0];
rz(-3.047909) q[0];
rz(0.64812135) q[2];
sx q[2];
rz(-2.0378049) q[2];
sx q[2];
rz(1.311122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30444333) q[1];
sx q[1];
rz(-1.7177561) q[1];
sx q[1];
rz(-1.8069581) q[1];
x q[2];
rz(1.3917543) q[3];
sx q[3];
rz(-2.2125707) q[3];
sx q[3];
rz(0.87953506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79167241) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(0.50393528) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(-2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85916096) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(-0.082745634) q[0];
rz(2.4619608) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(-0.98714978) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50901978) q[0];
sx q[0];
rz(-0.97826725) q[0];
sx q[0];
rz(-0.49125262) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5613902) q[2];
sx q[2];
rz(-0.9270037) q[2];
sx q[2];
rz(0.26088342) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6189177) q[1];
sx q[1];
rz(-2.4213311) q[1];
sx q[1];
rz(-2.656225) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6721252) q[3];
sx q[3];
rz(-2.7469027) q[3];
sx q[3];
rz(2.8916388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2510898) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(2.9303072) q[2];
rz(-0.42090297) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(-3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(0.791839) q[0];
rz(2.1461398) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(1.8621559) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148221) q[0];
sx q[0];
rz(-1.5855256) q[0];
sx q[0];
rz(1.8526088) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.860414) q[2];
sx q[2];
rz(-2.1589303) q[2];
sx q[2];
rz(-0.16298018) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93950677) q[1];
sx q[1];
rz(-1.0491976) q[1];
sx q[1];
rz(0.034394666) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0002954) q[3];
sx q[3];
rz(-1.391489) q[3];
sx q[3];
rz(-2.3395204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8906158) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(-0.77077579) q[2];
rz(-1.6714913) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8188748) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(3.0859257) q[0];
rz(0.21559134) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(3.1380222) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63821793) q[0];
sx q[0];
rz(-2.5711381) q[0];
sx q[0];
rz(2.4128782) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9622757) q[2];
sx q[2];
rz(-0.85748312) q[2];
sx q[2];
rz(1.5228524) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0842972) q[1];
sx q[1];
rz(-0.97505403) q[1];
sx q[1];
rz(0.36235313) q[1];
x q[2];
rz(-3.0655541) q[3];
sx q[3];
rz(-1.7685316) q[3];
sx q[3];
rz(-0.34464371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5443762) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(2.9428633) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(-2.7238817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6253117) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(1.7238808) q[0];
rz(-2.7334546) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(-0.67869854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40020254) q[0];
sx q[0];
rz(-0.42660248) q[0];
sx q[0];
rz(-2.461117) q[0];
rz(-pi) q[1];
rz(-2.0754201) q[2];
sx q[2];
rz(-0.80427158) q[2];
sx q[2];
rz(-1.7240766) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8933967) q[1];
sx q[1];
rz(-2.0940603) q[1];
sx q[1];
rz(0.072695331) q[1];
rz(2.7160866) q[3];
sx q[3];
rz(-1.7153499) q[3];
sx q[3];
rz(-2.1394465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(0.91782451) q[2];
rz(1.9994036) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98638242) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(0.25979364) q[0];
rz(-0.70867509) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(-2.6146467) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0748358) q[0];
sx q[0];
rz(-1.2196676) q[0];
sx q[0];
rz(-1.8021276) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91295816) q[2];
sx q[2];
rz(-1.7857988) q[2];
sx q[2];
rz(-2.1195597) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37104169) q[1];
sx q[1];
rz(-0.41401225) q[1];
sx q[1];
rz(-0.71055926) q[1];
x q[2];
rz(-1.280904) q[3];
sx q[3];
rz(-1.1416417) q[3];
sx q[3];
rz(-3.1115301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8156585) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(-2.3507067) q[2];
rz(2.7549426) q[3];
sx q[3];
rz(-2.1270042) q[3];
sx q[3];
rz(2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40795046) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(2.1561484) q[0];
rz(0.5685637) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-0.3607761) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4496778) q[0];
sx q[0];
rz(-1.5458108) q[0];
sx q[0];
rz(-1.7382394) q[0];
x q[1];
rz(-3.0214494) q[2];
sx q[2];
rz(-1.3573682) q[2];
sx q[2];
rz(-0.10211589) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39721397) q[1];
sx q[1];
rz(-0.75512868) q[1];
sx q[1];
rz(-2.932764) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1694069) q[3];
sx q[3];
rz(-1.7161887) q[3];
sx q[3];
rz(-1.8190847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(2.1981751) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5744793) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(1.7760989) q[2];
sx q[2];
rz(-1.4751954) q[2];
sx q[2];
rz(1.8196646) q[2];
rz(0.75136649) q[3];
sx q[3];
rz(-1.02117) q[3];
sx q[3];
rz(1.1985967) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
