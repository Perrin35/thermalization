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
rz(3.1130528) q[0];
sx q[0];
rz(-0.8172577) q[0];
sx q[0];
rz(0.28491268) q[0];
rz(-2.3218396) q[1];
sx q[1];
rz(-2.2567891) q[1];
sx q[1];
rz(1.1362145) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65596555) q[0];
sx q[0];
rz(-2.1871217) q[0];
sx q[0];
rz(0.095353145) q[0];
rz(-pi) q[1];
rz(-0.66464632) q[2];
sx q[2];
rz(-0.79865361) q[2];
sx q[2];
rz(-2.2448774) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6416671) q[1];
sx q[1];
rz(-1.8537775) q[1];
sx q[1];
rz(2.6989515) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6993027) q[3];
sx q[3];
rz(-1.7768404) q[3];
sx q[3];
rz(-1.1161436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7457486) q[2];
sx q[2];
rz(-2.6555847) q[2];
sx q[2];
rz(-0.29224545) q[2];
rz(0.62119836) q[3];
sx q[3];
rz(-1.0750333) q[3];
sx q[3];
rz(-0.21875374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.578823) q[0];
sx q[0];
rz(-1.1172453) q[0];
sx q[0];
rz(-0.086409464) q[0];
rz(-2.7482765) q[1];
sx q[1];
rz(-1.4540648) q[1];
sx q[1];
rz(-1.6578065) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4230712) q[0];
sx q[0];
rz(-2.4212061) q[0];
sx q[0];
rz(-1.9299401) q[0];
rz(3.0034054) q[2];
sx q[2];
rz(-1.8805247) q[2];
sx q[2];
rz(-0.85339862) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68951974) q[1];
sx q[1];
rz(-1.7240925) q[1];
sx q[1];
rz(0.91233715) q[1];
rz(0.57043907) q[3];
sx q[3];
rz(-1.3469056) q[3];
sx q[3];
rz(2.4233832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94574404) q[2];
sx q[2];
rz(-1.5986779) q[2];
sx q[2];
rz(0.0010541218) q[2];
rz(-0.70683181) q[3];
sx q[3];
rz(-2.2471434) q[3];
sx q[3];
rz(0.27757525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23866776) q[0];
sx q[0];
rz(-0.9280197) q[0];
sx q[0];
rz(2.6631885) q[0];
rz(-2.2967285) q[1];
sx q[1];
rz(-1.076315) q[1];
sx q[1];
rz(-2.1873651) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060392875) q[0];
sx q[0];
rz(-2.0234479) q[0];
sx q[0];
rz(1.9863818) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61561959) q[2];
sx q[2];
rz(-0.93608817) q[2];
sx q[2];
rz(2.8826432) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.79151691) q[1];
sx q[1];
rz(-1.9261987) q[1];
sx q[1];
rz(2.3018964) q[1];
rz(-pi) q[2];
rz(1.30458) q[3];
sx q[3];
rz(-2.2053218) q[3];
sx q[3];
rz(-0.75283829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93762952) q[2];
sx q[2];
rz(-2.8014247) q[2];
sx q[2];
rz(3.0703239) q[2];
rz(-2.1313306) q[3];
sx q[3];
rz(-1.9805757) q[3];
sx q[3];
rz(-2.8205813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6846652) q[0];
sx q[0];
rz(-1.2926956) q[0];
sx q[0];
rz(2.7780823) q[0];
rz(0.89206308) q[1];
sx q[1];
rz(-1.0328247) q[1];
sx q[1];
rz(1.2537289) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1680856) q[0];
sx q[0];
rz(-2.0225443) q[0];
sx q[0];
rz(-2.7049881) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2066098) q[2];
sx q[2];
rz(-1.4763583) q[2];
sx q[2];
rz(2.0409453) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2334785) q[1];
sx q[1];
rz(-1.6122979) q[1];
sx q[1];
rz(0.25021489) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8785973) q[3];
sx q[3];
rz(-1.1914763) q[3];
sx q[3];
rz(2.5130481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7518647) q[2];
sx q[2];
rz(-2.3675297) q[2];
sx q[2];
rz(1.5285726) q[2];
rz(-2.6914237) q[3];
sx q[3];
rz(-1.3521786) q[3];
sx q[3];
rz(-1.2990797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12440974) q[0];
sx q[0];
rz(-0.15394177) q[0];
sx q[0];
rz(0.85087878) q[0];
rz(1.7393913) q[1];
sx q[1];
rz(-1.5910999) q[1];
sx q[1];
rz(0.15241399) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2172594) q[0];
sx q[0];
rz(-1.3314864) q[0];
sx q[0];
rz(-3.13525) q[0];
rz(-pi) q[1];
rz(-1.1609116) q[2];
sx q[2];
rz(-2.4992547) q[2];
sx q[2];
rz(-2.040524) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20402699) q[1];
sx q[1];
rz(-1.7364677) q[1];
sx q[1];
rz(-1.284164) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7431454) q[3];
sx q[3];
rz(-0.51880793) q[3];
sx q[3];
rz(1.102953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5695345) q[2];
sx q[2];
rz(-2.2172838) q[2];
sx q[2];
rz(-2.7657236) q[2];
rz(2.87319) q[3];
sx q[3];
rz(-1.0207876) q[3];
sx q[3];
rz(1.0387897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1035476) q[0];
sx q[0];
rz(-2.2232942) q[0];
sx q[0];
rz(-0.12575664) q[0];
rz(3.0962931) q[1];
sx q[1];
rz(-0.89729512) q[1];
sx q[1];
rz(-0.50517857) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7369469) q[0];
sx q[0];
rz(-2.8789799) q[0];
sx q[0];
rz(0.24554952) q[0];
rz(-pi) q[1];
rz(1.7749106) q[2];
sx q[2];
rz(-2.6417126) q[2];
sx q[2];
rz(1.7244888) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69174978) q[1];
sx q[1];
rz(-2.3556752) q[1];
sx q[1];
rz(-2.9346026) q[1];
rz(-pi) q[2];
rz(-1.1041548) q[3];
sx q[3];
rz(-1.8710051) q[3];
sx q[3];
rz(1.2128346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3230285) q[2];
sx q[2];
rz(-2.8105141) q[2];
sx q[2];
rz(0.25311145) q[2];
rz(-2.3969635) q[3];
sx q[3];
rz(-1.5537477) q[3];
sx q[3];
rz(0.26960069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44148663) q[0];
sx q[0];
rz(-1.8550669) q[0];
sx q[0];
rz(0.043524608) q[0];
rz(-1.9684017) q[1];
sx q[1];
rz(-1.2930608) q[1];
sx q[1];
rz(2.3399369) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64555135) q[0];
sx q[0];
rz(-1.6694549) q[0];
sx q[0];
rz(2.990574) q[0];
rz(-pi) q[1];
x q[1];
rz(2.413676) q[2];
sx q[2];
rz(-0.9093098) q[2];
sx q[2];
rz(-0.9582516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3139627) q[1];
sx q[1];
rz(-0.68601841) q[1];
sx q[1];
rz(1.3242701) q[1];
x q[2];
rz(-1.3255865) q[3];
sx q[3];
rz(-1.86487) q[3];
sx q[3];
rz(-0.40286703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7511071) q[2];
sx q[2];
rz(-1.0932873) q[2];
sx q[2];
rz(0.92457479) q[2];
rz(-0.081192406) q[3];
sx q[3];
rz(-1.743318) q[3];
sx q[3];
rz(2.413522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8586332) q[0];
sx q[0];
rz(-1.2456243) q[0];
sx q[0];
rz(1.0754732) q[0];
rz(1.7643499) q[1];
sx q[1];
rz(-1.7694764) q[1];
sx q[1];
rz(-0.65518728) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7542962) q[0];
sx q[0];
rz(-1.6248024) q[0];
sx q[0];
rz(-1.61458) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2776112) q[2];
sx q[2];
rz(-0.80481968) q[2];
sx q[2];
rz(-0.41560666) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.97444087) q[1];
sx q[1];
rz(-2.2658684) q[1];
sx q[1];
rz(1.6593169) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.029441802) q[3];
sx q[3];
rz(-1.2323151) q[3];
sx q[3];
rz(-1.4127617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36771557) q[2];
sx q[2];
rz(-2.750062) q[2];
sx q[2];
rz(-1.9850622) q[2];
rz(-1.997442) q[3];
sx q[3];
rz(-0.63026989) q[3];
sx q[3];
rz(1.8088079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3211806) q[0];
sx q[0];
rz(-0.92083609) q[0];
sx q[0];
rz(-1.17571) q[0];
rz(-1.3395478) q[1];
sx q[1];
rz(-2.1422155) q[1];
sx q[1];
rz(2.7340926) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9391607) q[0];
sx q[0];
rz(-2.6337886) q[0];
sx q[0];
rz(-2.5277471) q[0];
rz(-pi) q[1];
rz(2.6596498) q[2];
sx q[2];
rz(-0.53496581) q[2];
sx q[2];
rz(2.9819161) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4983771) q[1];
sx q[1];
rz(-0.98340323) q[1];
sx q[1];
rz(0.024557928) q[1];
rz(2.8698233) q[3];
sx q[3];
rz(-1.5021694) q[3];
sx q[3];
rz(2.8420699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0672062) q[2];
sx q[2];
rz(-1.9830474) q[2];
sx q[2];
rz(-2.1161533) q[2];
rz(-2.6440559) q[3];
sx q[3];
rz(-2.189744) q[3];
sx q[3];
rz(2.735125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936754) q[0];
sx q[0];
rz(-1.1578639) q[0];
sx q[0];
rz(-1.7561308) q[0];
rz(-1.0476184) q[1];
sx q[1];
rz(-1.7966929) q[1];
sx q[1];
rz(-2.1687543) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0467918) q[0];
sx q[0];
rz(-2.5908906) q[0];
sx q[0];
rz(1.7355376) q[0];
rz(-2.3054584) q[2];
sx q[2];
rz(-2.4651066) q[2];
sx q[2];
rz(-1.033178) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7501513) q[1];
sx q[1];
rz(-1.4285993) q[1];
sx q[1];
rz(2.7846257) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7601401) q[3];
sx q[3];
rz(-0.90899639) q[3];
sx q[3];
rz(1.0009777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4086548) q[2];
sx q[2];
rz(-2.813377) q[2];
sx q[2];
rz(1.1019863) q[2];
rz(1.7105626) q[3];
sx q[3];
rz(-2.3165063) q[3];
sx q[3];
rz(-1.2380606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2796606) q[0];
sx q[0];
rz(-1.4623549) q[0];
sx q[0];
rz(-2.097492) q[0];
rz(-2.4109789) q[1];
sx q[1];
rz(-2.5082671) q[1];
sx q[1];
rz(-2.3309753) q[1];
rz(-1.8865449) q[2];
sx q[2];
rz(-0.98636711) q[2];
sx q[2];
rz(-1.2489088) q[2];
rz(-1.295119) q[3];
sx q[3];
rz(-1.2296005) q[3];
sx q[3];
rz(-2.8212764) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
