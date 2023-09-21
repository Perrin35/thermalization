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
rz(-2.8924195) q[0];
rz(-0.063440032) q[1];
sx q[1];
rz(-2.1698706) q[1];
sx q[1];
rz(0.5501774) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7101947) q[0];
sx q[0];
rz(-2.9203127) q[0];
sx q[0];
rz(1.6943323) q[0];
rz(0.43843856) q[2];
sx q[2];
rz(-2.0619832) q[2];
sx q[2];
rz(-3.0868798) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3783592) q[1];
sx q[1];
rz(-2.2765056) q[1];
sx q[1];
rz(-2.6325429) q[1];
rz(-pi) q[2];
rz(2.1816741) q[3];
sx q[3];
rz(-1.5023408) q[3];
sx q[3];
rz(1.6144891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7258519) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(-2.501781) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8752276) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(0.27045989) q[0];
rz(-0.71331435) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.6289904) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0727901) q[0];
sx q[0];
rz(-1.8913942) q[0];
sx q[0];
rz(-2.1945303) q[0];
rz(-pi) q[1];
rz(-1.5674044) q[2];
sx q[2];
rz(-1.8545824) q[2];
sx q[2];
rz(-2.0397182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0257033) q[1];
sx q[1];
rz(-1.70417) q[1];
sx q[1];
rz(0.89110903) q[1];
rz(-pi) q[2];
rz(-0.87725957) q[3];
sx q[3];
rz(-0.3650529) q[3];
sx q[3];
rz(-0.045543268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9397395) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(-2.3382323) q[2];
rz(-2.0837636) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(-2.846068) q[0];
rz(-0.23513901) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(-2.3957516) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5416109) q[0];
sx q[0];
rz(-2.9020502) q[0];
sx q[0];
rz(1.177686) q[0];
rz(-pi) q[1];
rz(1.7350115) q[2];
sx q[2];
rz(-2.0914372) q[2];
sx q[2];
rz(2.8525713) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44887603) q[1];
sx q[1];
rz(-1.8489031) q[1];
sx q[1];
rz(-2.6443308) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3219222) q[3];
sx q[3];
rz(-2.0119152) q[3];
sx q[3];
rz(0.53019023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.088034257) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(-2.4528465) q[2];
rz(3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.922309) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(-3.02137) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(-2.8682958) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8768339) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(-1.8434974) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69584537) q[2];
sx q[2];
rz(-0.77866422) q[2];
sx q[2];
rz(2.3455182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8371493) q[1];
sx q[1];
rz(-1.4238365) q[1];
sx q[1];
rz(1.8069581) q[1];
x q[2];
rz(2.9076505) q[3];
sx q[3];
rz(-2.4787239) q[3];
sx q[3];
rz(2.5556504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(-0.50393528) q[2];
rz(0.079581633) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85916096) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(3.058847) q[0];
rz(2.4619608) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(-0.98714978) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3521096) q[0];
sx q[0];
rz(-1.1687359) q[0];
sx q[0];
rz(-0.91870086) q[0];
rz(2.4977788) q[2];
sx q[2];
rz(-1.5632731) q[2];
sx q[2];
rz(1.8260337) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.134569) q[1];
sx q[1];
rz(-2.1937074) q[1];
sx q[1];
rz(1.9593777) q[1];
x q[2];
rz(-1.6721252) q[3];
sx q[3];
rz(-2.7469027) q[3];
sx q[3];
rz(-0.24995382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2510898) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(-2.9303072) q[2];
rz(2.7206897) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(-0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.485065) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(0.791839) q[0];
rz(-0.99545288) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(-1.8621559) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6364481) q[0];
sx q[0];
rz(-0.2821869) q[0];
sx q[0];
rz(1.6237153) q[0];
x q[1];
rz(2.5336669) q[2];
sx q[2];
rz(-1.330901) q[2];
sx q[2];
rz(1.243967) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8705604) q[1];
sx q[1];
rz(-2.6189657) q[1];
sx q[1];
rz(1.51103) q[1];
rz(0.1412973) q[3];
sx q[3];
rz(-1.391489) q[3];
sx q[3];
rz(-0.80207223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8906158) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(2.3708169) q[2];
rz(-1.4701014) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8188748) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(0.055667002) q[0];
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
rz(-2.9650426) q[0];
sx q[0];
rz(-1.1561484) q[0];
sx q[0];
rz(1.9745757) q[0];
rz(-pi) q[1];
rz(2.3891719) q[2];
sx q[2];
rz(-1.2781029) q[2];
sx q[2];
rz(2.9257286) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0842972) q[1];
sx q[1];
rz(-0.97505403) q[1];
sx q[1];
rz(0.36235313) q[1];
x q[2];
rz(0.076038578) q[3];
sx q[3];
rz(-1.3730611) q[3];
sx q[3];
rz(0.34464371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59721649) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(-2.9428633) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(-0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6253117) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(1.7238808) q[0];
rz(2.7334546) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(-2.4628941) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40020254) q[0];
sx q[0];
rz(-2.7149902) q[0];
sx q[0];
rz(0.68047561) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0754201) q[2];
sx q[2];
rz(-2.3373211) q[2];
sx q[2];
rz(1.417516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2862257) q[1];
sx q[1];
rz(-1.6337506) q[1];
sx q[1];
rz(1.0463868) q[1];
rz(-pi) q[2];
rz(-1.4123165) q[3];
sx q[3];
rz(-1.9915808) q[3];
sx q[3];
rz(2.5077523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9644908) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(0.91782451) q[2];
rz(-1.142189) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1552102) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(2.881799) q[0];
rz(0.70867509) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(-0.52694595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6078867) q[0];
sx q[0];
rz(-2.7237646) q[0];
sx q[0];
rz(2.582344) q[0];
rz(0.2692659) q[2];
sx q[2];
rz(-2.2109647) q[2];
sx q[2];
rz(0.38538853) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2745797) q[1];
sx q[1];
rz(-1.8362987) q[1];
sx q[1];
rz(2.8200779) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55836375) q[3];
sx q[3];
rz(-0.51279587) q[3];
sx q[3];
rz(0.59190291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32593411) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(-2.3507067) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-2.1270042) q[3];
sx q[3];
rz(0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336422) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(0.98544425) q[0];
rz(-0.5685637) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(2.7808166) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4496778) q[0];
sx q[0];
rz(-1.5957818) q[0];
sx q[0];
rz(-1.7382394) q[0];
rz(-3.0214494) q[2];
sx q[2];
rz(-1.7842245) q[2];
sx q[2];
rz(-3.0394768) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4611778) q[1];
sx q[1];
rz(-0.83592452) q[1];
sx q[1];
rz(1.3780891) q[1];
rz(1.4233227) q[3];
sx q[3];
rz(-1.7383988) q[3];
sx q[3];
rz(2.9180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(2.1981751) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(-1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5744793) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.3600596) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(2.0104682) q[2];
sx q[2];
rz(-0.22618539) q[2];
sx q[2];
rz(-2.4629081) q[2];
rz(0.73137024) q[3];
sx q[3];
rz(-0.89805713) q[3];
sx q[3];
rz(-3.0039207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
