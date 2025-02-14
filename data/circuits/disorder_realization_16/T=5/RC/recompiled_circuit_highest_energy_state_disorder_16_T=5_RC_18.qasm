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
rz(1.3311812) q[0];
sx q[0];
rz(-1.5067195) q[0];
sx q[0];
rz(0.17568406) q[0];
rz(1.0653347) q[1];
sx q[1];
rz(4.4720286) q[1];
sx q[1];
rz(8.7695697) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53090838) q[0];
sx q[0];
rz(-1.887868) q[0];
sx q[0];
rz(-1.0560587) q[0];
x q[1];
rz(2.6309825) q[2];
sx q[2];
rz(-3.0282927) q[2];
sx q[2];
rz(-2.4093546) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6933999) q[1];
sx q[1];
rz(-1.498135) q[1];
sx q[1];
rz(3.0334515) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8484905) q[3];
sx q[3];
rz(-2.1610073) q[3];
sx q[3];
rz(2.9410998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6931927) q[2];
sx q[2];
rz(-3.0312263) q[2];
sx q[2];
rz(-0.34090647) q[2];
rz(-2.36813) q[3];
sx q[3];
rz(-1.0186467) q[3];
sx q[3];
rz(-0.10507467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4648436) q[0];
sx q[0];
rz(-2.8571547) q[0];
sx q[0];
rz(-0.25962096) q[0];
rz(-2.3530841) q[1];
sx q[1];
rz(-0.85616833) q[1];
sx q[1];
rz(-1.2451008) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9282824) q[0];
sx q[0];
rz(-1.4699555) q[0];
sx q[0];
rz(-2.1891602) q[0];
x q[1];
rz(-1.4711667) q[2];
sx q[2];
rz(-1.2567695) q[2];
sx q[2];
rz(3.0977566) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.40971995) q[1];
sx q[1];
rz(-1.9202246) q[1];
sx q[1];
rz(-0.81373377) q[1];
rz(-1.7379846) q[3];
sx q[3];
rz(-2.072635) q[3];
sx q[3];
rz(1.5305647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9022687) q[2];
sx q[2];
rz(-2.6332492) q[2];
sx q[2];
rz(0.90947214) q[2];
rz(2.4746312) q[3];
sx q[3];
rz(-1.704155) q[3];
sx q[3];
rz(-0.10233574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7143836) q[0];
sx q[0];
rz(-0.41218555) q[0];
sx q[0];
rz(-1.3877731) q[0];
rz(-0.99675238) q[1];
sx q[1];
rz(-2.4501188) q[1];
sx q[1];
rz(-2.6549285) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5285641) q[0];
sx q[0];
rz(-1.7232304) q[0];
sx q[0];
rz(1.9804762) q[0];
x q[1];
rz(-0.85993729) q[2];
sx q[2];
rz(-0.87274466) q[2];
sx q[2];
rz(-2.6280478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0906252) q[1];
sx q[1];
rz(-1.0542467) q[1];
sx q[1];
rz(-1.4314639) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14778845) q[3];
sx q[3];
rz(-0.94919616) q[3];
sx q[3];
rz(0.10336598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77728689) q[2];
sx q[2];
rz(-1.2591668) q[2];
sx q[2];
rz(0.49059179) q[2];
rz(0.81654882) q[3];
sx q[3];
rz(-0.74211636) q[3];
sx q[3];
rz(-1.0402927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8051324) q[0];
sx q[0];
rz(-1.8959624) q[0];
sx q[0];
rz(-0.98139393) q[0];
rz(-2.7384752) q[1];
sx q[1];
rz(-2.6336481) q[1];
sx q[1];
rz(-0.53781646) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4889398) q[0];
sx q[0];
rz(-1.8171628) q[0];
sx q[0];
rz(0.15423488) q[0];
rz(-pi) q[1];
rz(-2.2914406) q[2];
sx q[2];
rz(-1.4701029) q[2];
sx q[2];
rz(-3.0792011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.28603468) q[1];
sx q[1];
rz(-1.0893309) q[1];
sx q[1];
rz(-0.43586327) q[1];
rz(-pi) q[2];
rz(1.9547525) q[3];
sx q[3];
rz(-0.71184413) q[3];
sx q[3];
rz(1.8432338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8889403) q[2];
sx q[2];
rz(-2.6344968) q[2];
sx q[2];
rz(-1.0520774) q[2];
rz(-0.55225152) q[3];
sx q[3];
rz(-1.735732) q[3];
sx q[3];
rz(-1.9210057) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1110765) q[0];
sx q[0];
rz(-1.3330326) q[0];
sx q[0];
rz(0.037574969) q[0];
rz(-0.75282085) q[1];
sx q[1];
rz(-1.564874) q[1];
sx q[1];
rz(0.084223025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60606979) q[0];
sx q[0];
rz(-1.3227934) q[0];
sx q[0];
rz(2.1169784) q[0];
rz(2.7780121) q[2];
sx q[2];
rz(-1.8474794) q[2];
sx q[2];
rz(-0.83438736) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3854691) q[1];
sx q[1];
rz(-1.9944571) q[1];
sx q[1];
rz(-3.0961354) q[1];
rz(-pi) q[2];
rz(1.7315167) q[3];
sx q[3];
rz(-1.525021) q[3];
sx q[3];
rz(1.9327527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.61421824) q[2];
sx q[2];
rz(-0.97443333) q[2];
sx q[2];
rz(2.7091889) q[2];
rz(1.6480986) q[3];
sx q[3];
rz(-1.3727539) q[3];
sx q[3];
rz(2.4026332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37345165) q[0];
sx q[0];
rz(-1.6141163) q[0];
sx q[0];
rz(-0.76171184) q[0];
rz(-1.0234045) q[1];
sx q[1];
rz(-2.6506212) q[1];
sx q[1];
rz(-0.34058079) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7135431) q[0];
sx q[0];
rz(-0.78589702) q[0];
sx q[0];
rz(0.12554306) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7271785) q[2];
sx q[2];
rz(-2.1264653) q[2];
sx q[2];
rz(-1.8506236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5300776) q[1];
sx q[1];
rz(-2.7002091) q[1];
sx q[1];
rz(-2.2468021) q[1];
rz(-pi) q[2];
rz(-1.4177001) q[3];
sx q[3];
rz(-2.3153437) q[3];
sx q[3];
rz(1.5797256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9684101) q[2];
sx q[2];
rz(-1.1334271) q[2];
sx q[2];
rz(-2.3859712) q[2];
rz(-0.38718811) q[3];
sx q[3];
rz(-1.3914934) q[3];
sx q[3];
rz(-0.79425991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25283915) q[0];
sx q[0];
rz(-0.08789739) q[0];
sx q[0];
rz(1.2095691) q[0];
rz(1.6190489) q[1];
sx q[1];
rz(-2.3698898) q[1];
sx q[1];
rz(1.7542155) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6722058) q[0];
sx q[0];
rz(-1.9090561) q[0];
sx q[0];
rz(-2.8126008) q[0];
rz(-pi) q[1];
rz(2.1437239) q[2];
sx q[2];
rz(-1.838316) q[2];
sx q[2];
rz(-0.47779412) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55303326) q[1];
sx q[1];
rz(-0.84399763) q[1];
sx q[1];
rz(-0.56808205) q[1];
x q[2];
rz(-2.6764836) q[3];
sx q[3];
rz(-1.8583663) q[3];
sx q[3];
rz(0.63048922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1986971) q[2];
sx q[2];
rz(-2.8479452) q[2];
sx q[2];
rz(-2.7315308) q[2];
rz(2.6041218) q[3];
sx q[3];
rz(-1.042807) q[3];
sx q[3];
rz(-2.0736096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.033567) q[0];
sx q[0];
rz(-1.9060059) q[0];
sx q[0];
rz(1.0066475) q[0];
rz(2.0058696) q[1];
sx q[1];
rz(-0.4907116) q[1];
sx q[1];
rz(2.8699285) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.131529) q[0];
sx q[0];
rz(-1.9926785) q[0];
sx q[0];
rz(-0.58229851) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2211248) q[2];
sx q[2];
rz(-1.0050541) q[2];
sx q[2];
rz(2.5793902) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9342088) q[1];
sx q[1];
rz(-1.0107688) q[1];
sx q[1];
rz(-2.4337342) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24231301) q[3];
sx q[3];
rz(-2.3829082) q[3];
sx q[3];
rz(3.098947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2037619) q[2];
sx q[2];
rz(-0.93272847) q[2];
sx q[2];
rz(2.0358098) q[2];
rz(-1.3984937) q[3];
sx q[3];
rz(-2.1574557) q[3];
sx q[3];
rz(0.85247803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6451013) q[0];
sx q[0];
rz(-0.95106769) q[0];
sx q[0];
rz(2.7880461) q[0];
rz(-1.4191267) q[1];
sx q[1];
rz(-1.5357693) q[1];
sx q[1];
rz(3.0790192) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69819151) q[0];
sx q[0];
rz(-1.0834604) q[0];
sx q[0];
rz(-0.62701012) q[0];
rz(2.1252421) q[2];
sx q[2];
rz(-0.89134634) q[2];
sx q[2];
rz(-1.045595) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3741604) q[1];
sx q[1];
rz(-1.2632152) q[1];
sx q[1];
rz(2.6454283) q[1];
rz(0.66257368) q[3];
sx q[3];
rz(-1.6443313) q[3];
sx q[3];
rz(-0.92360332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.320437) q[2];
sx q[2];
rz(-1.374036) q[2];
sx q[2];
rz(2.721411) q[2];
rz(-0.34475103) q[3];
sx q[3];
rz(-0.77778608) q[3];
sx q[3];
rz(-2.2043118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.47278136) q[0];
sx q[0];
rz(-2.9720699) q[0];
sx q[0];
rz(1.2813168) q[0];
rz(-1.5538813) q[1];
sx q[1];
rz(-0.80931598) q[1];
sx q[1];
rz(0.70714998) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90072322) q[0];
sx q[0];
rz(-1.5828195) q[0];
sx q[0];
rz(-2.2012086) q[0];
rz(-2.1936839) q[2];
sx q[2];
rz(-2.4595408) q[2];
sx q[2];
rz(-2.1932604) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36701074) q[1];
sx q[1];
rz(-0.38644192) q[1];
sx q[1];
rz(-2.5143753) q[1];
rz(-pi) q[2];
rz(-0.78451802) q[3];
sx q[3];
rz(-0.79567474) q[3];
sx q[3];
rz(-2.8592542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.82004929) q[2];
sx q[2];
rz(-2.0115435) q[2];
sx q[2];
rz(-0.496544) q[2];
rz(1.8494362) q[3];
sx q[3];
rz(-0.29785952) q[3];
sx q[3];
rz(-2.7636102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-0.15688607) q[0];
sx q[0];
rz(-0.26837415) q[0];
sx q[0];
rz(-1.0347086) q[0];
rz(0.59255076) q[1];
sx q[1];
rz(-2.4609346) q[1];
sx q[1];
rz(1.5998283) q[1];
rz(-1.8027923) q[2];
sx q[2];
rz(-1.6101888) q[2];
sx q[2];
rz(-2.0801795) q[2];
rz(2.5377688) q[3];
sx q[3];
rz(-1.3841938) q[3];
sx q[3];
rz(-2.9740372) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
