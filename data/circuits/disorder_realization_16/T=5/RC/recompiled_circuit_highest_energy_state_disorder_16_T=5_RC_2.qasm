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
rz(1.6348732) q[0];
sx q[0];
rz(9.2490939) q[0];
rz(-2.076258) q[1];
sx q[1];
rz(-1.330436) q[1];
sx q[1];
rz(-2.4863844) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86533812) q[0];
sx q[0];
rz(-1.0840345) q[0];
sx q[0];
rz(-2.7810762) q[0];
x q[1];
rz(1.5152446) q[2];
sx q[2];
rz(-1.4719988) q[2];
sx q[2];
rz(-2.9227119) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0111078) q[1];
sx q[1];
rz(-1.678651) q[1];
sx q[1];
rz(1.6438831) q[1];
x q[2];
rz(1.9780065) q[3];
sx q[3];
rz(-2.4904479) q[3];
sx q[3];
rz(2.4442087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4483999) q[2];
sx q[2];
rz(-3.0312263) q[2];
sx q[2];
rz(2.8006862) q[2];
rz(0.77346268) q[3];
sx q[3];
rz(-2.1229459) q[3];
sx q[3];
rz(0.10507467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4648436) q[0];
sx q[0];
rz(-0.28443795) q[0];
sx q[0];
rz(-0.25962096) q[0];
rz(-2.3530841) q[1];
sx q[1];
rz(-0.85616833) q[1];
sx q[1];
rz(1.8964918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7126078) q[0];
sx q[0];
rz(-2.1855506) q[0];
sx q[0];
rz(3.0180467) q[0];
rz(-pi) q[1];
rz(-2.8261024) q[2];
sx q[2];
rz(-1.6655388) q[2];
sx q[2];
rz(1.4960932) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81345612) q[1];
sx q[1];
rz(-0.81902796) q[1];
sx q[1];
rz(2.0585895) q[1];
x q[2];
rz(0.29446843) q[3];
sx q[3];
rz(-0.52669062) q[3];
sx q[3];
rz(1.9484438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9022687) q[2];
sx q[2];
rz(-2.6332492) q[2];
sx q[2];
rz(-2.2321205) q[2];
rz(0.66696143) q[3];
sx q[3];
rz(-1.704155) q[3];
sx q[3];
rz(-3.0392569) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42720902) q[0];
sx q[0];
rz(-0.41218555) q[0];
sx q[0];
rz(-1.7538196) q[0];
rz(-0.99675238) q[1];
sx q[1];
rz(-0.69147384) q[1];
sx q[1];
rz(2.6549285) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62138689) q[0];
sx q[0];
rz(-2.7059816) q[0];
sx q[0];
rz(1.9389047) q[0];
rz(2.2816554) q[2];
sx q[2];
rz(-2.268848) q[2];
sx q[2];
rz(2.6280478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4506766) q[1];
sx q[1];
rz(-1.6918536) q[1];
sx q[1];
rz(-0.52074213) q[1];
x q[2];
rz(-2.1975945) q[3];
sx q[3];
rz(-1.6907915) q[3];
sx q[3];
rz(1.7606408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3643058) q[2];
sx q[2];
rz(-1.8824258) q[2];
sx q[2];
rz(-2.6510009) q[2];
rz(-2.3250438) q[3];
sx q[3];
rz(-2.3994763) q[3];
sx q[3];
rz(-2.1013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33646026) q[0];
sx q[0];
rz(-1.2456303) q[0];
sx q[0];
rz(-2.1601987) q[0];
rz(0.40311748) q[1];
sx q[1];
rz(-0.50794452) q[1];
sx q[1];
rz(0.53781646) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6526529) q[0];
sx q[0];
rz(-1.8171628) q[0];
sx q[0];
rz(0.15423488) q[0];
rz(-0.13366416) q[2];
sx q[2];
rz(-2.2870009) q[2];
sx q[2];
rz(1.72124) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0672978) q[1];
sx q[1];
rz(-2.5038669) q[1];
sx q[1];
rz(2.2504351) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8290121) q[3];
sx q[3];
rz(-2.2214032) q[3];
sx q[3];
rz(-1.788511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2526523) q[2];
sx q[2];
rz(-2.6344968) q[2];
sx q[2];
rz(1.0520774) q[2];
rz(2.5893411) q[3];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1110765) q[0];
sx q[0];
rz(-1.3330326) q[0];
sx q[0];
rz(-3.1040177) q[0];
rz(-2.3887718) q[1];
sx q[1];
rz(-1.564874) q[1];
sx q[1];
rz(-0.084223025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5606623) q[0];
sx q[0];
rz(-0.59460527) q[0];
sx q[0];
rz(-2.0243852) q[0];
x q[1];
rz(-1.8657612) q[2];
sx q[2];
rz(-1.2216481) q[2];
sx q[2];
rz(-2.5087506) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4956717) q[1];
sx q[1];
rz(-0.42594566) q[1];
sx q[1];
rz(1.6712214) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4100759) q[3];
sx q[3];
rz(-1.525021) q[3];
sx q[3];
rz(1.9327527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.61421824) q[2];
sx q[2];
rz(-0.97443333) q[2];
sx q[2];
rz(-2.7091889) q[2];
rz(1.6480986) q[3];
sx q[3];
rz(-1.7688388) q[3];
sx q[3];
rz(0.73895946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.768141) q[0];
sx q[0];
rz(-1.5274763) q[0];
sx q[0];
rz(-2.3798808) q[0];
rz(-1.0234045) q[1];
sx q[1];
rz(-0.49097148) q[1];
sx q[1];
rz(0.34058079) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7135431) q[0];
sx q[0];
rz(-2.3556956) q[0];
sx q[0];
rz(-0.12554306) q[0];
rz(-pi) q[1];
rz(0.24575524) q[2];
sx q[2];
rz(-0.57502103) q[2];
sx q[2];
rz(-1.0005282) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5300776) q[1];
sx q[1];
rz(-2.7002091) q[1];
sx q[1];
rz(2.2468021) q[1];
rz(-2.3911796) q[3];
sx q[3];
rz(-1.6831796) q[3];
sx q[3];
rz(3.0284798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1731825) q[2];
sx q[2];
rz(-2.0081655) q[2];
sx q[2];
rz(-0.75562149) q[2];
rz(0.38718811) q[3];
sx q[3];
rz(-1.3914934) q[3];
sx q[3];
rz(-2.3473327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8887535) q[0];
sx q[0];
rz(-0.08789739) q[0];
sx q[0];
rz(1.2095691) q[0];
rz(1.6190489) q[1];
sx q[1];
rz(-0.77170283) q[1];
sx q[1];
rz(1.3873772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98859859) q[0];
sx q[0];
rz(-1.261088) q[0];
sx q[0];
rz(1.214908) q[0];
rz(1.102669) q[2];
sx q[2];
rz(-2.515677) q[2];
sx q[2];
rz(2.4374838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.61664192) q[1];
sx q[1];
rz(-1.9845647) q[1];
sx q[1];
rz(-2.3828808) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8903743) q[3];
sx q[3];
rz(-1.1261903) q[3];
sx q[3];
rz(0.79892677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9428955) q[2];
sx q[2];
rz(-0.2936475) q[2];
sx q[2];
rz(-0.4100619) q[2];
rz(-2.6041218) q[3];
sx q[3];
rz(-1.042807) q[3];
sx q[3];
rz(2.0736096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.033567) q[0];
sx q[0];
rz(-1.2355868) q[0];
sx q[0];
rz(-1.0066475) q[0];
rz(2.0058696) q[1];
sx q[1];
rz(-2.6508811) q[1];
sx q[1];
rz(-2.8699285) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.131529) q[0];
sx q[0];
rz(-1.9926785) q[0];
sx q[0];
rz(-0.58229851) q[0];
rz(-pi) q[1];
rz(-2.2211248) q[2];
sx q[2];
rz(-1.0050541) q[2];
sx q[2];
rz(2.5793902) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9342088) q[1];
sx q[1];
rz(-1.0107688) q[1];
sx q[1];
rz(0.70785849) q[1];
x q[2];
rz(-1.3471421) q[3];
sx q[3];
rz(-2.3021379) q[3];
sx q[3];
rz(2.770693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9378308) q[2];
sx q[2];
rz(-0.93272847) q[2];
sx q[2];
rz(2.0358098) q[2];
rz(1.743099) q[3];
sx q[3];
rz(-0.98413697) q[3];
sx q[3];
rz(2.2891146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6451013) q[0];
sx q[0];
rz(-0.95106769) q[0];
sx q[0];
rz(2.7880461) q[0];
rz(1.722466) q[1];
sx q[1];
rz(-1.6058233) q[1];
sx q[1];
rz(0.062573418) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4434011) q[0];
sx q[0];
rz(-1.0834604) q[0];
sx q[0];
rz(0.62701012) q[0];
rz(-2.563971) q[2];
sx q[2];
rz(-2.2934539) q[2];
sx q[2];
rz(-1.8236782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.313301) q[1];
sx q[1];
rz(-2.5646659) q[1];
sx q[1];
rz(2.5531656) q[1];
rz(-0.11918862) q[3];
sx q[3];
rz(-2.475563) q[3];
sx q[3];
rz(2.5882849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.320437) q[2];
sx q[2];
rz(-1.374036) q[2];
sx q[2];
rz(0.42018166) q[2];
rz(0.34475103) q[3];
sx q[3];
rz(-2.3638066) q[3];
sx q[3];
rz(0.93728089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47278136) q[0];
sx q[0];
rz(-2.9720699) q[0];
sx q[0];
rz(-1.8602759) q[0];
rz(-1.5877113) q[1];
sx q[1];
rz(-0.80931598) q[1];
sx q[1];
rz(-0.70714998) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90072322) q[0];
sx q[0];
rz(-1.5828195) q[0];
sx q[0];
rz(2.2012086) q[0];
x q[1];
rz(-2.1538582) q[2];
sx q[2];
rz(-1.1942004) q[2];
sx q[2];
rz(1.131112) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3465298) q[1];
sx q[1];
rz(-1.3477541) q[1];
sx q[1];
rz(2.8233379) q[1];
rz(2.1955802) q[3];
sx q[3];
rz(-2.1008228) q[3];
sx q[3];
rz(-1.241714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3215434) q[2];
sx q[2];
rz(-1.1300491) q[2];
sx q[2];
rz(-2.6450487) q[2];
rz(1.2921565) q[3];
sx q[3];
rz(-0.29785952) q[3];
sx q[3];
rz(-0.37798247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9847066) q[0];
sx q[0];
rz(-0.26837415) q[0];
sx q[0];
rz(-1.0347086) q[0];
rz(-0.59255076) q[1];
sx q[1];
rz(-0.68065803) q[1];
sx q[1];
rz(-1.5417644) q[1];
rz(1.401027) q[2];
sx q[2];
rz(-2.9063354) q[2];
sx q[2];
rz(-0.34420446) q[2];
rz(2.8205806) q[3];
sx q[3];
rz(-0.62855084) q[3];
sx q[3];
rz(1.4756065) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
