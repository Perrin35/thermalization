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
rz(-2.076258) q[1];
sx q[1];
rz(-1.330436) q[1];
sx q[1];
rz(0.65520823) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86533812) q[0];
sx q[0];
rz(-2.0575581) q[0];
sx q[0];
rz(-2.7810762) q[0];
rz(-pi) q[1];
rz(1.6263481) q[2];
sx q[2];
rz(-1.6695938) q[2];
sx q[2];
rz(0.21888079) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.44819278) q[1];
sx q[1];
rz(-1.6434577) q[1];
sx q[1];
rz(-3.0334515) q[1];
x q[2];
rz(-1.1635861) q[3];
sx q[3];
rz(-2.4904479) q[3];
sx q[3];
rz(2.4442087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6931927) q[2];
sx q[2];
rz(-0.11036631) q[2];
sx q[2];
rz(-2.8006862) q[2];
rz(-2.36813) q[3];
sx q[3];
rz(-1.0186467) q[3];
sx q[3];
rz(-0.10507467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(2.3530841) q[1];
sx q[1];
rz(-0.85616833) q[1];
sx q[1];
rz(-1.8964918) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9282824) q[0];
sx q[0];
rz(-1.6716372) q[0];
sx q[0];
rz(0.95243246) q[0];
x q[1];
rz(0.31549021) q[2];
sx q[2];
rz(-1.6655388) q[2];
sx q[2];
rz(1.4960932) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7318727) q[1];
sx q[1];
rz(-1.9202246) q[1];
sx q[1];
rz(-2.3278589) q[1];
rz(-0.29446843) q[3];
sx q[3];
rz(-2.614902) q[3];
sx q[3];
rz(-1.1931488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9022687) q[2];
sx q[2];
rz(-2.6332492) q[2];
sx q[2];
rz(-0.90947214) q[2];
rz(0.66696143) q[3];
sx q[3];
rz(-1.4374377) q[3];
sx q[3];
rz(3.0392569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7143836) q[0];
sx q[0];
rz(-0.41218555) q[0];
sx q[0];
rz(1.7538196) q[0];
rz(0.99675238) q[1];
sx q[1];
rz(-2.4501188) q[1];
sx q[1];
rz(-0.48666418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1179817) q[0];
sx q[0];
rz(-1.9754462) q[0];
sx q[0];
rz(0.16594529) q[0];
x q[1];
rz(-2.2816554) q[2];
sx q[2];
rz(-2.268848) q[2];
sx q[2];
rz(0.51354487) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0906252) q[1];
sx q[1];
rz(-1.0542467) q[1];
sx q[1];
rz(-1.7101287) q[1];
x q[2];
rz(0.9439982) q[3];
sx q[3];
rz(-1.4508012) q[3];
sx q[3];
rz(-1.7606408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77728689) q[2];
sx q[2];
rz(-1.8824258) q[2];
sx q[2];
rz(-2.6510009) q[2];
rz(2.3250438) q[3];
sx q[3];
rz(-2.3994763) q[3];
sx q[3];
rz(-1.0402927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33646026) q[0];
sx q[0];
rz(-1.2456303) q[0];
sx q[0];
rz(2.1601987) q[0];
rz(0.40311748) q[1];
sx q[1];
rz(-0.50794452) q[1];
sx q[1];
rz(0.53781646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6526529) q[0];
sx q[0];
rz(-1.8171628) q[0];
sx q[0];
rz(-0.15423488) q[0];
rz(-pi) q[1];
rz(0.13366416) q[2];
sx q[2];
rz(-2.2870009) q[2];
sx q[2];
rz(-1.72124) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.855558) q[1];
sx q[1];
rz(-2.0522617) q[1];
sx q[1];
rz(0.43586327) q[1];
x q[2];
rz(0.31258055) q[3];
sx q[3];
rz(-0.9201895) q[3];
sx q[3];
rz(1.788511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2526523) q[2];
sx q[2];
rz(-2.6344968) q[2];
sx q[2];
rz(-2.0895152) q[2];
rz(-2.5893411) q[3];
sx q[3];
rz(-1.4058607) q[3];
sx q[3];
rz(-1.9210057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1110765) q[0];
sx q[0];
rz(-1.80856) q[0];
sx q[0];
rz(-3.1040177) q[0];
rz(-0.75282085) q[1];
sx q[1];
rz(-1.5767187) q[1];
sx q[1];
rz(3.0573696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5606623) q[0];
sx q[0];
rz(-2.5469874) q[0];
sx q[0];
rz(2.0243852) q[0];
rz(-pi) q[1];
rz(1.2758314) q[2];
sx q[2];
rz(-1.9199445) q[2];
sx q[2];
rz(2.5087506) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3854691) q[1];
sx q[1];
rz(-1.9944571) q[1];
sx q[1];
rz(0.045457289) q[1];
x q[2];
rz(-1.2920078) q[3];
sx q[3];
rz(-2.9745347) q[3];
sx q[3];
rz(-0.63705772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5273744) q[2];
sx q[2];
rz(-0.97443333) q[2];
sx q[2];
rz(-2.7091889) q[2];
rz(-1.6480986) q[3];
sx q[3];
rz(-1.3727539) q[3];
sx q[3];
rz(0.73895946) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.768141) q[0];
sx q[0];
rz(-1.5274763) q[0];
sx q[0];
rz(2.3798808) q[0];
rz(2.1181882) q[1];
sx q[1];
rz(-0.49097148) q[1];
sx q[1];
rz(0.34058079) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5368333) q[0];
sx q[0];
rz(-0.7927466) q[0];
sx q[0];
rz(-1.446108) q[0];
rz(-pi) q[1];
rz(2.8958374) q[2];
sx q[2];
rz(-0.57502103) q[2];
sx q[2];
rz(1.0005282) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88587609) q[1];
sx q[1];
rz(-1.9105366) q[1];
sx q[1];
rz(-2.8541616) q[1];
rz(-pi) q[2];
rz(-0.75041308) q[3];
sx q[3];
rz(-1.6831796) q[3];
sx q[3];
rz(0.11311287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9684101) q[2];
sx q[2];
rz(-2.0081655) q[2];
sx q[2];
rz(2.3859712) q[2];
rz(-0.38718811) q[3];
sx q[3];
rz(-1.3914934) q[3];
sx q[3];
rz(-0.79425991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8887535) q[0];
sx q[0];
rz(-3.0536953) q[0];
sx q[0];
rz(1.9320236) q[0];
rz(-1.6190489) q[1];
sx q[1];
rz(-2.3698898) q[1];
sx q[1];
rz(1.3873772) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4693869) q[0];
sx q[0];
rz(-1.9090561) q[0];
sx q[0];
rz(-2.8126008) q[0];
rz(-pi) q[1];
rz(0.99786873) q[2];
sx q[2];
rz(-1.838316) q[2];
sx q[2];
rz(0.47779412) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5249507) q[1];
sx q[1];
rz(-1.9845647) q[1];
sx q[1];
rz(-2.3828808) q[1];
rz(-2.6764836) q[3];
sx q[3];
rz(-1.2832264) q[3];
sx q[3];
rz(-0.63048922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9428955) q[2];
sx q[2];
rz(-2.8479452) q[2];
sx q[2];
rz(2.7315308) q[2];
rz(0.53747082) q[3];
sx q[3];
rz(-1.042807) q[3];
sx q[3];
rz(-1.0679831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10802565) q[0];
sx q[0];
rz(-1.2355868) q[0];
sx q[0];
rz(-2.1349452) q[0];
rz(1.135723) q[1];
sx q[1];
rz(-2.6508811) q[1];
sx q[1];
rz(2.8699285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.131529) q[0];
sx q[0];
rz(-1.9926785) q[0];
sx q[0];
rz(2.5592941) q[0];
rz(-pi) q[1];
rz(2.2211248) q[2];
sx q[2];
rz(-2.1365385) q[2];
sx q[2];
rz(-0.56220245) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2073839) q[1];
sx q[1];
rz(-2.1308239) q[1];
sx q[1];
rz(0.70785849) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74389462) q[3];
sx q[3];
rz(-1.7366323) q[3];
sx q[3];
rz(1.3506571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2037619) q[2];
sx q[2];
rz(-2.2088642) q[2];
sx q[2];
rz(1.1057828) q[2];
rz(-1.3984937) q[3];
sx q[3];
rz(-0.98413697) q[3];
sx q[3];
rz(2.2891146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49649134) q[0];
sx q[0];
rz(-2.190525) q[0];
sx q[0];
rz(2.7880461) q[0];
rz(-1.722466) q[1];
sx q[1];
rz(-1.6058233) q[1];
sx q[1];
rz(3.0790192) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8427575) q[0];
sx q[0];
rz(-2.3681545) q[0];
sx q[0];
rz(-2.406975) q[0];
rz(-pi) q[1];
rz(-2.381787) q[2];
sx q[2];
rz(-1.9927597) q[2];
sx q[2];
rz(-2.9874731) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.500675) q[1];
sx q[1];
rz(-1.0998678) q[1];
sx q[1];
rz(-1.2241609) q[1];
rz(-pi) q[2];
rz(-1.4776286) q[3];
sx q[3];
rz(-0.9103295) q[3];
sx q[3];
rz(0.70445337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8211557) q[2];
sx q[2];
rz(-1.374036) q[2];
sx q[2];
rz(-0.42018166) q[2];
rz(2.7968416) q[3];
sx q[3];
rz(-2.3638066) q[3];
sx q[3];
rz(-0.93728089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47278136) q[0];
sx q[0];
rz(-0.16952276) q[0];
sx q[0];
rz(1.8602759) q[0];
rz(1.5538813) q[1];
sx q[1];
rz(-2.3322767) q[1];
sx q[1];
rz(0.70714998) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67884655) q[0];
sx q[0];
rz(-2.2011559) q[0];
sx q[0];
rz(-0.014883777) q[0];
x q[1];
rz(-2.699171) q[2];
sx q[2];
rz(-2.1083064) q[2];
sx q[2];
rz(0.20172449) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7745819) q[1];
sx q[1];
rz(-0.38644192) q[1];
sx q[1];
rz(2.5143753) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62561927) q[3];
sx q[3];
rz(-1.0418001) q[3];
sx q[3];
rz(2.4629018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3215434) q[2];
sx q[2];
rz(-2.0115435) q[2];
sx q[2];
rz(0.496544) q[2];
rz(-1.2921565) q[3];
sx q[3];
rz(-0.29785952) q[3];
sx q[3];
rz(-2.7636102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15688607) q[0];
sx q[0];
rz(-2.8732185) q[0];
sx q[0];
rz(2.1068841) q[0];
rz(-2.5490419) q[1];
sx q[1];
rz(-2.4609346) q[1];
sx q[1];
rz(1.5998283) q[1];
rz(3.101117) q[2];
sx q[2];
rz(-1.3389836) q[2];
sx q[2];
rz(-0.51868696) q[2];
rz(-2.5377688) q[3];
sx q[3];
rz(-1.7573988) q[3];
sx q[3];
rz(0.16755541) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
