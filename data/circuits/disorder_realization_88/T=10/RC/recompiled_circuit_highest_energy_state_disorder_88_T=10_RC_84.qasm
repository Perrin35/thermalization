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
rz(2.9945381) q[0];
sx q[0];
rz(-1.7400063) q[0];
sx q[0];
rz(2.2214878) q[0];
rz(-0.080634557) q[1];
sx q[1];
rz(3.7205003) q[1];
sx q[1];
rz(10.401934) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1950705) q[0];
sx q[0];
rz(-1.7818639) q[0];
sx q[0];
rz(0.37977438) q[0];
x q[1];
rz(0.10213587) q[2];
sx q[2];
rz(-0.13291026) q[2];
sx q[2];
rz(-1.6249715) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5861533) q[1];
sx q[1];
rz(-1.7057014) q[1];
sx q[1];
rz(0.6026661) q[1];
rz(-pi) q[2];
rz(2.4529934) q[3];
sx q[3];
rz(-1.3312201) q[3];
sx q[3];
rz(-1.857615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5795634) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(-0.62082949) q[2];
rz(-0.51554716) q[3];
sx q[3];
rz(-1.7190869) q[3];
sx q[3];
rz(2.2990885) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2466549) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(2.5625693) q[0];
rz(0.82194263) q[1];
sx q[1];
rz(-1.5162946) q[1];
sx q[1];
rz(2.6878405) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5229086) q[0];
sx q[0];
rz(-2.4755619) q[0];
sx q[0];
rz(2.7722107) q[0];
rz(2.6751509) q[2];
sx q[2];
rz(-1.6714408) q[2];
sx q[2];
rz(-0.60324861) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.0099185506) q[1];
sx q[1];
rz(-1.6955396) q[1];
sx q[1];
rz(2.591955) q[1];
rz(-pi) q[2];
rz(1.0250799) q[3];
sx q[3];
rz(-1.1835872) q[3];
sx q[3];
rz(-0.51147616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.73432505) q[2];
sx q[2];
rz(-0.13886034) q[2];
sx q[2];
rz(-0.56488758) q[2];
rz(-0.35483739) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(-1.423665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9722209) q[0];
sx q[0];
rz(-0.2114978) q[0];
sx q[0];
rz(2.5900904) q[0];
rz(-2.7768199) q[1];
sx q[1];
rz(-0.33939895) q[1];
sx q[1];
rz(-0.50484467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0009815) q[0];
sx q[0];
rz(-1.4112844) q[0];
sx q[0];
rz(-0.41054748) q[0];
x q[1];
rz(2.9302858) q[2];
sx q[2];
rz(-0.42102764) q[2];
sx q[2];
rz(0.94079921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80124679) q[1];
sx q[1];
rz(-2.1643504) q[1];
sx q[1];
rz(2.1562063) q[1];
rz(-pi) q[2];
rz(-2.6608493) q[3];
sx q[3];
rz(-1.2337483) q[3];
sx q[3];
rz(-2.2729006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.43368936) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(-3.0943387) q[2];
rz(-2.3047678) q[3];
sx q[3];
rz(-2.4182726) q[3];
sx q[3];
rz(-1.0251454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4243917) q[0];
sx q[0];
rz(-1.0294788) q[0];
sx q[0];
rz(-1.3264054) q[0];
rz(-1.6560417) q[1];
sx q[1];
rz(-2.0573719) q[1];
sx q[1];
rz(0.63492376) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.062092) q[0];
sx q[0];
rz(-1.5720815) q[0];
sx q[0];
rz(1.5888402) q[0];
rz(-pi) q[1];
rz(2.1702386) q[2];
sx q[2];
rz(-0.58341372) q[2];
sx q[2];
rz(-2.3177878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.288347) q[1];
sx q[1];
rz(-1.4985871) q[1];
sx q[1];
rz(1.9068524) q[1];
rz(2.1011971) q[3];
sx q[3];
rz(-2.1302855) q[3];
sx q[3];
rz(-0.44236576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2403468) q[2];
sx q[2];
rz(-2.2525658) q[2];
sx q[2];
rz(-1.7690313) q[2];
rz(-1.9339804) q[3];
sx q[3];
rz(-0.6438846) q[3];
sx q[3];
rz(-0.27455583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32886252) q[0];
sx q[0];
rz(-0.48506081) q[0];
sx q[0];
rz(3.0364756) q[0];
rz(0.33991995) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(-2.0786659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859131) q[0];
sx q[0];
rz(-2.4300008) q[0];
sx q[0];
rz(3.1200073) q[0];
rz(-pi) q[1];
x q[1];
rz(2.379871) q[2];
sx q[2];
rz(-1.3889379) q[2];
sx q[2];
rz(-1.3052502) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8430994) q[1];
sx q[1];
rz(-1.5122422) q[1];
sx q[1];
rz(1.5565518) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5751198) q[3];
sx q[3];
rz(-2.6427442) q[3];
sx q[3];
rz(-2.3259142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0792599) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(1.3366535) q[2];
rz(-0.054232728) q[3];
sx q[3];
rz(-1.9510061) q[3];
sx q[3];
rz(0.59468734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17276758) q[0];
sx q[0];
rz(-0.69197881) q[0];
sx q[0];
rz(-1.2139976) q[0];
rz(-1.0460151) q[1];
sx q[1];
rz(-1.0449301) q[1];
sx q[1];
rz(0.85375839) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3873685) q[0];
sx q[0];
rz(-1.7048536) q[0];
sx q[0];
rz(-3.0813974) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5327024) q[2];
sx q[2];
rz(-0.85431803) q[2];
sx q[2];
rz(-2.4520055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56928192) q[1];
sx q[1];
rz(-1.0195059) q[1];
sx q[1];
rz(-0.88030073) q[1];
x q[2];
rz(2.7925909) q[3];
sx q[3];
rz(-1.7141388) q[3];
sx q[3];
rz(2.3239612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13009109) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(0.74756527) q[2];
rz(-2.2085564) q[3];
sx q[3];
rz(-1.3676164) q[3];
sx q[3];
rz(-2.8396377) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0316684) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(-0.64144301) q[0];
rz(2.172442) q[1];
sx q[1];
rz(-2.0717924) q[1];
sx q[1];
rz(-2.0279121) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71834356) q[0];
sx q[0];
rz(-1.166543) q[0];
sx q[0];
rz(0.1057616) q[0];
x q[1];
rz(0.53345726) q[2];
sx q[2];
rz(-1.0048303) q[2];
sx q[2];
rz(-0.12803687) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8889474) q[1];
sx q[1];
rz(-1.749649) q[1];
sx q[1];
rz(-2.0986544) q[1];
x q[2];
rz(-1.6051172) q[3];
sx q[3];
rz(-2.0908818) q[3];
sx q[3];
rz(-1.8994562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5668737) q[2];
sx q[2];
rz(-0.69872624) q[2];
sx q[2];
rz(0.75801545) q[2];
rz(-2.8356683) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(-2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4161943) q[0];
sx q[0];
rz(-0.60878009) q[0];
sx q[0];
rz(2.388227) q[0];
rz(-0.92195177) q[1];
sx q[1];
rz(-1.4267068) q[1];
sx q[1];
rz(-3.0070378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3938569) q[0];
sx q[0];
rz(-1.6574291) q[0];
sx q[0];
rz(1.2133693) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33498165) q[2];
sx q[2];
rz(-2.4726395) q[2];
sx q[2];
rz(1.9425336) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2462897) q[1];
sx q[1];
rz(-0.091769204) q[1];
sx q[1];
rz(-0.20905881) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0287168) q[3];
sx q[3];
rz(-0.94967128) q[3];
sx q[3];
rz(0.21645138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77384633) q[2];
sx q[2];
rz(-1.6951122) q[2];
sx q[2];
rz(-1.194225) q[2];
rz(1.9959244) q[3];
sx q[3];
rz(-2.001389) q[3];
sx q[3];
rz(-2.0755419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1016178) q[0];
sx q[0];
rz(-0.45926738) q[0];
sx q[0];
rz(-0.38247821) q[0];
rz(-2.3756012) q[1];
sx q[1];
rz(-1.6440369) q[1];
sx q[1];
rz(1.3409748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4744953) q[0];
sx q[0];
rz(-1.4112368) q[0];
sx q[0];
rz(1.6612396) q[0];
rz(2.3843896) q[2];
sx q[2];
rz(-1.6348038) q[2];
sx q[2];
rz(0.369095) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7786918) q[1];
sx q[1];
rz(-2.3039989) q[1];
sx q[1];
rz(-1.2630839) q[1];
x q[2];
rz(-1.2148592) q[3];
sx q[3];
rz(-1.4934469) q[3];
sx q[3];
rz(-2.429395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0535023) q[2];
sx q[2];
rz(-2.221205) q[2];
sx q[2];
rz(0.11605334) q[2];
rz(-1.8807489) q[3];
sx q[3];
rz(-2.0033658) q[3];
sx q[3];
rz(1.457823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55105227) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(-0.1846479) q[0];
rz(0.91839904) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(1.437423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1383863) q[0];
sx q[0];
rz(-2.0420364) q[0];
sx q[0];
rz(2.9127761) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.060861) q[2];
sx q[2];
rz(-2.417832) q[2];
sx q[2];
rz(-0.91659509) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.559222) q[1];
sx q[1];
rz(-2.5032024) q[1];
sx q[1];
rz(1.2450144) q[1];
x q[2];
rz(2.4653696) q[3];
sx q[3];
rz(-1.8340602) q[3];
sx q[3];
rz(-2.3674008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47716466) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(-2.5028382) q[2];
rz(-2.3264558) q[3];
sx q[3];
rz(-0.4709979) q[3];
sx q[3];
rz(0.42069978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7246134) q[0];
sx q[0];
rz(-1.8920349) q[0];
sx q[0];
rz(-2.98988) q[0];
rz(0.78085113) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(-1.8506321) q[2];
sx q[2];
rz(-1.4536413) q[2];
sx q[2];
rz(-2.8870256) q[2];
rz(-1.5315957) q[3];
sx q[3];
rz(-1.3269674) q[3];
sx q[3];
rz(1.5404601) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
