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
rz(0.55573207) q[0];
sx q[0];
rz(4.421173) q[0];
sx q[0];
rz(9.0968994) q[0];
rz(0.15283395) q[1];
sx q[1];
rz(-0.489355) q[1];
sx q[1];
rz(-2.1305003) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3837102) q[0];
sx q[0];
rz(-0.61311904) q[0];
sx q[0];
rz(-2.8611373) q[0];
x q[1];
rz(-0.4083486) q[2];
sx q[2];
rz(-1.4065925) q[2];
sx q[2];
rz(2.4105031) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3950306) q[1];
sx q[1];
rz(-0.85828188) q[1];
sx q[1];
rz(2.0601574) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1527083) q[3];
sx q[3];
rz(-1.7888513) q[3];
sx q[3];
rz(-0.59948987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8622387) q[2];
sx q[2];
rz(-0.78247672) q[2];
sx q[2];
rz(1.5834825) q[2];
rz(2.8067348) q[3];
sx q[3];
rz(-2.0575276) q[3];
sx q[3];
rz(0.28086942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8849477) q[0];
sx q[0];
rz(-2.5769951) q[0];
sx q[0];
rz(2.8096492) q[0];
rz(0.360082) q[1];
sx q[1];
rz(-1.8372476) q[1];
sx q[1];
rz(-0.28775451) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5482564) q[0];
sx q[0];
rz(-1.7711444) q[0];
sx q[0];
rz(-0.25357539) q[0];
x q[1];
rz(2.4096978) q[2];
sx q[2];
rz(-1.8127155) q[2];
sx q[2];
rz(0.11920028) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1102241) q[1];
sx q[1];
rz(-1.9166295) q[1];
sx q[1];
rz(2.7909173) q[1];
x q[2];
rz(0.23116206) q[3];
sx q[3];
rz(-1.729768) q[3];
sx q[3];
rz(2.8249225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.687279) q[2];
sx q[2];
rz(-1.5328898) q[2];
sx q[2];
rz(-1.3909371) q[2];
rz(-2.2823997) q[3];
sx q[3];
rz(-1.2733368) q[3];
sx q[3];
rz(2.9505762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39621064) q[0];
sx q[0];
rz(-1.6062382) q[0];
sx q[0];
rz(0.23042738) q[0];
rz(-0.51741171) q[1];
sx q[1];
rz(-2.0473862) q[1];
sx q[1];
rz(0.28883019) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5342425) q[0];
sx q[0];
rz(-1.6301883) q[0];
sx q[0];
rz(1.1618105) q[0];
x q[1];
rz(0.46858139) q[2];
sx q[2];
rz(-1.6372674) q[2];
sx q[2];
rz(-2.2867212) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7728665) q[1];
sx q[1];
rz(-1.9802367) q[1];
sx q[1];
rz(-0.057842908) q[1];
x q[2];
rz(2.2346441) q[3];
sx q[3];
rz(-0.28093279) q[3];
sx q[3];
rz(2.8317766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1032224) q[2];
sx q[2];
rz(-0.73870814) q[2];
sx q[2];
rz(3.1008516) q[2];
rz(3.0401958) q[3];
sx q[3];
rz(-1.3359759) q[3];
sx q[3];
rz(-1.066677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3876225) q[0];
sx q[0];
rz(-3.1356223) q[0];
sx q[0];
rz(-0.23400865) q[0];
rz(-0.19800828) q[1];
sx q[1];
rz(-2.104069) q[1];
sx q[1];
rz(2.1580946) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9390181) q[0];
sx q[0];
rz(-1.6133623) q[0];
sx q[0];
rz(2.0545511) q[0];
rz(-0.16093851) q[2];
sx q[2];
rz(-2.6223824) q[2];
sx q[2];
rz(0.96566654) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8810087) q[1];
sx q[1];
rz(-2.1420292) q[1];
sx q[1];
rz(-0.84239475) q[1];
rz(-1.6487021) q[3];
sx q[3];
rz(-0.39038218) q[3];
sx q[3];
rz(-1.4161033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46821758) q[2];
sx q[2];
rz(-0.28926352) q[2];
sx q[2];
rz(-2.3667228) q[2];
rz(1.3261999) q[3];
sx q[3];
rz(-1.9522791) q[3];
sx q[3];
rz(-2.6302122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4020017) q[0];
sx q[0];
rz(-2.2690161) q[0];
sx q[0];
rz(-3.1091029) q[0];
rz(1.2292817) q[1];
sx q[1];
rz(-0.928343) q[1];
sx q[1];
rz(-0.083018735) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56606021) q[0];
sx q[0];
rz(-2.0216746) q[0];
sx q[0];
rz(1.321248) q[0];
rz(-pi) q[1];
rz(-0.26814383) q[2];
sx q[2];
rz(-1.4073155) q[2];
sx q[2];
rz(-3.1250866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8889035) q[1];
sx q[1];
rz(-1.1071883) q[1];
sx q[1];
rz(-0.75485922) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1260499) q[3];
sx q[3];
rz(-2.2076026) q[3];
sx q[3];
rz(-1.6890845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3518389) q[2];
sx q[2];
rz(-2.3408076) q[2];
sx q[2];
rz(-2.9808673) q[2];
rz(-1.1927346) q[3];
sx q[3];
rz(-2.9294117) q[3];
sx q[3];
rz(0.76782697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47650325) q[0];
sx q[0];
rz(-2.022321) q[0];
sx q[0];
rz(-0.29092586) q[0];
rz(1.7581958) q[1];
sx q[1];
rz(-1.29888) q[1];
sx q[1];
rz(0.49984041) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2384278) q[0];
sx q[0];
rz(-0.060783371) q[0];
sx q[0];
rz(-0.77068909) q[0];
rz(-2.1792214) q[2];
sx q[2];
rz(-2.314724) q[2];
sx q[2];
rz(0.27368557) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0306727) q[1];
sx q[1];
rz(-2.3706712) q[1];
sx q[1];
rz(-2.4708807) q[1];
rz(2.6074516) q[3];
sx q[3];
rz(-0.93399601) q[3];
sx q[3];
rz(-2.4780688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2073888) q[2];
sx q[2];
rz(-2.5298205) q[2];
sx q[2];
rz(-2.440051) q[2];
rz(-0.97405854) q[3];
sx q[3];
rz(-2.3249224) q[3];
sx q[3];
rz(-0.90132236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3153673) q[0];
sx q[0];
rz(-2.3235445) q[0];
sx q[0];
rz(-2.4819964) q[0];
rz(1.2391799) q[1];
sx q[1];
rz(-2.0678803) q[1];
sx q[1];
rz(2.0674131) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.488193) q[0];
sx q[0];
rz(-1.5563413) q[0];
sx q[0];
rz(3.1281592) q[0];
x q[1];
rz(2.7137382) q[2];
sx q[2];
rz(-1.9249467) q[2];
sx q[2];
rz(0.71989518) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82141187) q[1];
sx q[1];
rz(-1.3548676) q[1];
sx q[1];
rz(1.126299) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76403615) q[3];
sx q[3];
rz(-1.2437399) q[3];
sx q[3];
rz(-2.5491109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85019511) q[2];
sx q[2];
rz(-2.3829298) q[2];
sx q[2];
rz(2.5568753) q[2];
rz(2.0753453) q[3];
sx q[3];
rz(-1.9138347) q[3];
sx q[3];
rz(-0.40759459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19121118) q[0];
sx q[0];
rz(-1.1703015) q[0];
sx q[0];
rz(2.6608652) q[0];
rz(-1.9302543) q[1];
sx q[1];
rz(-1.6119266) q[1];
sx q[1];
rz(2.5239677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0636198) q[0];
sx q[0];
rz(-1.5750066) q[0];
sx q[0];
rz(1.8959037) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5343393) q[2];
sx q[2];
rz(-1.1285845) q[2];
sx q[2];
rz(-0.84458015) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27965266) q[1];
sx q[1];
rz(-1.8649001) q[1];
sx q[1];
rz(0.22464629) q[1];
x q[2];
rz(1.5354352) q[3];
sx q[3];
rz(-1.3786982) q[3];
sx q[3];
rz(2.0336322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9628613) q[2];
sx q[2];
rz(-1.3946673) q[2];
sx q[2];
rz(2.2507131) q[2];
rz(2.1122872) q[3];
sx q[3];
rz(-1.3903214) q[3];
sx q[3];
rz(-1.5911969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7903098) q[0];
sx q[0];
rz(-0.93515486) q[0];
sx q[0];
rz(-1.9679605) q[0];
rz(-1.952518) q[1];
sx q[1];
rz(-2.3937841) q[1];
sx q[1];
rz(2.660451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8136776) q[0];
sx q[0];
rz(-2.5338533) q[0];
sx q[0];
rz(-0.39791664) q[0];
rz(1.4699773) q[2];
sx q[2];
rz(-2.1283796) q[2];
sx q[2];
rz(-2.4876311) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.18080841) q[1];
sx q[1];
rz(-2.167302) q[1];
sx q[1];
rz(-0.1779314) q[1];
x q[2];
rz(-2.6974996) q[3];
sx q[3];
rz(-1.0837348) q[3];
sx q[3];
rz(2.5626078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.22566191) q[2];
sx q[2];
rz(-1.2292726) q[2];
sx q[2];
rz(-1.4813598) q[2];
rz(0.046317421) q[3];
sx q[3];
rz(-1.9527083) q[3];
sx q[3];
rz(-1.2272629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.34761053) q[0];
sx q[0];
rz(-2.1078258) q[0];
sx q[0];
rz(3.0718497) q[0];
rz(-0.28930411) q[1];
sx q[1];
rz(-1.8569088) q[1];
sx q[1];
rz(-2.0893673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7523868) q[0];
sx q[0];
rz(-1.2121887) q[0];
sx q[0];
rz(1.2239271) q[0];
rz(-1.7830816) q[2];
sx q[2];
rz(-0.67321482) q[2];
sx q[2];
rz(-2.5517983) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6695646) q[1];
sx q[1];
rz(-1.9306364) q[1];
sx q[1];
rz(-1.9739975) q[1];
x q[2];
rz(-1.1484365) q[3];
sx q[3];
rz(-1.1760654) q[3];
sx q[3];
rz(-2.0572061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59794402) q[2];
sx q[2];
rz(-1.1755377) q[2];
sx q[2];
rz(0.0607461) q[2];
rz(-2.8525823) q[3];
sx q[3];
rz(-1.3649536) q[3];
sx q[3];
rz(1.2750767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.48988265) q[0];
sx q[0];
rz(-1.6429506) q[0];
sx q[0];
rz(2.0605675) q[0];
rz(-1.0501077) q[1];
sx q[1];
rz(-2.0790015) q[1];
sx q[1];
rz(1.9008295) q[1];
rz(2.0403258) q[2];
sx q[2];
rz(-1.7055837) q[2];
sx q[2];
rz(-2.3619426) q[2];
rz(-1.1393094) q[3];
sx q[3];
rz(-1.3130975) q[3];
sx q[3];
rz(1.3102875) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
