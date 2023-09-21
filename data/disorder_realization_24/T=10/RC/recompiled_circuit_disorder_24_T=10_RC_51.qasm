OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(2.0547325) q[0];
sx q[0];
rz(7.6261043) q[0];
rz(-1.5919332) q[1];
sx q[1];
rz(-3.0631493) q[1];
sx q[1];
rz(0.6426386) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6519878) q[0];
sx q[0];
rz(-1.542924) q[0];
sx q[0];
rz(-2.6023988) q[0];
rz(-pi) q[1];
rz(2.6334555) q[2];
sx q[2];
rz(-1.4027486) q[2];
sx q[2];
rz(2.7009168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.027015162) q[1];
sx q[1];
rz(-1.4511961) q[1];
sx q[1];
rz(-3.0199865) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23407614) q[3];
sx q[3];
rz(-1.5353068) q[3];
sx q[3];
rz(1.2897829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6671483) q[2];
sx q[2];
rz(-2.4953304) q[2];
sx q[2];
rz(-1.8623964) q[2];
rz(2.4228418) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(2.8180502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6779125) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(-2.1610778) q[0];
rz(-0.15788831) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(0.78871361) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4855027) q[0];
sx q[0];
rz(-0.1682818) q[0];
sx q[0];
rz(-0.82505723) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5390225) q[2];
sx q[2];
rz(-1.3147768) q[2];
sx q[2];
rz(-0.94871828) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84856725) q[1];
sx q[1];
rz(-1.8752521) q[1];
sx q[1];
rz(0.57363631) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85487811) q[3];
sx q[3];
rz(-2.2023871) q[3];
sx q[3];
rz(-2.7090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.366189) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(-2.9555087) q[2];
rz(-0.6535334) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(0.11793605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2114975) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(2.511456) q[0];
rz(0.12763003) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(-0.72174597) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16746178) q[0];
sx q[0];
rz(-2.2142017) q[0];
sx q[0];
rz(-1.2533623) q[0];
rz(-1.1100936) q[2];
sx q[2];
rz(-2.2113872) q[2];
sx q[2];
rz(-3.1415423) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0112146) q[1];
sx q[1];
rz(-1.6092669) q[1];
sx q[1];
rz(-2.6973666) q[1];
rz(1.1214764) q[3];
sx q[3];
rz(-2.46393) q[3];
sx q[3];
rz(-0.47798702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7599941) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-0.90908137) q[2];
rz(-2.6233853) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(2.853493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58406126) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(-0.042536143) q[0];
rz(2.361239) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(-2.3775878) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9873136) q[0];
sx q[0];
rz(-2.2494082) q[0];
sx q[0];
rz(-0.016088386) q[0];
rz(-pi) q[1];
rz(2.9742083) q[2];
sx q[2];
rz(-0.65132729) q[2];
sx q[2];
rz(-1.960388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.15061513) q[1];
sx q[1];
rz(-0.58826485) q[1];
sx q[1];
rz(2.2554643) q[1];
rz(-pi) q[2];
rz(-2.571991) q[3];
sx q[3];
rz(-2.3929425) q[3];
sx q[3];
rz(-0.99018712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2970695) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(2.6085473) q[2];
rz(2.879203) q[3];
sx q[3];
rz(-0.636594) q[3];
sx q[3];
rz(-2.6791402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2247291) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(-1.2438783) q[0];
rz(0.22661701) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(2.7382543) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4606844) q[0];
sx q[0];
rz(-1.4049238) q[0];
sx q[0];
rz(2.9537863) q[0];
rz(-pi) q[1];
rz(-1.2071768) q[2];
sx q[2];
rz(-2.223613) q[2];
sx q[2];
rz(1.2448685) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12521872) q[1];
sx q[1];
rz(-2.5902936) q[1];
sx q[1];
rz(0.2972879) q[1];
rz(-2.1020528) q[3];
sx q[3];
rz(-0.96635339) q[3];
sx q[3];
rz(1.8464551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.32101813) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(0.78249758) q[2];
rz(-1.1123505) q[3];
sx q[3];
rz(-0.4370884) q[3];
sx q[3];
rz(-2.1330244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7111506) q[0];
sx q[0];
rz(-2.7395881) q[0];
sx q[0];
rz(-0.45561403) q[0];
rz(0.09952155) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(-0.0064370357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71953668) q[0];
sx q[0];
rz(-0.7932084) q[0];
sx q[0];
rz(2.7842245) q[0];
x q[1];
rz(-0.10613425) q[2];
sx q[2];
rz(-1.1960293) q[2];
sx q[2];
rz(2.8449164) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41758075) q[1];
sx q[1];
rz(-1.9829826) q[1];
sx q[1];
rz(-2.8860303) q[1];
x q[2];
rz(2.3859343) q[3];
sx q[3];
rz(-1.0010127) q[3];
sx q[3];
rz(2.8017686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7727938) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(-0.84645611) q[2];
rz(0.99772292) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(-1.6830106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0434175) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(3.085882) q[0];
rz(-0.78272351) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.9810716) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46471483) q[0];
sx q[0];
rz(-2.5589295) q[0];
sx q[0];
rz(2.23784) q[0];
x q[1];
rz(-0.84476446) q[2];
sx q[2];
rz(-0.65462199) q[2];
sx q[2];
rz(2.1824333) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5934505) q[1];
sx q[1];
rz(-1.1957809) q[1];
sx q[1];
rz(-2.1795991) q[1];
x q[2];
rz(-0.37961827) q[3];
sx q[3];
rz(-1.2109204) q[3];
sx q[3];
rz(0.66756638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.065585) q[2];
sx q[2];
rz(-0.92131725) q[2];
sx q[2];
rz(0.78545061) q[2];
rz(0.75585946) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(-2.7261962) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8287559) q[0];
sx q[0];
rz(-3.0506595) q[0];
sx q[0];
rz(-1.1428517) q[0];
rz(1.8354592) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(-0.41608861) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0932255) q[0];
sx q[0];
rz(-1.0972293) q[0];
sx q[0];
rz(2.918539) q[0];
rz(2.9990254) q[2];
sx q[2];
rz(-1.1066184) q[2];
sx q[2];
rz(-0.67205059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9302952) q[1];
sx q[1];
rz(-2.7080309) q[1];
sx q[1];
rz(-0.3890721) q[1];
x q[2];
rz(-1.7851402) q[3];
sx q[3];
rz(-2.6594901) q[3];
sx q[3];
rz(1.6899504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0351506) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(-2.8184334) q[2];
rz(0.20251003) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37339661) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(-1.9966104) q[0];
rz(1.9454983) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(-0.51913613) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3648758) q[0];
sx q[0];
rz(-1.4380699) q[0];
sx q[0];
rz(1.9507292) q[0];
rz(-pi) q[1];
rz(0.81705117) q[2];
sx q[2];
rz(-2.8231986) q[2];
sx q[2];
rz(1.5621834) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12431006) q[1];
sx q[1];
rz(-0.78622422) q[1];
sx q[1];
rz(3.0148274) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9037644) q[3];
sx q[3];
rz(-1.7416218) q[3];
sx q[3];
rz(-1.886614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3141994) q[2];
sx q[2];
rz(-1.2131571) q[2];
sx q[2];
rz(0.040977565) q[2];
rz(0.86769062) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39559078) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(1.6145153) q[0];
rz(1.7136259) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(-1.4987) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.98793) q[0];
sx q[0];
rz(-1.4953574) q[0];
sx q[0];
rz(3.1213785) q[0];
rz(-pi) q[1];
rz(-1.4156878) q[2];
sx q[2];
rz(-1.6181706) q[2];
sx q[2];
rz(-0.93946379) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0000671) q[1];
sx q[1];
rz(-2.2728517) q[1];
sx q[1];
rz(-0.4483923) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4478217) q[3];
sx q[3];
rz(-1.7719367) q[3];
sx q[3];
rz(0.5826544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3830118) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(0.14653462) q[2];
rz(-2.3274029) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(-0.41671419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9183337) q[0];
sx q[0];
rz(-1.2467361) q[0];
sx q[0];
rz(0.99714739) q[0];
rz(1.8756443) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(-2.0064034) q[2];
sx q[2];
rz(-0.26293593) q[2];
sx q[2];
rz(1.5756366) q[2];
rz(2.4155865) q[3];
sx q[3];
rz(-0.86961679) q[3];
sx q[3];
rz(1.1179954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];