OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8799514) q[0];
sx q[0];
rz(-2.9961442) q[0];
sx q[0];
rz(2.8737336) q[0];
rz(1.574006) q[1];
sx q[1];
rz(-2.9730453) q[1];
sx q[1];
rz(-0.57810098) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12500873) q[0];
sx q[0];
rz(-1.8468231) q[0];
sx q[0];
rz(0.87645032) q[0];
rz(-pi) q[1];
rz(0.16139754) q[2];
sx q[2];
rz(-0.70490743) q[2];
sx q[2];
rz(-2.7053331) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0904304) q[1];
sx q[1];
rz(-1.1641065) q[1];
sx q[1];
rz(2.6654408) q[1];
rz(-pi) q[2];
rz(0.86712305) q[3];
sx q[3];
rz(-1.4438585) q[3];
sx q[3];
rz(-2.1849887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9428923) q[2];
sx q[2];
rz(-2.7304724) q[2];
sx q[2];
rz(-1.4760419) q[2];
rz(2.5623411) q[3];
sx q[3];
rz(-1.9871291) q[3];
sx q[3];
rz(0.66592413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85775527) q[0];
sx q[0];
rz(-2.0960161) q[0];
sx q[0];
rz(-3.0864518) q[0];
rz(-2.227123) q[1];
sx q[1];
rz(-1.6998467) q[1];
sx q[1];
rz(1.9257911) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8912635) q[0];
sx q[0];
rz(-0.065527409) q[0];
sx q[0];
rz(2.854268) q[0];
rz(-pi) q[1];
rz(2.6138805) q[2];
sx q[2];
rz(-2.0047422) q[2];
sx q[2];
rz(-2.3134856) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3275274) q[1];
sx q[1];
rz(-1.1720997) q[1];
sx q[1];
rz(1.5908123) q[1];
rz(-1.7082105) q[3];
sx q[3];
rz(-2.0169037) q[3];
sx q[3];
rz(-0.29665369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9332283) q[2];
sx q[2];
rz(-2.6231982) q[2];
sx q[2];
rz(-0.62057692) q[2];
rz(1.4536475) q[3];
sx q[3];
rz(-1.0317289) q[3];
sx q[3];
rz(-1.0769963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8715912) q[0];
sx q[0];
rz(-2.8102165) q[0];
sx q[0];
rz(1.510386) q[0];
rz(-2.467678) q[1];
sx q[1];
rz(-1.8464512) q[1];
sx q[1];
rz(-1.5162226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5249859) q[0];
sx q[0];
rz(-1.5956889) q[0];
sx q[0];
rz(-2.1014433) q[0];
x q[1];
rz(-2.2095334) q[2];
sx q[2];
rz(-1.8839508) q[2];
sx q[2];
rz(1.0822288) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.060884692) q[1];
sx q[1];
rz(-1.9906033) q[1];
sx q[1];
rz(-1.3221198) q[1];
rz(2.5862891) q[3];
sx q[3];
rz(-1.318299) q[3];
sx q[3];
rz(-0.50336526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7565833) q[2];
sx q[2];
rz(-1.5417121) q[2];
sx q[2];
rz(-0.6655244) q[2];
rz(-0.74337983) q[3];
sx q[3];
rz(-2.0205108) q[3];
sx q[3];
rz(-2.9140748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-1.7706364) q[0];
sx q[0];
rz(-2.0891068) q[0];
sx q[0];
rz(1.3002522) q[0];
rz(3.0216253) q[1];
sx q[1];
rz(-2.0610466) q[1];
sx q[1];
rz(-1.0467451) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20489731) q[0];
sx q[0];
rz(-0.012210695) q[0];
sx q[0];
rz(-1.8050675) q[0];
rz(-pi) q[1];
rz(0.36645269) q[2];
sx q[2];
rz(-0.88694807) q[2];
sx q[2];
rz(2.2506335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1018104) q[1];
sx q[1];
rz(-2.6521195) q[1];
sx q[1];
rz(1.4533978) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98174121) q[3];
sx q[3];
rz(-1.1906173) q[3];
sx q[3];
rz(-1.5721934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7946502) q[2];
sx q[2];
rz(-2.0796937) q[2];
sx q[2];
rz(2.3700355) q[2];
rz(-2.0906585) q[3];
sx q[3];
rz(-2.1080878) q[3];
sx q[3];
rz(-1.9482025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13450384) q[0];
sx q[0];
rz(-1.4480042) q[0];
sx q[0];
rz(-0.80528468) q[0];
rz(1.6979506) q[1];
sx q[1];
rz(-0.81595683) q[1];
sx q[1];
rz(3.1051292) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6556824) q[0];
sx q[0];
rz(-0.56942372) q[0];
sx q[0];
rz(-1.6605366) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1389937) q[2];
sx q[2];
rz(-1.5130965) q[2];
sx q[2];
rz(2.1926243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7005706) q[1];
sx q[1];
rz(-2.6614958) q[1];
sx q[1];
rz(1.455485) q[1];
x q[2];
rz(2.8920435) q[3];
sx q[3];
rz(-0.59246906) q[3];
sx q[3];
rz(-1.0147926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5263851) q[2];
sx q[2];
rz(-2.0168763) q[2];
sx q[2];
rz(-2.5649694) q[2];
rz(-1.5455101) q[3];
sx q[3];
rz(-0.85448623) q[3];
sx q[3];
rz(0.57687783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1103519) q[0];
sx q[0];
rz(-1.4875655) q[0];
sx q[0];
rz(-2.5417969) q[0];
rz(2.339263) q[1];
sx q[1];
rz(-2.0498514) q[1];
sx q[1];
rz(2.4937627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8819274) q[0];
sx q[0];
rz(-0.81810942) q[0];
sx q[0];
rz(2.1955793) q[0];
rz(-pi) q[1];
rz(0.0070096891) q[2];
sx q[2];
rz(-1.3607549) q[2];
sx q[2];
rz(-1.5633068) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.105567) q[1];
sx q[1];
rz(-2.470825) q[1];
sx q[1];
rz(0.48700602) q[1];
x q[2];
rz(2.4100001) q[3];
sx q[3];
rz(-0.84699455) q[3];
sx q[3];
rz(-0.38672894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1442147) q[2];
sx q[2];
rz(-2.7684559) q[2];
sx q[2];
rz(1.0972265) q[2];
rz(-1.0002452) q[3];
sx q[3];
rz(-0.92366832) q[3];
sx q[3];
rz(-1.3057115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10709396) q[0];
sx q[0];
rz(-1.9871563) q[0];
sx q[0];
rz(0.19919285) q[0];
rz(-1.8432603) q[1];
sx q[1];
rz(-0.36111626) q[1];
sx q[1];
rz(-2.4995506) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10810971) q[0];
sx q[0];
rz(-0.75844049) q[0];
sx q[0];
rz(1.2940501) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0429616) q[2];
sx q[2];
rz(-1.2073852) q[2];
sx q[2];
rz(0.052415457) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.047028001) q[1];
sx q[1];
rz(-2.92647) q[1];
sx q[1];
rz(0.38472979) q[1];
x q[2];
rz(0.20241356) q[3];
sx q[3];
rz(-2.0127986) q[3];
sx q[3];
rz(-1.3998264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9127427) q[2];
sx q[2];
rz(-2.0759089) q[2];
sx q[2];
rz(1.9672811) q[2];
rz(-2.2187388) q[3];
sx q[3];
rz(-0.43761161) q[3];
sx q[3];
rz(-0.16930425) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4787503) q[0];
sx q[0];
rz(-1.6599382) q[0];
sx q[0];
rz(-2.1424275) q[0];
rz(2.303458) q[1];
sx q[1];
rz(-0.90845388) q[1];
sx q[1];
rz(1.025544) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1555754) q[0];
sx q[0];
rz(-0.83283573) q[0];
sx q[0];
rz(-2.2283594) q[0];
rz(2.3633358) q[2];
sx q[2];
rz(-0.35936752) q[2];
sx q[2];
rz(0.12675135) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.57069381) q[1];
sx q[1];
rz(-2.0624196) q[1];
sx q[1];
rz(-0.44906434) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5426264) q[3];
sx q[3];
rz(-2.3271797) q[3];
sx q[3];
rz(1.6391476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.020236882) q[2];
sx q[2];
rz(-1.7842224) q[2];
sx q[2];
rz(-2.9244002) q[2];
rz(1.1413261) q[3];
sx q[3];
rz(-0.3749899) q[3];
sx q[3];
rz(0.91514897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7104257) q[0];
sx q[0];
rz(-1.2942261) q[0];
sx q[0];
rz(-3.1173832) q[0];
rz(-1.0912033) q[1];
sx q[1];
rz(-1.1187436) q[1];
sx q[1];
rz(0.81407636) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5494848) q[0];
sx q[0];
rz(-2.4257437) q[0];
sx q[0];
rz(-3.0107193) q[0];
x q[1];
rz(0.8701088) q[2];
sx q[2];
rz(-2.3860156) q[2];
sx q[2];
rz(2.8062964) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7634552) q[1];
sx q[1];
rz(-2.3296851) q[1];
sx q[1];
rz(-0.0017536963) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2132675) q[3];
sx q[3];
rz(-2.8983064) q[3];
sx q[3];
rz(-1.4668087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39390627) q[2];
sx q[2];
rz(-0.86204356) q[2];
sx q[2];
rz(1.6115335) q[2];
rz(-2.090442) q[3];
sx q[3];
rz(-1.3408778) q[3];
sx q[3];
rz(3.0401201) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53973389) q[0];
sx q[0];
rz(-2.1506385) q[0];
sx q[0];
rz(0.47716004) q[0];
rz(-1.5513783) q[1];
sx q[1];
rz(-1.4789707) q[1];
sx q[1];
rz(1.8345376) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73471776) q[0];
sx q[0];
rz(-2.2860323) q[0];
sx q[0];
rz(0.060013219) q[0];
x q[1];
rz(-2.3315694) q[2];
sx q[2];
rz(-2.5794583) q[2];
sx q[2];
rz(0.13371828) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.160678) q[1];
sx q[1];
rz(-1.0458993) q[1];
sx q[1];
rz(0.69634931) q[1];
rz(-pi) q[2];
rz(-2.1784276) q[3];
sx q[3];
rz(-1.9279216) q[3];
sx q[3];
rz(-2.4357093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9755379) q[2];
sx q[2];
rz(-2.1259978) q[2];
sx q[2];
rz(0.75237742) q[2];
rz(3.0289529) q[3];
sx q[3];
rz(-0.17387667) q[3];
sx q[3];
rz(-1.9788474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.024121506) q[0];
sx q[0];
rz(-1.7565256) q[0];
sx q[0];
rz(1.1275445) q[0];
rz(-2.7680001) q[1];
sx q[1];
rz(-1.6289381) q[1];
sx q[1];
rz(-2.1388114) q[1];
rz(-2.3898771) q[2];
sx q[2];
rz(-0.57705078) q[2];
sx q[2];
rz(-2.5129872) q[2];
rz(0.74736377) q[3];
sx q[3];
rz(-1.7186281) q[3];
sx q[3];
rz(-2.4627198) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
