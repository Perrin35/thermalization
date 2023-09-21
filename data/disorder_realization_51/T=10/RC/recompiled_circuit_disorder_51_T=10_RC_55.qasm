OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.59453073) q[0];
sx q[0];
rz(-1.1214331) q[0];
sx q[0];
rz(-2.9601331) q[0];
rz(-1.0815066) q[1];
sx q[1];
rz(-2.4681611) q[1];
sx q[1];
rz(-2.0884617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7979413) q[0];
sx q[0];
rz(-2.6290253) q[0];
sx q[0];
rz(1.1287862) q[0];
x q[1];
rz(1.2744781) q[2];
sx q[2];
rz(-0.9689435) q[2];
sx q[2];
rz(-1.559343) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2145572) q[1];
sx q[1];
rz(-1.3753478) q[1];
sx q[1];
rz(0.77597015) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0041204) q[3];
sx q[3];
rz(-2.7702799) q[3];
sx q[3];
rz(2.5759047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9782605) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(-1.367761) q[2];
rz(1.0129499) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0435836) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(2.5464771) q[0];
rz(-1.0455421) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(1.6275303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9428974) q[0];
sx q[0];
rz(-1.9778344) q[0];
sx q[0];
rz(-1.0242978) q[0];
rz(-pi) q[1];
rz(-2.3081231) q[2];
sx q[2];
rz(-1.5916628) q[2];
sx q[2];
rz(1.8881063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6914312) q[1];
sx q[1];
rz(-1.7427923) q[1];
sx q[1];
rz(-1.7683692) q[1];
rz(-1.1176923) q[3];
sx q[3];
rz(-3.0278904) q[3];
sx q[3];
rz(2.5642455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4864768) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(0.79616037) q[2];
rz(-2.1697309) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(0.17175737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765091) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(0.064095108) q[0];
rz(2.8308716) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(1.4583189) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1761988) q[0];
sx q[0];
rz(-1.5237234) q[0];
sx q[0];
rz(1.4233627) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4086191) q[2];
sx q[2];
rz(-0.79967116) q[2];
sx q[2];
rz(-0.016499585) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8628143) q[1];
sx q[1];
rz(-1.1618975) q[1];
sx q[1];
rz(-1.6294953) q[1];
rz(-pi) q[2];
x q[2];
rz(0.052555368) q[3];
sx q[3];
rz(-2.8263546) q[3];
sx q[3];
rz(-2.0656352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1594499) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(1.3726161) q[3];
sx q[3];
rz(-1.8434098) q[3];
sx q[3];
rz(2.1931271) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452334) q[0];
sx q[0];
rz(-2.5722752) q[0];
sx q[0];
rz(1.1244208) q[0];
rz(-1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(1.6569998) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9254018) q[0];
sx q[0];
rz(-1.6214217) q[0];
sx q[0];
rz(-0.030406818) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9494483) q[2];
sx q[2];
rz(-1.0612744) q[2];
sx q[2];
rz(-1.1332878) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9323737) q[1];
sx q[1];
rz(-2.0062431) q[1];
sx q[1];
rz(2.2030764) q[1];
rz(-pi) q[2];
rz(-0.72317601) q[3];
sx q[3];
rz(-1.3190862) q[3];
sx q[3];
rz(-2.6808006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1241887) q[2];
sx q[2];
rz(-1.2839395) q[2];
sx q[2];
rz(1.0162639) q[2];
rz(-1.7381564) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(-1.0884292) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6338585) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(-2.741709) q[0];
rz(1.9790861) q[1];
sx q[1];
rz(-1.3299273) q[1];
sx q[1];
rz(0.17366017) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50786103) q[0];
sx q[0];
rz(-1.5426239) q[0];
sx q[0];
rz(-3.0690166) q[0];
x q[1];
rz(1.3599672) q[2];
sx q[2];
rz(-1.1846917) q[2];
sx q[2];
rz(-2.2112276) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.083399051) q[1];
sx q[1];
rz(-0.10903437) q[1];
sx q[1];
rz(-1.8341679) q[1];
rz(-pi) q[2];
rz(-0.55032702) q[3];
sx q[3];
rz(-0.384207) q[3];
sx q[3];
rz(1.3889544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6614723) q[2];
sx q[2];
rz(-1.1602594) q[2];
sx q[2];
rz(-2.373467) q[2];
rz(-2.2875732) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(-0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1059234) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(-1.6712028) q[0];
rz(2.629783) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(1.1434198) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2238732) q[0];
sx q[0];
rz(-1.4871162) q[0];
sx q[0];
rz(2.9907945) q[0];
rz(-2.6799455) q[2];
sx q[2];
rz(-0.6558154) q[2];
sx q[2];
rz(-1.1193502) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1868362) q[1];
sx q[1];
rz(-1.7446767) q[1];
sx q[1];
rz(1.1980921) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0295168) q[3];
sx q[3];
rz(-1.8347077) q[3];
sx q[3];
rz(1.2332066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.64289552) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(-1.6112304) q[2];
rz(1.6879843) q[3];
sx q[3];
rz(-1.0691103) q[3];
sx q[3];
rz(-0.24916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0221508) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(-1.4402333) q[0];
rz(0.72921905) q[1];
sx q[1];
rz(-1.9897285) q[1];
sx q[1];
rz(2.008332) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5232552) q[0];
sx q[0];
rz(-1.6132857) q[0];
sx q[0];
rz(2.8190814) q[0];
rz(-pi) q[1];
rz(1.5845756) q[2];
sx q[2];
rz(-0.85859495) q[2];
sx q[2];
rz(-2.6028002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82517351) q[1];
sx q[1];
rz(-2.2146041) q[1];
sx q[1];
rz(2.7979147) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4134737) q[3];
sx q[3];
rz(-1.3120578) q[3];
sx q[3];
rz(-1.7595921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7363654) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(0.48842946) q[2];
rz(-1.8296261) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(0.64129889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0768123) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(-3.0704165) q[0];
rz(-3.1094303) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(-1.2088998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4874408) q[0];
sx q[0];
rz(-1.7983266) q[0];
sx q[0];
rz(0.80765101) q[0];
x q[1];
rz(2.9572763) q[2];
sx q[2];
rz(-1.3375998) q[2];
sx q[2];
rz(1.2554982) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2541131) q[1];
sx q[1];
rz(-1.5867609) q[1];
sx q[1];
rz(1.70114) q[1];
rz(-1.204793) q[3];
sx q[3];
rz(-2.1605948) q[3];
sx q[3];
rz(0.78192151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9341087) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(-2.2682244) q[3];
sx q[3];
rz(-1.7375172) q[3];
sx q[3];
rz(-0.35203716) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97380012) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(2.8299676) q[0];
rz(-2.3198126) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(1.6315546) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85743839) q[0];
sx q[0];
rz(-0.3785924) q[0];
sx q[0];
rz(-0.5666825) q[0];
x q[1];
rz(2.0779583) q[2];
sx q[2];
rz(-0.38804752) q[2];
sx q[2];
rz(-2.0575112) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65731243) q[1];
sx q[1];
rz(-1.3687951) q[1];
sx q[1];
rz(0.61985086) q[1];
rz(-pi) q[2];
rz(0.66442499) q[3];
sx q[3];
rz(-2.8492152) q[3];
sx q[3];
rz(-1.2053306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(2.3256425) q[2];
rz(-0.50968918) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508535) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(-0.62581217) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(0.65840107) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9394768) q[0];
sx q[0];
rz(-0.80047031) q[0];
sx q[0];
rz(2.4329484) q[0];
rz(2.7294331) q[2];
sx q[2];
rz(-2.4240652) q[2];
sx q[2];
rz(-1.621643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8734332) q[1];
sx q[1];
rz(-0.50926103) q[1];
sx q[1];
rz(1.9164273) q[1];
rz(0.53604605) q[3];
sx q[3];
rz(-1.2671766) q[3];
sx q[3];
rz(-2.924078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6752424) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(-0.33774439) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(-1.9256928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6170549) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(1.9032003) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
rz(-1.8160309) q[2];
sx q[2];
rz(-0.64074466) q[2];
sx q[2];
rz(-1.7111039) q[2];
rz(-0.57237207) q[3];
sx q[3];
rz(-1.6245681) q[3];
sx q[3];
rz(1.4369042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];