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
rz(0.18145951) q[0];
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
rz(-1.6183137) q[0];
sx q[0];
rz(-1.7821454) q[0];
sx q[0];
rz(1.1002512) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2744781) q[2];
sx q[2];
rz(-0.9689435) q[2];
sx q[2];
rz(1.5822496) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83208132) q[1];
sx q[1];
rz(-2.32825) q[1];
sx q[1];
rz(-1.3002212) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.231108) q[3];
sx q[3];
rz(-1.7237444) q[3];
sx q[3];
rz(-1.7294401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9782605) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(1.367761) q[2];
rz(1.0129499) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(-3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.098009) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(2.5464771) q[0];
rz(2.0960506) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(-1.5140623) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1986952) q[0];
sx q[0];
rz(-1.1637582) q[0];
sx q[0];
rz(-2.1172949) q[0];
rz(-pi) q[1];
rz(-1.5397649) q[2];
sx q[2];
rz(-2.4040262) q[2];
sx q[2];
rz(2.8013128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6914312) q[1];
sx q[1];
rz(-1.7427923) q[1];
sx q[1];
rz(1.7683692) q[1];
rz(-pi) q[2];
rz(1.1176923) q[3];
sx q[3];
rz(-3.0278904) q[3];
sx q[3];
rz(-2.5642455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4864768) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(-0.79616037) q[2];
rz(2.1697309) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(0.17175737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765091) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(-3.0774975) q[0];
rz(-0.31072101) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(-1.6832738) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9653939) q[0];
sx q[0];
rz(-1.5237234) q[0];
sx q[0];
rz(-1.7182299) q[0];
rz(-pi) q[1];
rz(-0.65285039) q[2];
sx q[2];
rz(-2.0712426) q[2];
sx q[2];
rz(0.99393883) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4255193) q[1];
sx q[1];
rz(-2.7287373) q[1];
sx q[1];
rz(0.13456657) q[1];
rz(-pi) q[2];
x q[2];
rz(0.052555368) q[3];
sx q[3];
rz(-2.8263546) q[3];
sx q[3];
rz(1.0759575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1594499) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(1.3726161) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(-2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2963592) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(-2.0171719) q[0];
rz(1.2202948) q[1];
sx q[1];
rz(-1.9492457) q[1];
sx q[1];
rz(-1.4845928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9254018) q[0];
sx q[0];
rz(-1.6214217) q[0];
sx q[0];
rz(3.1111858) q[0];
rz(-pi) q[1];
x q[1];
rz(1.053327) q[2];
sx q[2];
rz(-1.4033068) q[2];
sx q[2];
rz(-0.53211624) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.209219) q[1];
sx q[1];
rz(-2.0062431) q[1];
sx q[1];
rz(0.93851628) q[1];
rz(0.37064151) q[3];
sx q[3];
rz(-0.7581884) q[3];
sx q[3];
rz(-1.7565808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1241887) q[2];
sx q[2];
rz(-1.2839395) q[2];
sx q[2];
rz(1.0162639) q[2];
rz(1.7381564) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(1.0884292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50773412) q[0];
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
rz(1.0608873) q[0];
sx q[0];
rz(-1.6433435) q[0];
sx q[0];
rz(-1.5425496) q[0];
rz(2.6661751) q[2];
sx q[2];
rz(-2.7042275) q[2];
sx q[2];
rz(1.4471444) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.793321) q[1];
sx q[1];
rz(-1.6760567) q[1];
sx q[1];
rz(-3.1131016) q[1];
rz(-pi) q[2];
rz(1.3624304) q[3];
sx q[3];
rz(-1.8959798) q[3];
sx q[3];
rz(-1.1680101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6614723) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(-0.76812569) q[2];
rz(-2.2875732) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(2.1693726) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1059234) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(-1.6712028) q[0];
rz(2.629783) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(1.9981729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91771942) q[0];
sx q[0];
rz(-1.6544764) q[0];
sx q[0];
rz(-0.15079817) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6799455) q[2];
sx q[2];
rz(-0.6558154) q[2];
sx q[2];
rz(1.1193502) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6900942) q[1];
sx q[1];
rz(-1.203981) q[1];
sx q[1];
rz(0.18641285) q[1];
rz(-2.8361736) q[3];
sx q[3];
rz(-1.050204) q[3];
sx q[3];
rz(-2.6484495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4986971) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(1.5303622) q[2];
rz(1.4536084) q[3];
sx q[3];
rz(-1.0691103) q[3];
sx q[3];
rz(0.24916515) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11944184) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(1.7013593) q[0];
rz(2.4123736) q[1];
sx q[1];
rz(-1.9897285) q[1];
sx q[1];
rz(-2.008332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92111174) q[0];
sx q[0];
rz(-0.32520121) q[0];
sx q[0];
rz(3.0082506) q[0];
x q[1];
rz(1.557017) q[2];
sx q[2];
rz(-0.85859495) q[2];
sx q[2];
rz(2.6028002) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1843695) q[1];
sx q[1];
rz(-1.2979227) q[1];
sx q[1];
rz(2.2437614) q[1];
rz(-pi) q[2];
rz(-0.53448581) q[3];
sx q[3];
rz(-2.839698) q[3];
sx q[3];
rz(-1.2045977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7363654) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(0.48842946) q[2];
rz(1.3119665) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(-0.64129889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064780386) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(3.0704165) q[0];
rz(3.1094303) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(1.9326928) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9933388) q[0];
sx q[0];
rz(-2.3518666) q[0];
sx q[0];
rz(1.8940311) q[0];
x q[1];
rz(-0.9135984) q[2];
sx q[2];
rz(-0.29619869) q[2];
sx q[2];
rz(-1.2072472) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4603685) q[1];
sx q[1];
rz(-1.7011233) q[1];
sx q[1];
rz(-0.016101109) q[1];
rz(-pi) q[2];
x q[2];
rz(1.204793) q[3];
sx q[3];
rz(-2.1605948) q[3];
sx q[3];
rz(2.3596711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20748392) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(-1.570638) q[2];
rz(2.2682244) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(-0.35203716) q[3];
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
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1677925) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(-0.31162509) q[0];
rz(0.82178003) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(-1.5100381) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8845997) q[0];
sx q[0];
rz(-1.25367) q[0];
sx q[0];
rz(1.7811799) q[0];
rz(-pi) q[1];
rz(-2.9456003) q[2];
sx q[2];
rz(-1.2336944) q[2];
sx q[2];
rz(-0.54346426) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4842802) q[1];
sx q[1];
rz(-1.7727976) q[1];
sx q[1];
rz(-0.61985086) q[1];
x q[2];
rz(0.66442499) q[3];
sx q[3];
rz(-2.8492152) q[3];
sx q[3];
rz(1.9362621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(2.3256425) q[2];
rz(2.6319035) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(-1.9036487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8907392) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(0.2302641) q[0];
rz(2.5157805) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(-2.4831916) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9068245) q[0];
sx q[0];
rz(-1.0848197) q[0];
sx q[0];
rz(-0.66396873) q[0];
x q[1];
rz(-1.9071104) q[2];
sx q[2];
rz(-0.92421495) q[2];
sx q[2];
rz(0.99415776) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8734332) q[1];
sx q[1];
rz(-2.6323316) q[1];
sx q[1];
rz(1.9164273) q[1];
rz(-1.2213311) q[3];
sx q[3];
rz(-1.06171) q[3];
sx q[3];
rz(-1.6125319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6752424) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(-2.8038483) q[2];
rz(1.0234458) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(-1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6170549) q[0];
sx q[0];
rz(-2.239997) q[0];
sx q[0];
rz(-0.064185113) q[0];
rz(1.2383923) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(-2.1970489) q[2];
sx q[2];
rz(-1.4251475) q[2];
sx q[2];
rz(0.057694358) q[2];
rz(0.099048793) q[3];
sx q[3];
rz(-0.5746114) q[3];
sx q[3];
rz(-0.050669908) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
