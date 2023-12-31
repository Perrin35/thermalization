OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5470619) q[0];
sx q[0];
rz(4.2630258) q[0];
sx q[0];
rz(6.1017258) q[0];
rz(-1.0815066) q[1];
sx q[1];
rz(-2.4681611) q[1];
sx q[1];
rz(-2.0884617) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5232789) q[0];
sx q[0];
rz(-1.7821454) q[0];
sx q[0];
rz(2.0413415) q[0];
rz(-pi) q[1];
rz(-2.7396169) q[2];
sx q[2];
rz(-0.66265124) q[2];
sx q[2];
rz(-1.0647917) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3095113) q[1];
sx q[1];
rz(-2.32825) q[1];
sx q[1];
rz(-1.8413715) q[1];
rz(-pi) q[2];
rz(-1.231108) q[3];
sx q[3];
rz(-1.7237444) q[3];
sx q[3];
rz(1.4121526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16333214) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(1.7738316) q[2];
rz(1.0129499) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.098009) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(-0.59511551) q[0];
rz(-1.0455421) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(-1.5140623) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13578116) q[0];
sx q[0];
rz(-2.0683056) q[0];
sx q[0];
rz(2.6742629) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1134084) q[2];
sx q[2];
rz(-2.3079254) q[2];
sx q[2];
rz(2.84323) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6914312) q[1];
sx q[1];
rz(-1.3988004) q[1];
sx q[1];
rz(1.7683692) q[1];
x q[2];
rz(-1.1176923) q[3];
sx q[3];
rz(-3.0278904) q[3];
sx q[3];
rz(2.5642455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4864768) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(-0.79616037) q[2];
rz(-2.1697309) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3765091) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(-3.0774975) q[0];
rz(2.8308716) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(-1.6832738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70143914) q[0];
sx q[0];
rz(-0.15471409) q[0];
sx q[0];
rz(-1.8811149) q[0];
x q[1];
rz(2.4086191) q[2];
sx q[2];
rz(-0.79967116) q[2];
sx q[2];
rz(0.016499585) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7160733) q[1];
sx q[1];
rz(-0.41285535) q[1];
sx q[1];
rz(3.0070261) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5536669) q[3];
sx q[3];
rz(-1.8855842) q[3];
sx q[3];
rz(-1.1312315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1594499) q[2];
sx q[2];
rz(-2.2950256) q[2];
sx q[2];
rz(0.67908755) q[2];
rz(-1.3726161) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452334) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(-1.1244208) q[0];
rz(-1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(1.6569998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9254018) q[0];
sx q[0];
rz(-1.5201709) q[0];
sx q[0];
rz(3.1111858) q[0];
x q[1];
rz(-2.9494483) q[2];
sx q[2];
rz(-2.0803183) q[2];
sx q[2];
rz(-1.1332878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.209219) q[1];
sx q[1];
rz(-2.0062431) q[1];
sx q[1];
rz(-0.93851628) q[1];
rz(-0.72317601) q[3];
sx q[3];
rz(-1.3190862) q[3];
sx q[3];
rz(0.46079208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0174039) q[2];
sx q[2];
rz(-1.2839395) q[2];
sx q[2];
rz(1.0162639) q[2];
rz(-1.7381564) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50773412) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(-2.741709) q[0];
rz(-1.1625066) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(-0.17366017) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0608873) q[0];
sx q[0];
rz(-1.6433435) q[0];
sx q[0];
rz(-1.5990431) q[0];
x q[1];
rz(2.747614) q[2];
sx q[2];
rz(-1.3756961) q[2];
sx q[2];
rz(0.56001284) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.083399051) q[1];
sx q[1];
rz(-3.0325583) q[1];
sx q[1];
rz(-1.3074247) q[1];
rz(1.3624304) q[3];
sx q[3];
rz(-1.2456129) q[3];
sx q[3];
rz(-1.9735826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48012039) q[2];
sx q[2];
rz(-1.1602594) q[2];
sx q[2];
rz(-2.373467) q[2];
rz(-0.85401946) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(2.1693726) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0356692) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(-1.6712028) q[0];
rz(0.51180965) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(-1.9981729) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2238732) q[0];
sx q[0];
rz(-1.6544764) q[0];
sx q[0];
rz(-2.9907945) q[0];
x q[1];
rz(-2.5383699) q[2];
sx q[2];
rz(-1.8458741) q[2];
sx q[2];
rz(-0.075866931) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.032635078) q[1];
sx q[1];
rz(-2.7320478) q[1];
sx q[1];
rz(-2.0202548) q[1];
x q[2];
rz(-0.30541909) q[3];
sx q[3];
rz(-2.0913887) q[3];
sx q[3];
rz(0.49314317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64289552) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(-1.5303622) q[2];
rz(-1.6879843) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(-0.24916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11944184) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(1.4402333) q[0];
rz(-0.72921905) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(-1.1332606) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2204809) q[0];
sx q[0];
rz(-2.8163914) q[0];
sx q[0];
rz(0.13334206) q[0];
x q[1];
rz(2.4293443) q[2];
sx q[2];
rz(-1.5603666) q[2];
sx q[2];
rz(-1.041009) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3164191) q[1];
sx q[1];
rz(-2.2146041) q[1];
sx q[1];
rz(-2.7979147) q[1];
rz(0.26182884) q[3];
sx q[3];
rz(-1.7228408) q[3];
sx q[3];
rz(0.1482299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40522727) q[2];
sx q[2];
rz(-1.1232802) q[2];
sx q[2];
rz(2.6531632) q[2];
rz(1.8296261) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(0.64129889) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064780386) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(0.07117614) q[0];
rz(-0.03216234) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(1.9326928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29584822) q[0];
sx q[0];
rz(-2.309572) q[0];
sx q[0];
rz(-2.8315298) q[0];
rz(-pi) q[1];
rz(-1.8078631) q[2];
sx q[2];
rz(-1.7500688) q[2];
sx q[2];
rz(2.8693503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4603685) q[1];
sx q[1];
rz(-1.4404693) q[1];
sx q[1];
rz(3.1254915) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6505359) q[3];
sx q[3];
rz(-0.68248442) q[3];
sx q[3];
rz(-1.3852937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.20748392) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(2.2682244) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(-0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1677925) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(2.8299676) q[0];
rz(0.82178003) q[1];
sx q[1];
rz(-2.5506134) q[1];
sx q[1];
rz(1.5100381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8845997) q[0];
sx q[0];
rz(-1.25367) q[0];
sx q[0];
rz(1.3604128) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2276149) q[2];
sx q[2];
rz(-1.7556264) q[2];
sx q[2];
rz(0.96175805) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.65731243) q[1];
sx q[1];
rz(-1.7727976) q[1];
sx q[1];
rz(2.5217418) q[1];
rz(-pi) q[2];
rz(-0.23267965) q[3];
sx q[3];
rz(-1.7494697) q[3];
sx q[3];
rz(2.132638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(2.3256425) q[2];
rz(0.50968918) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(-1.9036487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508535) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(2.9113286) q[0];
rz(-0.62581217) q[1];
sx q[1];
rz(-0.9451378) q[1];
sx q[1];
rz(2.4831916) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9068245) q[0];
sx q[0];
rz(-2.0567729) q[0];
sx q[0];
rz(0.66396873) q[0];
rz(-pi) q[1];
rz(-2.4670829) q[2];
sx q[2];
rz(-1.837338) q[2];
sx q[2];
rz(0.36905497) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2681594) q[1];
sx q[1];
rz(-2.6323316) q[1];
sx q[1];
rz(1.2251653) q[1];
rz(0.55024054) q[3];
sx q[3];
rz(-2.5329258) q[3];
sx q[3];
rz(-2.254571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6752424) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(-0.33774439) q[2];
rz(-2.1181469) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(1.9256928) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52453775) q[0];
sx q[0];
rz(-2.239997) q[0];
sx q[0];
rz(-0.064185113) q[0];
rz(-1.2383923) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
rz(1.8160309) q[2];
sx q[2];
rz(-2.500848) q[2];
sx q[2];
rz(1.4304888) q[2];
rz(-1.5068549) q[3];
sx q[3];
rz(-2.1422374) q[3];
sx q[3];
rz(-0.16850785) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
