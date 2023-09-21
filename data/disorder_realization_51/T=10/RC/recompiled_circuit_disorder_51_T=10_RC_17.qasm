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
rz(2.060086) q[1];
sx q[1];
rz(-0.67343155) q[1];
sx q[1];
rz(-1.0531309) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9877729) q[0];
sx q[0];
rz(-2.0300555) q[0];
sx q[0];
rz(0.23621969) q[0];
rz(-pi) q[1];
rz(-2.7396169) q[2];
sx q[2];
rz(-2.4789414) q[2];
sx q[2];
rz(1.0647917) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83208132) q[1];
sx q[1];
rz(-2.32825) q[1];
sx q[1];
rz(1.3002212) q[1];
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
rz(-1.1085199) q[2];
sx q[2];
rz(-1.7738316) q[2];
rz(-2.1286428) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(3.0701385) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.098009) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(-2.5464771) q[0];
rz(1.0455421) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(1.5140623) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1986952) q[0];
sx q[0];
rz(-1.9778344) q[0];
sx q[0];
rz(-1.0242978) q[0];
rz(-pi) q[1];
rz(-3.1134084) q[2];
sx q[2];
rz(-0.83366725) q[2];
sx q[2];
rz(2.84323) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4501614) q[1];
sx q[1];
rz(-1.7427923) q[1];
sx q[1];
rz(-1.7683692) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1176923) q[3];
sx q[3];
rz(-0.11370224) q[3];
sx q[3];
rz(-2.5642455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4864768) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(-0.79616037) q[2];
rz(-0.97186175) q[3];
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
rz(-2.7650836) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(-3.0774975) q[0];
rz(-2.8308716) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(-1.6832738) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70143914) q[0];
sx q[0];
rz(-2.9868786) q[0];
sx q[0];
rz(-1.2604777) q[0];
x q[1];
rz(2.4086191) q[2];
sx q[2];
rz(-2.3419215) q[2];
sx q[2];
rz(-0.016499585) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.872936) q[1];
sx q[1];
rz(-1.6246512) q[1];
sx q[1];
rz(-0.40952803) q[1];
rz(-pi) q[2];
rz(1.5536669) q[3];
sx q[3];
rz(-1.2560085) q[3];
sx q[3];
rz(-1.1312315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.98214275) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(1.7689765) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452334) q[0];
sx q[0];
rz(-2.5722752) q[0];
sx q[0];
rz(2.0171719) q[0];
rz(1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(-1.6569998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2161908) q[0];
sx q[0];
rz(-1.5201709) q[0];
sx q[0];
rz(0.030406818) q[0];
x q[1];
rz(-0.19214432) q[2];
sx q[2];
rz(-2.0803183) q[2];
sx q[2];
rz(-2.0083049) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.93814072) q[1];
sx q[1];
rz(-2.1362937) q[1];
sx q[1];
rz(-0.52312619) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2403537) q[3];
sx q[3];
rz(-2.2664865) q[3];
sx q[3];
rz(-0.89360039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1241887) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(2.1253288) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50773412) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(0.39988363) q[0];
rz(-1.9790861) q[1];
sx q[1];
rz(-1.3299273) q[1];
sx q[1];
rz(2.9679325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6337316) q[0];
sx q[0];
rz(-1.5989688) q[0];
sx q[0];
rz(0.072576056) q[0];
rz(-pi) q[1];
rz(2.747614) q[2];
sx q[2];
rz(-1.7658965) q[2];
sx q[2];
rz(-0.56001284) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.793321) q[1];
sx q[1];
rz(-1.4655359) q[1];
sx q[1];
rz(0.028491032) q[1];
x q[2];
rz(1.7791622) q[3];
sx q[3];
rz(-1.2456129) q[3];
sx q[3];
rz(1.9735826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6614723) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(0.76812569) q[2];
rz(-2.2875732) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(-2.1693726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0356692) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(-1.4703898) q[0];
rz(0.51180965) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(-1.1434198) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66577673) q[0];
sx q[0];
rz(-1.4205298) q[0];
sx q[0];
rz(1.6554324) q[0];
rz(-pi) q[1];
rz(-2.6799455) q[2];
sx q[2];
rz(-0.6558154) q[2];
sx q[2];
rz(-1.1193502) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45149849) q[1];
sx q[1];
rz(-1.9376117) q[1];
sx q[1];
rz(2.9551798) q[1];
rz(2.1120758) q[3];
sx q[3];
rz(-1.8347077) q[3];
sx q[3];
rz(-1.2332066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64289552) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(-1.6112304) q[2];
rz(1.4536084) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(2.8924275) q[3];
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
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6183375) q[0];
sx q[0];
rz(-1.528307) q[0];
sx q[0];
rz(-0.32251127) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.015958162) q[2];
sx q[2];
rz(-2.4292813) q[2];
sx q[2];
rz(-0.51770891) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1843695) q[1];
sx q[1];
rz(-1.8436699) q[1];
sx q[1];
rz(2.2437614) q[1];
rz(-pi) q[2];
rz(-1.4134737) q[3];
sx q[3];
rz(-1.3120578) q[3];
sx q[3];
rz(1.7595921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7363654) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(-0.48842946) q[2];
rz(1.8296261) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(-2.5002938) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064780386) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(-0.07117614) q[0];
rz(-3.1094303) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(1.9326928) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6541518) q[0];
sx q[0];
rz(-1.7983266) q[0];
sx q[0];
rz(2.3339416) q[0];
x q[1];
rz(2.9572763) q[2];
sx q[2];
rz(-1.3375998) q[2];
sx q[2];
rz(1.2554982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4603685) q[1];
sx q[1];
rz(-1.7011233) q[1];
sx q[1];
rz(3.1254915) q[1];
rz(-0.49105673) q[3];
sx q[3];
rz(-0.68248442) q[3];
sx q[3];
rz(-1.3852937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.20748392) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(1.570638) q[2];
rz(-0.87336826) q[3];
sx q[3];
rz(-1.7375172) q[3];
sx q[3];
rz(0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97380012) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(-2.8299676) q[0];
rz(-2.3198126) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(1.6315546) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25699297) q[0];
sx q[0];
rz(-1.8879226) q[0];
sx q[0];
rz(1.7811799) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9139778) q[2];
sx q[2];
rz(-1.7556264) q[2];
sx q[2];
rz(-2.1798346) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9540625) q[1];
sx q[1];
rz(-0.64779753) q[1];
sx q[1];
rz(2.8026583) q[1];
rz(-pi) q[2];
rz(-2.4771677) q[3];
sx q[3];
rz(-0.29237745) q[3];
sx q[3];
rz(-1.9362621) q[3];
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
rz(-0.81595016) q[2];
rz(0.50968918) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(1.9036487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508535) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(-0.2302641) q[0];
rz(-2.5157805) q[1];
sx q[1];
rz(-0.9451378) q[1];
sx q[1];
rz(-2.4831916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31357665) q[0];
sx q[0];
rz(-0.99452924) q[0];
sx q[0];
rz(-2.1616031) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9071104) q[2];
sx q[2];
rz(-2.2173777) q[2];
sx q[2];
rz(0.99415776) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.876993) q[1];
sx q[1];
rz(-2.0473285) q[1];
sx q[1];
rz(-0.18696733) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55024054) q[3];
sx q[3];
rz(-2.5329258) q[3];
sx q[3];
rz(2.254571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.6752424) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(0.33774439) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.7815536) q[3];
sx q[3];
rz(-1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52453775) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(1.2383923) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(1.8160309) q[2];
sx q[2];
rz(-2.500848) q[2];
sx q[2];
rz(1.4304888) q[2];
rz(0.57237207) q[3];
sx q[3];
rz(-1.5170245) q[3];
sx q[3];
rz(-1.7046884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
