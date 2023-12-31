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
rz(1.0531309) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9877729) q[0];
sx q[0];
rz(-1.1115371) q[0];
sx q[0];
rz(0.23621969) q[0];
rz(-pi) q[1];
rz(0.40197576) q[2];
sx q[2];
rz(-2.4789414) q[2];
sx q[2];
rz(-2.0768009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44838201) q[1];
sx q[1];
rz(-2.3464077) q[1];
sx q[1];
rz(2.8661212) q[1];
rz(-pi) q[2];
rz(-2.0041204) q[3];
sx q[3];
rz(-0.37131272) q[3];
sx q[3];
rz(-0.565688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.16333214) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(1.7738316) q[2];
rz(-2.1286428) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(-0.071454123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.098009) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(-2.5464771) q[0];
rz(1.0455421) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(1.5140623) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0058115) q[0];
sx q[0];
rz(-1.0732871) q[0];
sx q[0];
rz(0.46732975) q[0];
x q[1];
rz(2.3081231) q[2];
sx q[2];
rz(-1.5499299) q[2];
sx q[2];
rz(-1.2534864) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6914312) q[1];
sx q[1];
rz(-1.3988004) q[1];
sx q[1];
rz(1.7683692) q[1];
rz(-pi) q[2];
rz(2.0239003) q[3];
sx q[3];
rz(-3.0278904) q[3];
sx q[3];
rz(2.5642455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.65511584) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(2.3454323) q[2];
rz(-2.1697309) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(0.17175737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-0.3765091) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(-3.0774975) q[0];
rz(0.31072101) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(1.4583189) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1761988) q[0];
sx q[0];
rz(-1.5237234) q[0];
sx q[0];
rz(1.4233627) q[0];
rz(-pi) q[1];
rz(-2.4086191) q[2];
sx q[2];
rz(-0.79967116) q[2];
sx q[2];
rz(-0.016499585) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4255193) q[1];
sx q[1];
rz(-0.41285535) q[1];
sx q[1];
rz(3.0070261) q[1];
x q[2];
rz(-0.31483105) q[3];
sx q[3];
rz(-1.587084) q[3];
sx q[3];
rz(-2.6967238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98214275) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-2.4625051) q[2];
rz(-1.3726161) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452334) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(-1.1244208) q[0];
rz(-1.9212978) q[1];
sx q[1];
rz(-1.9492457) q[1];
sx q[1];
rz(1.6569998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9254018) q[0];
sx q[0];
rz(-1.5201709) q[0];
sx q[0];
rz(3.1111858) q[0];
rz(1.2414615) q[2];
sx q[2];
rz(-0.54154684) q[2];
sx q[2];
rz(-2.3878218) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2034519) q[1];
sx q[1];
rz(-1.005299) q[1];
sx q[1];
rz(-2.6184665) q[1];
rz(0.37064151) q[3];
sx q[3];
rz(-2.3834043) q[3];
sx q[3];
rz(1.7565808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1241887) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(-2.1253288) q[2];
rz(-1.7381564) q[3];
sx q[3];
rz(-1.4973463) q[3];
sx q[3];
rz(1.0884292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6338585) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(0.39988363) q[0];
rz(1.1625066) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(0.17366017) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7090209) q[0];
sx q[0];
rz(-3.0637494) q[0];
sx q[0];
rz(2.7709333) q[0];
rz(-pi) q[1];
rz(2.747614) q[2];
sx q[2];
rz(-1.7658965) q[2];
sx q[2];
rz(2.5815798) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9160737) q[1];
sx q[1];
rz(-1.5991296) q[1];
sx q[1];
rz(1.4654935) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3624304) q[3];
sx q[3];
rz(-1.8959798) q[3];
sx q[3];
rz(-1.1680101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6614723) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(-2.373467) q[2];
rz(-2.2875732) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1059234) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(-1.6712028) q[0];
rz(-0.51180965) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(1.1434198) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66577673) q[0];
sx q[0];
rz(-1.4205298) q[0];
sx q[0];
rz(1.6554324) q[0];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6900942) q[1];
sx q[1];
rz(-1.9376117) q[1];
sx q[1];
rz(0.18641285) q[1];
rz(-1.0295168) q[3];
sx q[3];
rz(-1.3068849) q[3];
sx q[3];
rz(-1.908386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4986971) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(-1.5303622) q[2];
rz(-1.4536084) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(-2.8924275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0221508) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(-1.4402333) q[0];
rz(0.72921905) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(-2.008332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0798577) q[0];
sx q[0];
rz(-1.2485866) q[0];
sx q[0];
rz(1.6155924) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5845756) q[2];
sx q[2];
rz(-2.2829977) q[2];
sx q[2];
rz(-2.6028002) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95722317) q[1];
sx q[1];
rz(-1.2979227) q[1];
sx q[1];
rz(-2.2437614) q[1];
rz(2.8797638) q[3];
sx q[3];
rz(-1.7228408) q[3];
sx q[3];
rz(2.9933628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7363654) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(-0.48842946) q[2];
rz(-1.8296261) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(2.5002938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064780386) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(3.0704165) q[0];
rz(0.03216234) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(1.2088998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9933388) q[0];
sx q[0];
rz(-2.3518666) q[0];
sx q[0];
rz(-1.8940311) q[0];
x q[1];
rz(-0.9135984) q[2];
sx q[2];
rz(-2.845394) q[2];
sx q[2];
rz(-1.9343455) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3370918) q[1];
sx q[1];
rz(-3.0102804) q[1];
sx q[1];
rz(-1.4485703) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62187059) q[3];
sx q[3];
rz(-1.8727881) q[3];
sx q[3];
rz(2.5627476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.20748392) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(1.5709546) q[2];
rz(-2.2682244) q[3];
sx q[3];
rz(-1.7375172) q[3];
sx q[3];
rz(2.7895555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1677925) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(-2.8299676) q[0];
rz(0.82178003) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(-1.5100381) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85743839) q[0];
sx q[0];
rz(-0.3785924) q[0];
sx q[0];
rz(-2.5749102) q[0];
x q[1];
rz(-0.19599234) q[2];
sx q[2];
rz(-1.9078983) q[2];
sx q[2];
rz(-0.54346426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.65731243) q[1];
sx q[1];
rz(-1.7727976) q[1];
sx q[1];
rz(0.61985086) q[1];
rz(-pi) q[2];
rz(-1.3872836) q[3];
sx q[3];
rz(-1.7997051) q[3];
sx q[3];
rz(-2.6218417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(2.3256425) q[2];
rz(0.50968918) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(1.9036487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9394768) q[0];
sx q[0];
rz(-2.3411223) q[0];
sx q[0];
rz(-2.4329484) q[0];
rz(1.2344822) q[2];
sx q[2];
rz(-2.2173777) q[2];
sx q[2];
rz(-0.99415776) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60724466) q[1];
sx q[1];
rz(-1.4048647) q[1];
sx q[1];
rz(1.0870618) q[1];
x q[2];
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
rz(-pi/2) q[1];
rz(-2.4663503) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(-0.33774439) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.7815536) q[3];
sx q[3];
rz(-1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.9032003) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(-0.94454371) q[2];
sx q[2];
rz(-1.7164451) q[2];
sx q[2];
rz(-3.0838983) q[2];
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
