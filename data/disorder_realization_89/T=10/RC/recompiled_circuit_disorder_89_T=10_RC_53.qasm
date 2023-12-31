OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.75582957) q[0];
sx q[0];
rz(-1.4094149) q[0];
sx q[0];
rz(-0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(-2.7198305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5428107) q[0];
sx q[0];
rz(-0.069617696) q[0];
sx q[0];
rz(-1.1845864) q[0];
rz(3.0938901) q[2];
sx q[2];
rz(-0.68714266) q[2];
sx q[2];
rz(1.3840165) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4946343) q[1];
sx q[1];
rz(-1.1263493) q[1];
sx q[1];
rz(1.7723945) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20538879) q[3];
sx q[3];
rz(-2.5377512) q[3];
sx q[3];
rz(2.3103034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1224147) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(1.0100693) q[2];
rz(1.646237) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0881969) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(2.2825867) q[0];
rz(2.7711218) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.6765615) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37061781) q[0];
sx q[0];
rz(-0.3807225) q[0];
sx q[0];
rz(-2.5230663) q[0];
rz(-pi) q[1];
rz(-2.6071762) q[2];
sx q[2];
rz(-2.538531) q[2];
sx q[2];
rz(2.7433928) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.78039353) q[1];
sx q[1];
rz(-0.91460278) q[1];
sx q[1];
rz(1.111163) q[1];
rz(-0.46918418) q[3];
sx q[3];
rz(-0.57933925) q[3];
sx q[3];
rz(-0.34512305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4375962) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(1.0393418) q[3];
sx q[3];
rz(-1.9266409) q[3];
sx q[3];
rz(-0.59282747) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(0.16361374) q[0];
rz(0.71584654) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(-1.3735501) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93689504) q[0];
sx q[0];
rz(-2.097313) q[0];
sx q[0];
rz(1.599647) q[0];
rz(-0.85767304) q[2];
sx q[2];
rz(-1.8048986) q[2];
sx q[2];
rz(2.8785994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.19792412) q[1];
sx q[1];
rz(-1.0218044) q[1];
sx q[1];
rz(-2.3180069) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80188607) q[3];
sx q[3];
rz(-2.826414) q[3];
sx q[3];
rz(-1.8568045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(-2.1219357) q[2];
rz(2.1608458) q[3];
sx q[3];
rz(-2.5648983) q[3];
sx q[3];
rz(0.61292928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402128) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(2.3024978) q[0];
rz(1.6038731) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(0.61027169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2903039) q[0];
sx q[0];
rz(-1.5820869) q[0];
sx q[0];
rz(1.6800866) q[0];
rz(-pi) q[1];
rz(0.35595603) q[2];
sx q[2];
rz(-1.2095371) q[2];
sx q[2];
rz(-1.4796096) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1000881) q[1];
sx q[1];
rz(-1.3106292) q[1];
sx q[1];
rz(-0.86160223) q[1];
rz(-pi) q[2];
rz(-2.1498333) q[3];
sx q[3];
rz(-1.1513396) q[3];
sx q[3];
rz(-1.6384244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.443976) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(2.9042517) q[2];
rz(2.234263) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(-1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6326555) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(1.6089815) q[0];
rz(-0.73348796) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(-1.8283432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26966306) q[0];
sx q[0];
rz(-0.73686826) q[0];
sx q[0];
rz(-1.0760197) q[0];
rz(-pi) q[1];
rz(-1.6275431) q[2];
sx q[2];
rz(-1.5350255) q[2];
sx q[2];
rz(-2.2030731) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.207706) q[1];
sx q[1];
rz(-1.943214) q[1];
sx q[1];
rz(1.1697342) q[1];
rz(2.5562708) q[3];
sx q[3];
rz(-2.1748741) q[3];
sx q[3];
rz(-2.4620591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0017172) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(0.66429794) q[2];
rz(2.4364566) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(0.10372182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649081) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(-0.054811906) q[0];
rz(1.0143771) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(-2.6584113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85855243) q[0];
sx q[0];
rz(-2.9498093) q[0];
sx q[0];
rz(1.900308) q[0];
rz(-pi) q[1];
x q[1];
rz(0.032776253) q[2];
sx q[2];
rz(-1.1083229) q[2];
sx q[2];
rz(-2.0605161) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.253787) q[1];
sx q[1];
rz(-1.0711545) q[1];
sx q[1];
rz(-0.81412022) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0275434) q[3];
sx q[3];
rz(-1.9012791) q[3];
sx q[3];
rz(-1.198311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15423916) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(-0.085263578) q[2];
rz(-1.9412458) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7744556) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(-2.1869587) q[0];
rz(0.57149354) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(0.25209299) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.945767) q[0];
sx q[0];
rz(-2.6951163) q[0];
sx q[0];
rz(2.8004942) q[0];
rz(-0.78105314) q[2];
sx q[2];
rz(-0.99596802) q[2];
sx q[2];
rz(-0.99036723) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.568087) q[1];
sx q[1];
rz(-1.2099669) q[1];
sx q[1];
rz(-0.87330694) q[1];
rz(-pi) q[2];
rz(2.0310568) q[3];
sx q[3];
rz(-1.2292976) q[3];
sx q[3];
rz(-0.5865435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38348848) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(-0.90448109) q[2];
rz(1.0036184) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43338183) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(1.7277539) q[0];
rz(-0.66043234) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(-1.3716912) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6721281) q[0];
sx q[0];
rz(-2.1882595) q[0];
sx q[0];
rz(0.10829627) q[0];
rz(-0.22050627) q[2];
sx q[2];
rz(-0.77324152) q[2];
sx q[2];
rz(-0.79794937) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.129442) q[1];
sx q[1];
rz(-2.4852942) q[1];
sx q[1];
rz(0.49883573) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36088862) q[3];
sx q[3];
rz(-1.4286255) q[3];
sx q[3];
rz(1.7796381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82683212) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(0.39789847) q[2];
rz(0.49063101) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9234377) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(-1.461347) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(-0.27059069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4280677) q[0];
sx q[0];
rz(-1.7881835) q[0];
sx q[0];
rz(-0.2268976) q[0];
rz(2.9491896) q[2];
sx q[2];
rz(-1.9115598) q[2];
sx q[2];
rz(2.2601688) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8671659) q[1];
sx q[1];
rz(-1.567489) q[1];
sx q[1];
rz(-0.13722384) q[1];
rz(-pi) q[2];
rz(-1.7646592) q[3];
sx q[3];
rz(-2.6541775) q[3];
sx q[3];
rz(0.91538115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33971912) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(1.6516997) q[2];
rz(2.4387032) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(-1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(2.9872966) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(2.399209) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5267747) q[0];
sx q[0];
rz(-2.1702538) q[0];
sx q[0];
rz(-2.4012698) q[0];
rz(-pi) q[1];
rz(0.97147271) q[2];
sx q[2];
rz(-0.68442548) q[2];
sx q[2];
rz(1.3542092) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4780477) q[1];
sx q[1];
rz(-1.5819342) q[1];
sx q[1];
rz(-0.34578172) q[1];
rz(2.9418482) q[3];
sx q[3];
rz(-0.14530694) q[3];
sx q[3];
rz(2.6655281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8538889) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(-1.5819736) q[2];
rz(2.93086) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(-0.5334841) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29006526) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(-1.3416946) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(2.3201597) q[2];
sx q[2];
rz(-1.3315622) q[2];
sx q[2];
rz(-0.17377725) q[2];
rz(0.13989862) q[3];
sx q[3];
rz(-1.4362873) q[3];
sx q[3];
rz(0.39713058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
