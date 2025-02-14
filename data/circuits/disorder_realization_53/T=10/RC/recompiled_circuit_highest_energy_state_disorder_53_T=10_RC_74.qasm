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
rz(-0.80225575) q[0];
sx q[0];
rz(-1.7576317) q[0];
sx q[0];
rz(1.2686165) q[0];
rz(0.4624548) q[1];
sx q[1];
rz(6.3422536) q[1];
sx q[1];
rz(8.0191945) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.424762) q[0];
sx q[0];
rz(-0.66436902) q[0];
sx q[0];
rz(0.82358255) q[0];
x q[1];
rz(1.8008158) q[2];
sx q[2];
rz(-0.91311306) q[2];
sx q[2];
rz(0.44154007) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4753603) q[1];
sx q[1];
rz(-1.9812498) q[1];
sx q[1];
rz(-0.79373856) q[1];
rz(0.89449331) q[3];
sx q[3];
rz(-1.749012) q[3];
sx q[3];
rz(-1.3823929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9872226) q[2];
sx q[2];
rz(-1.9742249) q[2];
sx q[2];
rz(-0.78987375) q[2];
rz(-0.023898276) q[3];
sx q[3];
rz(-1.8022715) q[3];
sx q[3];
rz(0.53643119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8928878) q[0];
sx q[0];
rz(-1.4905812) q[0];
sx q[0];
rz(1.8875246) q[0];
rz(-1.6528543) q[1];
sx q[1];
rz(-2.2658927) q[1];
sx q[1];
rz(2.658433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7308759) q[0];
sx q[0];
rz(-0.19598254) q[0];
sx q[0];
rz(-1.8394801) q[0];
rz(-2.580609) q[2];
sx q[2];
rz(-0.81811726) q[2];
sx q[2];
rz(2.4647692) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70139685) q[1];
sx q[1];
rz(-0.45397511) q[1];
sx q[1];
rz(-0.63164901) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53584953) q[3];
sx q[3];
rz(-2.8070309) q[3];
sx q[3];
rz(1.1862792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.88863215) q[2];
sx q[2];
rz(-1.6855719) q[2];
sx q[2];
rz(1.6271094) q[2];
rz(-0.0557946) q[3];
sx q[3];
rz(-1.5443708) q[3];
sx q[3];
rz(-1.4896711) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2273939) q[0];
sx q[0];
rz(-0.44454235) q[0];
sx q[0];
rz(-3.0134873) q[0];
rz(-0.57614342) q[1];
sx q[1];
rz(-2.4100401) q[1];
sx q[1];
rz(2.3760956) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2554995) q[0];
sx q[0];
rz(-1.0168726) q[0];
sx q[0];
rz(1.163068) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1942991) q[2];
sx q[2];
rz(-2.2222179) q[2];
sx q[2];
rz(-2.5848856) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.24352267) q[1];
sx q[1];
rz(-2.0546034) q[1];
sx q[1];
rz(1.0462199) q[1];
rz(0.14344826) q[3];
sx q[3];
rz(-1.5960403) q[3];
sx q[3];
rz(-2.6526566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3196044) q[2];
sx q[2];
rz(-0.87279785) q[2];
sx q[2];
rz(2.2785211) q[2];
rz(-0.34267628) q[3];
sx q[3];
rz(-2.7211012) q[3];
sx q[3];
rz(1.6531403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89219013) q[0];
sx q[0];
rz(-0.96977314) q[0];
sx q[0];
rz(-0.48200193) q[0];
rz(0.30872289) q[1];
sx q[1];
rz(-0.32544193) q[1];
sx q[1];
rz(2.5293005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4929313) q[0];
sx q[0];
rz(-1.0928286) q[0];
sx q[0];
rz(-0.98832358) q[0];
rz(1.4050414) q[2];
sx q[2];
rz(-3.0566772) q[2];
sx q[2];
rz(1.9254534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4548851) q[1];
sx q[1];
rz(-0.60777449) q[1];
sx q[1];
rz(-1.4922423) q[1];
rz(-pi) q[2];
x q[2];
rz(0.044030146) q[3];
sx q[3];
rz(-0.14942871) q[3];
sx q[3];
rz(1.7946515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.65773949) q[2];
sx q[2];
rz(-0.95119363) q[2];
sx q[2];
rz(1.8335906) q[2];
rz(-0.43404239) q[3];
sx q[3];
rz(-1.3515892) q[3];
sx q[3];
rz(-1.8724117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9798715) q[0];
sx q[0];
rz(-1.5890108) q[0];
sx q[0];
rz(2.8671434) q[0];
rz(-1.5212003) q[1];
sx q[1];
rz(-1.043964) q[1];
sx q[1];
rz(1.257198) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8306106) q[0];
sx q[0];
rz(-0.3468967) q[0];
sx q[0];
rz(2.3603289) q[0];
rz(-2.1047701) q[2];
sx q[2];
rz(-1.2495012) q[2];
sx q[2];
rz(-2.9095813) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2272531) q[1];
sx q[1];
rz(-1.3554595) q[1];
sx q[1];
rz(-2.1244177) q[1];
rz(0.32088478) q[3];
sx q[3];
rz(-2.2498331) q[3];
sx q[3];
rz(1.3037552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7460798) q[2];
sx q[2];
rz(-1.505625) q[2];
sx q[2];
rz(2.5858322) q[2];
rz(0.79103509) q[3];
sx q[3];
rz(-1.0020703) q[3];
sx q[3];
rz(1.5033495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.898734) q[0];
sx q[0];
rz(-1.4396311) q[0];
sx q[0];
rz(-0.71204251) q[0];
rz(0.60246077) q[1];
sx q[1];
rz(-2.7280877) q[1];
sx q[1];
rz(-3.0625524) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9475381) q[0];
sx q[0];
rz(-1.506516) q[0];
sx q[0];
rz(0.044919515) q[0];
rz(-1.6557793) q[2];
sx q[2];
rz(-1.4956317) q[2];
sx q[2];
rz(0.70460457) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.95922071) q[1];
sx q[1];
rz(-2.3010265) q[1];
sx q[1];
rz(1.8664136) q[1];
x q[2];
rz(2.8435264) q[3];
sx q[3];
rz(-1.6692999) q[3];
sx q[3];
rz(-2.5677266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.041212335) q[2];
sx q[2];
rz(-0.93116394) q[2];
sx q[2];
rz(0.97990123) q[2];
rz(1.6117217) q[3];
sx q[3];
rz(-1.2414705) q[3];
sx q[3];
rz(-0.25947586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.74705446) q[0];
sx q[0];
rz(-1.8600445) q[0];
sx q[0];
rz(-0.10093149) q[0];
rz(-0.55502564) q[1];
sx q[1];
rz(-0.9340159) q[1];
sx q[1];
rz(-1.2947882) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0880343) q[0];
sx q[0];
rz(-0.71672601) q[0];
sx q[0];
rz(-1.929856) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1406222) q[2];
sx q[2];
rz(-1.2764659) q[2];
sx q[2];
rz(-0.4145588) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6375512) q[1];
sx q[1];
rz(-2.6890272) q[1];
sx q[1];
rz(0.33337793) q[1];
rz(-1.0247308) q[3];
sx q[3];
rz(-2.6203558) q[3];
sx q[3];
rz(1.2574399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9635705) q[2];
sx q[2];
rz(-2.2691085) q[2];
sx q[2];
rz(-2.8998609) q[2];
rz(-2.669615) q[3];
sx q[3];
rz(-0.21502544) q[3];
sx q[3];
rz(1.5836466) q[3];
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
rz(-pi/2) q[0];
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
rz(0.80411512) q[0];
sx q[0];
rz(-0.98144704) q[0];
sx q[0];
rz(-0.61019439) q[0];
rz(-0.21774165) q[1];
sx q[1];
rz(-0.95548958) q[1];
sx q[1];
rz(-2.311923) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1581081) q[0];
sx q[0];
rz(-2.2451441) q[0];
sx q[0];
rz(0.83263065) q[0];
x q[1];
rz(-2.81189) q[2];
sx q[2];
rz(-1.9118983) q[2];
sx q[2];
rz(-1.4280048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9431994) q[1];
sx q[1];
rz(-1.3855318) q[1];
sx q[1];
rz(2.7279502) q[1];
rz(-1.4368966) q[3];
sx q[3];
rz(-2.4256676) q[3];
sx q[3];
rz(-1.1435103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3024451) q[2];
sx q[2];
rz(-1.7586917) q[2];
sx q[2];
rz(-2.006532) q[2];
rz(2.6531687) q[3];
sx q[3];
rz(-1.4819744) q[3];
sx q[3];
rz(1.2357014) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9957073) q[0];
sx q[0];
rz(-1.1111525) q[0];
sx q[0];
rz(0.97440326) q[0];
rz(-2.3459332) q[1];
sx q[1];
rz(-0.37745044) q[1];
sx q[1];
rz(2.5239351) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97516135) q[0];
sx q[0];
rz(-1.5733061) q[0];
sx q[0];
rz(3.0581492) q[0];
rz(1.1716003) q[2];
sx q[2];
rz(-2.0621057) q[2];
sx q[2];
rz(1.0688865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1682081) q[1];
sx q[1];
rz(-2.1323556) q[1];
sx q[1];
rz(-1.8902768) q[1];
rz(-pi) q[2];
rz(-0.32989721) q[3];
sx q[3];
rz(-1.1217692) q[3];
sx q[3];
rz(-2.2638418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.016102942) q[2];
sx q[2];
rz(-1.2697271) q[2];
sx q[2];
rz(-1.7529091) q[2];
rz(-1.8359418) q[3];
sx q[3];
rz(-2.4180222) q[3];
sx q[3];
rz(1.2828264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1284803) q[0];
sx q[0];
rz(-2.4025669) q[0];
sx q[0];
rz(-0.59984961) q[0];
rz(0.57303095) q[1];
sx q[1];
rz(-2.532798) q[1];
sx q[1];
rz(2.1175687) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1896418) q[0];
sx q[0];
rz(-1.1770403) q[0];
sx q[0];
rz(-1.1313361) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5285792) q[2];
sx q[2];
rz(-2.2229872) q[2];
sx q[2];
rz(3.1162645) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3249859) q[1];
sx q[1];
rz(-1.4868951) q[1];
sx q[1];
rz(-2.744426) q[1];
x q[2];
rz(-2.7598466) q[3];
sx q[3];
rz(-1.4305528) q[3];
sx q[3];
rz(0.5776757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.92274252) q[2];
sx q[2];
rz(-1.7538193) q[2];
sx q[2];
rz(-2.3663523) q[2];
rz(1.5357337) q[3];
sx q[3];
rz(-1.9550867) q[3];
sx q[3];
rz(1.2576013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53104644) q[0];
sx q[0];
rz(-0.95120593) q[0];
sx q[0];
rz(-1.9134941) q[0];
rz(1.3950521) q[1];
sx q[1];
rz(-1.6680622) q[1];
sx q[1];
rz(-0.39573085) q[1];
rz(-2.5644518) q[2];
sx q[2];
rz(-0.52684035) q[2];
sx q[2];
rz(3.1160804) q[2];
rz(-0.80202924) q[3];
sx q[3];
rz(-1.5939972) q[3];
sx q[3];
rz(0.93929285) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
