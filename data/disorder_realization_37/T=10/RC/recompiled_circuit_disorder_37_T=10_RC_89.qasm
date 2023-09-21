OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0377169) q[0];
sx q[0];
rz(-1.2020943) q[0];
sx q[0];
rz(-1.9934959) q[0];
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(1.747945) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924661) q[0];
sx q[0];
rz(-1.2611715) q[0];
sx q[0];
rz(-1.5729088) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68140985) q[2];
sx q[2];
rz(-2.133956) q[2];
sx q[2];
rz(-1.6207221) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4268036) q[1];
sx q[1];
rz(-1.1183294) q[1];
sx q[1];
rz(1.2697551) q[1];
rz(-pi) q[2];
rz(-1.0127134) q[3];
sx q[3];
rz(-1.3289684) q[3];
sx q[3];
rz(-1.946132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48774886) q[2];
sx q[2];
rz(-1.2922492) q[2];
sx q[2];
rz(-0.1208819) q[2];
rz(0.17928784) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(2.9860935) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091846175) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(-0.13277408) q[0];
rz(-1.4615387) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(2.9002088) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3791083) q[0];
sx q[0];
rz(-2.603841) q[0];
sx q[0];
rz(0.82378973) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74626211) q[2];
sx q[2];
rz(-0.66821874) q[2];
sx q[2];
rz(-1.8909188) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.452255) q[1];
sx q[1];
rz(-0.66879767) q[1];
sx q[1];
rz(1.7178839) q[1];
rz(1.4199735) q[3];
sx q[3];
rz(-1.5449459) q[3];
sx q[3];
rz(-2.6005656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0138578) q[2];
sx q[2];
rz(-1.8227791) q[2];
sx q[2];
rz(-1.1068809) q[2];
rz(1.3876623) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36104193) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(0.74044359) q[0];
rz(0.45117798) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(2.0887451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6360639) q[0];
sx q[0];
rz(-1.9331421) q[0];
sx q[0];
rz(-0.58689582) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89715965) q[2];
sx q[2];
rz(-1.8067915) q[2];
sx q[2];
rz(1.2981851) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9991418) q[1];
sx q[1];
rz(-1.72292) q[1];
sx q[1];
rz(-1.692148) q[1];
x q[2];
rz(-2.0009082) q[3];
sx q[3];
rz(-1.579869) q[3];
sx q[3];
rz(0.7113925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.25049245) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(2.5946674) q[2];
rz(2.8524103) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.7893715) q[0];
rz(2.9317454) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(-2.9052177) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.198092) q[0];
sx q[0];
rz(-0.31790942) q[0];
sx q[0];
rz(0.99347465) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1430074) q[2];
sx q[2];
rz(-2.5479655) q[2];
sx q[2];
rz(-1.3015391) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35509767) q[1];
sx q[1];
rz(-1.7384643) q[1];
sx q[1];
rz(-3.1349036) q[1];
x q[2];
rz(-1.4218016) q[3];
sx q[3];
rz(-1.0862724) q[3];
sx q[3];
rz(-1.5414343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9005047) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(-0.55348712) q[2];
rz(-0.91529804) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(0.47167641) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(-0.70415235) q[0];
rz(-1.0559121) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(-0.59590894) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0376315) q[0];
sx q[0];
rz(-0.2346633) q[0];
sx q[0];
rz(1.1681359) q[0];
rz(2.9688409) q[2];
sx q[2];
rz(-2.4606332) q[2];
sx q[2];
rz(-1.3576042) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6223645) q[1];
sx q[1];
rz(-1.9740826) q[1];
sx q[1];
rz(-1.3695903) q[1];
rz(-1.3866053) q[3];
sx q[3];
rz(-1.2053688) q[3];
sx q[3];
rz(0.75152498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.492505) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(0.55111432) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(1.5392039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15774396) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(2.1898848) q[0];
rz(2.6668008) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(0.29528433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26074043) q[0];
sx q[0];
rz(-1.4319001) q[0];
sx q[0];
rz(1.9689346) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7252543) q[2];
sx q[2];
rz(-1.8784349) q[2];
sx q[2];
rz(2.6442106) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8413137) q[1];
sx q[1];
rz(-2.0804555) q[1];
sx q[1];
rz(-1.0213486) q[1];
x q[2];
rz(1.0601677) q[3];
sx q[3];
rz(-1.5139297) q[3];
sx q[3];
rz(2.1631654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39020145) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(-0.93377101) q[2];
rz(-1.6823403) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(-1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0657601) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(0.24205762) q[0];
rz(-2.4767955) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(-0.40329969) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4664073) q[0];
sx q[0];
rz(-1.9403606) q[0];
sx q[0];
rz(-2.7778366) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4629668) q[2];
sx q[2];
rz(-1.7947065) q[2];
sx q[2];
rz(-0.14397552) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2029018) q[1];
sx q[1];
rz(-2.9314329) q[1];
sx q[1];
rz(-1.5797257) q[1];
rz(0.025267406) q[3];
sx q[3];
rz(-1.380904) q[3];
sx q[3];
rz(0.59786284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.91281259) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(-2.7071803) q[2];
rz(1.0007535) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0077165724) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(2.6469321) q[0];
rz(-1.5085295) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(-0.60044926) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4279815) q[0];
sx q[0];
rz(-2.5117154) q[0];
sx q[0];
rz(-2.4226818) q[0];
rz(1.6161643) q[2];
sx q[2];
rz(-2.6160935) q[2];
sx q[2];
rz(-0.85277992) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8200092) q[1];
sx q[1];
rz(-1.6899741) q[1];
sx q[1];
rz(-0.18305852) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8377257) q[3];
sx q[3];
rz(-1.3715203) q[3];
sx q[3];
rz(2.5602333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0042469) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(0.11631575) q[2];
rz(-0.4256734) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(-1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15329926) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.4784038) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(-0.41752648) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077174515) q[0];
sx q[0];
rz(-1.6155956) q[0];
sx q[0];
rz(1.0391462) q[0];
rz(-pi) q[1];
rz(1.3067901) q[2];
sx q[2];
rz(-1.4679113) q[2];
sx q[2];
rz(-1.7057982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36996499) q[1];
sx q[1];
rz(-1.4039478) q[1];
sx q[1];
rz(-0.19373993) q[1];
rz(-pi) q[2];
rz(-2.0017654) q[3];
sx q[3];
rz(-1.3943854) q[3];
sx q[3];
rz(-2.924502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4801165) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(-0.49368668) q[2];
rz(0.72475973) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(-2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9676554) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(0.68558145) q[0];
rz(2.8441692) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(-1.1313653) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81998527) q[0];
sx q[0];
rz(-1.2258343) q[0];
sx q[0];
rz(-2.5489775) q[0];
rz(-2.4360043) q[2];
sx q[2];
rz(-1.0552647) q[2];
sx q[2];
rz(2.3956092) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86913925) q[1];
sx q[1];
rz(-1.1932045) q[1];
sx q[1];
rz(1.8121522) q[1];
rz(-2.0496619) q[3];
sx q[3];
rz(-2.1225956) q[3];
sx q[3];
rz(-2.5294876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3161105) q[2];
sx q[2];
rz(-1.9448514) q[2];
sx q[2];
rz(-2.4342009) q[2];
rz(-2.3317544) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(-2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.0512882) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(-2.9293625) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(-0.99954188) q[2];
sx q[2];
rz(-2.0959601) q[2];
sx q[2];
rz(-3.0053896) q[2];
rz(2.2255696) q[3];
sx q[3];
rz(-2.3229204) q[3];
sx q[3];
rz(-1.7182072) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
