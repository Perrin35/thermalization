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
rz(1.9011693) q[0];
sx q[0];
rz(-1.9440396) q[0];
sx q[0];
rz(-2.9252606) q[0];
rz(-2.1842015) q[1];
sx q[1];
rz(-0.67617813) q[1];
sx q[1];
rz(1.6300936) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0441598) q[0];
sx q[0];
rz(-2.5852647) q[0];
sx q[0];
rz(-3.1233643) q[0];
rz(-pi) q[1];
rz(-0.58403973) q[2];
sx q[2];
rz(-0.55070832) q[2];
sx q[2];
rz(-2.5322897) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18876702) q[1];
sx q[1];
rz(-2.5270542) q[1];
sx q[1];
rz(-0.60093083) q[1];
x q[2];
rz(2.0262296) q[3];
sx q[3];
rz(-0.75616992) q[3];
sx q[3];
rz(-2.2952473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.034885255) q[2];
sx q[2];
rz(-2.7896176) q[2];
sx q[2];
rz(2.0931639) q[2];
rz(-0.18167051) q[3];
sx q[3];
rz(-0.964966) q[3];
sx q[3];
rz(-0.65339965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6492017) q[0];
sx q[0];
rz(-0.94434706) q[0];
sx q[0];
rz(-2.6948068) q[0];
rz(-2.2611639) q[1];
sx q[1];
rz(-1.3648405) q[1];
sx q[1];
rz(-2.3588038) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.114417) q[0];
sx q[0];
rz(-1.8150629) q[0];
sx q[0];
rz(0.063112325) q[0];
rz(-pi) q[1];
rz(1.9982463) q[2];
sx q[2];
rz(-2.1437217) q[2];
sx q[2];
rz(1.7597511) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0228717) q[1];
sx q[1];
rz(-1.1069655) q[1];
sx q[1];
rz(1.9074341) q[1];
rz(-2.1332333) q[3];
sx q[3];
rz(-2.376412) q[3];
sx q[3];
rz(2.6443114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.030152628) q[2];
sx q[2];
rz(-1.4806662) q[2];
sx q[2];
rz(-0.97935575) q[2];
rz(0.3848981) q[3];
sx q[3];
rz(-1.9210457) q[3];
sx q[3];
rz(0.16429193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48591831) q[0];
sx q[0];
rz(-3.0874708) q[0];
sx q[0];
rz(0.78980494) q[0];
rz(-2.9549331) q[1];
sx q[1];
rz(-1.4013545) q[1];
sx q[1];
rz(-2.1479215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45897608) q[0];
sx q[0];
rz(-1.0280711) q[0];
sx q[0];
rz(1.8629462) q[0];
rz(0.31618677) q[2];
sx q[2];
rz(-2.5044166) q[2];
sx q[2];
rz(-3.1115465) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6216065) q[1];
sx q[1];
rz(-2.4944759) q[1];
sx q[1];
rz(2.3220329) q[1];
rz(-pi) q[2];
rz(0.72088269) q[3];
sx q[3];
rz(-1.7768806) q[3];
sx q[3];
rz(-0.64180798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32187244) q[2];
sx q[2];
rz(-1.1178144) q[2];
sx q[2];
rz(2.7703088) q[2];
rz(-2.7927981) q[3];
sx q[3];
rz(-1.0943509) q[3];
sx q[3];
rz(-2.546052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45469859) q[0];
sx q[0];
rz(-2.1214387) q[0];
sx q[0];
rz(-1.7373079) q[0];
rz(2.4644201) q[1];
sx q[1];
rz(-1.155747) q[1];
sx q[1];
rz(1.6489395) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20805222) q[0];
sx q[0];
rz(-1.786199) q[0];
sx q[0];
rz(-2.893704) q[0];
rz(1.1392966) q[2];
sx q[2];
rz(-1.8427263) q[2];
sx q[2];
rz(1.4828585) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6552393) q[1];
sx q[1];
rz(-2.7905373) q[1];
sx q[1];
rz(0.69610657) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8294213) q[3];
sx q[3];
rz(-2.2709284) q[3];
sx q[3];
rz(2.2362102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6877785) q[2];
sx q[2];
rz(-2.1731589) q[2];
sx q[2];
rz(2.7933534) q[2];
rz(-1.4549152) q[3];
sx q[3];
rz(-1.7125407) q[3];
sx q[3];
rz(-1.7961563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9652047) q[0];
sx q[0];
rz(-2.5768953) q[0];
sx q[0];
rz(2.2494466) q[0];
rz(2.6773894) q[1];
sx q[1];
rz(-1.2449539) q[1];
sx q[1];
rz(-1.8468599) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2703122) q[0];
sx q[0];
rz(-2.0593606) q[0];
sx q[0];
rz(1.0315007) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53769298) q[2];
sx q[2];
rz(-1.6811996) q[2];
sx q[2];
rz(-0.7796208) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6695822) q[1];
sx q[1];
rz(-1.8947487) q[1];
sx q[1];
rz(1.1060017) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.283461) q[3];
sx q[3];
rz(-0.81908021) q[3];
sx q[3];
rz(-1.302296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8190454) q[2];
sx q[2];
rz(-0.61083856) q[2];
sx q[2];
rz(2.5757705) q[2];
rz(0.081341751) q[3];
sx q[3];
rz(-0.9809202) q[3];
sx q[3];
rz(-2.5210023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6767947) q[0];
sx q[0];
rz(-0.066635266) q[0];
sx q[0];
rz(1.5860522) q[0];
rz(2.0809035) q[1];
sx q[1];
rz(-1.5763177) q[1];
sx q[1];
rz(2.5097844) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9264797) q[0];
sx q[0];
rz(-2.8936912) q[0];
sx q[0];
rz(2.8126024) q[0];
rz(2.1366227) q[2];
sx q[2];
rz(-0.7938876) q[2];
sx q[2];
rz(-0.53295202) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8579156) q[1];
sx q[1];
rz(-1.2322958) q[1];
sx q[1];
rz(-2.1382911) q[1];
x q[2];
rz(1.4238213) q[3];
sx q[3];
rz(-0.89255652) q[3];
sx q[3];
rz(3.105046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.98171988) q[2];
sx q[2];
rz(-1.1184511) q[2];
sx q[2];
rz(0.49883207) q[2];
rz(-1.8286797) q[3];
sx q[3];
rz(-0.73176089) q[3];
sx q[3];
rz(1.4341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53443921) q[0];
sx q[0];
rz(-0.3648912) q[0];
sx q[0];
rz(2.4420807) q[0];
rz(0.46547678) q[1];
sx q[1];
rz(-0.87124467) q[1];
sx q[1];
rz(-0.93719283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1904966) q[0];
sx q[0];
rz(-1.7825923) q[0];
sx q[0];
rz(-1.7429211) q[0];
rz(-pi) q[1];
x q[1];
rz(1.466339) q[2];
sx q[2];
rz(-1.3323931) q[2];
sx q[2];
rz(-1.4739715) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4027462) q[1];
sx q[1];
rz(-2.7838696) q[1];
sx q[1];
rz(0.92836942) q[1];
rz(-pi) q[2];
rz(2.3714957) q[3];
sx q[3];
rz(-0.40989629) q[3];
sx q[3];
rz(0.65576762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.297544) q[2];
sx q[2];
rz(-1.8397477) q[2];
sx q[2];
rz(-0.37332264) q[2];
rz(-1.1897872) q[3];
sx q[3];
rz(-2.6326284) q[3];
sx q[3];
rz(-2.6313307) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74476403) q[0];
sx q[0];
rz(-1.0823534) q[0];
sx q[0];
rz(-2.3936791) q[0];
rz(-0.76639908) q[1];
sx q[1];
rz(-0.26873573) q[1];
sx q[1];
rz(-3.1386197) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.317694) q[0];
sx q[0];
rz(-2.4269773) q[0];
sx q[0];
rz(1.2750285) q[0];
rz(0.019549088) q[2];
sx q[2];
rz(-1.6416993) q[2];
sx q[2];
rz(2.8343458) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.06719477) q[1];
sx q[1];
rz(-1.7247827) q[1];
sx q[1];
rz(-3.0976553) q[1];
x q[2];
rz(-0.19488867) q[3];
sx q[3];
rz(-1.3759383) q[3];
sx q[3];
rz(0.080435924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.030674) q[2];
sx q[2];
rz(-1.7577533) q[2];
sx q[2];
rz(-2.0745011) q[2];
rz(3.057632) q[3];
sx q[3];
rz(-2.6559918) q[3];
sx q[3];
rz(-0.77897227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.344051) q[0];
sx q[0];
rz(-0.9386971) q[0];
sx q[0];
rz(0.10636605) q[0];
rz(-0.97995177) q[1];
sx q[1];
rz(-1.4915024) q[1];
sx q[1];
rz(0.76593691) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2625807) q[0];
sx q[0];
rz(-2.0011138) q[0];
sx q[0];
rz(0.58027123) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9418403) q[2];
sx q[2];
rz(-1.6589266) q[2];
sx q[2];
rz(2.5823808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83726604) q[1];
sx q[1];
rz(-0.89511739) q[1];
sx q[1];
rz(-1.6153687) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5670789) q[3];
sx q[3];
rz(-1.5848985) q[3];
sx q[3];
rz(-0.64025793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.47436675) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(-2.4412947) q[2];
rz(2.4750366) q[3];
sx q[3];
rz(-1.3349814) q[3];
sx q[3];
rz(-2.0535645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.67548442) q[0];
sx q[0];
rz(-1.0270783) q[0];
sx q[0];
rz(-1.745537) q[0];
rz(1.667977) q[1];
sx q[1];
rz(-1.8831848) q[1];
sx q[1];
rz(2.4748763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4748866) q[0];
sx q[0];
rz(-2.1726554) q[0];
sx q[0];
rz(-0.78420297) q[0];
x q[1];
rz(-2.4820508) q[2];
sx q[2];
rz(-1.7983984) q[2];
sx q[2];
rz(1.890581) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0565225) q[1];
sx q[1];
rz(-1.4398972) q[1];
sx q[1];
rz(-2.8225949) q[1];
rz(-pi) q[2];
rz(-0.57228831) q[3];
sx q[3];
rz(-2.0438571) q[3];
sx q[3];
rz(0.285404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0746158) q[2];
sx q[2];
rz(-0.65698996) q[2];
sx q[2];
rz(2.2508049) q[2];
rz(-2.7095419) q[3];
sx q[3];
rz(-1.0703577) q[3];
sx q[3];
rz(2.3945358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-2.4074832) q[0];
sx q[0];
rz(-1.4984087) q[0];
sx q[0];
rz(1.8468504) q[0];
rz(-1.2474077) q[1];
sx q[1];
rz(-2.4220962) q[1];
sx q[1];
rz(-1.9069506) q[1];
rz(-1.9966077) q[2];
sx q[2];
rz(-2.7663284) q[2];
sx q[2];
rz(-0.13079499) q[2];
rz(-3.0790764) q[3];
sx q[3];
rz(-0.48988952) q[3];
sx q[3];
rz(2.9174138) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
