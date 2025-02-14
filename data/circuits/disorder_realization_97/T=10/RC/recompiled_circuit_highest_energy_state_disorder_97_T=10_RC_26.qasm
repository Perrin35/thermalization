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
rz(-0.47973862) q[0];
sx q[0];
rz(-2.0070183) q[0];
sx q[0];
rz(-1.467508) q[0];
rz(1.6788586) q[1];
sx q[1];
rz(-2.1958513) q[1];
sx q[1];
rz(0.50409627) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5228793) q[0];
sx q[0];
rz(-1.0600495) q[0];
sx q[0];
rz(-2.5070531) q[0];
x q[1];
rz(0.43457793) q[2];
sx q[2];
rz(-2.761026) q[2];
sx q[2];
rz(-1.5411045) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2254653) q[1];
sx q[1];
rz(-2.230495) q[1];
sx q[1];
rz(-0.74972357) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.619574) q[3];
sx q[3];
rz(-1.2360473) q[3];
sx q[3];
rz(2.3693645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41245875) q[2];
sx q[2];
rz(-1.2428186) q[2];
sx q[2];
rz(2.9553555) q[2];
rz(-1.3650182) q[3];
sx q[3];
rz(-2.5297207) q[3];
sx q[3];
rz(-1.2046825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3314826) q[0];
sx q[0];
rz(-0.45833603) q[0];
sx q[0];
rz(-3.0895184) q[0];
rz(-2.1049818) q[1];
sx q[1];
rz(-0.60984817) q[1];
sx q[1];
rz(-0.61526543) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4724169) q[0];
sx q[0];
rz(-0.76509464) q[0];
sx q[0];
rz(-1.0932176) q[0];
rz(-pi) q[1];
rz(-2.9040163) q[2];
sx q[2];
rz(-1.4427836) q[2];
sx q[2];
rz(1.0281171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8456455) q[1];
sx q[1];
rz(-1.6270279) q[1];
sx q[1];
rz(2.2241431) q[1];
rz(-pi) q[2];
rz(1.7795981) q[3];
sx q[3];
rz(-0.80392814) q[3];
sx q[3];
rz(-1.7957791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44876862) q[2];
sx q[2];
rz(-1.8961467) q[2];
sx q[2];
rz(-1.9057062) q[2];
rz(-1.1290733) q[3];
sx q[3];
rz(-0.90714199) q[3];
sx q[3];
rz(-0.83704078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6579984) q[0];
sx q[0];
rz(-0.75355419) q[0];
sx q[0];
rz(-1.3760706) q[0];
rz(-2.6013069) q[1];
sx q[1];
rz(-1.0397725) q[1];
sx q[1];
rz(-1.4556063) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3106011) q[0];
sx q[0];
rz(-1.6414223) q[0];
sx q[0];
rz(0.40343209) q[0];
rz(-0.20353453) q[2];
sx q[2];
rz(-1.6586721) q[2];
sx q[2];
rz(0.7754625) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1108577) q[1];
sx q[1];
rz(-2.1848715) q[1];
sx q[1];
rz(-0.43621896) q[1];
rz(-pi) q[2];
rz(-0.40850477) q[3];
sx q[3];
rz(-0.91445476) q[3];
sx q[3];
rz(-2.07043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49440631) q[2];
sx q[2];
rz(-2.8072085) q[2];
sx q[2];
rz(1.6953267) q[2];
rz(1.1930195) q[3];
sx q[3];
rz(-2.1665067) q[3];
sx q[3];
rz(-1.4421991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0874264) q[0];
sx q[0];
rz(-2.8717201) q[0];
sx q[0];
rz(0.42359459) q[0];
rz(0.95701796) q[1];
sx q[1];
rz(-1.3514163) q[1];
sx q[1];
rz(-0.097600309) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4159629) q[0];
sx q[0];
rz(-1.0495674) q[0];
sx q[0];
rz(2.002191) q[0];
x q[1];
rz(-1.3988858) q[2];
sx q[2];
rz(-0.70016501) q[2];
sx q[2];
rz(-1.8309878) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.017189055) q[1];
sx q[1];
rz(-1.2252955) q[1];
sx q[1];
rz(-2.8082737) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37577263) q[3];
sx q[3];
rz(-2.3762581) q[3];
sx q[3];
rz(-2.1735512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61509722) q[2];
sx q[2];
rz(-0.47250938) q[2];
sx q[2];
rz(0.16461593) q[2];
rz(0.047867157) q[3];
sx q[3];
rz(-1.7367626) q[3];
sx q[3];
rz(0.78711787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1614302) q[0];
sx q[0];
rz(-0.4466559) q[0];
sx q[0];
rz(1.584378) q[0];
rz(-2.7730675) q[1];
sx q[1];
rz(-1.3216113) q[1];
sx q[1];
rz(1.6065067) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7212352) q[0];
sx q[0];
rz(-2.0518655) q[0];
sx q[0];
rz(2.1522983) q[0];
rz(-pi) q[1];
rz(-2.755411) q[2];
sx q[2];
rz(-2.2791822) q[2];
sx q[2];
rz(-0.86175534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6905744) q[1];
sx q[1];
rz(-1.3214045) q[1];
sx q[1];
rz(1.0728157) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3103043) q[3];
sx q[3];
rz(-0.29803571) q[3];
sx q[3];
rz(-0.71825114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2313472) q[2];
sx q[2];
rz(-0.37800899) q[2];
sx q[2];
rz(-2.0965915) q[2];
rz(2.9735978) q[3];
sx q[3];
rz(-1.1863656) q[3];
sx q[3];
rz(0.35429889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0198233) q[0];
sx q[0];
rz(-2.1124117) q[0];
sx q[0];
rz(1.2044915) q[0];
rz(-0.80351859) q[1];
sx q[1];
rz(-0.81356994) q[1];
sx q[1];
rz(2.4258851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7835282) q[0];
sx q[0];
rz(-0.77371374) q[0];
sx q[0];
rz(2.2893002) q[0];
x q[1];
rz(-1.3722695) q[2];
sx q[2];
rz(-2.0110705) q[2];
sx q[2];
rz(-1.0834875) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83166158) q[1];
sx q[1];
rz(-0.67218053) q[1];
sx q[1];
rz(-2.9596555) q[1];
rz(2.1200772) q[3];
sx q[3];
rz(-2.6166354) q[3];
sx q[3];
rz(-2.9385645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7314926) q[2];
sx q[2];
rz(-2.4918753) q[2];
sx q[2];
rz(2.0984207) q[2];
rz(0.10797524) q[3];
sx q[3];
rz(-0.94111809) q[3];
sx q[3];
rz(2.030355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1825948) q[0];
sx q[0];
rz(-1.3866871) q[0];
sx q[0];
rz(-1.9166272) q[0];
rz(-2.909868) q[1];
sx q[1];
rz(-2.2544315) q[1];
sx q[1];
rz(-1.8909594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62364336) q[0];
sx q[0];
rz(-2.3139489) q[0];
sx q[0];
rz(1.4046937) q[0];
rz(3.1194889) q[2];
sx q[2];
rz(-1.7689509) q[2];
sx q[2];
rz(-2.1111859) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96737427) q[1];
sx q[1];
rz(-1.3506445) q[1];
sx q[1];
rz(-2.8323814) q[1];
rz(1.5579079) q[3];
sx q[3];
rz(-1.0143712) q[3];
sx q[3];
rz(-2.9146359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0687678) q[2];
sx q[2];
rz(-2.4710957) q[2];
sx q[2];
rz(-0.13216275) q[2];
rz(0.69563785) q[3];
sx q[3];
rz(-1.9591103) q[3];
sx q[3];
rz(-2.3343991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19145963) q[0];
sx q[0];
rz(-1.5605254) q[0];
sx q[0];
rz(2.4323442) q[0];
rz(-1.5059772) q[1];
sx q[1];
rz(-1.9874856) q[1];
sx q[1];
rz(1.0848612) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5800672) q[0];
sx q[0];
rz(-2.170772) q[0];
sx q[0];
rz(-1.5011494) q[0];
rz(-pi) q[1];
rz(-0.80471595) q[2];
sx q[2];
rz(-1.448481) q[2];
sx q[2];
rz(1.027077) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0604531) q[1];
sx q[1];
rz(-2.6067002) q[1];
sx q[1];
rz(0.51523955) q[1];
rz(1.510511) q[3];
sx q[3];
rz(-1.5593646) q[3];
sx q[3];
rz(2.722198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.52013493) q[2];
sx q[2];
rz(-2.1389213) q[2];
sx q[2];
rz(-1.3168859) q[2];
rz(-0.78553158) q[3];
sx q[3];
rz(-1.1979878) q[3];
sx q[3];
rz(0.42640105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.921628) q[0];
sx q[0];
rz(-1.6554609) q[0];
sx q[0];
rz(2.7332136) q[0];
rz(2.1829103) q[1];
sx q[1];
rz(-0.32902333) q[1];
sx q[1];
rz(-1.6201409) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0199838) q[0];
sx q[0];
rz(-0.43916288) q[0];
sx q[0];
rz(-2.1497221) q[0];
x q[1];
rz(-1.5217811) q[2];
sx q[2];
rz(-1.00178) q[2];
sx q[2];
rz(-2.4328277) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1679528) q[1];
sx q[1];
rz(-2.3427737) q[1];
sx q[1];
rz(-2.1583907) q[1];
rz(1.2580885) q[3];
sx q[3];
rz(-1.3888479) q[3];
sx q[3];
rz(3.0489717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7659605) q[2];
sx q[2];
rz(-1.1124632) q[2];
sx q[2];
rz(2.5755303) q[2];
rz(1.798897) q[3];
sx q[3];
rz(-0.415396) q[3];
sx q[3];
rz(0.28877637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.600243) q[0];
sx q[0];
rz(-0.039027795) q[0];
sx q[0];
rz(-1.716123) q[0];
rz(1.0936945) q[1];
sx q[1];
rz(-0.92233557) q[1];
sx q[1];
rz(2.7519382) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3411639) q[0];
sx q[0];
rz(-2.6376056) q[0];
sx q[0];
rz(2.1227073) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94691838) q[2];
sx q[2];
rz(-1.3585261) q[2];
sx q[2];
rz(-2.61643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0691955) q[1];
sx q[1];
rz(-2.715413) q[1];
sx q[1];
rz(1.6276433) q[1];
x q[2];
rz(-0.52441729) q[3];
sx q[3];
rz(-2.5811385) q[3];
sx q[3];
rz(2.2956212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23715544) q[2];
sx q[2];
rz(-2.6025786) q[2];
sx q[2];
rz(1.944444) q[2];
rz(0.14687471) q[3];
sx q[3];
rz(-2.4683888) q[3];
sx q[3];
rz(-0.98744121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7214397) q[0];
sx q[0];
rz(-1.0316105) q[0];
sx q[0];
rz(1.6461865) q[0];
rz(-2.738476) q[1];
sx q[1];
rz(-1.8164201) q[1];
sx q[1];
rz(1.4774189) q[1];
rz(2.5667122) q[2];
sx q[2];
rz(-1.4670062) q[2];
sx q[2];
rz(0.23679096) q[2];
rz(-2.7561989) q[3];
sx q[3];
rz(-1.3968395) q[3];
sx q[3];
rz(-2.2341961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
