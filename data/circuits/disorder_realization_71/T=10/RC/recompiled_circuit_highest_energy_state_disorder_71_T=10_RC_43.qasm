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
rz(3.0687301) q[0];
sx q[0];
rz(-1.0412359) q[0];
sx q[0];
rz(-0.6676724) q[0];
rz(3.0296037) q[1];
sx q[1];
rz(-1.4854687) q[1];
sx q[1];
rz(2.8356584) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95655694) q[0];
sx q[0];
rz(-1.4887012) q[0];
sx q[0];
rz(2.7285485) q[0];
x q[1];
rz(2.5247113) q[2];
sx q[2];
rz(-1.9494257) q[2];
sx q[2];
rz(-0.69064349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0623603) q[1];
sx q[1];
rz(-1.1992593) q[1];
sx q[1];
rz(1.6420235) q[1];
x q[2];
rz(-2.4351727) q[3];
sx q[3];
rz(-1.6230246) q[3];
sx q[3];
rz(2.4097648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2970807) q[2];
sx q[2];
rz(-1.7252555) q[2];
sx q[2];
rz(-2.1616006) q[2];
rz(-0.32198191) q[3];
sx q[3];
rz(-1.201509) q[3];
sx q[3];
rz(0.15596786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62671536) q[0];
sx q[0];
rz(-1.5276271) q[0];
sx q[0];
rz(-0.8901341) q[0];
rz(-2.0251677) q[1];
sx q[1];
rz(-2.2869488) q[1];
sx q[1];
rz(1.1847624) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6044004) q[0];
sx q[0];
rz(-2.064115) q[0];
sx q[0];
rz(1.7866578) q[0];
rz(-pi) q[1];
rz(1.0035788) q[2];
sx q[2];
rz(-0.62892709) q[2];
sx q[2];
rz(0.24404446) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5819232) q[1];
sx q[1];
rz(-1.6861177) q[1];
sx q[1];
rz(1.6827464) q[1];
rz(-pi) q[2];
rz(-3.1050097) q[3];
sx q[3];
rz(-0.42538339) q[3];
sx q[3];
rz(-0.98664397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6637791) q[2];
sx q[2];
rz(-2.9728643) q[2];
sx q[2];
rz(1.7185877) q[2];
rz(2.6723828) q[3];
sx q[3];
rz(-1.6459295) q[3];
sx q[3];
rz(2.7928228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94217268) q[0];
sx q[0];
rz(-1.1976137) q[0];
sx q[0];
rz(-2.0592101) q[0];
rz(-0.17467817) q[1];
sx q[1];
rz(-1.382261) q[1];
sx q[1];
rz(-2.6686525) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4648165) q[0];
sx q[0];
rz(-1.4107786) q[0];
sx q[0];
rz(-1.141044) q[0];
rz(2.7874417) q[2];
sx q[2];
rz(-2.2406497) q[2];
sx q[2];
rz(2.4669569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.096361209) q[1];
sx q[1];
rz(-1.0013818) q[1];
sx q[1];
rz(3.1147733) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93227902) q[3];
sx q[3];
rz(-2.4014056) q[3];
sx q[3];
rz(2.1719751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.48367721) q[2];
sx q[2];
rz(-1.4853073) q[2];
sx q[2];
rz(1.8243054) q[2];
rz(0.45051908) q[3];
sx q[3];
rz(-2.2299288) q[3];
sx q[3];
rz(1.7821504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4860185) q[0];
sx q[0];
rz(-2.0756606) q[0];
sx q[0];
rz(2.4701212) q[0];
rz(-0.94054049) q[1];
sx q[1];
rz(-2.7214919) q[1];
sx q[1];
rz(-0.95983517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67755015) q[0];
sx q[0];
rz(-1.3565956) q[0];
sx q[0];
rz(-1.1826128) q[0];
rz(-0.28428712) q[2];
sx q[2];
rz(-0.686336) q[2];
sx q[2];
rz(-0.057675408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18273045) q[1];
sx q[1];
rz(-2.7736543) q[1];
sx q[1];
rz(2.5065305) q[1];
x q[2];
rz(1.6790326) q[3];
sx q[3];
rz(-2.0104791) q[3];
sx q[3];
rz(0.81793438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3313716) q[2];
sx q[2];
rz(-1.081531) q[2];
sx q[2];
rz(-3.1312211) q[2];
rz(-1.2822019) q[3];
sx q[3];
rz(-2.5344262) q[3];
sx q[3];
rz(-0.31266323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4717167) q[0];
sx q[0];
rz(-2.3343021) q[0];
sx q[0];
rz(-2.5655991) q[0];
rz(0.55201375) q[1];
sx q[1];
rz(-1.3060528) q[1];
sx q[1];
rz(-0.920151) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7940038) q[0];
sx q[0];
rz(-2.2146642) q[0];
sx q[0];
rz(-2.3809303) q[0];
x q[1];
rz(1.5414562) q[2];
sx q[2];
rz(-1.27904) q[2];
sx q[2];
rz(-2.2911366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.52251053) q[1];
sx q[1];
rz(-0.42370161) q[1];
sx q[1];
rz(-0.9724618) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68814419) q[3];
sx q[3];
rz(-2.9811908) q[3];
sx q[3];
rz(0.80996603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0612023) q[2];
sx q[2];
rz(-1.3866501) q[2];
sx q[2];
rz(2.6913255) q[2];
rz(0.14084147) q[3];
sx q[3];
rz(-1.9167269) q[3];
sx q[3];
rz(2.5950477) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30215728) q[0];
sx q[0];
rz(-1.7385229) q[0];
sx q[0];
rz(-0.19254011) q[0];
rz(0.64388609) q[1];
sx q[1];
rz(-1.5637249) q[1];
sx q[1];
rz(3.0582757) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5697428) q[0];
sx q[0];
rz(-0.62456399) q[0];
sx q[0];
rz(1.3032105) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0599454) q[2];
sx q[2];
rz(-1.6892216) q[2];
sx q[2];
rz(-1.7497964) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6308438) q[1];
sx q[1];
rz(-2.7845195) q[1];
sx q[1];
rz(-2.0791173) q[1];
rz(1.313435) q[3];
sx q[3];
rz(-1.93522) q[3];
sx q[3];
rz(2.6217883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7328428) q[2];
sx q[2];
rz(-2.7538444) q[2];
sx q[2];
rz(2.5852618) q[2];
rz(2.2522816) q[3];
sx q[3];
rz(-1.5143062) q[3];
sx q[3];
rz(-2.6591532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0675875) q[0];
sx q[0];
rz(-2.6370625) q[0];
sx q[0];
rz(2.0489847) q[0];
rz(0.69674528) q[1];
sx q[1];
rz(-0.69843355) q[1];
sx q[1];
rz(2.4901857) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3676746) q[0];
sx q[0];
rz(-2.2332397) q[0];
sx q[0];
rz(1.7791616) q[0];
rz(-pi) q[1];
rz(-2.0938854) q[2];
sx q[2];
rz(-2.4019057) q[2];
sx q[2];
rz(-0.51442622) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5719381) q[1];
sx q[1];
rz(-1.2982315) q[1];
sx q[1];
rz(2.1361208) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6241996) q[3];
sx q[3];
rz(-1.5861476) q[3];
sx q[3];
rz(2.9563167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6394028) q[2];
sx q[2];
rz(-0.71523372) q[2];
sx q[2];
rz(-2.5779842) q[2];
rz(3.0230076) q[3];
sx q[3];
rz(-2.1309116) q[3];
sx q[3];
rz(-2.3732869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3350155) q[0];
sx q[0];
rz(-1.7105569) q[0];
sx q[0];
rz(3.0203876) q[0];
rz(-0.14713261) q[1];
sx q[1];
rz(-1.9857261) q[1];
sx q[1];
rz(2.3187231) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8504353) q[0];
sx q[0];
rz(-1.4802093) q[0];
sx q[0];
rz(-0.25048243) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0953832) q[2];
sx q[2];
rz(-2.5760057) q[2];
sx q[2];
rz(-2.2639745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40892021) q[1];
sx q[1];
rz(-0.80574811) q[1];
sx q[1];
rz(0.21269704) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9328281) q[3];
sx q[3];
rz(-2.9144807) q[3];
sx q[3];
rz(2.1696573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0269028) q[2];
sx q[2];
rz(-0.59640408) q[2];
sx q[2];
rz(2.5531595) q[2];
rz(2.9450997) q[3];
sx q[3];
rz(-1.4339707) q[3];
sx q[3];
rz(0.65832037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7115985) q[0];
sx q[0];
rz(-2.5007432) q[0];
sx q[0];
rz(0.95630056) q[0];
rz(-0.20005964) q[1];
sx q[1];
rz(-1.6911036) q[1];
sx q[1];
rz(-2.6106581) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5945608) q[0];
sx q[0];
rz(-2.6480125) q[0];
sx q[0];
rz(1.4303239) q[0];
rz(-2.7199581) q[2];
sx q[2];
rz(-0.68071013) q[2];
sx q[2];
rz(1.1816848) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0643474) q[1];
sx q[1];
rz(-2.5750011) q[1];
sx q[1];
rz(0.30736012) q[1];
rz(0.28778971) q[3];
sx q[3];
rz(-0.46817259) q[3];
sx q[3];
rz(-0.60540199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42813534) q[2];
sx q[2];
rz(-2.0811446) q[2];
sx q[2];
rz(-2.0088081) q[2];
rz(2.1652341) q[3];
sx q[3];
rz(-0.68621695) q[3];
sx q[3];
rz(1.8648225) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3499632) q[0];
sx q[0];
rz(-0.11194734) q[0];
sx q[0];
rz(-1.5891225) q[0];
rz(-2.7516229) q[1];
sx q[1];
rz(-1.5733122) q[1];
sx q[1];
rz(0.76036298) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50770201) q[0];
sx q[0];
rz(-1.639524) q[0];
sx q[0];
rz(0.46724369) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5993714) q[2];
sx q[2];
rz(-2.0325538) q[2];
sx q[2];
rz(-0.28959488) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8514357) q[1];
sx q[1];
rz(-0.37237114) q[1];
sx q[1];
rz(-1.6769767) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30837183) q[3];
sx q[3];
rz(-1.0532612) q[3];
sx q[3];
rz(-0.3299768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4752263) q[2];
sx q[2];
rz(-0.97576371) q[2];
sx q[2];
rz(-2.547612) q[2];
rz(2.5178759) q[3];
sx q[3];
rz(-1.2902322) q[3];
sx q[3];
rz(-2.311603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8754616) q[0];
sx q[0];
rz(-0.13127357) q[0];
sx q[0];
rz(1.776478) q[0];
rz(-3.1214556) q[1];
sx q[1];
rz(-0.29300856) q[1];
sx q[1];
rz(-1.551052) q[1];
rz(-1.5503442) q[2];
sx q[2];
rz(-2.7566789) q[2];
sx q[2];
rz(-2.821417) q[2];
rz(-1.5539215) q[3];
sx q[3];
rz(-1.3951835) q[3];
sx q[3];
rz(1.8512948) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
