OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(1.2020943) q[0];
sx q[0];
rz(10.572875) q[0];
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(1.747945) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924661) q[0];
sx q[0];
rz(-1.8804212) q[0];
sx q[0];
rz(-1.5686839) q[0];
rz(0.78656466) q[2];
sx q[2];
rz(-0.85430356) q[2];
sx q[2];
rz(0.63209817) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8092825) q[1];
sx q[1];
rz(-2.6039632) q[1];
sx q[1];
rz(2.5938631) q[1];
rz(-2.8586219) q[3];
sx q[3];
rz(-2.1108147) q[3];
sx q[3];
rz(2.6178544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.48774886) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(-0.1208819) q[2];
rz(0.17928784) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(0.15549913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091846175) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(-0.13277408) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(-0.24138385) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5854908) q[0];
sx q[0];
rz(-1.9560768) q[0];
sx q[0];
rz(0.38498621) q[0];
x q[1];
rz(2.0627459) q[2];
sx q[2];
rz(-2.0430806) q[2];
sx q[2];
rz(-2.7578596) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.0028210359) q[1];
sx q[1];
rz(-1.661794) q[1];
sx q[1];
rz(0.9072733) q[1];
rz(1.7412132) q[3];
sx q[3];
rz(-2.9885871) q[3];
sx q[3];
rz(-1.1982329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36104193) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(-2.4011491) q[0];
rz(2.6904147) q[1];
sx q[1];
rz(-1.2606882) q[1];
sx q[1];
rz(-1.0528475) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42441472) q[0];
sx q[0];
rz(-2.4632235) q[0];
sx q[0];
rz(2.5413187) q[0];
rz(2.8430976) q[2];
sx q[2];
rz(-2.2224991) q[2];
sx q[2];
rz(-0.088108206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14245089) q[1];
sx q[1];
rz(-1.4186727) q[1];
sx q[1];
rz(1.4494447) q[1];
rz(1.5490407) q[3];
sx q[3];
rz(-2.7113911) q[3];
sx q[3];
rz(2.3019626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.25049245) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(2.5946674) q[2];
rz(2.8524103) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(-3.1291936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(1.7893715) q[0];
rz(-2.9317454) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(2.9052177) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5445697) q[0];
sx q[0];
rz(-1.835808) q[0];
sx q[0];
rz(-0.17770627) q[0];
rz(-pi) q[1];
rz(2.1430074) q[2];
sx q[2];
rz(-0.59362715) q[2];
sx q[2];
rz(1.8400536) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35509767) q[1];
sx q[1];
rz(-1.7384643) q[1];
sx q[1];
rz(-3.1349036) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.719791) q[3];
sx q[3];
rz(-2.0553203) q[3];
sx q[3];
rz(-1.5414343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9005047) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(0.55348712) q[2];
rz(-0.91529804) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6699162) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(0.70415235) q[0];
rz(-1.0559121) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(2.5456837) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0674954) q[0];
sx q[0];
rz(-1.4795545) q[0];
sx q[0];
rz(1.7872966) q[0];
x q[1];
rz(1.4324097) q[2];
sx q[2];
rz(-0.90183898) q[2];
sx q[2];
rz(-1.578518) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.039918598) q[1];
sx q[1];
rz(-2.6933751) q[1];
sx q[1];
rz(0.43804534) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3866053) q[3];
sx q[3];
rz(-1.2053688) q[3];
sx q[3];
rz(2.3900677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.492505) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(0.55111432) q[2];
rz(0.20714949) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15774396) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(-2.1898848) q[0];
rz(-2.6668008) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(-0.29528433) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8808522) q[0];
sx q[0];
rz(-1.4319001) q[0];
sx q[0];
rz(1.9689346) q[0];
x q[1];
rz(-1.4163383) q[2];
sx q[2];
rz(-1.2631577) q[2];
sx q[2];
rz(2.6442106) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5438248) q[1];
sx q[1];
rz(-0.73111594) q[1];
sx q[1];
rz(2.3901229) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6867562) q[3];
sx q[3];
rz(-0.51350683) q[3];
sx q[3];
rz(-2.6503369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.39020145) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-2.2078216) q[2];
rz(1.4592524) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0657601) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(-2.899535) q[0];
rz(0.66479713) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(-2.738293) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1093275) q[0];
sx q[0];
rz(-1.9089713) q[0];
sx q[0];
rz(-1.1778675) q[0];
rz(2.7935739) q[2];
sx q[2];
rz(-0.70901477) q[2];
sx q[2];
rz(1.4460756) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.94782103) q[1];
sx q[1];
rz(-1.7809476) q[1];
sx q[1];
rz(-0.0019046849) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7607479) q[3];
sx q[3];
rz(-1.5459832) q[3];
sx q[3];
rz(2.1638889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.91281259) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(2.7071803) q[2];
rz(-2.1408391) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(2.6312857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0077165724) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(0.49466053) q[0];
rz(-1.6330632) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(-2.5411434) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531882) q[0];
sx q[0];
rz(-2.0300403) q[0];
sx q[0];
rz(-2.0183536) q[0];
rz(-pi) q[1];
rz(-1.6161643) q[2];
sx q[2];
rz(-2.6160935) q[2];
sx q[2];
rz(0.85277992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2272051) q[1];
sx q[1];
rz(-1.7525418) q[1];
sx q[1];
rz(-1.449613) q[1];
rz(-pi) q[2];
rz(2.8377257) q[3];
sx q[3];
rz(-1.3715203) q[3];
sx q[3];
rz(2.5602333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1373458) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(-3.0252769) q[2];
rz(-0.4256734) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.9882934) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(-1.4784038) q[0];
rz(2.2019745) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(2.7240662) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.674304) q[0];
sx q[0];
rz(-2.1018565) q[0];
sx q[0];
rz(-0.051961016) q[0];
x q[1];
rz(-0.1065503) q[2];
sx q[2];
rz(-1.8333734) q[2];
sx q[2];
rz(3.0343461) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6431943) q[1];
sx q[1];
rz(-0.25499757) q[1];
sx q[1];
rz(-2.4229089) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9741251) q[3];
sx q[3];
rz(-2.6780193) q[3];
sx q[3];
rz(-1.7183255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6614762) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(-0.49368668) q[2];
rz(0.72475973) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9676554) q[0];
sx q[0];
rz(-0.65615654) q[0];
sx q[0];
rz(0.68558145) q[0];
rz(0.29742345) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(-2.0102274) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5269276) q[0];
sx q[0];
rz(-2.1242495) q[0];
sx q[0];
rz(1.1620031) q[0];
rz(-pi) q[1];
rz(-2.4234424) q[2];
sx q[2];
rz(-0.8469204) q[2];
sx q[2];
rz(-1.7921599) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2724534) q[1];
sx q[1];
rz(-1.1932045) q[1];
sx q[1];
rz(-1.8121522) q[1];
rz(1.0919308) q[3];
sx q[3];
rz(-2.1225956) q[3];
sx q[3];
rz(0.61210504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3161105) q[2];
sx q[2];
rz(-1.9448514) q[2];
sx q[2];
rz(-2.4342009) q[2];
rz(-2.3317544) q[3];
sx q[3];
rz(-2.4707068) q[3];
sx q[3];
rz(2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903044) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(0.21223016) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(2.390776) q[2];
sx q[2];
rz(-2.3859947) q[2];
sx q[2];
rz(1.044556) q[2];
rz(0.91602305) q[3];
sx q[3];
rz(-0.81867221) q[3];
sx q[3];
rz(1.4233854) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];