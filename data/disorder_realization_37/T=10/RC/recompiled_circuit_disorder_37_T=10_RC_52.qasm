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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0855334) q[0];
sx q[0];
rz(-0.30963184) q[0];
sx q[0];
rz(-3.1349896) q[0];
rz(-0.78656466) q[2];
sx q[2];
rz(-2.2872891) q[2];
sx q[2];
rz(0.63209817) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8092825) q[1];
sx q[1];
rz(-2.6039632) q[1];
sx q[1];
rz(2.5938631) q[1];
rz(-1.1349036) q[3];
sx q[3];
rz(-0.60308686) q[3];
sx q[3];
rz(3.1325504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6538438) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(3.0207108) q[2];
rz(0.17928784) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(-0.15549913) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.091846175) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(0.13277408) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(2.9002088) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55610181) q[0];
sx q[0];
rz(-1.1855159) q[0];
sx q[0];
rz(-2.7566064) q[0];
rz(-pi) q[1];
rz(0.52526371) q[2];
sx q[2];
rz(-1.1366476) q[2];
sx q[2];
rz(2.1936552) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6893377) q[1];
sx q[1];
rz(-2.472795) q[1];
sx q[1];
rz(-1.4237088) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7216191) q[3];
sx q[3];
rz(-1.5449459) q[3];
sx q[3];
rz(-0.54102708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1277348) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(1.1068809) q[2];
rz(1.7539304) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36104193) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(0.74044359) q[0];
rz(2.6904147) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(1.0528475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2967865) q[0];
sx q[0];
rz(-1.0264945) q[0];
sx q[0];
rz(1.9980206) q[0];
rz(-pi) q[1];
rz(-0.89715965) q[2];
sx q[2];
rz(-1.8067915) q[2];
sx q[2];
rz(1.2981851) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82008541) q[1];
sx q[1];
rz(-0.19430375) q[1];
sx q[1];
rz(-0.66837515) q[1];
rz(-pi) q[2];
rz(0.0099817688) q[3];
sx q[3];
rz(-1.1407033) q[3];
sx q[3];
rz(0.86356589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8911002) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(0.54692522) q[2];
rz(-2.8524103) q[3];
sx q[3];
rz(-2.6047891) q[3];
sx q[3];
rz(-3.1291936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.96994394) q[0];
sx q[0];
rz(-0.26370731) q[0];
sx q[0];
rz(-1.3522211) q[0];
rz(0.2098473) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(2.9052177) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94350067) q[0];
sx q[0];
rz(-0.31790942) q[0];
sx q[0];
rz(0.99347465) q[0];
rz(2.0868446) q[2];
sx q[2];
rz(-1.8785254) q[2];
sx q[2];
rz(-2.9204521) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.826556) q[1];
sx q[1];
rz(-2.9737925) q[1];
sx q[1];
rz(1.6102953) q[1];
x q[2];
rz(-2.6524622) q[3];
sx q[3];
rz(-1.7025347) q[3];
sx q[3];
rz(3.1011503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2410879) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(0.55348712) q[2];
rz(-0.91529804) q[3];
sx q[3];
rz(-1.2585636) q[3];
sx q[3];
rz(0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6699162) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(2.4374403) q[0];
rz(1.0559121) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(-0.59590894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51673698) q[0];
sx q[0];
rz(-1.3552109) q[0];
sx q[0];
rz(-3.0481824) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4679568) q[2];
sx q[2];
rz(-1.4623702) q[2];
sx q[2];
rz(-3.0631531) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0101498) q[1];
sx q[1];
rz(-1.3859268) q[1];
sx q[1];
rz(-0.4106945) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3866053) q[3];
sx q[3];
rz(-1.2053688) q[3];
sx q[3];
rz(-2.3900677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6490877) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(-0.55111432) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(1.6023887) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838487) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(0.95170784) q[0];
rz(-2.6668008) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(-0.29528433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5135358) q[0];
sx q[0];
rz(-0.4204458) q[0];
sx q[0];
rz(1.9168617) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45093243) q[2];
sx q[2];
rz(-2.7984598) q[2];
sx q[2];
rz(-2.1692838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5608279) q[1];
sx q[1];
rz(-1.0974713) q[1];
sx q[1];
rz(0.58014371) q[1];
rz(-pi) q[2];
rz(-1.6867562) q[3];
sx q[3];
rz(-2.6280858) q[3];
sx q[3];
rz(-0.49125571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7513912) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-0.93377101) q[2];
rz(1.6823403) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.075832531) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(-0.24205762) q[0];
rz(2.4767955) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(2.738293) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032265183) q[0];
sx q[0];
rz(-1.9089713) q[0];
sx q[0];
rz(1.9637252) q[0];
x q[1];
rz(-2.4629668) q[2];
sx q[2];
rz(-1.3468862) q[2];
sx q[2];
rz(-0.14397552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.93869081) q[1];
sx q[1];
rz(-2.9314329) q[1];
sx q[1];
rz(-1.5797257) q[1];
rz(1.7014916) q[3];
sx q[3];
rz(-2.9500467) q[3];
sx q[3];
rz(0.46476118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91281259) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(0.43441233) q[2];
rz(1.0007535) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(-2.6312857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1338761) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(-2.6469321) q[0];
rz(-1.6330632) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(-0.60044926) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7136112) q[0];
sx q[0];
rz(-0.62987721) q[0];
sx q[0];
rz(0.71891086) q[0];
rz(-0.026293228) q[2];
sx q[2];
rz(-2.0956989) q[2];
sx q[2];
rz(-0.80034791) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9143876) q[1];
sx q[1];
rz(-1.3890508) q[1];
sx q[1];
rz(-1.6919796) q[1];
rz(-pi) q[2];
rz(-0.59372254) q[3];
sx q[3];
rz(-2.7799118) q[3];
sx q[3];
rz(-0.42632521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0042469) q[2];
sx q[2];
rz(-2.1477951) q[2];
sx q[2];
rz(3.0252769) q[2];
rz(0.4256734) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(-11*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15329926) q[0];
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
rz(-3.0644181) q[0];
sx q[0];
rz(-1.525997) q[0];
sx q[0];
rz(-2.1024465) q[0];
rz(-pi) q[1];
rz(-3.0350424) q[2];
sx q[2];
rz(-1.3082192) q[2];
sx q[2];
rz(3.0343461) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.49839834) q[1];
sx q[1];
rz(-2.8865951) q[1];
sx q[1];
rz(0.71868371) q[1];
rz(-pi) q[2];
rz(-1.9741251) q[3];
sx q[3];
rz(-0.46357337) q[3];
sx q[3];
rz(1.4232672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4801165) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(0.49368668) q[2];
rz(-2.4168329) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(0.31989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(0.17393728) q[0];
sx q[0];
rz(-0.65615654) q[0];
sx q[0];
rz(2.4560112) q[0];
rz(-0.29742345) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(2.0102274) q[1];
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
x q[1];
rz(2.4234424) q[2];
sx q[2];
rz(-0.8469204) q[2];
sx q[2];
rz(-1.3494327) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.28045052) q[1];
sx q[1];
rz(-0.44499731) q[1];
sx q[1];
rz(2.599237) q[1];
rz(-pi) q[2];
rz(-2.0496619) q[3];
sx q[3];
rz(-1.0189971) q[3];
sx q[3];
rz(2.5294876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82548213) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(-0.70739174) q[2];
rz(-0.80983821) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512882) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(0.21223016) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(-0.99954188) q[2];
sx q[2];
rz(-2.0959601) q[2];
sx q[2];
rz(-3.0053896) q[2];
rz(2.5645732) q[3];
sx q[3];
rz(-2.1885625) q[3];
sx q[3];
rz(2.267005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
