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
rz(-1.3936477) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05605927) q[0];
sx q[0];
rz(-2.8319608) q[0];
sx q[0];
rz(-0.0066030533) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2533478) q[2];
sx q[2];
rz(-2.1324854) q[2];
sx q[2];
rz(-2.782926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71478903) q[1];
sx q[1];
rz(-2.0232632) q[1];
sx q[1];
rz(-1.2697551) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2829708) q[3];
sx q[3];
rz(-1.030778) q[3];
sx q[3];
rz(2.6178544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.48774886) q[2];
sx q[2];
rz(-1.2922492) q[2];
sx q[2];
rz(0.1208819) q[2];
rz(0.17928784) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(0.15549913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091846175) q[0];
sx q[0];
rz(-0.76773983) q[0];
sx q[0];
rz(0.13277408) q[0];
rz(-1.6800539) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(-0.24138385) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7624843) q[0];
sx q[0];
rz(-2.603841) q[0];
sx q[0];
rz(2.3178029) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0788467) q[2];
sx q[2];
rz(-2.0430806) q[2];
sx q[2];
rz(2.7578596) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5026958) q[1];
sx q[1];
rz(-2.2310889) q[1];
sx q[1];
rz(3.026282) q[1];
rz(-pi) q[2];
rz(-1.7412132) q[3];
sx q[3];
rz(-2.9885871) q[3];
sx q[3];
rz(1.1982329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1277348) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(1.1068809) q[2];
rz(-1.3876623) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.36104193) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(-2.4011491) q[0];
rz(-2.6904147) q[1];
sx q[1];
rz(-1.2606882) q[1];
sx q[1];
rz(1.0528475) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2967865) q[0];
sx q[0];
rz(-1.0264945) q[0];
sx q[0];
rz(1.9980206) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2028678) q[2];
sx q[2];
rz(-2.433948) q[2];
sx q[2];
rz(-0.55756535) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.82008541) q[1];
sx q[1];
rz(-0.19430375) q[1];
sx q[1];
rz(-2.4732175) q[1];
rz(1.5925519) q[3];
sx q[3];
rz(-2.7113911) q[3];
sx q[3];
rz(0.83963001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25049245) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(-0.54692522) q[2];
rz(0.28918239) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(3.1291936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1716487) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.3522211) q[0];
rz(2.9317454) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(2.9052177) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.597023) q[0];
sx q[0];
rz(-1.835808) q[0];
sx q[0];
rz(-0.17770627) q[0];
x q[1];
rz(-0.9985853) q[2];
sx q[2];
rz(-0.59362715) q[2];
sx q[2];
rz(-1.3015391) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.826556) q[1];
sx q[1];
rz(-2.9737925) q[1];
sx q[1];
rz(-1.5312974) q[1];
rz(-pi) q[2];
rz(-1.4218016) q[3];
sx q[3];
rz(-2.0553203) q[3];
sx q[3];
rz(-1.6001584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9005047) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(0.55348712) q[2];
rz(0.91529804) q[3];
sx q[3];
rz(-1.2585636) q[3];
sx q[3];
rz(-0.64490157) q[3];
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
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47167641) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(0.70415235) q[0];
rz(-2.0856805) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(-0.59590894) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51673698) q[0];
sx q[0];
rz(-1.3552109) q[0];
sx q[0];
rz(3.0481824) q[0];
rz(-pi) q[1];
rz(2.4679568) q[2];
sx q[2];
rz(-1.6792225) q[2];
sx q[2];
rz(0.078439586) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1314428) q[1];
sx q[1];
rz(-1.7556659) q[1];
sx q[1];
rz(-0.4106945) q[1];
rz(-pi) q[2];
rz(0.44645198) q[3];
sx q[3];
rz(-2.7342396) q[3];
sx q[3];
rz(0.27094597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6490877) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(-2.5904783) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15774396) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(2.1898848) q[0];
rz(0.47479182) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(-2.8463083) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2518894) q[0];
sx q[0];
rz(-1.9648874) q[0];
sx q[0];
rz(0.15051145) q[0];
x q[1];
rz(1.7252543) q[2];
sx q[2];
rz(-1.8784349) q[2];
sx q[2];
rz(0.49738202) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5438248) q[1];
sx q[1];
rz(-0.73111594) q[1];
sx q[1];
rz(-0.75146971) q[1];
rz(-pi) q[2];
rz(0.065159273) q[3];
sx q[3];
rz(-2.0805196) q[3];
sx q[3];
rz(2.5173957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39020145) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(2.2078216) q[2];
rz(1.6823403) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(-1.4565844) q[3];
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
rz(-0.075832531) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(-0.24205762) q[0];
rz(0.66479713) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(0.40329969) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1093275) q[0];
sx q[0];
rz(-1.9089713) q[0];
sx q[0];
rz(1.1778675) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4629668) q[2];
sx q[2];
rz(-1.7947065) q[2];
sx q[2];
rz(-2.9976171) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62337263) q[1];
sx q[1];
rz(-1.5726591) q[1];
sx q[1];
rz(1.3606447) q[1];
rz(-pi) q[2];
rz(1.4401011) q[3];
sx q[3];
rz(-0.19154597) q[3];
sx q[3];
rz(-2.6768315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91281259) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(2.7071803) q[2];
rz(1.0007535) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(-0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1338761) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(-2.6469321) q[0];
rz(1.5085295) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(0.60044926) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8884044) q[0];
sx q[0];
rz(-1.1115523) q[0];
sx q[0];
rz(2.0183536) q[0];
rz(1.5254283) q[2];
sx q[2];
rz(-2.6160935) q[2];
sx q[2];
rz(0.85277992) q[2];
x q[3];
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
rz(0.59372254) q[3];
sx q[3];
rz(-2.7799118) q[3];
sx q[3];
rz(-2.7152674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-2.1477951) q[2];
sx q[2];
rz(-3.0252769) q[2];
rz(-0.4256734) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(-11*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9882934) q[0];
sx q[0];
rz(-0.17803742) q[0];
sx q[0];
rz(-1.4784038) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(-0.41752648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4672887) q[0];
sx q[0];
rz(-2.1018565) q[0];
sx q[0];
rz(3.0896316) q[0];
x q[1];
rz(1.9475627) q[2];
sx q[2];
rz(-0.28290877) q[2];
sx q[2];
rz(0.49809581) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.233404) q[1];
sx q[1];
rz(-1.7618124) q[1];
sx q[1];
rz(1.4008271) q[1];
x q[2];
rz(-2.9478491) q[3];
sx q[3];
rz(-1.1469517) q[3];
sx q[3];
rz(1.8684052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6614762) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9676554) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(-2.4560112) q[0];
rz(0.29742345) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(1.1313653) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.92537) q[0];
sx q[0];
rz(-2.4664481) q[0];
sx q[0];
rz(0.57168369) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4360043) q[2];
sx q[2];
rz(-1.0552647) q[2];
sx q[2];
rz(2.3956092) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3494306) q[1];
sx q[1];
rz(-1.3467448) q[1];
sx q[1];
rz(-0.38778023) q[1];
rz(0.60641842) q[3];
sx q[3];
rz(-1.1675721) q[3];
sx q[3];
rz(1.2244146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3161105) q[2];
sx q[2];
rz(-1.9448514) q[2];
sx q[2];
rz(2.4342009) q[2];
rz(2.3317544) q[3];
sx q[3];
rz(-2.4707068) q[3];
sx q[3];
rz(0.18856089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903044) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(2.9293625) q[1];
sx q[1];
rz(-1.6925015) q[1];
sx q[1];
rz(-0.5136516) q[1];
rz(-2.390776) q[2];
sx q[2];
rz(-0.75559794) q[2];
sx q[2];
rz(-2.0970367) q[2];
rz(-0.86757579) q[3];
sx q[3];
rz(-2.0316364) q[3];
sx q[3];
rz(0.33566725) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
