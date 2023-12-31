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
rz(-1.9394983) q[0];
sx q[0];
rz(1.9934959) q[0];
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(1.747945) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5210261) q[0];
sx q[0];
rz(-1.5728083) q[0];
sx q[0];
rz(-2.8319671) q[0];
rz(-pi) q[1];
rz(0.68140985) q[2];
sx q[2];
rz(-1.0076367) q[2];
sx q[2];
rz(1.6207221) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.72109556) q[1];
sx q[1];
rz(-1.3008529) q[1];
sx q[1];
rz(-2.6707778) q[1];
x q[2];
rz(2.8586219) q[3];
sx q[3];
rz(-1.030778) q[3];
sx q[3];
rz(2.6178544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48774886) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(-3.0207108) q[2];
rz(0.17928784) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(-2.9860935) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.091846175) q[0];
sx q[0];
rz(-0.76773983) q[0];
sx q[0];
rz(-0.13277408) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(-0.24138385) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55610181) q[0];
sx q[0];
rz(-1.9560768) q[0];
sx q[0];
rz(-2.7566064) q[0];
rz(-1.0788467) q[2];
sx q[2];
rz(-1.098512) q[2];
sx q[2];
rz(2.7578596) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1387716) q[1];
sx q[1];
rz(-1.661794) q[1];
sx q[1];
rz(0.9072733) q[1];
rz(-pi) q[2];
x q[2];
rz(0.026147141) q[3];
sx q[3];
rz(-1.7215683) q[3];
sx q[3];
rz(-2.1157516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1277348) q[2];
sx q[2];
rz(-1.8227791) q[2];
sx q[2];
rz(-2.0347118) q[2];
rz(-1.3876623) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36104193) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(2.4011491) q[0];
rz(0.45117798) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(2.0887451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171779) q[0];
sx q[0];
rz(-2.4632235) q[0];
sx q[0];
rz(-2.5413187) q[0];
rz(1.2028678) q[2];
sx q[2];
rz(-2.433948) q[2];
sx q[2];
rz(0.55756535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3215072) q[1];
sx q[1];
rz(-2.9472889) q[1];
sx q[1];
rz(-2.4732175) q[1];
rz(-1.5490407) q[3];
sx q[3];
rz(-2.7113911) q[3];
sx q[3];
rz(0.83963001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8911002) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(0.54692522) q[2];
rz(0.28918239) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(-0.012399013) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.3522211) q[0];
rz(0.2098473) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(-0.23637493) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5445697) q[0];
sx q[0];
rz(-1.3057846) q[0];
sx q[0];
rz(0.17770627) q[0];
x q[1];
rz(-2.0868446) q[2];
sx q[2];
rz(-1.8785254) q[2];
sx q[2];
rz(2.9204521) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35509767) q[1];
sx q[1];
rz(-1.7384643) q[1];
sx q[1];
rz(3.1349036) q[1];
x q[2];
rz(1.4218016) q[3];
sx q[3];
rz(-2.0553203) q[3];
sx q[3];
rz(1.6001584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2410879) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(2.5881055) q[2];
rz(2.2262946) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.47167641) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(-0.70415235) q[0];
rz(-1.0559121) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(2.5456837) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0674954) q[0];
sx q[0];
rz(-1.4795545) q[0];
sx q[0];
rz(1.3542961) q[0];
x q[1];
rz(2.9688409) q[2];
sx q[2];
rz(-2.4606332) q[2];
sx q[2];
rz(1.7839884) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1016741) q[1];
sx q[1];
rz(-2.6933751) q[1];
sx q[1];
rz(0.43804534) q[1];
rz(-pi) q[2];
rz(-1.7549873) q[3];
sx q[3];
rz(-1.2053688) q[3];
sx q[3];
rz(2.3900677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.492505) q[2];
sx q[2];
rz(-2.0272144) q[2];
sx q[2];
rz(0.55111432) q[2];
rz(0.20714949) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(-1.5392039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15774396) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(-0.95170784) q[0];
rz(0.47479182) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(-2.8463083) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8897032) q[0];
sx q[0];
rz(-1.1767052) q[0];
sx q[0];
rz(-0.15051145) q[0];
rz(-0.45093243) q[2];
sx q[2];
rz(-0.34313289) q[2];
sx q[2];
rz(2.1692838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5807647) q[1];
sx q[1];
rz(-2.0441214) q[1];
sx q[1];
rz(-2.5614489) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6867562) q[3];
sx q[3];
rz(-0.51350683) q[3];
sx q[3];
rz(0.49125571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7513912) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-2.2078216) q[2];
rz(-1.6823403) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(-1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075832531) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(0.24205762) q[0];
rz(-0.66479713) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(0.40329969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86352578) q[0];
sx q[0];
rz(-2.6289872) q[0];
sx q[0];
rz(2.313732) q[0];
rz(-pi) q[1];
rz(0.67862582) q[2];
sx q[2];
rz(-1.3468862) q[2];
sx q[2];
rz(2.9976171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2029018) q[1];
sx q[1];
rz(-2.9314329) q[1];
sx q[1];
rz(-1.5797257) q[1];
rz(-1.7014916) q[3];
sx q[3];
rz(-2.9500467) q[3];
sx q[3];
rz(2.6768315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2287801) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(2.7071803) q[2];
rz(-1.0007535) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1338761) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(2.6469321) q[0];
rz(-1.6330632) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(-0.60044926) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6143215) q[0];
sx q[0];
rz(-1.9691756) q[0];
sx q[0];
rz(0.50171731) q[0];
x q[1];
rz(-1.5254283) q[2];
sx q[2];
rz(-2.6160935) q[2];
sx q[2];
rz(2.2888127) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2272051) q[1];
sx q[1];
rz(-1.7525418) q[1];
sx q[1];
rz(1.449613) q[1];
rz(1.3622215) q[3];
sx q[3];
rz(-1.8684636) q[3];
sx q[3];
rz(1.0514333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-2.1477951) q[2];
sx q[2];
rz(0.11631575) q[2];
rz(2.7159193) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(-1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15329926) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(-1.6631888) q[0];
rz(2.2019745) q[1];
sx q[1];
rz(-1.3213108) q[1];
sx q[1];
rz(0.41752648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077174515) q[0];
sx q[0];
rz(-1.6155956) q[0];
sx q[0];
rz(-1.0391462) q[0];
rz(-pi) q[1];
rz(3.0350424) q[2];
sx q[2];
rz(-1.8333734) q[2];
sx q[2];
rz(-0.1072466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6431943) q[1];
sx q[1];
rz(-0.25499757) q[1];
sx q[1];
rz(0.71868371) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19374356) q[3];
sx q[3];
rz(-1.994641) q[3];
sx q[3];
rz(1.2731874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6614762) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(0.49368668) q[2];
rz(-0.72475973) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(0.31989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9676554) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(2.4560112) q[0];
rz(-2.8441692) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(-2.0102274) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5269276) q[0];
sx q[0];
rz(-2.1242495) q[0];
sx q[0];
rz(-1.1620031) q[0];
rz(0.93090242) q[2];
sx q[2];
rz(-2.1702592) q[2];
sx q[2];
rz(2.714389) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2724534) q[1];
sx q[1];
rz(-1.9483882) q[1];
sx q[1];
rz(-1.8121522) q[1];
x q[2];
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
rz(-0.82548213) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(-2.4342009) q[2];
rz(-0.80983821) q[3];
sx q[3];
rz(-2.4707068) q[3];
sx q[3];
rz(0.18856089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512882) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(0.21223016) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(2.1420508) q[2];
sx q[2];
rz(-2.0959601) q[2];
sx q[2];
rz(-3.0053896) q[2];
rz(0.57701941) q[3];
sx q[3];
rz(-0.95303017) q[3];
sx q[3];
rz(-0.87458761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
