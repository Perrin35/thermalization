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
rz(-1.3936477) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924661) q[0];
sx q[0];
rz(-1.2611715) q[0];
sx q[0];
rz(-1.5729088) q[0];
x q[1];
rz(-0.68140985) q[2];
sx q[2];
rz(-1.0076367) q[2];
sx q[2];
rz(1.5208706) q[2];
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
rz(-0.54772954) q[1];
rz(-2.8586219) q[3];
sx q[3];
rz(-2.1108147) q[3];
sx q[3];
rz(-0.52373826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48774886) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(-3.0207108) q[2];
rz(0.17928784) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(0.15549913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(3.0497465) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(0.13277408) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(-2.9002088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2780212) q[0];
sx q[0];
rz(-1.9262505) q[0];
sx q[0];
rz(-1.9832719) q[0];
x q[1];
rz(-2.6163289) q[2];
sx q[2];
rz(-2.004945) q[2];
sx q[2];
rz(0.94793749) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1387716) q[1];
sx q[1];
rz(-1.4797987) q[1];
sx q[1];
rz(2.2343193) q[1];
x q[2];
rz(1.7412132) q[3];
sx q[3];
rz(-0.15300551) q[3];
sx q[3];
rz(-1.9433598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1277348) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(2.0347118) q[2];
rz(1.7539304) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(1.2319516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7805507) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(2.4011491) q[0];
rz(-0.45117798) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(1.0528475) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5055288) q[0];
sx q[0];
rz(-1.2084506) q[0];
sx q[0];
rz(-2.5546968) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9387248) q[2];
sx q[2];
rz(-2.433948) q[2];
sx q[2];
rz(0.55756535) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14245089) q[1];
sx q[1];
rz(-1.4186727) q[1];
sx q[1];
rz(1.692148) q[1];
rz(-pi) q[2];
rz(-3.1316109) q[3];
sx q[3];
rz(-2.0008893) q[3];
sx q[3];
rz(-0.86356589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8911002) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(0.54692522) q[2];
rz(2.8524103) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1716487) q[0];
sx q[0];
rz(-0.26370731) q[0];
sx q[0];
rz(-1.7893715) q[0];
rz(-2.9317454) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(2.9052177) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5445697) q[0];
sx q[0];
rz(-1.3057846) q[0];
sx q[0];
rz(-0.17770627) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1430074) q[2];
sx q[2];
rz(-0.59362715) q[2];
sx q[2];
rz(-1.3015391) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.216815) q[1];
sx q[1];
rz(-1.5642011) q[1];
sx q[1];
rz(-1.738468) q[1];
rz(-1.4218016) q[3];
sx q[3];
rz(-2.0553203) q[3];
sx q[3];
rz(1.5414343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9005047) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(0.55348712) q[2];
rz(-2.2262946) q[3];
sx q[3];
rz(-1.2585636) q[3];
sx q[3];
rz(2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47167641) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(-0.70415235) q[0];
rz(1.0559121) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(-2.5456837) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0376315) q[0];
sx q[0];
rz(-2.9069293) q[0];
sx q[0];
rz(-1.9734567) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9688409) q[2];
sx q[2];
rz(-2.4606332) q[2];
sx q[2];
rz(1.3576042) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1314428) q[1];
sx q[1];
rz(-1.3859268) q[1];
sx q[1];
rz(2.7308982) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7704352) q[3];
sx q[3];
rz(-1.3988929) q[3];
sx q[3];
rz(-0.88574823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6490877) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(0.55111432) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(-1.5392039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15774396) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(-0.95170784) q[0];
rz(-2.6668008) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(-0.29528433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8897032) q[0];
sx q[0];
rz(-1.1767052) q[0];
sx q[0];
rz(-2.9910812) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6906602) q[2];
sx q[2];
rz(-2.7984598) q[2];
sx q[2];
rz(0.97230881) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59776781) q[1];
sx q[1];
rz(-2.4104767) q[1];
sx q[1];
rz(-0.75146971) q[1];
rz(-pi) q[2];
rz(-0.065159273) q[3];
sx q[3];
rz(-2.0805196) q[3];
sx q[3];
rz(-2.5173957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7513912) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-2.2078216) q[2];
rz(1.4592524) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(-1.4565844) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075832531) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(-0.24205762) q[0];
rz(-0.66479713) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(-0.40329969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2780669) q[0];
sx q[0];
rz(-0.51260548) q[0];
sx q[0];
rz(0.82786064) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34801872) q[2];
sx q[2];
rz(-2.4325779) q[2];
sx q[2];
rz(-1.4460756) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94782103) q[1];
sx q[1];
rz(-1.7809476) q[1];
sx q[1];
rz(-3.139688) q[1];
rz(-0.025267406) q[3];
sx q[3];
rz(-1.380904) q[3];
sx q[3];
rz(-0.59786284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2287801) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(0.43441233) q[2];
rz(-1.0007535) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1338761) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(-0.49466053) q[0];
rz(-1.5085295) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(2.5411434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4279815) q[0];
sx q[0];
rz(-2.5117154) q[0];
sx q[0];
rz(0.71891086) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1152994) q[2];
sx q[2];
rz(-2.0956989) q[2];
sx q[2];
rz(-0.80034791) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2272051) q[1];
sx q[1];
rz(-1.3890508) q[1];
sx q[1];
rz(1.6919796) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7793711) q[3];
sx q[3];
rz(-1.8684636) q[3];
sx q[3];
rz(1.0514333) q[3];
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
rz(2.7159193) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(-1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15329926) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(-1.6631888) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(2.7240662) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5696213) q[0];
sx q[0];
rz(-0.53335359) q[0];
sx q[0];
rz(-1.4825975) q[0];
x q[1];
rz(-1.19403) q[2];
sx q[2];
rz(-2.8586839) q[2];
sx q[2];
rz(-0.49809581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7716277) q[1];
sx q[1];
rz(-1.7376448) q[1];
sx q[1];
rz(-0.19373993) q[1];
rz(-pi) q[2];
rz(-2.9478491) q[3];
sx q[3];
rz(-1.994641) q[3];
sx q[3];
rz(-1.8684052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4801165) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(-2.647906) q[2];
rz(-0.72475973) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(0.31989583) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9676554) q[0];
sx q[0];
rz(-0.65615654) q[0];
sx q[0];
rz(-0.68558145) q[0];
rz(-2.8441692) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(-2.0102274) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.92537) q[0];
sx q[0];
rz(-2.4664481) q[0];
sx q[0];
rz(2.569909) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70558833) q[2];
sx q[2];
rz(-1.0552647) q[2];
sx q[2];
rz(0.74598344) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3494306) q[1];
sx q[1];
rz(-1.7948479) q[1];
sx q[1];
rz(2.7538124) q[1];
x q[2];
rz(2.4990436) q[3];
sx q[3];
rz(-0.71392871) q[3];
sx q[3];
rz(-0.16845265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3161105) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(2.4342009) q[2];
rz(0.80983821) q[3];
sx q[3];
rz(-2.4707068) q[3];
sx q[3];
rz(-0.18856089) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512882) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(-0.21223016) q[1];
sx q[1];
rz(-1.6925015) q[1];
sx q[1];
rz(-0.5136516) q[1];
rz(-2.1420508) q[2];
sx q[2];
rz(-1.0456325) q[2];
sx q[2];
rz(0.13620302) q[2];
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