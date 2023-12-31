OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5503791) q[0];
sx q[0];
rz(3.1449218) q[0];
sx q[0];
rz(9.7700906) q[0];
rz(-1.3357063) q[1];
sx q[1];
rz(-0.3392646) q[1];
sx q[1];
rz(0.27944922) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1169259) q[0];
sx q[0];
rz(-1.3660087) q[0];
sx q[0];
rz(-1.8860399) q[0];
rz(0.34502132) q[2];
sx q[2];
rz(-0.81232386) q[2];
sx q[2];
rz(1.617384) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2742259) q[1];
sx q[1];
rz(-0.19051954) q[1];
sx q[1];
rz(1.6412559) q[1];
rz(-pi) q[2];
rz(2.4996098) q[3];
sx q[3];
rz(-1.6025474) q[3];
sx q[3];
rz(-0.25105219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7001069) q[2];
sx q[2];
rz(-1.4725279) q[2];
sx q[2];
rz(-2.1842365) q[2];
rz(-0.29933128) q[3];
sx q[3];
rz(-0.39756164) q[3];
sx q[3];
rz(2.7295952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1608202) q[0];
sx q[0];
rz(-2.1919724) q[0];
sx q[0];
rz(-0.41369307) q[0];
rz(1.7970239) q[1];
sx q[1];
rz(-0.78318703) q[1];
sx q[1];
rz(-2.5059674) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18938633) q[0];
sx q[0];
rz(-1.4385907) q[0];
sx q[0];
rz(0.044747523) q[0];
x q[1];
rz(-1.9102155) q[2];
sx q[2];
rz(-3.0017827) q[2];
sx q[2];
rz(-1.2466873) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0483413) q[1];
sx q[1];
rz(-1.3158568) q[1];
sx q[1];
rz(-1.3202207) q[1];
rz(-0.70062462) q[3];
sx q[3];
rz(-1.0430338) q[3];
sx q[3];
rz(0.81567314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3893434) q[2];
sx q[2];
rz(-1.8698591) q[2];
sx q[2];
rz(0.52865571) q[2];
rz(1.901249) q[3];
sx q[3];
rz(-0.35959187) q[3];
sx q[3];
rz(-2.8958029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5813331) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(2.8833959) q[0];
rz(1.6351581) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(-2.7071276) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4397944) q[0];
sx q[0];
rz(-2.5520303) q[0];
sx q[0];
rz(-2.2586285) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0412695) q[2];
sx q[2];
rz(-1.5142421) q[2];
sx q[2];
rz(2.2296485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79954051) q[1];
sx q[1];
rz(-0.5040579) q[1];
sx q[1];
rz(2.6964158) q[1];
x q[2];
rz(-0.86359777) q[3];
sx q[3];
rz(-2.2798385) q[3];
sx q[3];
rz(1.1532702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7210641) q[2];
sx q[2];
rz(-1.8177744) q[2];
sx q[2];
rz(-0.055796441) q[2];
rz(0.93786401) q[3];
sx q[3];
rz(-0.28356975) q[3];
sx q[3];
rz(-0.24308932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043512251) q[0];
sx q[0];
rz(-0.94399095) q[0];
sx q[0];
rz(-3.0010624) q[0];
rz(0.17164104) q[1];
sx q[1];
rz(-1.834603) q[1];
sx q[1];
rz(-0.26352873) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7292273) q[0];
sx q[0];
rz(-1.0050887) q[0];
sx q[0];
rz(1.6914781) q[0];
rz(-pi) q[1];
rz(0.66480555) q[2];
sx q[2];
rz(-1.1603328) q[2];
sx q[2];
rz(-0.38924205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4315223) q[1];
sx q[1];
rz(-0.78774161) q[1];
sx q[1];
rz(-0.4962994) q[1];
rz(-pi) q[2];
rz(-2.2534091) q[3];
sx q[3];
rz(-0.41999751) q[3];
sx q[3];
rz(0.94204599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(1.8919224) q[2];
rz(-1.9469056) q[3];
sx q[3];
rz(-1.0281111) q[3];
sx q[3];
rz(-2.4627114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2018305) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(-0.029065954) q[0];
rz(-1.4020231) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(0.32863858) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0120481) q[0];
sx q[0];
rz(-0.49783266) q[0];
sx q[0];
rz(-1.7853907) q[0];
x q[1];
rz(-1.3939912) q[2];
sx q[2];
rz(-1.8362152) q[2];
sx q[2];
rz(3.0826498) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29341104) q[1];
sx q[1];
rz(-1.1383346) q[1];
sx q[1];
rz(-1.5007988) q[1];
x q[2];
rz(-2.394649) q[3];
sx q[3];
rz(-1.1960104) q[3];
sx q[3];
rz(0.67583109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11792004) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(-2.9809791) q[2];
rz(-1.1462071) q[3];
sx q[3];
rz(-1.3137772) q[3];
sx q[3];
rz(2.5207991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6495431) q[0];
sx q[0];
rz(-2.2821125) q[0];
sx q[0];
rz(-2.3550526) q[0];
rz(-0.37711626) q[1];
sx q[1];
rz(-2.2943594) q[1];
sx q[1];
rz(-0.4424817) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9839358) q[0];
sx q[0];
rz(-1.929495) q[0];
sx q[0];
rz(1.1789807) q[0];
rz(-0.57420337) q[2];
sx q[2];
rz(-1.310537) q[2];
sx q[2];
rz(1.1299709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9128531) q[1];
sx q[1];
rz(-0.40816669) q[1];
sx q[1];
rz(2.2250882) q[1];
rz(-1.8507067) q[3];
sx q[3];
rz(-0.32752447) q[3];
sx q[3];
rz(2.1675828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4601712) q[2];
sx q[2];
rz(-2.5632016) q[2];
sx q[2];
rz(2.1703413) q[2];
rz(2.4462637) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(0.42738459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5044395) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(-1.0094281) q[0];
rz(-2.5908453) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(-1.1632464) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79272072) q[0];
sx q[0];
rz(-2.7808245) q[0];
sx q[0];
rz(-1.9804721) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86261729) q[2];
sx q[2];
rz(-1.7628981) q[2];
sx q[2];
rz(2.5196599) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0102331) q[1];
sx q[1];
rz(-1.7903923) q[1];
sx q[1];
rz(-1.4036199) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8323152) q[3];
sx q[3];
rz(-0.89865548) q[3];
sx q[3];
rz(2.6291763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.358868) q[2];
sx q[2];
rz(-1.9903368) q[2];
sx q[2];
rz(-0.7981832) q[2];
rz(-0.77945566) q[3];
sx q[3];
rz(-2.6051086) q[3];
sx q[3];
rz(-1.1727758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1786132) q[0];
sx q[0];
rz(-0.48291746) q[0];
sx q[0];
rz(3.0187507) q[0];
rz(0.12610647) q[1];
sx q[1];
rz(-1.5051196) q[1];
sx q[1];
rz(1.925148) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0053596) q[0];
sx q[0];
rz(-1.981712) q[0];
sx q[0];
rz(0.67816011) q[0];
rz(1.4037651) q[2];
sx q[2];
rz(-2.0207496) q[2];
sx q[2];
rz(0.049875967) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0876906) q[1];
sx q[1];
rz(-1.1864099) q[1];
sx q[1];
rz(-1.3575421) q[1];
rz(-pi) q[2];
rz(0.27463953) q[3];
sx q[3];
rz(-2.5317149) q[3];
sx q[3];
rz(-2.9230737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8806261) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(-0.043126062) q[2];
rz(0.17523266) q[3];
sx q[3];
rz(-0.8774811) q[3];
sx q[3];
rz(1.5847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50424987) q[0];
sx q[0];
rz(-0.054333996) q[0];
sx q[0];
rz(-0.20877008) q[0];
rz(1.6437221) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(-2.0170905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9546684) q[0];
sx q[0];
rz(-2.9483729) q[0];
sx q[0];
rz(-2.3528071) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8378216) q[2];
sx q[2];
rz(-0.73799947) q[2];
sx q[2];
rz(-2.8673025) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3885755) q[1];
sx q[1];
rz(-1.5696226) q[1];
sx q[1];
rz(-1.4374103) q[1];
rz(2.6526768) q[3];
sx q[3];
rz(-0.93547869) q[3];
sx q[3];
rz(-2.2152701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1407397) q[2];
sx q[2];
rz(-0.72379392) q[2];
sx q[2];
rz(-0.17803426) q[2];
rz(1.5464276) q[3];
sx q[3];
rz(-1.3083357) q[3];
sx q[3];
rz(-2.6509638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7889325) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(1.0349405) q[0];
rz(-0.79822284) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(-0.14990526) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7222544) q[0];
sx q[0];
rz(-2.1968578) q[0];
sx q[0];
rz(-1.1115848) q[0];
rz(-pi) q[1];
rz(-2.855004) q[2];
sx q[2];
rz(-1.5550685) q[2];
sx q[2];
rz(0.77428267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1356218) q[1];
sx q[1];
rz(-1.0711728) q[1];
sx q[1];
rz(0.47132229) q[1];
x q[2];
rz(2.4622915) q[3];
sx q[3];
rz(-1.7656109) q[3];
sx q[3];
rz(2.4904345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4118816) q[2];
sx q[2];
rz(-1.9103266) q[2];
sx q[2];
rz(-2.7330772) q[2];
rz(0.92489964) q[3];
sx q[3];
rz(-0.81245208) q[3];
sx q[3];
rz(-0.48172054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9020486) q[0];
sx q[0];
rz(-1.5477381) q[0];
sx q[0];
rz(1.5292194) q[0];
rz(-1.3760024) q[1];
sx q[1];
rz(-1.1767495) q[1];
sx q[1];
rz(-1.8935988) q[1];
rz(0.79509905) q[2];
sx q[2];
rz(-1.3394525) q[2];
sx q[2];
rz(0.48223334) q[2];
rz(0.33384791) q[3];
sx q[3];
rz(-2.4945989) q[3];
sx q[3];
rz(0.37024959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
