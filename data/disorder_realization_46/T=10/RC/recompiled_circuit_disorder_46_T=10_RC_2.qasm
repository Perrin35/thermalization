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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5292408) q[0];
sx q[0];
rz(-1.2623598) q[0];
sx q[0];
rz(0.2150857) q[0];
x q[1];
rz(0.34502132) q[2];
sx q[2];
rz(-2.3292688) q[2];
sx q[2];
rz(1.5242087) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5073519) q[1];
sx q[1];
rz(-1.5841286) q[1];
sx q[1];
rz(1.3807382) q[1];
rz(-pi) q[2];
rz(0.052993943) q[3];
sx q[3];
rz(-0.64265673) q[3];
sx q[3];
rz(1.3621804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44148579) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(0.95735615) q[2];
rz(2.8422614) q[3];
sx q[3];
rz(-2.744031) q[3];
sx q[3];
rz(0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1608202) q[0];
sx q[0];
rz(-0.94962025) q[0];
sx q[0];
rz(0.41369307) q[0];
rz(-1.3445688) q[1];
sx q[1];
rz(-2.3584056) q[1];
sx q[1];
rz(-0.63562524) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5168415) q[0];
sx q[0];
rz(-3.0020614) q[0];
sx q[0];
rz(-1.2463039) q[0];
rz(-pi) q[1];
rz(3.0947729) q[2];
sx q[2];
rz(-1.702582) q[2];
sx q[2];
rz(-1.552396) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0932514) q[1];
sx q[1];
rz(-1.3158568) q[1];
sx q[1];
rz(1.3202207) q[1];
x q[2];
rz(-0.73511519) q[3];
sx q[3];
rz(-0.8494091) q[3];
sx q[3];
rz(-1.2934367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.7522493) q[2];
sx q[2];
rz(-1.8698591) q[2];
sx q[2];
rz(0.52865571) q[2];
rz(1.901249) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(-0.24578978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56025958) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(2.8833959) q[0];
rz(-1.5064346) q[1];
sx q[1];
rz(-0.55748993) q[1];
sx q[1];
rz(-0.43446508) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6601352) q[0];
sx q[0];
rz(-2.0148206) q[0];
sx q[0];
rz(2.7399979) q[0];
rz(-0.065504727) q[2];
sx q[2];
rz(-1.0422049) q[2];
sx q[2];
rz(-2.5158109) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9745969) q[1];
sx q[1];
rz(-1.7803065) q[1];
sx q[1];
rz(-2.6796883) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4934019) q[3];
sx q[3];
rz(-2.18581) q[3];
sx q[3];
rz(2.0730413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7210641) q[2];
sx q[2];
rz(-1.8177744) q[2];
sx q[2];
rz(-0.055796441) q[2];
rz(-2.2037286) q[3];
sx q[3];
rz(-0.28356975) q[3];
sx q[3];
rz(2.8985033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0980804) q[0];
sx q[0];
rz(-0.94399095) q[0];
sx q[0];
rz(-3.0010624) q[0];
rz(-0.17164104) q[1];
sx q[1];
rz(-1.3069897) q[1];
sx q[1];
rz(-0.26352873) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4123654) q[0];
sx q[0];
rz(-2.1365039) q[0];
sx q[0];
rz(-1.6914781) q[0];
x q[1];
rz(-2.5272335) q[2];
sx q[2];
rz(-2.376997) q[2];
sx q[2];
rz(-1.6524397) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4315223) q[1];
sx q[1];
rz(-0.78774161) q[1];
sx q[1];
rz(2.6452933) q[1];
x q[2];
rz(1.9043546) q[3];
sx q[3];
rz(-1.3106489) q[3];
sx q[3];
rz(-0.0098269193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.043109743) q[2];
sx q[2];
rz(-2.7973599) q[2];
sx q[2];
rz(1.8919224) q[2];
rz(1.9469056) q[3];
sx q[3];
rz(-1.0281111) q[3];
sx q[3];
rz(2.4627114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2018305) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(0.029065954) q[0];
rz(-1.4020231) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(0.32863858) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12954457) q[0];
sx q[0];
rz(-2.64376) q[0];
sx q[0];
rz(-1.7853907) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5673037) q[2];
sx q[2];
rz(-0.31775489) q[2];
sx q[2];
rz(-0.53900915) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2480093) q[1];
sx q[1];
rz(-1.5072522) q[1];
sx q[1];
rz(-2.7081972) q[1];
x q[2];
rz(-0.5248431) q[3];
sx q[3];
rz(-2.3224324) q[3];
sx q[3];
rz(-0.51845779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0236726) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(-0.16061352) q[2];
rz(-1.1462071) q[3];
sx q[3];
rz(-1.3137772) q[3];
sx q[3];
rz(2.5207991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6495431) q[0];
sx q[0];
rz(-0.8594802) q[0];
sx q[0];
rz(2.3550526) q[0];
rz(-0.37711626) q[1];
sx q[1];
rz(-2.2943594) q[1];
sx q[1];
rz(-0.4424817) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8724976) q[0];
sx q[0];
rz(-1.9364534) q[0];
sx q[0];
rz(-0.38537607) q[0];
rz(-2.5673893) q[2];
sx q[2];
rz(-1.310537) q[2];
sx q[2];
rz(-1.1299709) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1861371) q[1];
sx q[1];
rz(-1.8147787) q[1];
sx q[1];
rz(-1.2402435) q[1];
x q[2];
rz(-3.0480012) q[3];
sx q[3];
rz(-1.2564661) q[3];
sx q[3];
rz(-0.67925727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4601712) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(-2.1703413) q[2];
rz(2.4462637) q[3];
sx q[3];
rz(-0.23614241) q[3];
sx q[3];
rz(-0.42738459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5044395) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(2.1321645) q[0];
rz(0.55074739) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(-1.1632464) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3488719) q[0];
sx q[0];
rz(-2.7808245) q[0];
sx q[0];
rz(1.1611206) q[0];
x q[1];
rz(1.2802358) q[2];
sx q[2];
rz(-0.72939789) q[2];
sx q[2];
rz(-1.1682208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0102331) q[1];
sx q[1];
rz(-1.3512003) q[1];
sx q[1];
rz(-1.4036199) q[1];
rz(-0.87485119) q[3];
sx q[3];
rz(-1.3303183) q[3];
sx q[3];
rz(-1.8868173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.78272468) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(0.7981832) q[2];
rz(2.362137) q[3];
sx q[3];
rz(-0.53648406) q[3];
sx q[3];
rz(1.1727758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1786132) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(-3.0187507) q[0];
rz(3.0154862) q[1];
sx q[1];
rz(-1.6364731) q[1];
sx q[1];
rz(-1.2164446) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.018369) q[0];
sx q[0];
rz(-2.1835678) q[0];
sx q[0];
rz(-1.0606674) q[0];
x q[1];
rz(-1.4037651) q[2];
sx q[2];
rz(-2.0207496) q[2];
sx q[2];
rz(-0.049875967) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0876906) q[1];
sx q[1];
rz(-1.1864099) q[1];
sx q[1];
rz(1.7840506) q[1];
x q[2];
rz(1.7580732) q[3];
sx q[3];
rz(-0.98687275) q[3];
sx q[3];
rz(-2.5919979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2609666) q[2];
sx q[2];
rz(-2.5239021) q[2];
sx q[2];
rz(-3.0984666) q[2];
rz(0.17523266) q[3];
sx q[3];
rz(-0.8774811) q[3];
sx q[3];
rz(1.5847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6373428) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(0.20877008) q[0];
rz(1.4978706) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(2.0170905) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9546684) q[0];
sx q[0];
rz(-2.9483729) q[0];
sx q[0];
rz(0.78878553) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3051946) q[2];
sx q[2];
rz(-2.2679066) q[2];
sx q[2];
rz(0.12649378) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75301718) q[1];
sx q[1];
rz(-1.5696226) q[1];
sx q[1];
rz(1.7041824) q[1];
x q[2];
rz(0.87499683) q[3];
sx q[3];
rz(-1.1831302) q[3];
sx q[3];
rz(0.33867237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1407397) q[2];
sx q[2];
rz(-0.72379392) q[2];
sx q[2];
rz(-2.9635584) q[2];
rz(1.595165) q[3];
sx q[3];
rz(-1.3083357) q[3];
sx q[3];
rz(2.6509638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.35266018) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(2.1066522) q[0];
rz(0.79822284) q[1];
sx q[1];
rz(-2.1612576) q[1];
sx q[1];
rz(2.9916874) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5750654) q[0];
sx q[0];
rz(-1.2034104) q[0];
sx q[0];
rz(-0.67879403) q[0];
rz(-pi) q[1];
x q[1];
rz(2.855004) q[2];
sx q[2];
rz(-1.5865241) q[2];
sx q[2];
rz(-2.36731) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3372911) q[1];
sx q[1];
rz(-1.9807439) q[1];
sx q[1];
rz(1.0211584) q[1];
rz(2.4622915) q[3];
sx q[3];
rz(-1.3759817) q[3];
sx q[3];
rz(0.65115813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4118816) q[2];
sx q[2];
rz(-1.9103266) q[2];
sx q[2];
rz(0.40851545) q[2];
rz(-0.92489964) q[3];
sx q[3];
rz(-0.81245208) q[3];
sx q[3];
rz(0.48172054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9020486) q[0];
sx q[0];
rz(-1.5938546) q[0];
sx q[0];
rz(-1.6123733) q[0];
rz(-1.3760024) q[1];
sx q[1];
rz(-1.1767495) q[1];
sx q[1];
rz(-1.8935988) q[1];
rz(-0.79509905) q[2];
sx q[2];
rz(-1.8021402) q[2];
sx q[2];
rz(-2.6593593) q[2];
rz(1.8134712) q[3];
sx q[3];
rz(-0.96488733) q[3];
sx q[3];
rz(0.78028954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
