OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5912136) q[0];
sx q[0];
rz(-0.0033291078) q[0];
sx q[0];
rz(2.79628) q[0];
rz(1.8058864) q[1];
sx q[1];
rz(-2.8023281) q[1];
sx q[1];
rz(-0.27944922) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011443519) q[0];
sx q[0];
rz(-0.37405095) q[0];
sx q[0];
rz(0.98056294) q[0];
rz(-pi) q[1];
rz(-1.9136393) q[2];
sx q[2];
rz(-0.81878412) q[2];
sx q[2];
rz(-2.0057099) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2024723) q[1];
sx q[1];
rz(-1.7608374) q[1];
sx q[1];
rz(-0.01357667) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.052993943) q[3];
sx q[3];
rz(-0.64265673) q[3];
sx q[3];
rz(-1.3621804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7001069) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(-0.95735615) q[2];
rz(-2.8422614) q[3];
sx q[3];
rz(-2.744031) q[3];
sx q[3];
rz(-0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
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
rz(0.63562524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3755075) q[0];
sx q[0];
rz(-1.6151531) q[0];
sx q[0];
rz(1.7031329) q[0];
rz(1.4388678) q[2];
sx q[2];
rz(-1.6172098) q[2];
sx q[2];
rz(3.1293491) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8860652) q[1];
sx q[1];
rz(-0.35554245) q[1];
sx q[1];
rz(-2.38106) q[1];
rz(-0.91931822) q[3];
sx q[3];
rz(-2.161536) q[3];
sx q[3];
rz(0.35349333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3893434) q[2];
sx q[2];
rz(-1.8698591) q[2];
sx q[2];
rz(2.6129369) q[2];
rz(1.2403437) q[3];
sx q[3];
rz(-0.35959187) q[3];
sx q[3];
rz(2.8958029) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56025958) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(-2.8833959) q[0];
rz(-1.5064346) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(-2.7071276) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4397944) q[0];
sx q[0];
rz(-2.5520303) q[0];
sx q[0];
rz(-0.88296417) q[0];
rz(-pi) q[1];
rz(-1.6824109) q[2];
sx q[2];
rz(-0.53225213) q[2];
sx q[2];
rz(-0.75512952) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30063054) q[1];
sx q[1];
rz(-1.1197487) q[1];
sx q[1];
rz(-1.3375963) q[1];
x q[2];
rz(2.2779949) q[3];
sx q[3];
rz(-2.2798385) q[3];
sx q[3];
rz(1.1532702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7210641) q[2];
sx q[2];
rz(-1.8177744) q[2];
sx q[2];
rz(0.055796441) q[2];
rz(0.93786401) q[3];
sx q[3];
rz(-2.8580229) q[3];
sx q[3];
rz(-2.8985033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0980804) q[0];
sx q[0];
rz(-2.1976017) q[0];
sx q[0];
rz(3.0010624) q[0];
rz(-0.17164104) q[1];
sx q[1];
rz(-1.3069897) q[1];
sx q[1];
rz(-0.26352873) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9182501) q[0];
sx q[0];
rz(-1.6726057) q[0];
sx q[0];
rz(-0.56901594) q[0];
rz(-pi) q[1];
rz(2.5272335) q[2];
sx q[2];
rz(-2.376997) q[2];
sx q[2];
rz(-1.489153) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0553953) q[1];
sx q[1];
rz(-2.2436884) q[1];
sx q[1];
rz(1.124568) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2534091) q[3];
sx q[3];
rz(-0.41999751) q[3];
sx q[3];
rz(-2.1995467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(1.2496703) q[2];
rz(-1.9469056) q[3];
sx q[3];
rz(-1.0281111) q[3];
sx q[3];
rz(0.67888129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93976218) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(0.029065954) q[0];
rz(1.4020231) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(2.8129541) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0279854) q[0];
sx q[0];
rz(-1.0853882) q[0];
sx q[0];
rz(0.11522449) q[0];
rz(-pi) q[1];
rz(2.5673037) q[2];
sx q[2];
rz(-0.31775489) q[2];
sx q[2];
rz(2.6025835) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6824274) q[1];
sx q[1];
rz(-2.703856) q[1];
sx q[1];
rz(-0.15037219) q[1];
rz(-0.7469437) q[3];
sx q[3];
rz(-1.1960104) q[3];
sx q[3];
rz(2.4657616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0236726) q[2];
sx q[2];
rz(-0.4549883) q[2];
sx q[2];
rz(-2.9809791) q[2];
rz(-1.9953856) q[3];
sx q[3];
rz(-1.8278154) q[3];
sx q[3];
rz(-0.62079352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49204957) q[0];
sx q[0];
rz(-0.8594802) q[0];
sx q[0];
rz(-0.78654003) q[0];
rz(-0.37711626) q[1];
sx q[1];
rz(-2.2943594) q[1];
sx q[1];
rz(2.699111) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2690951) q[0];
sx q[0];
rz(-1.9364534) q[0];
sx q[0];
rz(-0.38537607) q[0];
rz(2.5673893) q[2];
sx q[2];
rz(-1.8310556) q[2];
sx q[2];
rz(-1.1299709) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1861371) q[1];
sx q[1];
rz(-1.8147787) q[1];
sx q[1];
rz(1.9013491) q[1];
x q[2];
rz(0.093591452) q[3];
sx q[3];
rz(-1.8851265) q[3];
sx q[3];
rz(0.67925727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4601712) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(2.1703413) q[2];
rz(0.69532895) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(-0.42738459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-0.63715315) q[0];
sx q[0];
rz(-1.9313066) q[0];
sx q[0];
rz(2.1321645) q[0];
rz(2.5908453) q[1];
sx q[1];
rz(-1.2754722) q[1];
sx q[1];
rz(-1.1632464) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7834085) q[0];
sx q[0];
rz(-1.9005214) q[0];
sx q[0];
rz(0.14915906) q[0];
x q[1];
rz(0.86261729) q[2];
sx q[2];
rz(-1.3786945) q[2];
sx q[2];
rz(0.62193279) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7904661) q[1];
sx q[1];
rz(-0.27517056) q[1];
sx q[1];
rz(2.5009584) q[1];
rz(-pi) q[2];
x q[2];
rz(1.936124) q[3];
sx q[3];
rz(-2.4118773) q[3];
sx q[3];
rz(-3.1033033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.358868) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(0.7981832) q[2];
rz(-2.362137) q[3];
sx q[3];
rz(-2.6051086) q[3];
sx q[3];
rz(-1.9688169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9629795) q[0];
sx q[0];
rz(-0.48291746) q[0];
sx q[0];
rz(-3.0187507) q[0];
rz(0.12610647) q[1];
sx q[1];
rz(-1.5051196) q[1];
sx q[1];
rz(-1.2164446) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.018369) q[0];
sx q[0];
rz(-0.95802486) q[0];
sx q[0];
rz(-2.0809253) q[0];
rz(-pi) q[1];
rz(-0.33151303) q[2];
sx q[2];
rz(-2.6636332) q[2];
sx q[2];
rz(0.31994672) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.053902102) q[1];
sx q[1];
rz(-1.9551827) q[1];
sx q[1];
rz(1.7840506) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8669531) q[3];
sx q[3];
rz(-2.5317149) q[3];
sx q[3];
rz(0.21851893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.2609666) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(-0.043126062) q[2];
rz(0.17523266) q[3];
sx q[3];
rz(-0.8774811) q[3];
sx q[3];
rz(-1.5568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50424987) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(2.9328226) q[0];
rz(-1.4978706) q[1];
sx q[1];
rz(-1.4893963) q[1];
sx q[1];
rz(-1.1245022) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61125206) q[0];
sx q[0];
rz(-1.4350622) q[0];
sx q[0];
rz(1.432857) q[0];
x q[1];
rz(0.71473177) q[2];
sx q[2];
rz(-1.3681612) q[2];
sx q[2];
rz(-1.5243901) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3325602) q[1];
sx q[1];
rz(-0.13339116) q[1];
sx q[1];
rz(1.5619713) q[1];
rz(-pi) q[2];
rz(-0.87499683) q[3];
sx q[3];
rz(-1.9584624) q[3];
sx q[3];
rz(0.33867237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1407397) q[2];
sx q[2];
rz(-2.4177987) q[2];
sx q[2];
rz(-0.17803426) q[2];
rz(1.595165) q[3];
sx q[3];
rz(-1.3083357) q[3];
sx q[3];
rz(2.6509638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7889325) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(-2.1066522) q[0];
rz(-2.3433698) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(-2.9916874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7222544) q[0];
sx q[0];
rz(-2.1968578) q[0];
sx q[0];
rz(-2.0300079) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5543999) q[2];
sx q[2];
rz(-1.2842442) q[2];
sx q[2];
rz(0.80114844) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0059709) q[1];
sx q[1];
rz(-2.0704198) q[1];
sx q[1];
rz(-2.6702704) q[1];
rz(-pi) q[2];
rz(1.8191765) q[3];
sx q[3];
rz(-2.2349149) q[3];
sx q[3];
rz(0.76457232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7297111) q[2];
sx q[2];
rz(-1.9103266) q[2];
sx q[2];
rz(0.40851545) q[2];
rz(-2.216693) q[3];
sx q[3];
rz(-0.81245208) q[3];
sx q[3];
rz(-0.48172054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2395441) q[0];
sx q[0];
rz(-1.5938546) q[0];
sx q[0];
rz(-1.6123733) q[0];
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
rz(-0.61990191) q[3];
sx q[3];
rz(-1.3719659) q[3];
sx q[3];
rz(-0.93056783) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
