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
rz(-0.69775692) q[0];
sx q[0];
rz(-1.193576) q[0];
sx q[0];
rz(2.9370263) q[0];
rz(2.3808631) q[1];
sx q[1];
rz(-1.7525571) q[1];
sx q[1];
rz(-1.3995481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3239878) q[0];
sx q[0];
rz(-1.3496552) q[0];
sx q[0];
rz(-3.0482685) q[0];
rz(0.84487652) q[2];
sx q[2];
rz(-0.57327081) q[2];
sx q[2];
rz(0.63001652) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5198361) q[1];
sx q[1];
rz(-2.7108828) q[1];
sx q[1];
rz(0.21060305) q[1];
x q[2];
rz(-2.4321626) q[3];
sx q[3];
rz(-1.7755847) q[3];
sx q[3];
rz(0.11573175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.43510398) q[2];
sx q[2];
rz(-0.032328345) q[2];
sx q[2];
rz(0.48933634) q[2];
rz(-1.0232183) q[3];
sx q[3];
rz(-0.018298572) q[3];
sx q[3];
rz(-2.0160915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0959051) q[0];
sx q[0];
rz(-0.65653312) q[0];
sx q[0];
rz(-2.3387961) q[0];
rz(0.071391694) q[1];
sx q[1];
rz(-0.26669058) q[1];
sx q[1];
rz(0.057770483) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5481541) q[0];
sx q[0];
rz(-1.8280941) q[0];
sx q[0];
rz(2.1307039) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8550496) q[2];
sx q[2];
rz(-1.4689969) q[2];
sx q[2];
rz(1.941178) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.80177414) q[1];
sx q[1];
rz(-1.8027824) q[1];
sx q[1];
rz(1.07527) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9936419) q[3];
sx q[3];
rz(-0.86455621) q[3];
sx q[3];
rz(1.0195635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1959261) q[2];
sx q[2];
rz(-2.0464996) q[2];
sx q[2];
rz(-1.2954953) q[2];
rz(-2.1748491) q[3];
sx q[3];
rz(-0.77015489) q[3];
sx q[3];
rz(-2.3841592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1118689) q[0];
sx q[0];
rz(-1.3628549) q[0];
sx q[0];
rz(1.4335853) q[0];
rz(0.068610527) q[1];
sx q[1];
rz(-1.5674633) q[1];
sx q[1];
rz(0.57919085) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1807208) q[0];
sx q[0];
rz(-1.6287043) q[0];
sx q[0];
rz(-1.6625893) q[0];
x q[1];
rz(0.4757232) q[2];
sx q[2];
rz(-1.5724725) q[2];
sx q[2];
rz(0.36538183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.80651307) q[1];
sx q[1];
rz(-1.5349704) q[1];
sx q[1];
rz(2.9153009) q[1];
x q[2];
rz(-0.45626943) q[3];
sx q[3];
rz(-1.1351403) q[3];
sx q[3];
rz(1.6317451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27089831) q[2];
sx q[2];
rz(-1.8879994) q[2];
sx q[2];
rz(-2.9581621) q[2];
rz(2.3373248) q[3];
sx q[3];
rz(-0.99884123) q[3];
sx q[3];
rz(-2.6075294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19876984) q[0];
sx q[0];
rz(-0.11480055) q[0];
sx q[0];
rz(2.542069) q[0];
rz(-2.7291258) q[1];
sx q[1];
rz(-3.1209374) q[1];
sx q[1];
rz(-2.1309158) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3314047) q[0];
sx q[0];
rz(-1.4443026) q[0];
sx q[0];
rz(0.35270113) q[0];
rz(-pi) q[1];
rz(-0.93866556) q[2];
sx q[2];
rz(-1.5297339) q[2];
sx q[2];
rz(1.9384042) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4356723) q[1];
sx q[1];
rz(-1.3002035) q[1];
sx q[1];
rz(-2.1822004) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14673503) q[3];
sx q[3];
rz(-0.76125604) q[3];
sx q[3];
rz(-0.066731922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3707054) q[2];
sx q[2];
rz(-0.34660307) q[2];
sx q[2];
rz(-0.31751219) q[2];
rz(-0.62234771) q[3];
sx q[3];
rz(-0.93572664) q[3];
sx q[3];
rz(0.58421016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36251003) q[0];
sx q[0];
rz(-2.2267987) q[0];
sx q[0];
rz(2.1146178) q[0];
rz(2.5912071) q[1];
sx q[1];
rz(-3.0774979) q[1];
sx q[1];
rz(1.9245573) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69637042) q[0];
sx q[0];
rz(-2.1289941) q[0];
sx q[0];
rz(2.3675015) q[0];
rz(-pi) q[1];
rz(0.12995692) q[2];
sx q[2];
rz(-2.5995147) q[2];
sx q[2];
rz(-0.96058577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40081208) q[1];
sx q[1];
rz(-2.2601839) q[1];
sx q[1];
rz(-2.9184398) q[1];
rz(-pi) q[2];
rz(2.2812165) q[3];
sx q[3];
rz(-1.851322) q[3];
sx q[3];
rz(2.4376412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7050742) q[2];
sx q[2];
rz(-1.7784092) q[2];
sx q[2];
rz(2.3684033) q[2];
rz(0.1117205) q[3];
sx q[3];
rz(-1.819928) q[3];
sx q[3];
rz(2.2022061) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47950995) q[0];
sx q[0];
rz(-0.22308068) q[0];
sx q[0];
rz(-0.42698419) q[0];
rz(-2.2110979) q[1];
sx q[1];
rz(-0.016914802) q[1];
sx q[1];
rz(-0.46447909) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94356774) q[0];
sx q[0];
rz(-2.7901398) q[0];
sx q[0];
rz(0.60011421) q[0];
rz(-pi) q[1];
rz(2.0468087) q[2];
sx q[2];
rz(-0.65611984) q[2];
sx q[2];
rz(0.5364272) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.61349866) q[1];
sx q[1];
rz(-1.534675) q[1];
sx q[1];
rz(-2.7725459) q[1];
rz(-pi) q[2];
rz(0.031207009) q[3];
sx q[3];
rz(-1.1466807) q[3];
sx q[3];
rz(1.1348789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53133196) q[2];
sx q[2];
rz(-1.320763) q[2];
sx q[2];
rz(-0.28826928) q[2];
rz(-2.098295) q[3];
sx q[3];
rz(-0.61560029) q[3];
sx q[3];
rz(-0.73879761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4814602) q[0];
sx q[0];
rz(-1.8539424) q[0];
sx q[0];
rz(2.3840391) q[0];
rz(0.07269147) q[1];
sx q[1];
rz(-3.1158267) q[1];
sx q[1];
rz(0.048197897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5023247) q[0];
sx q[0];
rz(-1.6002965) q[0];
sx q[0];
rz(-0.98639368) q[0];
rz(-2.169007) q[2];
sx q[2];
rz(-1.3763104) q[2];
sx q[2];
rz(-2.011428) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0695677) q[1];
sx q[1];
rz(-1.1434907) q[1];
sx q[1];
rz(2.2263889) q[1];
rz(-0.052560135) q[3];
sx q[3];
rz(-1.1014928) q[3];
sx q[3];
rz(1.8183501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4797719) q[2];
sx q[2];
rz(-1.4996108) q[2];
sx q[2];
rz(-3.0721967) q[2];
rz(1.5875459) q[3];
sx q[3];
rz(-0.78726751) q[3];
sx q[3];
rz(-2.9230996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3875535) q[0];
sx q[0];
rz(-2.0856922) q[0];
sx q[0];
rz(1.4132502) q[0];
rz(-2.3130401) q[1];
sx q[1];
rz(-3.0998402) q[1];
sx q[1];
rz(-0.54263306) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4601645) q[0];
sx q[0];
rz(-1.6896473) q[0];
sx q[0];
rz(2.9959841) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75730308) q[2];
sx q[2];
rz(-2.6432163) q[2];
sx q[2];
rz(1.9973081) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5396351) q[1];
sx q[1];
rz(-1.6865206) q[1];
sx q[1];
rz(1.5166111) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7067634) q[3];
sx q[3];
rz(-1.9169589) q[3];
sx q[3];
rz(2.1157672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1594306) q[2];
sx q[2];
rz(-0.40951481) q[2];
sx q[2];
rz(-0.25992599) q[2];
rz(1.0204756) q[3];
sx q[3];
rz(-2.8838938) q[3];
sx q[3];
rz(0.79403383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6771773) q[0];
sx q[0];
rz(-2.9554415) q[0];
sx q[0];
rz(1.6625241) q[0];
rz(-1.5326477) q[1];
sx q[1];
rz(-2.1023991) q[1];
sx q[1];
rz(-2.398568) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31194132) q[0];
sx q[0];
rz(-0.52642979) q[0];
sx q[0];
rz(-1.3408324) q[0];
x q[1];
rz(2.114562) q[2];
sx q[2];
rz(-1.9025505) q[2];
sx q[2];
rz(-0.55032718) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2496532) q[1];
sx q[1];
rz(-1.5158476) q[1];
sx q[1];
rz(-1.5049388) q[1];
x q[2];
rz(-1.5398272) q[3];
sx q[3];
rz(-1.4659229) q[3];
sx q[3];
rz(-0.34377835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6741901) q[2];
sx q[2];
rz(-2.3154066) q[2];
sx q[2];
rz(2.3361333) q[2];
rz(1.6921267) q[3];
sx q[3];
rz(-1.9129246) q[3];
sx q[3];
rz(0.8031556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8886803) q[0];
sx q[0];
rz(-2.5997933) q[0];
sx q[0];
rz(-2.3895277) q[0];
rz(1.9883142) q[1];
sx q[1];
rz(-0.88429943) q[1];
sx q[1];
rz(2.8582252) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3842938) q[0];
sx q[0];
rz(-1.481825) q[0];
sx q[0];
rz(2.6020537) q[0];
rz(-pi) q[1];
rz(-0.14788515) q[2];
sx q[2];
rz(-2.3171632) q[2];
sx q[2];
rz(-1.3651207) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16285208) q[1];
sx q[1];
rz(-1.2417485) q[1];
sx q[1];
rz(0.63905893) q[1];
rz(-pi) q[2];
rz(-0.65431904) q[3];
sx q[3];
rz(-1.2551184) q[3];
sx q[3];
rz(-1.1512427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.89455426) q[2];
sx q[2];
rz(-0.082823195) q[2];
sx q[2];
rz(-1.438633) q[2];
rz(-0.29397193) q[3];
sx q[3];
rz(-0.014466244) q[3];
sx q[3];
rz(1.0283874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033584874) q[0];
sx q[0];
rz(-1.4172194) q[0];
sx q[0];
rz(-1.5237756) q[0];
rz(-0.53957466) q[1];
sx q[1];
rz(-2.3537666) q[1];
sx q[1];
rz(-2.981577) q[1];
rz(1.930936) q[2];
sx q[2];
rz(-1.4317206) q[2];
sx q[2];
rz(1.7862873) q[2];
rz(1.2992819) q[3];
sx q[3];
rz(-1.2012902) q[3];
sx q[3];
rz(1.741878) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
