OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.44242087) q[0];
sx q[0];
rz(-2.3306263) q[0];
sx q[0];
rz(-0.45642689) q[0];
rz(2.2189848) q[1];
sx q[1];
rz(-0.95093095) q[1];
sx q[1];
rz(2.9589597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1862926) q[0];
sx q[0];
rz(-2.4432123) q[0];
sx q[0];
rz(1.125505) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6668197) q[2];
sx q[2];
rz(-2.2150196) q[2];
sx q[2];
rz(0.34851532) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4306972) q[1];
sx q[1];
rz(-0.17696807) q[1];
sx q[1];
rz(-2.9269993) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6218518) q[3];
sx q[3];
rz(-1.9155353) q[3];
sx q[3];
rz(2.910579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7543588) q[2];
sx q[2];
rz(-2.0626455) q[2];
sx q[2];
rz(-2.8299502) q[2];
rz(-1.8841057) q[3];
sx q[3];
rz(-0.25185549) q[3];
sx q[3];
rz(-2.9529412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99609128) q[0];
sx q[0];
rz(-1.6956734) q[0];
sx q[0];
rz(-1.2392932) q[0];
rz(-0.39341012) q[1];
sx q[1];
rz(-1.0478123) q[1];
sx q[1];
rz(-1.0308824) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9728972) q[0];
sx q[0];
rz(-2.4542232) q[0];
sx q[0];
rz(2.1363791) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3412881) q[2];
sx q[2];
rz(-1.82845) q[2];
sx q[2];
rz(0.15783707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.70192897) q[1];
sx q[1];
rz(-1.2234634) q[1];
sx q[1];
rz(2.2442978) q[1];
rz(2.2575284) q[3];
sx q[3];
rz(-0.86694781) q[3];
sx q[3];
rz(2.2313909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3139412) q[2];
sx q[2];
rz(-2.2590019) q[2];
sx q[2];
rz(2.7871056) q[2];
rz(-2.3101824) q[3];
sx q[3];
rz(-0.76787132) q[3];
sx q[3];
rz(1.6262511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77446929) q[0];
sx q[0];
rz(-2.9614083) q[0];
sx q[0];
rz(-2.6572976) q[0];
rz(-2.2593185) q[1];
sx q[1];
rz(-0.87119281) q[1];
sx q[1];
rz(0.54214111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.079744) q[0];
sx q[0];
rz(-1.925191) q[0];
sx q[0];
rz(-2.8097879) q[0];
rz(-pi) q[1];
rz(-2.8201879) q[2];
sx q[2];
rz(-1.5084477) q[2];
sx q[2];
rz(-0.37918175) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0391991) q[1];
sx q[1];
rz(-2.3619283) q[1];
sx q[1];
rz(1.8236266) q[1];
rz(-pi) q[2];
x q[2];
rz(0.063850689) q[3];
sx q[3];
rz(-1.8076767) q[3];
sx q[3];
rz(-0.76343918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.813628) q[2];
sx q[2];
rz(-2.4317604) q[2];
sx q[2];
rz(2.1072809) q[2];
rz(1.3616925) q[3];
sx q[3];
rz(-1.1797649) q[3];
sx q[3];
rz(-0.55083197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66083241) q[0];
sx q[0];
rz(-1.3897422) q[0];
sx q[0];
rz(-2.0939636) q[0];
rz(1.818559) q[1];
sx q[1];
rz(-2.5593457) q[1];
sx q[1];
rz(0.38527647) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91222969) q[0];
sx q[0];
rz(-1.9989387) q[0];
sx q[0];
rz(2.1332801) q[0];
rz(-pi) q[1];
rz(-1.6462506) q[2];
sx q[2];
rz(-2.2860043) q[2];
sx q[2];
rz(-2.070141) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.90159081) q[1];
sx q[1];
rz(-1.4996254) q[1];
sx q[1];
rz(-2.932697) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.436211) q[3];
sx q[3];
rz(-2.453605) q[3];
sx q[3];
rz(2.1768895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4177527) q[2];
sx q[2];
rz(-2.1989792) q[2];
sx q[2];
rz(0.27457944) q[2];
rz(-2.5802021) q[3];
sx q[3];
rz(-1.0654819) q[3];
sx q[3];
rz(1.6339462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0059589) q[0];
sx q[0];
rz(-0.57666403) q[0];
sx q[0];
rz(-2.2744001) q[0];
rz(2.7658956) q[1];
sx q[1];
rz(-0.58940327) q[1];
sx q[1];
rz(-1.2295178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6331785) q[0];
sx q[0];
rz(-1.408218) q[0];
sx q[0];
rz(2.6988217) q[0];
rz(-pi) q[1];
rz(2.0784573) q[2];
sx q[2];
rz(-2.5637321) q[2];
sx q[2];
rz(-2.7323728) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7244959) q[1];
sx q[1];
rz(-2.1388106) q[1];
sx q[1];
rz(0.54996164) q[1];
x q[2];
rz(2.8035012) q[3];
sx q[3];
rz(-1.9623358) q[3];
sx q[3];
rz(-0.028926802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.025297252) q[2];
sx q[2];
rz(-1.9090434) q[2];
sx q[2];
rz(3.0787025) q[2];
rz(-2.7316015) q[3];
sx q[3];
rz(-0.72628179) q[3];
sx q[3];
rz(-2.9819152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.9922239) q[0];
sx q[0];
rz(-0.069644444) q[0];
sx q[0];
rz(-0.83576354) q[0];
rz(1.1831076) q[1];
sx q[1];
rz(-1.8712021) q[1];
sx q[1];
rz(-2.3979208) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90172529) q[0];
sx q[0];
rz(-2.5081303) q[0];
sx q[0];
rz(2.6153436) q[0];
x q[1];
rz(-1.7149107) q[2];
sx q[2];
rz(-1.6246572) q[2];
sx q[2];
rz(-0.28648057) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3144296) q[1];
sx q[1];
rz(-1.7531698) q[1];
sx q[1];
rz(-0.86079396) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8858769) q[3];
sx q[3];
rz(-1.5195432) q[3];
sx q[3];
rz(-2.0699786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30367294) q[2];
sx q[2];
rz(-0.49771365) q[2];
sx q[2];
rz(2.3480603) q[2];
rz(-2.7225336) q[3];
sx q[3];
rz(-1.6520809) q[3];
sx q[3];
rz(0.52217531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1917052) q[0];
sx q[0];
rz(-2.1623623) q[0];
sx q[0];
rz(1.1631843) q[0];
rz(0.78041068) q[1];
sx q[1];
rz(-0.33640948) q[1];
sx q[1];
rz(1.6132678) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1799729) q[0];
sx q[0];
rz(-2.442217) q[0];
sx q[0];
rz(1.4699303) q[0];
rz(-pi) q[1];
rz(-2.0397908) q[2];
sx q[2];
rz(-2.7534979) q[2];
sx q[2];
rz(2.3562252) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6032721) q[1];
sx q[1];
rz(-1.5919935) q[1];
sx q[1];
rz(1.3726329) q[1];
rz(2.8655474) q[3];
sx q[3];
rz(-1.1624884) q[3];
sx q[3];
rz(-2.9272872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0344326) q[2];
sx q[2];
rz(-1.1792504) q[2];
sx q[2];
rz(0.068664702) q[2];
rz(0.56985235) q[3];
sx q[3];
rz(-0.48044258) q[3];
sx q[3];
rz(-0.3705875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625075) q[0];
sx q[0];
rz(-2.469049) q[0];
sx q[0];
rz(-1.3154718) q[0];
rz(-1.8585809) q[1];
sx q[1];
rz(-0.43671572) q[1];
sx q[1];
rz(3.0924996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4955208) q[0];
sx q[0];
rz(-1.6611413) q[0];
sx q[0];
rz(3.0642088) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95048381) q[2];
sx q[2];
rz(-1.1178218) q[2];
sx q[2];
rz(-2.9633303) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.093542592) q[1];
sx q[1];
rz(-2.258) q[1];
sx q[1];
rz(-0.50205135) q[1];
rz(-0.92691874) q[3];
sx q[3];
rz(-1.4689494) q[3];
sx q[3];
rz(-0.97130126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7580238) q[2];
sx q[2];
rz(-0.87791666) q[2];
sx q[2];
rz(2.6007268) q[2];
rz(2.0567549) q[3];
sx q[3];
rz(-2.4386051) q[3];
sx q[3];
rz(-1.3802403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2329907) q[0];
sx q[0];
rz(-0.99642307) q[0];
sx q[0];
rz(-0.10502271) q[0];
rz(-2.5999293) q[1];
sx q[1];
rz(-0.88625208) q[1];
sx q[1];
rz(0.35596102) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5585238) q[0];
sx q[0];
rz(-1.6083058) q[0];
sx q[0];
rz(-1.7100699) q[0];
rz(-pi) q[1];
rz(-0.84635205) q[2];
sx q[2];
rz(-2.0931819) q[2];
sx q[2];
rz(-1.236793) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7193422) q[1];
sx q[1];
rz(-1.0021035) q[1];
sx q[1];
rz(-2.689792) q[1];
rz(-pi) q[2];
rz(1.8907053) q[3];
sx q[3];
rz(-1.6875132) q[3];
sx q[3];
rz(0.12091431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9291222) q[2];
sx q[2];
rz(-1.6419623) q[2];
sx q[2];
rz(-1.3193725) q[2];
rz(1.8170554) q[3];
sx q[3];
rz(-1.5683441) q[3];
sx q[3];
rz(-0.96778473) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9318555) q[0];
sx q[0];
rz(-0.034448817) q[0];
sx q[0];
rz(-1.6784278) q[0];
rz(-1.6819008) q[1];
sx q[1];
rz(-1.9821143) q[1];
sx q[1];
rz(-2.7957338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1347156) q[0];
sx q[0];
rz(-1.6505989) q[0];
sx q[0];
rz(0.047264506) q[0];
rz(2.9458617) q[2];
sx q[2];
rz(-0.32265857) q[2];
sx q[2];
rz(0.1041854) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9649557) q[1];
sx q[1];
rz(-1.1208911) q[1];
sx q[1];
rz(-0.69590203) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92697619) q[3];
sx q[3];
rz(-1.2937574) q[3];
sx q[3];
rz(-1.2848971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2807002) q[2];
sx q[2];
rz(-1.8701376) q[2];
sx q[2];
rz(-0.26091179) q[2];
rz(1.4704618) q[3];
sx q[3];
rz(-1.4025531) q[3];
sx q[3];
rz(-1.4298593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83229257) q[0];
sx q[0];
rz(-2.0335048) q[0];
sx q[0];
rz(-0.19620398) q[0];
rz(2.590754) q[1];
sx q[1];
rz(-1.486634) q[1];
sx q[1];
rz(-2.6015729) q[1];
rz(2.9020799) q[2];
sx q[2];
rz(-0.85776599) q[2];
sx q[2];
rz(0.89606482) q[2];
rz(-0.73055406) q[3];
sx q[3];
rz(-1.7777705) q[3];
sx q[3];
rz(-0.70845778) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
