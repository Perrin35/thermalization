OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9899848) q[0];
sx q[0];
rz(4.0816981) q[0];
sx q[0];
rz(8.8844086) q[0];
rz(-2.6481533) q[1];
sx q[1];
rz(-0.72055888) q[1];
sx q[1];
rz(0.61385733) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9073528) q[0];
sx q[0];
rz(-0.94830238) q[0];
sx q[0];
rz(1.3433775) q[0];
x q[1];
rz(-1.5905321) q[2];
sx q[2];
rz(-2.8497549) q[2];
sx q[2];
rz(2.8449051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.39345523) q[1];
sx q[1];
rz(-0.84987133) q[1];
sx q[1];
rz(1.9490521) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5924388) q[3];
sx q[3];
rz(-1.3236681) q[3];
sx q[3];
rz(-0.46268845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8453688) q[2];
sx q[2];
rz(-1.2426528) q[2];
sx q[2];
rz(-2.5073012) q[2];
rz(-2.2058709) q[3];
sx q[3];
rz(-2.8959385) q[3];
sx q[3];
rz(-2.3989357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.3212386) q[0];
sx q[0];
rz(-1.6004434) q[0];
sx q[0];
rz(2.6974005) q[0];
rz(2.3815637) q[1];
sx q[1];
rz(-2.0030231) q[1];
sx q[1];
rz(-0.98145032) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32641706) q[0];
sx q[0];
rz(-1.8627286) q[0];
sx q[0];
rz(-2.6089704) q[0];
rz(-pi) q[1];
x q[1];
rz(1.74238) q[2];
sx q[2];
rz(-1.7497471) q[2];
sx q[2];
rz(-1.5092261) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.98636282) q[1];
sx q[1];
rz(-1.6587509) q[1];
sx q[1];
rz(2.7025239) q[1];
rz(2.7090453) q[3];
sx q[3];
rz(-1.3290231) q[3];
sx q[3];
rz(1.6530619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0280219) q[2];
sx q[2];
rz(-2.5271723) q[2];
sx q[2];
rz(-1.2051955) q[2];
rz(-0.53660721) q[3];
sx q[3];
rz(-1.2380995) q[3];
sx q[3];
rz(-1.7369778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2236915) q[0];
sx q[0];
rz(-1.2811998) q[0];
sx q[0];
rz(0.83876383) q[0];
rz(-0.63610786) q[1];
sx q[1];
rz(-1.6517703) q[1];
sx q[1];
rz(-3.1210693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8021401) q[0];
sx q[0];
rz(-1.0217371) q[0];
sx q[0];
rz(0.86717506) q[0];
rz(2.0383561) q[2];
sx q[2];
rz(-2.1439432) q[2];
sx q[2];
rz(1.3064885) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1908326) q[1];
sx q[1];
rz(-1.5632707) q[1];
sx q[1];
rz(0.41295596) q[1];
x q[2];
rz(-0.15744029) q[3];
sx q[3];
rz(-1.2886815) q[3];
sx q[3];
rz(2.0705786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3402349) q[2];
sx q[2];
rz(-1.5284208) q[2];
sx q[2];
rz(0.94432962) q[2];
rz(1.0229735) q[3];
sx q[3];
rz(-1.2228271) q[3];
sx q[3];
rz(-2.003722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2706547) q[0];
sx q[0];
rz(-1.9671257) q[0];
sx q[0];
rz(-1.1573855) q[0];
rz(-1.6429139) q[1];
sx q[1];
rz(-1.6694371) q[1];
sx q[1];
rz(-1.9532983) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0152215) q[0];
sx q[0];
rz(-1.0980716) q[0];
sx q[0];
rz(-2.0822099) q[0];
rz(-pi) q[1];
rz(-1.6158197) q[2];
sx q[2];
rz(-2.1189711) q[2];
sx q[2];
rz(2.5803103) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3664361) q[1];
sx q[1];
rz(-0.85298733) q[1];
sx q[1];
rz(-2.3597673) q[1];
x q[2];
rz(1.9096987) q[3];
sx q[3];
rz(-1.7247685) q[3];
sx q[3];
rz(-3.0018501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1184065) q[2];
sx q[2];
rz(-1.1772757) q[2];
sx q[2];
rz(-2.4510621) q[2];
rz(-1.1150507) q[3];
sx q[3];
rz(-0.77459049) q[3];
sx q[3];
rz(-1.640865) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0951776) q[0];
sx q[0];
rz(-1.4354118) q[0];
sx q[0];
rz(2.114356) q[0];
rz(-1.3014303) q[1];
sx q[1];
rz(-2.4899028) q[1];
sx q[1];
rz(0.32381907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9607317) q[0];
sx q[0];
rz(-1.1723851) q[0];
sx q[0];
rz(2.2965527) q[0];
x q[1];
rz(-1.5478163) q[2];
sx q[2];
rz(-0.55292623) q[2];
sx q[2];
rz(-1.9000017) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.370812) q[1];
sx q[1];
rz(-0.64677927) q[1];
sx q[1];
rz(0.5495407) q[1];
rz(-pi) q[2];
rz(2.0824329) q[3];
sx q[3];
rz(-1.4338974) q[3];
sx q[3];
rz(2.9149027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7258437) q[2];
sx q[2];
rz(-2.9314633) q[2];
sx q[2];
rz(-2.6591163) q[2];
rz(-1.1680565) q[3];
sx q[3];
rz(-1.2169633) q[3];
sx q[3];
rz(2.4647958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93550682) q[0];
sx q[0];
rz(-0.85010234) q[0];
sx q[0];
rz(0.8859984) q[0];
rz(-0.63367263) q[1];
sx q[1];
rz(-0.88992563) q[1];
sx q[1];
rz(2.450313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46477485) q[0];
sx q[0];
rz(-2.4949843) q[0];
sx q[0];
rz(0.060641373) q[0];
rz(1.5070314) q[2];
sx q[2];
rz(-0.87997961) q[2];
sx q[2];
rz(-1.7036167) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4659652) q[1];
sx q[1];
rz(-1.4838929) q[1];
sx q[1];
rz(-0.81710941) q[1];
x q[2];
rz(-2.8915358) q[3];
sx q[3];
rz(-2.3365006) q[3];
sx q[3];
rz(-0.4522194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1401691) q[2];
sx q[2];
rz(-0.836335) q[2];
sx q[2];
rz(-2.8720065) q[2];
rz(-0.94830281) q[3];
sx q[3];
rz(-1.6092665) q[3];
sx q[3];
rz(-1.3986826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7798994) q[0];
sx q[0];
rz(-0.43710709) q[0];
sx q[0];
rz(1.0070739) q[0];
rz(2.7359447) q[1];
sx q[1];
rz(-2.5463153) q[1];
sx q[1];
rz(0.55353037) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6104975) q[0];
sx q[0];
rz(-2.880504) q[0];
sx q[0];
rz(1.7946662) q[0];
rz(-pi) q[1];
rz(-2.1184475) q[2];
sx q[2];
rz(-2.0549462) q[2];
sx q[2];
rz(2.1439056) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86846734) q[1];
sx q[1];
rz(-0.9573862) q[1];
sx q[1];
rz(-2.1347743) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5457492) q[3];
sx q[3];
rz(-1.1908997) q[3];
sx q[3];
rz(-2.8632426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3125399) q[2];
sx q[2];
rz(-2.0487787) q[2];
sx q[2];
rz(-1.2139758) q[2];
rz(-2.35516) q[3];
sx q[3];
rz(-1.7460456) q[3];
sx q[3];
rz(-3.1184375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19602747) q[0];
sx q[0];
rz(-0.29226154) q[0];
sx q[0];
rz(-1.8096402) q[0];
rz(-0.61344433) q[1];
sx q[1];
rz(-2.1211801) q[1];
sx q[1];
rz(-0.44949284) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4411104) q[0];
sx q[0];
rz(-1.4876517) q[0];
sx q[0];
rz(1.0821728) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9067326) q[2];
sx q[2];
rz(-1.4892231) q[2];
sx q[2];
rz(2.5334266) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4787138) q[1];
sx q[1];
rz(-0.7906853) q[1];
sx q[1];
rz(1.3776758) q[1];
rz(-pi) q[2];
rz(2.883955) q[3];
sx q[3];
rz(-1.6272568) q[3];
sx q[3];
rz(-0.27654058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.31676644) q[2];
sx q[2];
rz(-2.4591441) q[2];
sx q[2];
rz(2.5332434) q[2];
rz(-1.9781205) q[3];
sx q[3];
rz(-1.3694265) q[3];
sx q[3];
rz(-2.5782862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0332396) q[0];
sx q[0];
rz(-2.2305363) q[0];
sx q[0];
rz(0.60229993) q[0];
rz(2.1741518) q[1];
sx q[1];
rz(-0.93092218) q[1];
sx q[1];
rz(-2.0379351) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3924343) q[0];
sx q[0];
rz(-1.1811678) q[0];
sx q[0];
rz(0.97831877) q[0];
rz(-pi) q[1];
rz(0.63656143) q[2];
sx q[2];
rz(-1.5590073) q[2];
sx q[2];
rz(-1.7280886) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1608354) q[1];
sx q[1];
rz(-2.6904562) q[1];
sx q[1];
rz(1.5734929) q[1];
rz(-pi) q[2];
rz(0.067632631) q[3];
sx q[3];
rz(-0.99376947) q[3];
sx q[3];
rz(0.91995507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2063107) q[2];
sx q[2];
rz(-1.3056359) q[2];
sx q[2];
rz(-0.3375816) q[2];
rz(-2.857699) q[3];
sx q[3];
rz(-0.68113911) q[3];
sx q[3];
rz(-2.9547227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5934481) q[0];
sx q[0];
rz(-1.633506) q[0];
sx q[0];
rz(-0.067807587) q[0];
rz(1.1495122) q[1];
sx q[1];
rz(-1.5510473) q[1];
sx q[1];
rz(0.4745208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7260925) q[0];
sx q[0];
rz(-1.584238) q[0];
sx q[0];
rz(-1.8095762) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6340747) q[2];
sx q[2];
rz(-2.8675277) q[2];
sx q[2];
rz(0.12390359) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2429035) q[1];
sx q[1];
rz(-1.6575282) q[1];
sx q[1];
rz(-3.134722) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81409295) q[3];
sx q[3];
rz(-2.11907) q[3];
sx q[3];
rz(0.89422885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6600251) q[2];
sx q[2];
rz(-2.2109172) q[2];
sx q[2];
rz(2.068326) q[2];
rz(-2.977071) q[3];
sx q[3];
rz(-1.3273032) q[3];
sx q[3];
rz(-0.32065121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5545223) q[0];
sx q[0];
rz(-1.5656492) q[0];
sx q[0];
rz(1.5026305) q[0];
rz(-1.5837689) q[1];
sx q[1];
rz(-2.0638034) q[1];
sx q[1];
rz(2.6149909) q[1];
rz(2.4921992) q[2];
sx q[2];
rz(-0.55772256) q[2];
sx q[2];
rz(-1.100308) q[2];
rz(-2.6216636) q[3];
sx q[3];
rz(-2.9779696) q[3];
sx q[3];
rz(-0.61421052) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
