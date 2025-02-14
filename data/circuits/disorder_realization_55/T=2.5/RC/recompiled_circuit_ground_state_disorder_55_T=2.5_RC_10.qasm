OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80630535) q[0];
sx q[0];
rz(4.181554) q[0];
sx q[0];
rz(9.7747533) q[0];
rz(1.5179874) q[1];
sx q[1];
rz(-2.1972158) q[1];
sx q[1];
rz(1.1854393) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1040504) q[0];
sx q[0];
rz(-0.54706956) q[0];
sx q[0];
rz(0.68010917) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41857206) q[2];
sx q[2];
rz(-1.3579733) q[2];
sx q[2];
rz(-0.71453071) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2054641) q[1];
sx q[1];
rz(-1.3359114) q[1];
sx q[1];
rz(-2.111819) q[1];
x q[2];
rz(-1.8270742) q[3];
sx q[3];
rz(-2.5102237) q[3];
sx q[3];
rz(2.6169962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27652937) q[2];
sx q[2];
rz(-0.41215602) q[2];
sx q[2];
rz(0.40974799) q[2];
rz(2.3661738) q[3];
sx q[3];
rz(-1.5284458) q[3];
sx q[3];
rz(0.30737901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85768098) q[0];
sx q[0];
rz(-2.6429521) q[0];
sx q[0];
rz(0.45251244) q[0];
rz(2.9127938) q[1];
sx q[1];
rz(-2.6401873) q[1];
sx q[1];
rz(-0.14678821) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89167021) q[0];
sx q[0];
rz(-1.2091214) q[0];
sx q[0];
rz(-2.937243) q[0];
x q[1];
rz(0.54988523) q[2];
sx q[2];
rz(-2.7061979) q[2];
sx q[2];
rz(-0.76744902) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70554206) q[1];
sx q[1];
rz(-1.5311993) q[1];
sx q[1];
rz(3.0071053) q[1];
x q[2];
rz(0.66007075) q[3];
sx q[3];
rz(-1.1997031) q[3];
sx q[3];
rz(2.0160272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9734834) q[2];
sx q[2];
rz(-2.8066469) q[2];
sx q[2];
rz(0.87146634) q[2];
rz(-0.35428366) q[3];
sx q[3];
rz(-0.93577093) q[3];
sx q[3];
rz(1.833545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63534921) q[0];
sx q[0];
rz(-0.62017089) q[0];
sx q[0];
rz(1.1059906) q[0];
rz(-2.1982819) q[1];
sx q[1];
rz(-2.3425075) q[1];
sx q[1];
rz(-1.6518637) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36216736) q[0];
sx q[0];
rz(-0.12813103) q[0];
sx q[0];
rz(1.6511028) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9331384) q[2];
sx q[2];
rz(-1.4770419) q[2];
sx q[2];
rz(-1.1445647) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.024002687) q[1];
sx q[1];
rz(-1.6086714) q[1];
sx q[1];
rz(0.59549673) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2977887) q[3];
sx q[3];
rz(-2.892916) q[3];
sx q[3];
rz(-0.77588785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6279383) q[2];
sx q[2];
rz(-1.1955806) q[2];
sx q[2];
rz(0.99620831) q[2];
rz(-0.57322383) q[3];
sx q[3];
rz(-2.6300391) q[3];
sx q[3];
rz(2.5073124) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298249) q[0];
sx q[0];
rz(-1.8695762) q[0];
sx q[0];
rz(1.5856532) q[0];
rz(2.8171483) q[1];
sx q[1];
rz(-2.3093846) q[1];
sx q[1];
rz(1.197804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8666337) q[0];
sx q[0];
rz(-0.64104762) q[0];
sx q[0];
rz(-0.43305555) q[0];
x q[1];
rz(-1.7928855) q[2];
sx q[2];
rz(-2.54455) q[2];
sx q[2];
rz(0.82840234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50655925) q[1];
sx q[1];
rz(-2.6406619) q[1];
sx q[1];
rz(-1.2975724) q[1];
rz(-pi) q[2];
rz(-2.1230202) q[3];
sx q[3];
rz(-2.8365144) q[3];
sx q[3];
rz(-0.48459241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.045362554) q[2];
sx q[2];
rz(-2.5924293) q[2];
sx q[2];
rz(1.9722923) q[2];
rz(-0.62396389) q[3];
sx q[3];
rz(-1.4493161) q[3];
sx q[3];
rz(0.72871488) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28155926) q[0];
sx q[0];
rz(-0.39880729) q[0];
sx q[0];
rz(-1.9418203) q[0];
rz(1.8343605) q[1];
sx q[1];
rz(-2.2016826) q[1];
sx q[1];
rz(1.7050381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4685637) q[0];
sx q[0];
rz(-1.9583681) q[0];
sx q[0];
rz(3.0762071) q[0];
rz(0.13135037) q[2];
sx q[2];
rz(-1.3418014) q[2];
sx q[2];
rz(0.74818767) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87203997) q[1];
sx q[1];
rz(-1.8694573) q[1];
sx q[1];
rz(2.7725817) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1025401) q[3];
sx q[3];
rz(-2.4210586) q[3];
sx q[3];
rz(1.6417491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3522219) q[2];
sx q[2];
rz(-1.7033966) q[2];
sx q[2];
rz(0.5729202) q[2];
rz(-0.88417435) q[3];
sx q[3];
rz(-2.1638162) q[3];
sx q[3];
rz(2.6463215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5578515) q[0];
sx q[0];
rz(-1.252625) q[0];
sx q[0];
rz(-2.0649233) q[0];
rz(2.63511) q[1];
sx q[1];
rz(-2.5264637) q[1];
sx q[1];
rz(-1.3004251) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5426555) q[0];
sx q[0];
rz(-3.0910203) q[0];
sx q[0];
rz(-1.845832) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7188273) q[2];
sx q[2];
rz(-0.99642053) q[2];
sx q[2];
rz(0.21166052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4757931) q[1];
sx q[1];
rz(-0.25036004) q[1];
sx q[1];
rz(-1.0151035) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9229284) q[3];
sx q[3];
rz(-1.3826958) q[3];
sx q[3];
rz(-1.355483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41636813) q[2];
sx q[2];
rz(-2.29795) q[2];
sx q[2];
rz(-1.9492524) q[2];
rz(0.013817712) q[3];
sx q[3];
rz(-2.4255224) q[3];
sx q[3];
rz(-2.1684087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8927807) q[0];
sx q[0];
rz(-2.2154494) q[0];
sx q[0];
rz(-0.37175757) q[0];
rz(-2.2573076) q[1];
sx q[1];
rz(-0.50768948) q[1];
sx q[1];
rz(-2.5977792) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0293074) q[0];
sx q[0];
rz(-2.4196631) q[0];
sx q[0];
rz(0.63755905) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85141723) q[2];
sx q[2];
rz(-1.9236375) q[2];
sx q[2];
rz(-2.1746496) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8616756) q[1];
sx q[1];
rz(-1.1959198) q[1];
sx q[1];
rz(1.5779297) q[1];
rz(-pi) q[2];
rz(-1.722808) q[3];
sx q[3];
rz(-1.4626039) q[3];
sx q[3];
rz(-0.044807981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1112993) q[2];
sx q[2];
rz(-2.0798422) q[2];
sx q[2];
rz(-0.57478762) q[2];
rz(0.8542257) q[3];
sx q[3];
rz(-1.6653019) q[3];
sx q[3];
rz(2.0265719) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92693555) q[0];
sx q[0];
rz(-0.68100005) q[0];
sx q[0];
rz(2.8171203) q[0];
rz(-2.532161) q[1];
sx q[1];
rz(-2.29988) q[1];
sx q[1];
rz(0.51469523) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97159779) q[0];
sx q[0];
rz(-0.43525243) q[0];
sx q[0];
rz(-1.5488141) q[0];
rz(-pi) q[1];
rz(-2.2836907) q[2];
sx q[2];
rz(-1.2944572) q[2];
sx q[2];
rz(-0.954962) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.534651) q[1];
sx q[1];
rz(-1.80503) q[1];
sx q[1];
rz(-0.57778426) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10917337) q[3];
sx q[3];
rz(-1.2965373) q[3];
sx q[3];
rz(1.4281987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1631761) q[2];
sx q[2];
rz(-1.6360444) q[2];
sx q[2];
rz(0.27463883) q[2];
rz(-2.2091673) q[3];
sx q[3];
rz(-2.7115188) q[3];
sx q[3];
rz(-2.8453804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4579929) q[0];
sx q[0];
rz(-1.4608811) q[0];
sx q[0];
rz(-3.006111) q[0];
rz(-2.1645891) q[1];
sx q[1];
rz(-2.3833387) q[1];
sx q[1];
rz(0.0010842222) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1019395) q[0];
sx q[0];
rz(-0.51991594) q[0];
sx q[0];
rz(-2.824409) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3668225) q[2];
sx q[2];
rz(-0.53278953) q[2];
sx q[2];
rz(-2.8486657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26027706) q[1];
sx q[1];
rz(-2.1268357) q[1];
sx q[1];
rz(2.4350321) q[1];
rz(-0.75653177) q[3];
sx q[3];
rz(-2.0755141) q[3];
sx q[3];
rz(-2.0677419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.5929426) q[2];
sx q[2];
rz(-2.103001) q[2];
sx q[2];
rz(2.9366117) q[2];
rz(2.958988) q[3];
sx q[3];
rz(-2.9349116) q[3];
sx q[3];
rz(-2.4038147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7788727) q[0];
sx q[0];
rz(-0.39411476) q[0];
sx q[0];
rz(-0.47738281) q[0];
rz(-0.1167156) q[1];
sx q[1];
rz(-1.5201897) q[1];
sx q[1];
rz(0.83888549) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29541839) q[0];
sx q[0];
rz(-1.8356165) q[0];
sx q[0];
rz(-2.2841262) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9220542) q[2];
sx q[2];
rz(-0.87703401) q[2];
sx q[2];
rz(2.0036245) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15167639) q[1];
sx q[1];
rz(-0.71571678) q[1];
sx q[1];
rz(-2.8552961) q[1];
rz(-pi) q[2];
rz(1.3319043) q[3];
sx q[3];
rz(-1.5419067) q[3];
sx q[3];
rz(-0.98945557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3597151) q[2];
sx q[2];
rz(-0.30656591) q[2];
sx q[2];
rz(0.26620418) q[2];
rz(-2.6340458) q[3];
sx q[3];
rz(-0.72269732) q[3];
sx q[3];
rz(-2.7005196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1626641) q[0];
sx q[0];
rz(-1.6829818) q[0];
sx q[0];
rz(-0.80243954) q[0];
rz(-2.2294527) q[1];
sx q[1];
rz(-1.9504539) q[1];
sx q[1];
rz(0.84074195) q[1];
rz(0.24813454) q[2];
sx q[2];
rz(-1.4149181) q[2];
sx q[2];
rz(0.96542035) q[2];
rz(-1.9515362) q[3];
sx q[3];
rz(-1.1968812) q[3];
sx q[3];
rz(1.4637917) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
