OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1516079) q[0];
sx q[0];
rz(-0.94010544) q[0];
sx q[0];
rz(-2.6012233) q[0];
rz(-2.6481533) q[1];
sx q[1];
rz(-0.72055888) q[1];
sx q[1];
rz(0.61385733) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9073528) q[0];
sx q[0];
rz(-2.1932903) q[0];
sx q[0];
rz(1.3433775) q[0];
x q[1];
rz(-1.8625804) q[2];
sx q[2];
rz(-1.5651184) q[2];
sx q[2];
rz(1.8863854) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39345523) q[1];
sx q[1];
rz(-2.2917213) q[1];
sx q[1];
rz(1.9490521) q[1];
rz(0.085569445) q[3];
sx q[3];
rz(-2.8935379) q[3];
sx q[3];
rz(-0.55094592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2962239) q[2];
sx q[2];
rz(-1.8989398) q[2];
sx q[2];
rz(-0.63429147) q[2];
rz(-2.2058709) q[3];
sx q[3];
rz(-2.8959385) q[3];
sx q[3];
rz(0.74265695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.3212386) q[0];
sx q[0];
rz(-1.6004434) q[0];
sx q[0];
rz(-2.6974005) q[0];
rz(2.3815637) q[1];
sx q[1];
rz(-2.0030231) q[1];
sx q[1];
rz(-0.98145032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8151756) q[0];
sx q[0];
rz(-1.278864) q[0];
sx q[0];
rz(0.53262226) q[0];
x q[1];
rz(-1.3992127) q[2];
sx q[2];
rz(-1.7497471) q[2];
sx q[2];
rz(1.6323665) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5983898) q[1];
sx q[1];
rz(-2.0080505) q[1];
sx q[1];
rz(1.6679126) q[1];
rz(-pi) q[2];
rz(1.8360102) q[3];
sx q[3];
rz(-1.1516368) q[3];
sx q[3];
rz(-2.9492299) q[3];
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
rz(1.4046148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9179012) q[0];
sx q[0];
rz(-1.8603928) q[0];
sx q[0];
rz(0.83876383) q[0];
rz(0.63610786) q[1];
sx q[1];
rz(-1.4898224) q[1];
sx q[1];
rz(0.020523358) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9560709) q[0];
sx q[0];
rz(-0.98617109) q[0];
sx q[0];
rz(-2.4654075) q[0];
rz(-pi) q[1];
rz(-0.62600796) q[2];
sx q[2];
rz(-1.9591718) q[2];
sx q[2];
rz(0.0029390732) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.95076001) q[1];
sx q[1];
rz(-1.5783219) q[1];
sx q[1];
rz(-0.41295596) q[1];
rz(-pi) q[2];
rz(1.2853363) q[3];
sx q[3];
rz(-1.721964) q[3];
sx q[3];
rz(-2.5976439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8013578) q[2];
sx q[2];
rz(-1.5284208) q[2];
sx q[2];
rz(2.197263) q[2];
rz(1.0229735) q[3];
sx q[3];
rz(-1.9187656) q[3];
sx q[3];
rz(2.003722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.870938) q[0];
sx q[0];
rz(-1.174467) q[0];
sx q[0];
rz(-1.1573855) q[0];
rz(-1.4986787) q[1];
sx q[1];
rz(-1.6694371) q[1];
sx q[1];
rz(-1.1882943) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4469915) q[0];
sx q[0];
rz(-1.1199391) q[0];
sx q[0];
rz(-0.53038181) q[0];
x q[1];
rz(-3.068014) q[2];
sx q[2];
rz(-2.5917604) q[2];
sx q[2];
rz(2.4940707) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3664361) q[1];
sx q[1];
rz(-0.85298733) q[1];
sx q[1];
rz(-2.3597673) q[1];
x q[2];
rz(1.9096987) q[3];
sx q[3];
rz(-1.4168242) q[3];
sx q[3];
rz(-0.13974259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0231861) q[2];
sx q[2];
rz(-1.964317) q[2];
sx q[2];
rz(0.6905306) q[2];
rz(-2.0265419) q[3];
sx q[3];
rz(-0.77459049) q[3];
sx q[3];
rz(1.640865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0951776) q[0];
sx q[0];
rz(-1.4354118) q[0];
sx q[0];
rz(-1.0272367) q[0];
rz(1.8401624) q[1];
sx q[1];
rz(-2.4899028) q[1];
sx q[1];
rz(0.32381907) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9607317) q[0];
sx q[0];
rz(-1.9692076) q[0];
sx q[0];
rz(0.84503998) q[0];
rz(-2.1236046) q[2];
sx q[2];
rz(-1.5587285) q[2];
sx q[2];
rz(0.34876212) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4869838) q[1];
sx q[1];
rz(-1.2506079) q[1];
sx q[1];
rz(0.57211188) q[1];
rz(-pi) q[2];
rz(2.9848954) q[3];
sx q[3];
rz(-1.0644056) q[3];
sx q[3];
rz(1.2676257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41574898) q[2];
sx q[2];
rz(-0.21012935) q[2];
sx q[2];
rz(0.48247639) q[2];
rz(-1.1680565) q[3];
sx q[3];
rz(-1.9246293) q[3];
sx q[3];
rz(0.67679685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2060858) q[0];
sx q[0];
rz(-0.85010234) q[0];
sx q[0];
rz(-0.8859984) q[0];
rz(-0.63367263) q[1];
sx q[1];
rz(-0.88992563) q[1];
sx q[1];
rz(2.450313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9871498) q[0];
sx q[0];
rz(-1.5342752) q[0];
sx q[0];
rz(0.64572389) q[0];
rz(-pi) q[1];
rz(-3.0646661) q[2];
sx q[2];
rz(-2.4483213) q[2];
sx q[2];
rz(1.6037387) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.18610854) q[1];
sx q[1];
rz(-2.3209474) q[1];
sx q[1];
rz(-0.11891831) q[1];
x q[2];
rz(0.25005682) q[3];
sx q[3];
rz(-0.80509201) q[3];
sx q[3];
rz(0.4522194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1401691) q[2];
sx q[2];
rz(-2.3052577) q[2];
sx q[2];
rz(2.8720065) q[2];
rz(0.94830281) q[3];
sx q[3];
rz(-1.6092665) q[3];
sx q[3];
rz(1.3986826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36169323) q[0];
sx q[0];
rz(-2.7044856) q[0];
sx q[0];
rz(1.0070739) q[0];
rz(2.7359447) q[1];
sx q[1];
rz(-0.59527731) q[1];
sx q[1];
rz(-0.55353037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8419475) q[0];
sx q[0];
rz(-1.3163693) q[0];
sx q[0];
rz(3.0823452) q[0];
x q[1];
rz(2.361176) q[2];
sx q[2];
rz(-2.4274106) q[2];
sx q[2];
rz(-1.2250021) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.036052536) q[1];
sx q[1];
rz(-0.80784384) q[1];
sx q[1];
rz(2.4921662) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38000472) q[3];
sx q[3];
rz(-1.5475353) q[3];
sx q[3];
rz(-1.301736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8290528) q[2];
sx q[2];
rz(-2.0487787) q[2];
sx q[2];
rz(-1.9276169) q[2];
rz(0.78643262) q[3];
sx q[3];
rz(-1.7460456) q[3];
sx q[3];
rz(0.023155183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19602747) q[0];
sx q[0];
rz(-2.8493311) q[0];
sx q[0];
rz(-1.3319525) q[0];
rz(-0.61344433) q[1];
sx q[1];
rz(-2.1211801) q[1];
sx q[1];
rz(-0.44949284) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7004823) q[0];
sx q[0];
rz(-1.6539409) q[0];
sx q[0];
rz(-2.0594199) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.055213) q[2];
sx q[2];
rz(-1.9055718) q[2];
sx q[2];
rz(0.93418649) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7498988) q[1];
sx q[1];
rz(-0.79872455) q[1];
sx q[1];
rz(-2.9500089) q[1];
rz(0.21829101) q[3];
sx q[3];
rz(-2.8779753) q[3];
sx q[3];
rz(-1.5052312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31676644) q[2];
sx q[2];
rz(-2.4591441) q[2];
sx q[2];
rz(-2.5332434) q[2];
rz(-1.1634722) q[3];
sx q[3];
rz(-1.3694265) q[3];
sx q[3];
rz(-0.56330645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1083531) q[0];
sx q[0];
rz(-2.2305363) q[0];
sx q[0];
rz(-0.60229993) q[0];
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
rz(-2.8062144) q[0];
sx q[0];
rz(-0.69604448) q[0];
sx q[0];
rz(-0.93675128) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63656143) q[2];
sx q[2];
rz(-1.5590073) q[2];
sx q[2];
rz(1.7280886) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1608354) q[1];
sx q[1];
rz(-0.45113647) q[1];
sx q[1];
rz(-1.5734929) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6742485) q[3];
sx q[3];
rz(-2.561063) q[3];
sx q[3];
rz(-0.79642297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.935282) q[2];
sx q[2];
rz(-1.3056359) q[2];
sx q[2];
rz(-2.804011) q[2];
rz(-0.28389367) q[3];
sx q[3];
rz(-2.4604535) q[3];
sx q[3];
rz(-2.9547227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5481446) q[0];
sx q[0];
rz(-1.633506) q[0];
sx q[0];
rz(-0.067807587) q[0];
rz(1.1495122) q[1];
sx q[1];
rz(-1.5905453) q[1];
sx q[1];
rz(-0.4745208) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10013469) q[0];
sx q[0];
rz(-0.23915072) q[0];
sx q[0];
rz(1.5140223) q[0];
x q[1];
rz(-1.7065918) q[2];
sx q[2];
rz(-1.8095922) q[2];
sx q[2];
rz(0.39993024) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2429035) q[1];
sx q[1];
rz(-1.4840645) q[1];
sx q[1];
rz(3.134722) q[1];
x q[2];
rz(-2.2978333) q[3];
sx q[3];
rz(-0.90137989) q[3];
sx q[3];
rz(-2.9694174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48156753) q[2];
sx q[2];
rz(-0.93067545) q[2];
sx q[2];
rz(1.0732667) q[2];
rz(2.977071) q[3];
sx q[3];
rz(-1.8142895) q[3];
sx q[3];
rz(2.8209414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5545223) q[0];
sx q[0];
rz(-1.5656492) q[0];
sx q[0];
rz(1.5026305) q[0];
rz(1.5578237) q[1];
sx q[1];
rz(-2.0638034) q[1];
sx q[1];
rz(2.6149909) q[1];
rz(-1.9314968) q[2];
sx q[2];
rz(-2.0059235) q[2];
sx q[2];
rz(-0.37034482) q[2];
rz(-2.9992793) q[3];
sx q[3];
rz(-1.6518136) q[3];
sx q[3];
rz(1.470737) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
