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
rz(-2.5277353) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9073528) q[0];
sx q[0];
rz(-0.94830238) q[0];
sx q[0];
rz(1.3433775) q[0];
rz(-pi) q[1];
rz(-1.5905321) q[2];
sx q[2];
rz(-0.29183772) q[2];
sx q[2];
rz(-2.8449051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.93537027) q[1];
sx q[1];
rz(-2.3435278) q[1];
sx q[1];
rz(0.3978637) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5491539) q[3];
sx q[3];
rz(-1.3236681) q[3];
sx q[3];
rz(-0.46268845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2962239) q[2];
sx q[2];
rz(-1.2426528) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82035404) q[0];
sx q[0];
rz(-1.5411493) q[0];
sx q[0];
rz(-2.6974005) q[0];
rz(0.76002899) q[1];
sx q[1];
rz(-2.0030231) q[1];
sx q[1];
rz(0.98145032) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0652577) q[0];
sx q[0];
rz(-2.0786571) q[0];
sx q[0];
rz(1.2351551) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18155988) q[2];
sx q[2];
rz(-1.4019792) q[2];
sx q[2];
rz(0.09240514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.54320286) q[1];
sx q[1];
rz(-1.1335422) q[1];
sx q[1];
rz(1.6679126) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4325474) q[3];
sx q[3];
rz(-1.8125696) q[3];
sx q[3];
rz(1.4885308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0280219) q[2];
sx q[2];
rz(-0.61442033) q[2];
sx q[2];
rz(1.2051955) q[2];
rz(-2.6049854) q[3];
sx q[3];
rz(-1.2380995) q[3];
sx q[3];
rz(-1.4046148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2236915) q[0];
sx q[0];
rz(-1.2811998) q[0];
sx q[0];
rz(-0.83876383) q[0];
rz(-2.5054848) q[1];
sx q[1];
rz(-1.4898224) q[1];
sx q[1];
rz(-3.1210693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9560709) q[0];
sx q[0];
rz(-2.1554216) q[0];
sx q[0];
rz(0.67618518) q[0];
rz(-pi) q[1];
rz(-2.0383561) q[2];
sx q[2];
rz(-0.9976495) q[2];
sx q[2];
rz(1.3064885) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1908326) q[1];
sx q[1];
rz(-1.5783219) q[1];
sx q[1];
rz(2.7286367) q[1];
rz(-pi) q[2];
rz(2.9841524) q[3];
sx q[3];
rz(-1.2886815) q[3];
sx q[3];
rz(-1.071014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8013578) q[2];
sx q[2];
rz(-1.6131718) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.2706547) q[0];
sx q[0];
rz(-1.174467) q[0];
sx q[0];
rz(-1.1573855) q[0];
rz(1.6429139) q[1];
sx q[1];
rz(-1.4721556) q[1];
sx q[1];
rz(1.1882943) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.37876) q[0];
sx q[0];
rz(-0.68183696) q[0];
sx q[0];
rz(-2.3781611) q[0];
rz(-pi) q[1];
rz(1.525773) q[2];
sx q[2];
rz(-2.1189711) q[2];
sx q[2];
rz(2.5803103) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.374158) q[1];
sx q[1];
rz(-1.0113455) q[1];
sx q[1];
rz(-2.4591695) q[1];
rz(-pi) q[2];
rz(-2.0075624) q[3];
sx q[3];
rz(-2.7705857) q[3];
sx q[3];
rz(-1.3001668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1184065) q[2];
sx q[2];
rz(-1.1772757) q[2];
sx q[2];
rz(-2.4510621) q[2];
rz(1.1150507) q[3];
sx q[3];
rz(-0.77459049) q[3];
sx q[3];
rz(1.640865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0951776) q[0];
sx q[0];
rz(-1.4354118) q[0];
sx q[0];
rz(-1.0272367) q[0];
rz(1.3014303) q[1];
sx q[1];
rz(-2.4899028) q[1];
sx q[1];
rz(-0.32381907) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80211783) q[0];
sx q[0];
rz(-2.3314752) q[0];
sx q[0];
rz(-1.0056061) q[0];
x q[1];
rz(-1.5937763) q[2];
sx q[2];
rz(-0.55292623) q[2];
sx q[2];
rz(-1.241591) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4253772) q[1];
sx q[1];
rz(-1.0310804) q[1];
sx q[1];
rz(-1.9464689) q[1];
rz(-pi) q[2];
rz(-2.9848954) q[3];
sx q[3];
rz(-2.0771871) q[3];
sx q[3];
rz(1.2676257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7258437) q[2];
sx q[2];
rz(-0.21012935) q[2];
sx q[2];
rz(-2.6591163) q[2];
rz(1.9735362) q[3];
sx q[3];
rz(-1.2169633) q[3];
sx q[3];
rz(-0.67679685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93550682) q[0];
sx q[0];
rz(-2.2914903) q[0];
sx q[0];
rz(-0.8859984) q[0];
rz(-0.63367263) q[1];
sx q[1];
rz(-2.251667) q[1];
sx q[1];
rz(-2.450313) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9871498) q[0];
sx q[0];
rz(-1.5342752) q[0];
sx q[0];
rz(2.4958688) q[0];
rz(-pi) q[1];
rz(0.076926546) q[2];
sx q[2];
rz(-0.69327136) q[2];
sx q[2];
rz(1.537854) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1289802) q[1];
sx q[1];
rz(-0.75769934) q[1];
sx q[1];
rz(-1.6974259) q[1];
x q[2];
rz(0.25005682) q[3];
sx q[3];
rz(-2.3365006) q[3];
sx q[3];
rz(2.6893733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0014235) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36169323) q[0];
sx q[0];
rz(-2.7044856) q[0];
sx q[0];
rz(2.1345188) q[0];
rz(-0.40564793) q[1];
sx q[1];
rz(-0.59527731) q[1];
sx q[1];
rz(2.5880623) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2562228) q[0];
sx q[0];
rz(-1.6281343) q[0];
sx q[0];
rz(-1.3159412) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1184475) q[2];
sx q[2];
rz(-2.0549462) q[2];
sx q[2];
rz(2.1439056) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0515159) q[1];
sx q[1];
rz(-1.1184268) q[1];
sx q[1];
rz(0.69454792) q[1];
rz(-1.5457492) q[3];
sx q[3];
rz(-1.1908997) q[3];
sx q[3];
rz(0.27835007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8290528) q[2];
sx q[2];
rz(-1.092814) q[2];
sx q[2];
rz(-1.2139758) q[2];
rz(0.78643262) q[3];
sx q[3];
rz(-1.395547) q[3];
sx q[3];
rz(-0.023155183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19602747) q[0];
sx q[0];
rz(-0.29226154) q[0];
sx q[0];
rz(-1.3319525) q[0];
rz(2.5281483) q[1];
sx q[1];
rz(-1.0204126) q[1];
sx q[1];
rz(-2.6920998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7004823) q[0];
sx q[0];
rz(-1.6539409) q[0];
sx q[0];
rz(-1.0821728) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.055213) q[2];
sx q[2];
rz(-1.2360209) q[2];
sx q[2];
rz(2.2074062) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39169381) q[1];
sx q[1];
rz(-2.3428681) q[1];
sx q[1];
rz(-0.19158371) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2576377) q[3];
sx q[3];
rz(-1.6272568) q[3];
sx q[3];
rz(0.27654058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31676644) q[2];
sx q[2];
rz(-2.4591441) q[2];
sx q[2];
rz(-2.5332434) q[2];
rz(-1.1634722) q[3];
sx q[3];
rz(-1.3694265) q[3];
sx q[3];
rz(2.5782862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.1083531) q[0];
sx q[0];
rz(-0.9110564) q[0];
sx q[0];
rz(-2.5392927) q[0];
rz(-2.1741518) q[1];
sx q[1];
rz(-0.93092218) q[1];
sx q[1];
rz(-1.1036576) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7491584) q[0];
sx q[0];
rz(-1.1811678) q[0];
sx q[0];
rz(-0.97831877) q[0];
x q[1];
rz(3.121762) q[2];
sx q[2];
rz(-0.63665542) q[2];
sx q[2];
rz(-3.0002468) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.98075726) q[1];
sx q[1];
rz(-0.45113647) q[1];
sx q[1];
rz(-1.5680997) q[1];
rz(-pi) q[2];
rz(1.6742485) q[3];
sx q[3];
rz(-0.58052968) q[3];
sx q[3];
rz(2.3451697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2063107) q[2];
sx q[2];
rz(-1.8359567) q[2];
sx q[2];
rz(-2.804011) q[2];
rz(-0.28389367) q[3];
sx q[3];
rz(-0.68113911) q[3];
sx q[3];
rz(2.9547227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5934481) q[0];
sx q[0];
rz(-1.5080867) q[0];
sx q[0];
rz(-3.0737851) q[0];
rz(-1.1495122) q[1];
sx q[1];
rz(-1.5905453) q[1];
sx q[1];
rz(-2.6670719) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.041458) q[0];
sx q[0];
rz(-0.23915072) q[0];
sx q[0];
rz(-1.6275703) q[0];
x q[1];
rz(-1.7065918) q[2];
sx q[2];
rz(-1.3320005) q[2];
sx q[2];
rz(-0.39993024) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2429035) q[1];
sx q[1];
rz(-1.6575282) q[1];
sx q[1];
rz(-3.134722) q[1];
rz(0.84375937) q[3];
sx q[3];
rz(-2.2402128) q[3];
sx q[3];
rz(2.9694174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.48156753) q[2];
sx q[2];
rz(-0.93067545) q[2];
sx q[2];
rz(-1.0732667) q[2];
rz(-2.977071) q[3];
sx q[3];
rz(-1.8142895) q[3];
sx q[3];
rz(0.32065121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58707033) q[0];
sx q[0];
rz(-1.5759435) q[0];
sx q[0];
rz(-1.6389621) q[0];
rz(1.5837689) q[1];
sx q[1];
rz(-1.0777892) q[1];
sx q[1];
rz(-0.52660175) q[1];
rz(0.64939349) q[2];
sx q[2];
rz(-2.5838701) q[2];
sx q[2];
rz(2.0412847) q[2];
rz(2.6216636) q[3];
sx q[3];
rz(-0.16362301) q[3];
sx q[3];
rz(2.5273821) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
