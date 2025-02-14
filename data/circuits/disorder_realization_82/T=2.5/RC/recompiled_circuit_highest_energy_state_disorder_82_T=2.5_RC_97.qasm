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
rz(-0.43871969) q[0];
sx q[0];
rz(-2.5320142) q[0];
sx q[0];
rz(0.69019812) q[0];
rz(-1.8742427) q[1];
sx q[1];
rz(6.7725362) q[1];
sx q[1];
rz(9.0803774) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9716838) q[0];
sx q[0];
rz(-0.2580041) q[0];
sx q[0];
rz(2.1841316) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3252668) q[2];
sx q[2];
rz(-1.9787162) q[2];
sx q[2];
rz(1.7657042) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92235451) q[1];
sx q[1];
rz(-1.7686756) q[1];
sx q[1];
rz(-1.8333652) q[1];
rz(1.0788513) q[3];
sx q[3];
rz(-1.4189093) q[3];
sx q[3];
rz(1.9263751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0979746) q[2];
sx q[2];
rz(-2.2100885) q[2];
sx q[2];
rz(2.0841058) q[2];
rz(2.1438694) q[3];
sx q[3];
rz(-1.3910553) q[3];
sx q[3];
rz(-1.1725496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69720307) q[0];
sx q[0];
rz(-2.3219705) q[0];
sx q[0];
rz(2.3890553) q[0];
rz(1.2731816) q[1];
sx q[1];
rz(-1.9877142) q[1];
sx q[1];
rz(-2.2862327) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6763362) q[0];
sx q[0];
rz(-1.6688866) q[0];
sx q[0];
rz(1.7855745) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1925621) q[2];
sx q[2];
rz(-1.3675642) q[2];
sx q[2];
rz(-2.0469613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7400186) q[1];
sx q[1];
rz(-0.42045004) q[1];
sx q[1];
rz(0.40989532) q[1];
rz(3.0166808) q[3];
sx q[3];
rz(-0.20044573) q[3];
sx q[3];
rz(0.87504609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7684043) q[2];
sx q[2];
rz(-1.5311925) q[2];
sx q[2];
rz(-1.9739523) q[2];
rz(3.043637) q[3];
sx q[3];
rz(-2.8601213) q[3];
sx q[3];
rz(1.9907985) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34404594) q[0];
sx q[0];
rz(-0.56586376) q[0];
sx q[0];
rz(1.4669482) q[0];
rz(0.037874669) q[1];
sx q[1];
rz(-1.2118309) q[1];
sx q[1];
rz(2.0272592) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6303725) q[0];
sx q[0];
rz(-1.1603705) q[0];
sx q[0];
rz(-0.12271304) q[0];
rz(-pi) q[1];
rz(-0.75598209) q[2];
sx q[2];
rz(-2.4376483) q[2];
sx q[2];
rz(-2.9005425) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0689286) q[1];
sx q[1];
rz(-2.0004099) q[1];
sx q[1];
rz(1.0547897) q[1];
rz(2.15184) q[3];
sx q[3];
rz(-1.8268181) q[3];
sx q[3];
rz(2.6953002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7797839) q[2];
sx q[2];
rz(-0.59541687) q[2];
sx q[2];
rz(2.6105866) q[2];
rz(-0.94046721) q[3];
sx q[3];
rz(-2.2898424) q[3];
sx q[3];
rz(1.6865591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4130212) q[0];
sx q[0];
rz(-2.72609) q[0];
sx q[0];
rz(-0.79237932) q[0];
rz(-0.13262311) q[1];
sx q[1];
rz(-2.4515929) q[1];
sx q[1];
rz(-2.2161868) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26811582) q[0];
sx q[0];
rz(-0.24304427) q[0];
sx q[0];
rz(-2.2683709) q[0];
rz(1.4508574) q[2];
sx q[2];
rz(-2.687157) q[2];
sx q[2];
rz(1.1392405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1455411) q[1];
sx q[1];
rz(-2.5153603) q[1];
sx q[1];
rz(-1.224451) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6033591) q[3];
sx q[3];
rz(-0.58659121) q[3];
sx q[3];
rz(-0.80249062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35599071) q[2];
sx q[2];
rz(-1.2987368) q[2];
sx q[2];
rz(-1.7388657) q[2];
rz(1.2518903) q[3];
sx q[3];
rz(-1.8391049) q[3];
sx q[3];
rz(2.8437264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4414325) q[0];
sx q[0];
rz(-2.1692363) q[0];
sx q[0];
rz(-2.1527619) q[0];
rz(-1.5160457) q[1];
sx q[1];
rz(-1.4715618) q[1];
sx q[1];
rz(0.54738799) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47212663) q[0];
sx q[0];
rz(-1.6030114) q[0];
sx q[0];
rz(-2.6526582) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21018942) q[2];
sx q[2];
rz(-0.37156216) q[2];
sx q[2];
rz(-0.7084058) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2965387) q[1];
sx q[1];
rz(-0.88857874) q[1];
sx q[1];
rz(-2.6107772) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8655254) q[3];
sx q[3];
rz(-1.4763586) q[3];
sx q[3];
rz(0.51968473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0872515) q[2];
sx q[2];
rz(-2.6332899) q[2];
sx q[2];
rz(0.15092078) q[2];
rz(-1.5058676) q[3];
sx q[3];
rz(-1.8412291) q[3];
sx q[3];
rz(1.6747564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48581377) q[0];
sx q[0];
rz(-1.1968311) q[0];
sx q[0];
rz(2.6659513) q[0];
rz(2.7173243) q[1];
sx q[1];
rz(-1.2356267) q[1];
sx q[1];
rz(-1.7686527) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17775336) q[0];
sx q[0];
rz(-2.8256157) q[0];
sx q[0];
rz(-2.55992) q[0];
rz(-2.0654997) q[2];
sx q[2];
rz(-1.037179) q[2];
sx q[2];
rz(2.1442206) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2306421) q[1];
sx q[1];
rz(-1.2394718) q[1];
sx q[1];
rz(2.1403007) q[1];
rz(-pi) q[2];
rz(-1.3282534) q[3];
sx q[3];
rz(-0.45363344) q[3];
sx q[3];
rz(-0.88874528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0515392) q[2];
sx q[2];
rz(-2.1681483) q[2];
sx q[2];
rz(0.092183979) q[2];
rz(1.9696382) q[3];
sx q[3];
rz(-0.53297526) q[3];
sx q[3];
rz(0.91641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5343269) q[0];
sx q[0];
rz(-2.3054275) q[0];
sx q[0];
rz(-1.0323866) q[0];
rz(-3.0119925) q[1];
sx q[1];
rz(-1.383129) q[1];
sx q[1];
rz(1.3547156) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2808896) q[0];
sx q[0];
rz(-0.96183813) q[0];
sx q[0];
rz(-0.97228284) q[0];
rz(-pi) q[1];
rz(-1.2298035) q[2];
sx q[2];
rz(-1.2330867) q[2];
sx q[2];
rz(2.4203398) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70020218) q[1];
sx q[1];
rz(-2.1245777) q[1];
sx q[1];
rz(-0.9662083) q[1];
rz(1.4192102) q[3];
sx q[3];
rz(-2.6612284) q[3];
sx q[3];
rz(2.9472873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20595343) q[2];
sx q[2];
rz(-2.8796068) q[2];
sx q[2];
rz(-1.4865173) q[2];
rz(0.67241159) q[3];
sx q[3];
rz(-1.0786062) q[3];
sx q[3];
rz(-3.0730754) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9392149) q[0];
sx q[0];
rz(-3.0164533) q[0];
sx q[0];
rz(2.8676721) q[0];
rz(-2.9778453) q[1];
sx q[1];
rz(-1.3122908) q[1];
sx q[1];
rz(-0.40072498) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.816283) q[0];
sx q[0];
rz(-1.0772711) q[0];
sx q[0];
rz(1.8614344) q[0];
x q[1];
rz(3.0178304) q[2];
sx q[2];
rz(-1.1699737) q[2];
sx q[2];
rz(1.7289897) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5536062) q[1];
sx q[1];
rz(-2.4692236) q[1];
sx q[1];
rz(1.5410627) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2473281) q[3];
sx q[3];
rz(-2.6163963) q[3];
sx q[3];
rz(2.0390455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87216941) q[2];
sx q[2];
rz(-2.9672406) q[2];
sx q[2];
rz(-1.6677469) q[2];
rz(-0.10146865) q[3];
sx q[3];
rz(-1.9188107) q[3];
sx q[3];
rz(-1.3946704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2673016) q[0];
sx q[0];
rz(-2.3912781) q[0];
sx q[0];
rz(3.1399723) q[0];
rz(-0.56456176) q[1];
sx q[1];
rz(-1.8058585) q[1];
sx q[1];
rz(0.50216215) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82715568) q[0];
sx q[0];
rz(-1.0183882) q[0];
sx q[0];
rz(-2.3236426) q[0];
rz(-1.9829168) q[2];
sx q[2];
rz(-1.0555763) q[2];
sx q[2];
rz(2.5784022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3920711) q[1];
sx q[1];
rz(-1.3765613) q[1];
sx q[1];
rz(0.51699395) q[1];
rz(2.929964) q[3];
sx q[3];
rz(-2.0787079) q[3];
sx q[3];
rz(0.25903364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.042171176) q[2];
sx q[2];
rz(-1.1974988) q[2];
sx q[2];
rz(-1.2449123) q[2];
rz(-2.9391607) q[3];
sx q[3];
rz(-1.7499685) q[3];
sx q[3];
rz(-2.1421053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56364432) q[0];
sx q[0];
rz(-2.7751594) q[0];
sx q[0];
rz(1.4165437) q[0];
rz(-1.1997403) q[1];
sx q[1];
rz(-1.1734633) q[1];
sx q[1];
rz(0.27552584) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4352132) q[0];
sx q[0];
rz(-1.9322898) q[0];
sx q[0];
rz(2.778591) q[0];
rz(-pi) q[1];
rz(-0.24015719) q[2];
sx q[2];
rz(-2.7809745) q[2];
sx q[2];
rz(2.220253) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4135189) q[1];
sx q[1];
rz(-0.46183837) q[1];
sx q[1];
rz(1.5240662) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4472792) q[3];
sx q[3];
rz(-2.165861) q[3];
sx q[3];
rz(1.7473011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3755017) q[2];
sx q[2];
rz(-1.4644863) q[2];
sx q[2];
rz(2.7601833) q[2];
rz(-0.31960791) q[3];
sx q[3];
rz(-0.8173129) q[3];
sx q[3];
rz(-2.4735425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9850591) q[0];
sx q[0];
rz(-1.1548797) q[0];
sx q[0];
rz(-0.67697939) q[0];
rz(2.1398687) q[1];
sx q[1];
rz(-1.6698508) q[1];
sx q[1];
rz(2.225266) q[1];
rz(0.66522404) q[2];
sx q[2];
rz(-1.6679674) q[2];
sx q[2];
rz(-1.3178772) q[2];
rz(0.14409625) q[3];
sx q[3];
rz(-2.477705) q[3];
sx q[3];
rz(-2.4841819) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
