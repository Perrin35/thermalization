OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(0.16790976) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(-2.8462703) q[1];
sx q[1];
rz(0.056161031) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6457155) q[0];
sx q[0];
rz(-2.0464532) q[0];
sx q[0];
rz(-2.9021184) q[0];
rz(0.21284717) q[2];
sx q[2];
rz(-2.2058862) q[2];
sx q[2];
rz(1.1320621) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1311156) q[1];
sx q[1];
rz(-1.3033723) q[1];
sx q[1];
rz(-1.0130151) q[1];
rz(-pi) q[2];
rz(1.5264741) q[3];
sx q[3];
rz(-1.3111776) q[3];
sx q[3];
rz(-1.1322024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37796676) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(2.7089233) q[2];
rz(1.9487322) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6137961) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(1.8288076) q[0];
rz(2.9361172) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(-1.1516494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5702471) q[0];
sx q[0];
rz(-0.76128188) q[0];
sx q[0];
rz(-2.0061357) q[0];
rz(2.5580514) q[2];
sx q[2];
rz(-1.1569287) q[2];
sx q[2];
rz(-1.5577424) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9559905) q[1];
sx q[1];
rz(-0.62528175) q[1];
sx q[1];
rz(2.37466) q[1];
rz(-1.9563975) q[3];
sx q[3];
rz(-2.0413627) q[3];
sx q[3];
rz(3.1009931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(2.9197664) q[2];
rz(-2.7644073) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(-0.77243531) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8310228) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(3.1047399) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(3.085014) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73233561) q[0];
sx q[0];
rz(-1.1261254) q[0];
sx q[0];
rz(2.1952941) q[0];
x q[1];
rz(0.25643202) q[2];
sx q[2];
rz(-1.4415381) q[2];
sx q[2];
rz(-1.3916707) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0128855) q[1];
sx q[1];
rz(-1.371908) q[1];
sx q[1];
rz(2.835564) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49719663) q[3];
sx q[3];
rz(-0.70780863) q[3];
sx q[3];
rz(1.7266112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.27292192) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-0.92612129) q[2];
rz(-0.55666322) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(-2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2264003) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-2.3994989) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-0.46359584) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.050042) q[0];
sx q[0];
rz(-2.3766962) q[0];
sx q[0];
rz(-2.0043623) q[0];
x q[1];
rz(0.33275231) q[2];
sx q[2];
rz(-1.6289662) q[2];
sx q[2];
rz(-1.0156877) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3568748) q[1];
sx q[1];
rz(-0.80662913) q[1];
sx q[1];
rz(-2.9033317) q[1];
rz(-pi) q[2];
rz(-1.0158402) q[3];
sx q[3];
rz(-1.98588) q[3];
sx q[3];
rz(-0.89158981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7745557) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(0.049499361) q[2];
rz(0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0054935) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(2.8438925) q[0];
rz(-0.4822576) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(2.1972426) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0194861) q[0];
sx q[0];
rz(-1.6158316) q[0];
sx q[0];
rz(1.3615863) q[0];
rz(-0.77563939) q[2];
sx q[2];
rz(-1.8619985) q[2];
sx q[2];
rz(1.8679384) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7625092) q[1];
sx q[1];
rz(-0.06113872) q[1];
sx q[1];
rz(1.5361384) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3145507) q[3];
sx q[3];
rz(-1.991193) q[3];
sx q[3];
rz(2.8375569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8828316) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(-0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(0.32430696) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7047983) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-2.5571402) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(3.086673) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14558218) q[0];
sx q[0];
rz(-1.8390363) q[0];
sx q[0];
rz(3.0560188) q[0];
x q[1];
rz(1.247585) q[2];
sx q[2];
rz(-1.5885457) q[2];
sx q[2];
rz(-0.36349597) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.37327787) q[1];
sx q[1];
rz(-2.1149966) q[1];
sx q[1];
rz(2.0139704) q[1];
rz(0.46338007) q[3];
sx q[3];
rz(-1.0372835) q[3];
sx q[3];
rz(2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(1.0167271) q[2];
rz(-2.5975442) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(-1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(-2.8570535) q[0];
rz(0.94447213) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(-2.231266) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1658265) q[0];
sx q[0];
rz(-1.5162139) q[0];
sx q[0];
rz(1.635701) q[0];
rz(-pi) q[1];
rz(-1.1548642) q[2];
sx q[2];
rz(-1.2565194) q[2];
sx q[2];
rz(-1.2736125) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0525166) q[1];
sx q[1];
rz(-2.0564582) q[1];
sx q[1];
rz(-1.3658701) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4207553) q[3];
sx q[3];
rz(-1.1158873) q[3];
sx q[3];
rz(2.7618559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7515144) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(0.92203036) q[2];
rz(-0.56728029) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(-2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89408016) q[0];
sx q[0];
rz(-0.67665726) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(-0.63240504) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(2.8410889) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104615) q[0];
sx q[0];
rz(-2.3346402) q[0];
sx q[0];
rz(2.0956844) q[0];
x q[1];
rz(2.5007162) q[2];
sx q[2];
rz(-0.87101988) q[2];
sx q[2];
rz(-1.7388294) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.064425163) q[1];
sx q[1];
rz(-2.7556813) q[1];
sx q[1];
rz(-0.47972958) q[1];
x q[2];
rz(-0.31605966) q[3];
sx q[3];
rz(-2.3028767) q[3];
sx q[3];
rz(-0.72533208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5552716) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(-0.55109763) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5732116) q[0];
sx q[0];
rz(-1.9848354) q[0];
sx q[0];
rz(-3.0138299) q[0];
rz(-0.54221517) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(-0.75884563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7153873) q[0];
sx q[0];
rz(-1.9949159) q[0];
sx q[0];
rz(-3.12294) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60894572) q[2];
sx q[2];
rz(-1.7296089) q[2];
sx q[2];
rz(1.9132523) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.54770494) q[1];
sx q[1];
rz(-2.2624359) q[1];
sx q[1];
rz(-1.3681075) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5893448) q[3];
sx q[3];
rz(-2.4759001) q[3];
sx q[3];
rz(-1.1902283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(-2.6515567) q[2];
rz(1.7193433) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(-2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35995099) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-0.075335659) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(0.60992253) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9655351) q[0];
sx q[0];
rz(-0.81747222) q[0];
sx q[0];
rz(2.7479991) q[0];
rz(-pi) q[1];
rz(-0.37656017) q[2];
sx q[2];
rz(-1.3137523) q[2];
sx q[2];
rz(0.10274796) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8049106) q[1];
sx q[1];
rz(-1.8098127) q[1];
sx q[1];
rz(0.75093021) q[1];
rz(-2.0718594) q[3];
sx q[3];
rz(-2.2483629) q[3];
sx q[3];
rz(2.5354405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.909409) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(-2.4278736) q[2];
rz(-0.37832007) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338035) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(0.65080416) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(-1.4177633) q[2];
sx q[2];
rz(-0.19822181) q[2];
sx q[2];
rz(-0.76186686) q[2];
rz(2.0402504) q[3];
sx q[3];
rz(-2.6228842) q[3];
sx q[3];
rz(1.7026671) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
