OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.13088432) q[0];
sx q[0];
rz(-0.67357981) q[0];
sx q[0];
rz(-0.82290617) q[0];
rz(1.4505439) q[1];
sx q[1];
rz(2.4658642) q[1];
sx q[1];
rz(9.2171062) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7340644) q[0];
sx q[0];
rz(-2.3210566) q[0];
sx q[0];
rz(-2.0506917) q[0];
x q[1];
rz(-1.0124341) q[2];
sx q[2];
rz(-1.867138) q[2];
sx q[2];
rz(1.3989965) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77601469) q[1];
sx q[1];
rz(-1.509359) q[1];
sx q[1];
rz(-0.12428026) q[1];
rz(-pi) q[2];
x q[2];
rz(0.071524044) q[3];
sx q[3];
rz(-1.9404963) q[3];
sx q[3];
rz(-2.2184851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1755918) q[2];
sx q[2];
rz(-1.5755743) q[2];
sx q[2];
rz(-1.9123745) q[2];
rz(1.843533) q[3];
sx q[3];
rz(-1.7652054) q[3];
sx q[3];
rz(1.1840597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9640279) q[0];
sx q[0];
rz(-2.8819045) q[0];
sx q[0];
rz(2.2143256) q[0];
rz(0.19993965) q[1];
sx q[1];
rz(-1.4727458) q[1];
sx q[1];
rz(-2.1520069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7829076) q[0];
sx q[0];
rz(-1.9316422) q[0];
sx q[0];
rz(0.11519687) q[0];
rz(-3.0485905) q[2];
sx q[2];
rz(-2.7055397) q[2];
sx q[2];
rz(-2.2047037) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.67827077) q[1];
sx q[1];
rz(-3.0509147) q[1];
sx q[1];
rz(2.1342127) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58122509) q[3];
sx q[3];
rz(-0.94454256) q[3];
sx q[3];
rz(-0.083409781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9280615) q[2];
sx q[2];
rz(-1.7240883) q[2];
sx q[2];
rz(-2.7878063) q[2];
rz(0.95156041) q[3];
sx q[3];
rz(-2.3944201) q[3];
sx q[3];
rz(2.814754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78615171) q[0];
sx q[0];
rz(-2.1936301) q[0];
sx q[0];
rz(2.1308664) q[0];
rz(-0.058723681) q[1];
sx q[1];
rz(-1.6207691) q[1];
sx q[1];
rz(-0.40930632) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8008645) q[0];
sx q[0];
rz(-2.0052064) q[0];
sx q[0];
rz(0.77462642) q[0];
rz(0.1654314) q[2];
sx q[2];
rz(-1.1647875) q[2];
sx q[2];
rz(0.55961404) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5262437) q[1];
sx q[1];
rz(-1.0305291) q[1];
sx q[1];
rz(1.8002285) q[1];
x q[2];
rz(-2.6326551) q[3];
sx q[3];
rz(-1.5889611) q[3];
sx q[3];
rz(-0.048232676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.056444082) q[2];
sx q[2];
rz(-1.8935545) q[2];
sx q[2];
rz(2.9819581) q[2];
rz(1.2498445) q[3];
sx q[3];
rz(-1.7227453) q[3];
sx q[3];
rz(-2.0554525) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98642629) q[0];
sx q[0];
rz(-2.7041589) q[0];
sx q[0];
rz(1.0945818) q[0];
rz(-1.9704341) q[1];
sx q[1];
rz(-2.0234225) q[1];
sx q[1];
rz(1.0135244) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0243624) q[0];
sx q[0];
rz(-1.6147531) q[0];
sx q[0];
rz(1.5323601) q[0];
rz(-pi) q[1];
rz(3.1379478) q[2];
sx q[2];
rz(-2.0652152) q[2];
sx q[2];
rz(-2.078675) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5202433) q[1];
sx q[1];
rz(-1.831437) q[1];
sx q[1];
rz(-1.6127178) q[1];
rz(-pi) q[2];
rz(-1.7988241) q[3];
sx q[3];
rz(-2.6963391) q[3];
sx q[3];
rz(-2.5955615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1363498) q[2];
sx q[2];
rz(-1.3815657) q[2];
sx q[2];
rz(-1.3137777) q[2];
rz(0.93787307) q[3];
sx q[3];
rz(-1.8504986) q[3];
sx q[3];
rz(-3.042799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2336642) q[0];
sx q[0];
rz(-1.5282682) q[0];
sx q[0];
rz(1.044957) q[0];
rz(-2.9772671) q[1];
sx q[1];
rz(-1.3427837) q[1];
sx q[1];
rz(-0.16990653) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92238802) q[0];
sx q[0];
rz(-1.4996075) q[0];
sx q[0];
rz(1.7573331) q[0];
rz(-pi) q[1];
rz(-1.2163148) q[2];
sx q[2];
rz(-1.2574242) q[2];
sx q[2];
rz(-2.2872668) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3376209) q[1];
sx q[1];
rz(-0.70120431) q[1];
sx q[1];
rz(-1.9528725) q[1];
rz(-pi) q[2];
rz(1.702013) q[3];
sx q[3];
rz(-2.2926712) q[3];
sx q[3];
rz(3.1027681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.56875151) q[2];
sx q[2];
rz(-0.70256394) q[2];
sx q[2];
rz(2.11002) q[2];
rz(-2.0555563) q[3];
sx q[3];
rz(-2.3417754) q[3];
sx q[3];
rz(-1.4097376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.80424911) q[0];
sx q[0];
rz(-1.046109) q[0];
sx q[0];
rz(-1.401249) q[0];
rz(2.7119472) q[1];
sx q[1];
rz(-2.255217) q[1];
sx q[1];
rz(-1.8362129) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1032377) q[0];
sx q[0];
rz(-1.2631386) q[0];
sx q[0];
rz(0.75958138) q[0];
rz(-pi) q[1];
rz(0.0045417343) q[2];
sx q[2];
rz(-1.1318739) q[2];
sx q[2];
rz(-0.86290854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67638799) q[1];
sx q[1];
rz(-0.25034764) q[1];
sx q[1];
rz(-0.12488229) q[1];
rz(-pi) q[2];
rz(-2.1277818) q[3];
sx q[3];
rz(-0.28388043) q[3];
sx q[3];
rz(0.59040961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0171011) q[2];
sx q[2];
rz(-1.93511) q[2];
sx q[2];
rz(-0.98480946) q[2];
rz(-1.5752327) q[3];
sx q[3];
rz(-1.5519578) q[3];
sx q[3];
rz(-2.8563833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.259909) q[0];
sx q[0];
rz(-0.30170983) q[0];
sx q[0];
rz(-1.786422) q[0];
rz(3.1175218) q[1];
sx q[1];
rz(-1.5211952) q[1];
sx q[1];
rz(-2.7640061) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628869) q[0];
sx q[0];
rz(-2.2170904) q[0];
sx q[0];
rz(-1.3747526) q[0];
rz(-1.7837011) q[2];
sx q[2];
rz(-1.3981552) q[2];
sx q[2];
rz(-1.5162692) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7086975) q[1];
sx q[1];
rz(-0.55799498) q[1];
sx q[1];
rz(-0.79496986) q[1];
x q[2];
rz(2.2144775) q[3];
sx q[3];
rz(-0.52442951) q[3];
sx q[3];
rz(2.7853109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5923803) q[2];
sx q[2];
rz(-2.8040631) q[2];
sx q[2];
rz(-0.40360061) q[2];
rz(0.18925439) q[3];
sx q[3];
rz(-2.0892102) q[3];
sx q[3];
rz(-0.45690593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6117578) q[0];
sx q[0];
rz(-0.54946041) q[0];
sx q[0];
rz(-2.8305565) q[0];
rz(-2.1850736) q[1];
sx q[1];
rz(-1.4776968) q[1];
sx q[1];
rz(2.8172353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81203359) q[0];
sx q[0];
rz(-2.5987465) q[0];
sx q[0];
rz(-0.37522786) q[0];
x q[1];
rz(-1.6300522) q[2];
sx q[2];
rz(-1.0612592) q[2];
sx q[2];
rz(1.5006202) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.94453401) q[1];
sx q[1];
rz(-0.98577104) q[1];
sx q[1];
rz(0.50307517) q[1];
rz(2.0288386) q[3];
sx q[3];
rz(-0.73087091) q[3];
sx q[3];
rz(0.79938408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73584622) q[2];
sx q[2];
rz(-0.79068557) q[2];
sx q[2];
rz(2.9131367) q[2];
rz(-2.6089846) q[3];
sx q[3];
rz(-2.3753128) q[3];
sx q[3];
rz(-0.84958357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.547895) q[0];
sx q[0];
rz(-2.6314647) q[0];
sx q[0];
rz(-1.2702031) q[0];
rz(2.5529329) q[1];
sx q[1];
rz(-1.9344067) q[1];
sx q[1];
rz(0.38280815) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13148334) q[0];
sx q[0];
rz(-1.3416222) q[0];
sx q[0];
rz(0.025529941) q[0];
rz(-pi) q[1];
rz(-0.50254681) q[2];
sx q[2];
rz(-0.48641962) q[2];
sx q[2];
rz(-0.18196276) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8296903) q[1];
sx q[1];
rz(-1.4451348) q[1];
sx q[1];
rz(0.71255334) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0625913) q[3];
sx q[3];
rz(-1.8100836) q[3];
sx q[3];
rz(-0.21904473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9855014) q[2];
sx q[2];
rz(-1.568855) q[2];
sx q[2];
rz(-2.7413979) q[2];
rz(-2.8943446) q[3];
sx q[3];
rz(-1.6797545) q[3];
sx q[3];
rz(0.39321536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3006712) q[0];
sx q[0];
rz(-1.8074169) q[0];
sx q[0];
rz(-2.1260496) q[0];
rz(0.66889846) q[1];
sx q[1];
rz(-1.4543507) q[1];
sx q[1];
rz(-1.5060172) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43011623) q[0];
sx q[0];
rz(-1.073494) q[0];
sx q[0];
rz(0.74691746) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54803063) q[2];
sx q[2];
rz(-1.9491674) q[2];
sx q[2];
rz(3.1032094) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4989711) q[1];
sx q[1];
rz(-1.6097415) q[1];
sx q[1];
rz(1.0907111) q[1];
x q[2];
rz(-0.89370905) q[3];
sx q[3];
rz(-2.0314616) q[3];
sx q[3];
rz(-0.078848293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.16613913) q[2];
sx q[2];
rz(-1.0756805) q[2];
sx q[2];
rz(1.6157185) q[2];
rz(2.8499991) q[3];
sx q[3];
rz(-2.6988131) q[3];
sx q[3];
rz(0.35167882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6780554) q[0];
sx q[0];
rz(-2.1778477) q[0];
sx q[0];
rz(0.5974593) q[0];
rz(-2.8529104) q[1];
sx q[1];
rz(-0.79816993) q[1];
sx q[1];
rz(0.023963902) q[1];
rz(-0.409375) q[2];
sx q[2];
rz(-1.2843046) q[2];
sx q[2];
rz(1.6094111) q[2];
rz(2.6816159) q[3];
sx q[3];
rz(-1.3446076) q[3];
sx q[3];
rz(1.4190058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
