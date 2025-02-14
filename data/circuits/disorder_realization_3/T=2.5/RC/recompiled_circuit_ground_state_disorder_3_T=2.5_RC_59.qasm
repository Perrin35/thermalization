OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8752911) q[0];
sx q[0];
rz(-2.7346791) q[0];
sx q[0];
rz(-0.20242515) q[0];
rz(-1.3099194) q[1];
sx q[1];
rz(-1.1453495) q[1];
sx q[1];
rz(0.98869148) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3920445) q[0];
sx q[0];
rz(-2.364416) q[0];
sx q[0];
rz(-1.8424319) q[0];
x q[1];
rz(1.0568809) q[2];
sx q[2];
rz(-0.15794733) q[2];
sx q[2];
rz(-2.081209) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4943581) q[1];
sx q[1];
rz(-2.1796759) q[1];
sx q[1];
rz(-2.865164) q[1];
rz(0.12434953) q[3];
sx q[3];
rz(-1.1477587) q[3];
sx q[3];
rz(1.4165341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6053091) q[2];
sx q[2];
rz(-2.3399957) q[2];
sx q[2];
rz(1.0514222) q[2];
rz(1.5031523) q[3];
sx q[3];
rz(-1.7283311) q[3];
sx q[3];
rz(2.9063291) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94238344) q[0];
sx q[0];
rz(-1.0301882) q[0];
sx q[0];
rz(-2.8675766) q[0];
rz(0.38385299) q[1];
sx q[1];
rz(-2.5089788) q[1];
sx q[1];
rz(-1.8208108) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9649759) q[0];
sx q[0];
rz(-2.9701485) q[0];
sx q[0];
rz(1.9587396) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6422149) q[2];
sx q[2];
rz(-1.3670397) q[2];
sx q[2];
rz(0.2815385) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.49628851) q[1];
sx q[1];
rz(-0.72455614) q[1];
sx q[1];
rz(1.6425808) q[1];
rz(-pi) q[2];
rz(-2.3427547) q[3];
sx q[3];
rz(-0.96551408) q[3];
sx q[3];
rz(1.192361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5447834) q[2];
sx q[2];
rz(-1.7388672) q[2];
sx q[2];
rz(-1.7421494) q[2];
rz(-0.60747373) q[3];
sx q[3];
rz(-1.6739269) q[3];
sx q[3];
rz(2.7510711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6121599) q[0];
sx q[0];
rz(-2.2859892) q[0];
sx q[0];
rz(-2.8223619) q[0];
rz(-3.0901129) q[1];
sx q[1];
rz(-1.1795283) q[1];
sx q[1];
rz(1.0565588) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7565472) q[0];
sx q[0];
rz(-0.17160417) q[0];
sx q[0];
rz(-1.8770939) q[0];
rz(-pi) q[1];
rz(-0.86528565) q[2];
sx q[2];
rz(-1.533939) q[2];
sx q[2];
rz(1.3071905) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.33256862) q[1];
sx q[1];
rz(-1.7454552) q[1];
sx q[1];
rz(0.93775429) q[1];
rz(-pi) q[2];
rz(-1.7848739) q[3];
sx q[3];
rz(-2.7299462) q[3];
sx q[3];
rz(0.95142196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9487379) q[2];
sx q[2];
rz(-1.7098018) q[2];
sx q[2];
rz(1.7970435) q[2];
rz(0.91790849) q[3];
sx q[3];
rz(-2.5036006) q[3];
sx q[3];
rz(1.8857672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6767839) q[0];
sx q[0];
rz(-2.5510241) q[0];
sx q[0];
rz(-1.6326686) q[0];
rz(-0.50679874) q[1];
sx q[1];
rz(-1.2134039) q[1];
sx q[1];
rz(2.7510344) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9386425) q[0];
sx q[0];
rz(-0.0084059518) q[0];
sx q[0];
rz(0.52830835) q[0];
rz(-2.6919011) q[2];
sx q[2];
rz(-2.7241787) q[2];
sx q[2];
rz(-2.6455622) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4949172) q[1];
sx q[1];
rz(-0.61958379) q[1];
sx q[1];
rz(1.5219206) q[1];
x q[2];
rz(-2.9831577) q[3];
sx q[3];
rz(-1.6601175) q[3];
sx q[3];
rz(1.4591372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8916696) q[2];
sx q[2];
rz(-2.2107783) q[2];
sx q[2];
rz(2.2451952) q[2];
rz(2.2253288) q[3];
sx q[3];
rz(-0.93583411) q[3];
sx q[3];
rz(2.1808482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59178281) q[0];
sx q[0];
rz(-0.35583219) q[0];
sx q[0];
rz(-0.47019666) q[0];
rz(-1.4037941) q[1];
sx q[1];
rz(-2.4411968) q[1];
sx q[1];
rz(-1.9139003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0225343) q[0];
sx q[0];
rz(-1.5644363) q[0];
sx q[0];
rz(1.8409078) q[0];
x q[1];
rz(2.50493) q[2];
sx q[2];
rz(-2.4877254) q[2];
sx q[2];
rz(-0.9034397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.073138) q[1];
sx q[1];
rz(-0.92759354) q[1];
sx q[1];
rz(-3.0469608) q[1];
rz(2.6032734) q[3];
sx q[3];
rz(-1.3019058) q[3];
sx q[3];
rz(-1.5017832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1710743) q[2];
sx q[2];
rz(-0.71275622) q[2];
sx q[2];
rz(-1.7573382) q[2];
rz(0.61128831) q[3];
sx q[3];
rz(-1.7620112) q[3];
sx q[3];
rz(2.7454564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1432081) q[0];
sx q[0];
rz(-1.0118326) q[0];
sx q[0];
rz(0.55214733) q[0];
rz(-0.72969189) q[1];
sx q[1];
rz(-2.153986) q[1];
sx q[1];
rz(2.139835) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4491548) q[0];
sx q[0];
rz(-1.7040515) q[0];
sx q[0];
rz(-3.0603581) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92879734) q[2];
sx q[2];
rz(-1.1603519) q[2];
sx q[2];
rz(-2.9985119) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.6164847) q[1];
sx q[1];
rz(-1.387038) q[1];
sx q[1];
rz(0.77145578) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4625307) q[3];
sx q[3];
rz(-1.1521395) q[3];
sx q[3];
rz(-0.94840324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6843162) q[2];
sx q[2];
rz(-1.7924954) q[2];
sx q[2];
rz(1.3433733) q[2];
rz(0.7243048) q[3];
sx q[3];
rz(-1.2035921) q[3];
sx q[3];
rz(-2.8389285) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90172076) q[0];
sx q[0];
rz(-1.3322823) q[0];
sx q[0];
rz(2.4161762) q[0];
rz(1.8431009) q[1];
sx q[1];
rz(-2.6893171) q[1];
sx q[1];
rz(0.27870146) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.168047) q[0];
sx q[0];
rz(-0.5613882) q[0];
sx q[0];
rz(2.8532203) q[0];
rz(-pi) q[1];
rz(2.8093074) q[2];
sx q[2];
rz(-1.1934308) q[2];
sx q[2];
rz(-1.2720296) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78165752) q[1];
sx q[1];
rz(-1.8751448) q[1];
sx q[1];
rz(0.84142797) q[1];
rz(-1.0213357) q[3];
sx q[3];
rz(-2.3186893) q[3];
sx q[3];
rz(-1.8902035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7574888) q[2];
sx q[2];
rz(-2.8748547) q[2];
sx q[2];
rz(0.98814803) q[2];
rz(3.0366376) q[3];
sx q[3];
rz(-1.8433439) q[3];
sx q[3];
rz(-2.9102563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4526378) q[0];
sx q[0];
rz(-2.3404558) q[0];
sx q[0];
rz(-0.6193921) q[0];
rz(-3.1267005) q[1];
sx q[1];
rz(-0.89034021) q[1];
sx q[1];
rz(-0.68170086) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68530689) q[0];
sx q[0];
rz(-0.62236747) q[0];
sx q[0];
rz(-2.078889) q[0];
rz(0.9056717) q[2];
sx q[2];
rz(-1.5076414) q[2];
sx q[2];
rz(0.81005972) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8047844) q[1];
sx q[1];
rz(-0.35165916) q[1];
sx q[1];
rz(3.1218887) q[1];
rz(-2.3084749) q[3];
sx q[3];
rz(-1.4280115) q[3];
sx q[3];
rz(-0.24417434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2316078) q[2];
sx q[2];
rz(-2.3384428) q[2];
sx q[2];
rz(-2.8019359) q[2];
rz(0.61399442) q[3];
sx q[3];
rz(-1.3795779) q[3];
sx q[3];
rz(-1.5617255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15641077) q[0];
sx q[0];
rz(-2.7934533) q[0];
sx q[0];
rz(0.85439318) q[0];
rz(-2.5482381) q[1];
sx q[1];
rz(-1.0320458) q[1];
sx q[1];
rz(-0.29534435) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2345384) q[0];
sx q[0];
rz(-3.0772644) q[0];
sx q[0];
rz(0.5366908) q[0];
x q[1];
rz(1.6215084) q[2];
sx q[2];
rz(-2.7233633) q[2];
sx q[2];
rz(0.39722463) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1251245) q[1];
sx q[1];
rz(-2.0464026) q[1];
sx q[1];
rz(0.7866025) q[1];
x q[2];
rz(-1.909524) q[3];
sx q[3];
rz(-0.32707387) q[3];
sx q[3];
rz(-1.0504521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16193834) q[2];
sx q[2];
rz(-1.6831968) q[2];
sx q[2];
rz(2.4597607) q[2];
rz(-1.2006867) q[3];
sx q[3];
rz(-2.3537945) q[3];
sx q[3];
rz(2.70575) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1107263) q[0];
sx q[0];
rz(-1.2799355) q[0];
sx q[0];
rz(-2.6357292) q[0];
rz(-1.0181631) q[1];
sx q[1];
rz(-1.70645) q[1];
sx q[1];
rz(2.2811269) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2967401) q[0];
sx q[0];
rz(-2.2460947) q[0];
sx q[0];
rz(-1.5102415) q[0];
rz(-pi) q[1];
rz(-1.8637723) q[2];
sx q[2];
rz(-1.7070012) q[2];
sx q[2];
rz(-0.53060645) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.2194032) q[1];
sx q[1];
rz(-1.8199395) q[1];
sx q[1];
rz(-1.6046899) q[1];
x q[2];
rz(-0.56475909) q[3];
sx q[3];
rz(-2.107672) q[3];
sx q[3];
rz(0.91954654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0597509) q[2];
sx q[2];
rz(-1.1885208) q[2];
sx q[2];
rz(0.58319485) q[2];
rz(-3.1405295) q[3];
sx q[3];
rz(-0.93102264) q[3];
sx q[3];
rz(-0.49334905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4790333) q[0];
sx q[0];
rz(-1.4190577) q[0];
sx q[0];
rz(0.90202913) q[0];
rz(-2.5262911) q[1];
sx q[1];
rz(-1.6533783) q[1];
sx q[1];
rz(1.3396214) q[1];
rz(-1.8281219) q[2];
sx q[2];
rz(-2.0791292) q[2];
sx q[2];
rz(2.8544478) q[2];
rz(1.2801805) q[3];
sx q[3];
rz(-2.5463085) q[3];
sx q[3];
rz(0.13704722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
