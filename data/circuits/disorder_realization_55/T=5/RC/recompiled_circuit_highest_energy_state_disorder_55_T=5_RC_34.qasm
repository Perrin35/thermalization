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
rz(0.060184181) q[0];
sx q[0];
rz(4.2005778) q[0];
sx q[0];
rz(8.4146001) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(-0.29616907) q[1];
sx q[1];
rz(-2.8316166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84278569) q[0];
sx q[0];
rz(-1.1514542) q[0];
sx q[0];
rz(1.692665) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9363382) q[2];
sx q[2];
rz(-1.5488834) q[2];
sx q[2];
rz(-2.2147629) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.12116005) q[1];
sx q[1];
rz(-1.766664) q[1];
sx q[1];
rz(-1.905026) q[1];
x q[2];
rz(-1.7916075) q[3];
sx q[3];
rz(-0.75149262) q[3];
sx q[3];
rz(-1.3205075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6628722) q[2];
sx q[2];
rz(-0.24820776) q[2];
sx q[2];
rz(-0.09566801) q[2];
rz(-0.11224789) q[3];
sx q[3];
rz(-0.87533689) q[3];
sx q[3];
rz(1.7101425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2410759) q[0];
sx q[0];
rz(-0.88923419) q[0];
sx q[0];
rz(2.526793) q[0];
rz(2.2531033) q[1];
sx q[1];
rz(-1.588984) q[1];
sx q[1];
rz(-0.49957553) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1785806) q[0];
sx q[0];
rz(-2.7546429) q[0];
sx q[0];
rz(-1.011417) q[0];
rz(0.26878727) q[2];
sx q[2];
rz(-2.6013881) q[2];
sx q[2];
rz(-0.42551431) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.37812472) q[1];
sx q[1];
rz(-0.92510447) q[1];
sx q[1];
rz(-0.21685361) q[1];
x q[2];
rz(2.1533215) q[3];
sx q[3];
rz(-1.2856021) q[3];
sx q[3];
rz(2.5799283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2483612) q[2];
sx q[2];
rz(-1.2373368) q[2];
sx q[2];
rz(3.1207747) q[2];
rz(1.9013532) q[3];
sx q[3];
rz(-1.6720684) q[3];
sx q[3];
rz(1.1851236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5027387) q[0];
sx q[0];
rz(-0.26832142) q[0];
sx q[0];
rz(0.33367208) q[0];
rz(0.73792136) q[1];
sx q[1];
rz(-1.2260022) q[1];
sx q[1];
rz(-0.12942448) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5564011) q[0];
sx q[0];
rz(-1.4117774) q[0];
sx q[0];
rz(-0.88944737) q[0];
rz(-2.470181) q[2];
sx q[2];
rz(-2.0575626) q[2];
sx q[2];
rz(-1.1306131) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30363917) q[1];
sx q[1];
rz(-2.9898414) q[1];
sx q[1];
rz(2.2084153) q[1];
rz(0.78227104) q[3];
sx q[3];
rz(-1.9214298) q[3];
sx q[3];
rz(-0.80814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1239329) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(-1.7714436) q[2];
rz(0.74639368) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(0.082775041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133404) q[0];
sx q[0];
rz(-1.2097825) q[0];
sx q[0];
rz(-2.5764537) q[0];
rz(-0.42916974) q[1];
sx q[1];
rz(-2.4950835) q[1];
sx q[1];
rz(-2.3133004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9681184) q[0];
sx q[0];
rz(-2.3096363) q[0];
sx q[0];
rz(-1.4627187) q[0];
rz(2.056869) q[2];
sx q[2];
rz(-1.3157433) q[2];
sx q[2];
rz(-1.4866831) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.1431378) q[1];
sx q[1];
rz(-1.8363442) q[1];
sx q[1];
rz(-2.6396855) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24803646) q[3];
sx q[3];
rz(-1.0674849) q[3];
sx q[3];
rz(0.65500427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68507489) q[2];
sx q[2];
rz(-1.061331) q[2];
sx q[2];
rz(1.7100517) q[2];
rz(0.93959129) q[3];
sx q[3];
rz(-0.51274931) q[3];
sx q[3];
rz(-3.140894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0058896) q[0];
sx q[0];
rz(-2.8764184) q[0];
sx q[0];
rz(1.8121207) q[0];
rz(1.1520518) q[1];
sx q[1];
rz(-2.0538797) q[1];
sx q[1];
rz(-3.1413445) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5329929) q[0];
sx q[0];
rz(-1.3260576) q[0];
sx q[0];
rz(-0.18969638) q[0];
rz(0.087281422) q[2];
sx q[2];
rz(-1.0382663) q[2];
sx q[2];
rz(1.1171578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9226886) q[1];
sx q[1];
rz(-0.62059015) q[1];
sx q[1];
rz(2.1922078) q[1];
rz(-2.4100634) q[3];
sx q[3];
rz(-1.5168744) q[3];
sx q[3];
rz(1.9526854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0033215) q[2];
sx q[2];
rz(-1.2500117) q[2];
sx q[2];
rz(0.0058343466) q[2];
rz(0.88998574) q[3];
sx q[3];
rz(-2.4599288) q[3];
sx q[3];
rz(-2.0148923) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24790813) q[0];
sx q[0];
rz(-1.4434781) q[0];
sx q[0];
rz(1.4877315) q[0];
rz(0.030390175) q[1];
sx q[1];
rz(-2.0149714) q[1];
sx q[1];
rz(0.94246513) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1346325) q[0];
sx q[0];
rz(-0.31123268) q[0];
sx q[0];
rz(-0.65629058) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8428188) q[2];
sx q[2];
rz(-1.818383) q[2];
sx q[2];
rz(2.5556759) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2561431) q[1];
sx q[1];
rz(-2.0383308) q[1];
sx q[1];
rz(-0.3216775) q[1];
rz(1.0951772) q[3];
sx q[3];
rz(-1.3381357) q[3];
sx q[3];
rz(1.7847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5421062) q[2];
sx q[2];
rz(-0.78080559) q[2];
sx q[2];
rz(1.3055118) q[2];
rz(2.5944338) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(2.3356596) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18043537) q[0];
sx q[0];
rz(-1.0338217) q[0];
sx q[0];
rz(-0.10051522) q[0];
rz(-2.6590977) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(-0.85711342) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8554358) q[0];
sx q[0];
rz(-1.248994) q[0];
sx q[0];
rz(-2.4706173) q[0];
rz(-2.4947462) q[2];
sx q[2];
rz(-1.610272) q[2];
sx q[2];
rz(1.0919065) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4883721) q[1];
sx q[1];
rz(-1.3626507) q[1];
sx q[1];
rz(2.7339762) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1825652) q[3];
sx q[3];
rz(-0.60743466) q[3];
sx q[3];
rz(1.5608112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7414005) q[2];
sx q[2];
rz(-0.97101784) q[2];
sx q[2];
rz(0.3024438) q[2];
rz(0.97964573) q[3];
sx q[3];
rz(-2.1392348) q[3];
sx q[3];
rz(-1.8625331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.361146) q[0];
sx q[0];
rz(-0.13872153) q[0];
sx q[0];
rz(0.030990344) q[0];
rz(0.57394761) q[1];
sx q[1];
rz(-1.4731044) q[1];
sx q[1];
rz(-1.2773638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63420682) q[0];
sx q[0];
rz(-1.9045826) q[0];
sx q[0];
rz(2.7884363) q[0];
rz(1.0724284) q[2];
sx q[2];
rz(-1.785009) q[2];
sx q[2];
rz(1.5276078) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22945089) q[1];
sx q[1];
rz(-1.0849285) q[1];
sx q[1];
rz(-3.058601) q[1];
rz(-pi) q[2];
rz(1.2856917) q[3];
sx q[3];
rz(-1.630097) q[3];
sx q[3];
rz(1.7948732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0370827) q[2];
sx q[2];
rz(-3.0057378) q[2];
sx q[2];
rz(-2.6541397) q[2];
rz(2.4449352) q[3];
sx q[3];
rz(-0.928855) q[3];
sx q[3];
rz(3.0776183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071851991) q[0];
sx q[0];
rz(-2.6672279) q[0];
sx q[0];
rz(-2.790614) q[0];
rz(-0.17414302) q[1];
sx q[1];
rz(-1.5085647) q[1];
sx q[1];
rz(1.7399656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5183254) q[0];
sx q[0];
rz(-1.6931769) q[0];
sx q[0];
rz(0.24952475) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47408159) q[2];
sx q[2];
rz(-0.80898058) q[2];
sx q[2];
rz(-2.7681729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.40608866) q[1];
sx q[1];
rz(-1.5417409) q[1];
sx q[1];
rz(1.7543704) q[1];
rz(-2.4003567) q[3];
sx q[3];
rz(-2.0307856) q[3];
sx q[3];
rz(2.7335087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.031124) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(-0.6366716) q[2];
rz(-2.0311671) q[3];
sx q[3];
rz(-1.5444376) q[3];
sx q[3];
rz(-0.5184263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7745895) q[0];
sx q[0];
rz(-0.29357266) q[0];
sx q[0];
rz(-2.0830182) q[0];
rz(-3.1029347) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(-1.0640594) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.80262) q[0];
sx q[0];
rz(-1.5674599) q[0];
sx q[0];
rz(-3.1371069) q[0];
x q[1];
rz(0.27711192) q[2];
sx q[2];
rz(-1.849035) q[2];
sx q[2];
rz(0.81467512) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.56434332) q[1];
sx q[1];
rz(-1.5303622) q[1];
sx q[1];
rz(3.0811429) q[1];
rz(-0.45952103) q[3];
sx q[3];
rz(-1.6020613) q[3];
sx q[3];
rz(2.7783436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4260063) q[2];
sx q[2];
rz(-2.6435489) q[2];
sx q[2];
rz(-0.074020298) q[2];
rz(2.7600539) q[3];
sx q[3];
rz(-1.3149202) q[3];
sx q[3];
rz(2.7874302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1114125) q[0];
sx q[0];
rz(-1.3127865) q[0];
sx q[0];
rz(2.4319613) q[0];
rz(2.5297655) q[1];
sx q[1];
rz(-0.71129967) q[1];
sx q[1];
rz(1.5536972) q[1];
rz(2.5135573) q[2];
sx q[2];
rz(-2.5428094) q[2];
sx q[2];
rz(1.6369292) q[2];
rz(2.3930876) q[3];
sx q[3];
rz(-1.9288962) q[3];
sx q[3];
rz(0.68652912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
