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
rz(-2.7849164) q[0];
sx q[0];
rz(-1.4465605) q[0];
sx q[0];
rz(-2.2398228) q[0];
rz(-1.8094485) q[1];
sx q[1];
rz(-0.44337115) q[1];
sx q[1];
rz(0.765257) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5072696) q[0];
sx q[0];
rz(-2.1131086) q[0];
sx q[0];
rz(2.7141098) q[0];
x q[1];
rz(0.24306007) q[2];
sx q[2];
rz(-1.2420734) q[2];
sx q[2];
rz(2.6459733) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2697077) q[1];
sx q[1];
rz(-0.34075865) q[1];
sx q[1];
rz(-3.0213256) q[1];
rz(-0.13007836) q[3];
sx q[3];
rz(-1.8474694) q[3];
sx q[3];
rz(-0.94844669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29609933) q[2];
sx q[2];
rz(-2.0534434) q[2];
sx q[2];
rz(1.4720526) q[2];
rz(1.7146401) q[3];
sx q[3];
rz(-1.6441556) q[3];
sx q[3];
rz(0.53372395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.958309) q[0];
sx q[0];
rz(-1.5077718) q[0];
sx q[0];
rz(0.84877745) q[0];
rz(-0.035004184) q[1];
sx q[1];
rz(-1.1468381) q[1];
sx q[1];
rz(2.1880207) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5262819) q[0];
sx q[0];
rz(-1.0089951) q[0];
sx q[0];
rz(1.022382) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1814647) q[2];
sx q[2];
rz(-1.2876533) q[2];
sx q[2];
rz(-0.87923899) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.63605308) q[1];
sx q[1];
rz(-1.2314774) q[1];
sx q[1];
rz(2.9341594) q[1];
rz(-pi) q[2];
rz(0.93718174) q[3];
sx q[3];
rz(-0.92234334) q[3];
sx q[3];
rz(2.5536553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1635052) q[2];
sx q[2];
rz(-1.7814025) q[2];
sx q[2];
rz(-0.84954849) q[2];
rz(-0.16544011) q[3];
sx q[3];
rz(-1.4435507) q[3];
sx q[3];
rz(0.74976841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7716832) q[0];
sx q[0];
rz(-2.219438) q[0];
sx q[0];
rz(-0.40192303) q[0];
rz(2.9882714) q[1];
sx q[1];
rz(-0.17686495) q[1];
sx q[1];
rz(1.1361928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3055666) q[0];
sx q[0];
rz(-1.5587806) q[0];
sx q[0];
rz(-1.9703321) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49164518) q[2];
sx q[2];
rz(-0.53722135) q[2];
sx q[2];
rz(-1.1342837) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0979251) q[1];
sx q[1];
rz(-0.33722028) q[1];
sx q[1];
rz(-0.80101669) q[1];
rz(2.4506684) q[3];
sx q[3];
rz(-0.57583416) q[3];
sx q[3];
rz(1.3802647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.084126964) q[2];
sx q[2];
rz(-0.87856138) q[2];
sx q[2];
rz(-3.086536) q[2];
rz(1.347524) q[3];
sx q[3];
rz(-1.811458) q[3];
sx q[3];
rz(-2.1373035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7771626) q[0];
sx q[0];
rz(-2.936852) q[0];
sx q[0];
rz(1.5675911) q[0];
rz(-0.69152999) q[1];
sx q[1];
rz(-2.4722996) q[1];
sx q[1];
rz(0.18542586) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2381993) q[0];
sx q[0];
rz(-1.967634) q[0];
sx q[0];
rz(1.7522041) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71132423) q[2];
sx q[2];
rz(-2.5344116) q[2];
sx q[2];
rz(1.173623) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54900995) q[1];
sx q[1];
rz(-1.8235778) q[1];
sx q[1];
rz(1.2493253) q[1];
rz(1.3641649) q[3];
sx q[3];
rz(-2.1071599) q[3];
sx q[3];
rz(3.0230776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4039679) q[2];
sx q[2];
rz(-1.4289958) q[2];
sx q[2];
rz(3.1363623) q[2];
rz(-2.7541006) q[3];
sx q[3];
rz(-0.5046851) q[3];
sx q[3];
rz(-1.3030049) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4150647) q[0];
sx q[0];
rz(-2.1580577) q[0];
sx q[0];
rz(0.61597419) q[0];
rz(1.1888602) q[1];
sx q[1];
rz(-0.6183466) q[1];
sx q[1];
rz(2.628285) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5145288) q[0];
sx q[0];
rz(-1.5634057) q[0];
sx q[0];
rz(1.8523995) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72577621) q[2];
sx q[2];
rz(-1.3575524) q[2];
sx q[2];
rz(-2.2346614) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0541546) q[1];
sx q[1];
rz(-0.63870027) q[1];
sx q[1];
rz(3.020546) q[1];
rz(-pi) q[2];
rz(-3.0232459) q[3];
sx q[3];
rz(-1.1062628) q[3];
sx q[3];
rz(-1.1778414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1335699) q[2];
sx q[2];
rz(-1.2761389) q[2];
sx q[2];
rz(-0.7424773) q[2];
rz(-1.6978469) q[3];
sx q[3];
rz(-1.4092813) q[3];
sx q[3];
rz(-1.1013364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3105069) q[0];
sx q[0];
rz(-1.0467014) q[0];
sx q[0];
rz(-0.35980862) q[0];
rz(-2.2911435) q[1];
sx q[1];
rz(-1.1812187) q[1];
sx q[1];
rz(-0.46583072) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1547928) q[0];
sx q[0];
rz(-1.7868817) q[0];
sx q[0];
rz(-2.297202) q[0];
rz(-2.3783422) q[2];
sx q[2];
rz(-1.6692729) q[2];
sx q[2];
rz(-1.5710448) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.19037828) q[1];
sx q[1];
rz(-0.24628993) q[1];
sx q[1];
rz(-0.21065335) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5683018) q[3];
sx q[3];
rz(-2.0985641) q[3];
sx q[3];
rz(-0.19095574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3992074) q[2];
sx q[2];
rz(-1.1804322) q[2];
sx q[2];
rz(-0.36435374) q[2];
rz(-2.897701) q[3];
sx q[3];
rz(-0.049592169) q[3];
sx q[3];
rz(1.2172788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.57109443) q[0];
sx q[0];
rz(-0.97935337) q[0];
sx q[0];
rz(1.97557) q[0];
rz(-1.0150389) q[1];
sx q[1];
rz(-0.53703419) q[1];
sx q[1];
rz(2.7551415) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84595118) q[0];
sx q[0];
rz(-1.5401473) q[0];
sx q[0];
rz(1.3611193) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68604821) q[2];
sx q[2];
rz(-0.39271388) q[2];
sx q[2];
rz(0.60143747) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4402953) q[1];
sx q[1];
rz(-3.0616425) q[1];
sx q[1];
rz(-1.246212) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12549716) q[3];
sx q[3];
rz(-2.1096131) q[3];
sx q[3];
rz(-1.5244689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40547907) q[2];
sx q[2];
rz(-1.6454641) q[2];
sx q[2];
rz(-2.9138937) q[2];
rz(2.0685711) q[3];
sx q[3];
rz(-0.74879542) q[3];
sx q[3];
rz(-2.8480215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88184083) q[0];
sx q[0];
rz(-2.2417534) q[0];
sx q[0];
rz(-2.7135799) q[0];
rz(-1.0182861) q[1];
sx q[1];
rz(-2.7052453) q[1];
sx q[1];
rz(2.2705618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7334325) q[0];
sx q[0];
rz(-1.0116315) q[0];
sx q[0];
rz(2.5446822) q[0];
rz(1.2134654) q[2];
sx q[2];
rz(-1.7646731) q[2];
sx q[2];
rz(2.9952637) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5586413) q[1];
sx q[1];
rz(-1.4498596) q[1];
sx q[1];
rz(-1.4287097) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0139129) q[3];
sx q[3];
rz(-2.172951) q[3];
sx q[3];
rz(0.35246655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.009306) q[2];
sx q[2];
rz(-0.79557482) q[2];
sx q[2];
rz(1.2786678) q[2];
rz(2.5631185) q[3];
sx q[3];
rz(-2.5090802) q[3];
sx q[3];
rz(1.9875897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6744252) q[0];
sx q[0];
rz(-2.9245057) q[0];
sx q[0];
rz(2.5417852) q[0];
rz(-1.2990052) q[1];
sx q[1];
rz(-1.6103585) q[1];
sx q[1];
rz(0.64186796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80461795) q[0];
sx q[0];
rz(-1.0952767) q[0];
sx q[0];
rz(-1.9561951) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9974143) q[2];
sx q[2];
rz(-0.81071172) q[2];
sx q[2];
rz(1.415103) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1915726) q[1];
sx q[1];
rz(-2.2784462) q[1];
sx q[1];
rz(1.5465082) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5893552) q[3];
sx q[3];
rz(-2.304616) q[3];
sx q[3];
rz(0.53572922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.32144) q[2];
sx q[2];
rz(-1.7165311) q[2];
sx q[2];
rz(2.5313306) q[2];
rz(-2.6214456) q[3];
sx q[3];
rz(-1.0545694) q[3];
sx q[3];
rz(-0.49351969) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.801012) q[0];
sx q[0];
rz(-1.2435253) q[0];
sx q[0];
rz(-0.93801671) q[0];
rz(2.7998789) q[1];
sx q[1];
rz(-1.7594756) q[1];
sx q[1];
rz(-0.15728532) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2153335) q[0];
sx q[0];
rz(-2.4865827) q[0];
sx q[0];
rz(2.8192855) q[0];
rz(-0.75802676) q[2];
sx q[2];
rz(-1.6322513) q[2];
sx q[2];
rz(0.058493424) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5963635) q[1];
sx q[1];
rz(-2.5198054) q[1];
sx q[1];
rz(0.53485628) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0299133) q[3];
sx q[3];
rz(-0.4678886) q[3];
sx q[3];
rz(2.6998126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9226795) q[2];
sx q[2];
rz(-0.34505406) q[2];
sx q[2];
rz(1.867713) q[2];
rz(-0.0056565469) q[3];
sx q[3];
rz(-1.4414682) q[3];
sx q[3];
rz(-0.98102942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8511178) q[0];
sx q[0];
rz(-1.9708451) q[0];
sx q[0];
rz(-3.0921902) q[0];
rz(-2.6750917) q[1];
sx q[1];
rz(-1.3460881) q[1];
sx q[1];
rz(0.64429611) q[1];
rz(-2.1860547) q[2];
sx q[2];
rz(-2.370859) q[2];
sx q[2];
rz(-1.4683409) q[2];
rz(1.4884819) q[3];
sx q[3];
rz(-2.001279) q[3];
sx q[3];
rz(1.6829987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
