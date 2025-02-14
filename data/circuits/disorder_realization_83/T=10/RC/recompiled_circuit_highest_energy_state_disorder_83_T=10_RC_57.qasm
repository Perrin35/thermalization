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
rz(2.7163765) q[0];
sx q[0];
rz(-1.334231) q[0];
sx q[0];
rz(0.38036007) q[0];
rz(-1.3588139) q[1];
sx q[1];
rz(-2.9549197) q[1];
sx q[1];
rz(2.2169854) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7036669) q[0];
sx q[0];
rz(-2.3406918) q[0];
sx q[0];
rz(-0.4842224) q[0];
rz(-pi) q[1];
rz(-0.1670465) q[2];
sx q[2];
rz(-2.0809329) q[2];
sx q[2];
rz(-1.270592) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8262337) q[1];
sx q[1];
rz(-2.4118638) q[1];
sx q[1];
rz(-1.6297831) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0136631) q[3];
sx q[3];
rz(-1.5618213) q[3];
sx q[3];
rz(-0.25024807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1251462) q[2];
sx q[2];
rz(-2.1852198) q[2];
sx q[2];
rz(-2.1535786) q[2];
rz(0.2068578) q[3];
sx q[3];
rz(-1.4068973) q[3];
sx q[3];
rz(0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4619231) q[0];
sx q[0];
rz(-1.4689057) q[0];
sx q[0];
rz(-2.8894506) q[0];
rz(3.120046) q[1];
sx q[1];
rz(-0.69308678) q[1];
sx q[1];
rz(0.85539877) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1148672) q[0];
sx q[0];
rz(-1.572319) q[0];
sx q[0];
rz(1.5437838) q[0];
rz(2.5874675) q[2];
sx q[2];
rz(-1.9337855) q[2];
sx q[2];
rz(0.83152321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63355061) q[1];
sx q[1];
rz(-1.7910069) q[1];
sx q[1];
rz(-2.1160265) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7383201) q[3];
sx q[3];
rz(-2.4174815) q[3];
sx q[3];
rz(-1.4693174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97027332) q[2];
sx q[2];
rz(-3.1372941) q[2];
sx q[2];
rz(-2.8805736) q[2];
rz(-3.1018992) q[3];
sx q[3];
rz(-1.3986162) q[3];
sx q[3];
rz(-0.54350129) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4259341) q[0];
sx q[0];
rz(-1.4733227) q[0];
sx q[0];
rz(-2.4470827) q[0];
rz(2.014324) q[1];
sx q[1];
rz(-2.290461) q[1];
sx q[1];
rz(0.99004254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5607669) q[0];
sx q[0];
rz(-0.35835727) q[0];
sx q[0];
rz(0.23507765) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98857356) q[2];
sx q[2];
rz(-1.7372088) q[2];
sx q[2];
rz(-0.63936448) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5591919) q[1];
sx q[1];
rz(-1.2123931) q[1];
sx q[1];
rz(2.6646975) q[1];
rz(-pi) q[2];
rz(-0.86750908) q[3];
sx q[3];
rz(-1.9990307) q[3];
sx q[3];
rz(-3.0813062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.29828829) q[2];
sx q[2];
rz(-0.97323209) q[2];
sx q[2];
rz(-0.49015552) q[2];
rz(0.35689029) q[3];
sx q[3];
rz(-2.7481952) q[3];
sx q[3];
rz(1.2662158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0855584) q[0];
sx q[0];
rz(-1.1857251) q[0];
sx q[0];
rz(-2.2991142) q[0];
rz(0.8017686) q[1];
sx q[1];
rz(-2.8623878) q[1];
sx q[1];
rz(-2.3559949) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7467037) q[0];
sx q[0];
rz(-2.3869042) q[0];
sx q[0];
rz(-0.37209038) q[0];
rz(-pi) q[1];
rz(-1.4198205) q[2];
sx q[2];
rz(-2.2822126) q[2];
sx q[2];
rz(2.6670418) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.468049) q[1];
sx q[1];
rz(-2.1077413) q[1];
sx q[1];
rz(1.931744) q[1];
rz(-pi) q[2];
rz(2.9369041) q[3];
sx q[3];
rz(-0.83016268) q[3];
sx q[3];
rz(-1.9616753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.12513146) q[2];
sx q[2];
rz(-2.3931914) q[2];
sx q[2];
rz(1.008519) q[2];
rz(-0.85159167) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(2.0612702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1763024) q[0];
sx q[0];
rz(-2.1602614) q[0];
sx q[0];
rz(-2.0710131) q[0];
rz(-0.77752441) q[1];
sx q[1];
rz(-0.91064015) q[1];
sx q[1];
rz(1.0106962) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99791894) q[0];
sx q[0];
rz(-2.5287313) q[0];
sx q[0];
rz(0.65632239) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3959051) q[2];
sx q[2];
rz(-1.6351467) q[2];
sx q[2];
rz(-1.627587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3123828) q[1];
sx q[1];
rz(-2.0264033) q[1];
sx q[1];
rz(0.73789755) q[1];
rz(-pi) q[2];
rz(-1.4713953) q[3];
sx q[3];
rz(-2.1306921) q[3];
sx q[3];
rz(-2.3178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.30194482) q[2];
sx q[2];
rz(-0.10493111) q[2];
sx q[2];
rz(1.5555752) q[2];
rz(-2.1258866) q[3];
sx q[3];
rz(-1.1414707) q[3];
sx q[3];
rz(1.4122081) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509336) q[0];
sx q[0];
rz(-2.9582294) q[0];
sx q[0];
rz(0.80107981) q[0];
rz(3.0084897) q[1];
sx q[1];
rz(-1.8201273) q[1];
sx q[1];
rz(0.4020234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4426081) q[0];
sx q[0];
rz(-0.8359209) q[0];
sx q[0];
rz(-2.0280272) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0021462321) q[2];
sx q[2];
rz(-1.2800001) q[2];
sx q[2];
rz(0.69974297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.133144) q[1];
sx q[1];
rz(-0.56479543) q[1];
sx q[1];
rz(-1.0393591) q[1];
x q[2];
rz(-1.4457213) q[3];
sx q[3];
rz(-0.91642028) q[3];
sx q[3];
rz(2.666111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63783995) q[2];
sx q[2];
rz(-1.736234) q[2];
sx q[2];
rz(-0.77581882) q[2];
rz(1.4555629) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(2.1337401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0966454) q[0];
sx q[0];
rz(-0.89770397) q[0];
sx q[0];
rz(0.67895472) q[0];
rz(1.2848162) q[1];
sx q[1];
rz(-2.4285451) q[1];
sx q[1];
rz(-1.3791893) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5011936) q[0];
sx q[0];
rz(-1.7441347) q[0];
sx q[0];
rz(0.79357432) q[0];
rz(-pi) q[1];
rz(-1.2871615) q[2];
sx q[2];
rz(-0.29046187) q[2];
sx q[2];
rz(1.0303555) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0022565) q[1];
sx q[1];
rz(-0.69076194) q[1];
sx q[1];
rz(-0.73378566) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3706743) q[3];
sx q[3];
rz(-1.516023) q[3];
sx q[3];
rz(-1.9211662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3369559) q[2];
sx q[2];
rz(-0.75504428) q[2];
sx q[2];
rz(1.2389368) q[2];
rz(1.3321446) q[3];
sx q[3];
rz(-0.76056162) q[3];
sx q[3];
rz(-2.8968887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181344) q[0];
sx q[0];
rz(-1.9754388) q[0];
sx q[0];
rz(-3.0031257) q[0];
rz(2.9578517) q[1];
sx q[1];
rz(-2.467149) q[1];
sx q[1];
rz(2.8750681) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040695) q[0];
sx q[0];
rz(-0.36629391) q[0];
sx q[0];
rz(-1.1031325) q[0];
rz(-pi) q[1];
rz(1.9815923) q[2];
sx q[2];
rz(-0.60705429) q[2];
sx q[2];
rz(2.0701054) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.16921088) q[1];
sx q[1];
rz(-1.6270429) q[1];
sx q[1];
rz(-1.1509658) q[1];
rz(-1.5795535) q[3];
sx q[3];
rz(-2.2976795) q[3];
sx q[3];
rz(-1.383964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0900241) q[2];
sx q[2];
rz(-1.2568306) q[2];
sx q[2];
rz(2.9070692) q[2];
rz(-1.4759493) q[3];
sx q[3];
rz(-1.539295) q[3];
sx q[3];
rz(-2.3852111) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.804857) q[0];
sx q[0];
rz(-2.2754301) q[0];
sx q[0];
rz(-2.3418703) q[0];
rz(2.5526478) q[1];
sx q[1];
rz(-2.2942693) q[1];
sx q[1];
rz(2.934093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0711771) q[0];
sx q[0];
rz(-1.1780329) q[0];
sx q[0];
rz(-1.7123187) q[0];
x q[1];
rz(0.3090119) q[2];
sx q[2];
rz(-2.2591619) q[2];
sx q[2];
rz(2.7613044) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5821857) q[1];
sx q[1];
rz(-1.9684125) q[1];
sx q[1];
rz(0.84487265) q[1];
x q[2];
rz(0.028771632) q[3];
sx q[3];
rz(-2.040401) q[3];
sx q[3];
rz(-2.4503277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3346682) q[2];
sx q[2];
rz(-0.88403264) q[2];
sx q[2];
rz(2.7743288) q[2];
rz(1.7043097) q[3];
sx q[3];
rz(-1.9673012) q[3];
sx q[3];
rz(1.9678763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2631898) q[0];
sx q[0];
rz(-2.1387565) q[0];
sx q[0];
rz(0.81018418) q[0];
rz(1.1148249) q[1];
sx q[1];
rz(-0.38433847) q[1];
sx q[1];
rz(-2.5883163) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5636175) q[0];
sx q[0];
rz(-1.5153236) q[0];
sx q[0];
rz(0.93240057) q[0];
x q[1];
rz(0.79375879) q[2];
sx q[2];
rz(-2.0362034) q[2];
sx q[2];
rz(2.2528354) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3590433) q[1];
sx q[1];
rz(-1.7573253) q[1];
sx q[1];
rz(1.2839509) q[1];
x q[2];
rz(2.7891879) q[3];
sx q[3];
rz(-0.81768546) q[3];
sx q[3];
rz(-1.529983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1222003) q[2];
sx q[2];
rz(-2.4330752) q[2];
sx q[2];
rz(-2.9019287) q[2];
rz(1.8137118) q[3];
sx q[3];
rz(-1.0646822) q[3];
sx q[3];
rz(-2.6740668) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079530579) q[0];
sx q[0];
rz(-1.3539599) q[0];
sx q[0];
rz(0.39394105) q[0];
rz(-0.65555864) q[1];
sx q[1];
rz(-1.9656904) q[1];
sx q[1];
rz(0.94737731) q[1];
rz(2.4782933) q[2];
sx q[2];
rz(-2.6313783) q[2];
sx q[2];
rz(0.22102697) q[2];
rz(-2.7648457) q[3];
sx q[3];
rz(-0.6142817) q[3];
sx q[3];
rz(2.3913419) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
