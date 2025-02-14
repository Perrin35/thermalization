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
rz(-0.42521617) q[0];
sx q[0];
rz(-1.8073616) q[0];
sx q[0];
rz(-0.38036007) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(2.9549197) q[1];
sx q[1];
rz(8.5001707) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7036669) q[0];
sx q[0];
rz(-0.80090085) q[0];
sx q[0];
rz(-0.4842224) q[0];
rz(-pi) q[1];
rz(-2.9745462) q[2];
sx q[2];
rz(-2.0809329) q[2];
sx q[2];
rz(1.270592) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9052995) q[1];
sx q[1];
rz(-0.84262203) q[1];
sx q[1];
rz(0.052680101) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12792952) q[3];
sx q[3];
rz(-1.5797714) q[3];
sx q[3];
rz(-2.8913446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1251462) q[2];
sx q[2];
rz(-2.1852198) q[2];
sx q[2];
rz(-2.1535786) q[2];
rz(2.9347349) q[3];
sx q[3];
rz(-1.7346953) q[3];
sx q[3];
rz(0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4619231) q[0];
sx q[0];
rz(-1.4689057) q[0];
sx q[0];
rz(2.8894506) q[0];
rz(-3.120046) q[1];
sx q[1];
rz(-0.69308678) q[1];
sx q[1];
rz(-0.85539877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48777521) q[0];
sx q[0];
rz(-0.027055351) q[0];
sx q[0];
rz(-1.5144801) q[0];
rz(-pi) q[1];
rz(0.55412519) q[2];
sx q[2];
rz(-1.2078071) q[2];
sx q[2];
rz(0.83152321) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.63355061) q[1];
sx q[1];
rz(-1.7910069) q[1];
sx q[1];
rz(-1.0255662) q[1];
x q[2];
rz(-0.68287373) q[3];
sx q[3];
rz(-1.8337909) q[3];
sx q[3];
rz(2.9337286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1713193) q[2];
sx q[2];
rz(-3.1372941) q[2];
sx q[2];
rz(-2.8805736) q[2];
rz(-3.1018992) q[3];
sx q[3];
rz(-1.7429765) q[3];
sx q[3];
rz(-2.5980914) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71565851) q[0];
sx q[0];
rz(-1.4733227) q[0];
sx q[0];
rz(-2.4470827) q[0];
rz(2.014324) q[1];
sx q[1];
rz(-2.290461) q[1];
sx q[1];
rz(-2.1515501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8111525) q[0];
sx q[0];
rz(-1.9188723) q[0];
sx q[0];
rz(-1.483782) q[0];
rz(-pi) q[1];
rz(2.9431413) q[2];
sx q[2];
rz(-0.99764148) q[2];
sx q[2];
rz(2.1015374) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7496365) q[1];
sx q[1];
rz(-0.5881433) q[1];
sx q[1];
rz(-0.68444499) q[1];
rz(-0.95616266) q[3];
sx q[3];
rz(-0.80397881) q[3];
sx q[3];
rz(1.965919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29828829) q[2];
sx q[2];
rz(-2.1683606) q[2];
sx q[2];
rz(-2.6514371) q[2];
rz(2.7847024) q[3];
sx q[3];
rz(-2.7481952) q[3];
sx q[3];
rz(1.8753768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0855584) q[0];
sx q[0];
rz(-1.1857251) q[0];
sx q[0];
rz(-0.84247843) q[0];
rz(-0.8017686) q[1];
sx q[1];
rz(-0.27920488) q[1];
sx q[1];
rz(-2.3559949) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2548837) q[0];
sx q[0];
rz(-2.2629316) q[0];
sx q[0];
rz(-1.9002302) q[0];
rz(-pi) q[1];
rz(0.71707861) q[2];
sx q[2];
rz(-1.684965) q[2];
sx q[2];
rz(-1.1952497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.848104) q[1];
sx q[1];
rz(-1.2624718) q[1];
sx q[1];
rz(2.5749194) q[1];
rz(0.81964747) q[3];
sx q[3];
rz(-1.420212) q[3];
sx q[3];
rz(-2.8898847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0164612) q[2];
sx q[2];
rz(-0.74840122) q[2];
sx q[2];
rz(-1.008519) q[2];
rz(0.85159167) q[3];
sx q[3];
rz(-0.75751704) q[3];
sx q[3];
rz(2.0612702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9652902) q[0];
sx q[0];
rz(-0.98133123) q[0];
sx q[0];
rz(-1.0705795) q[0];
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
rz(0.010653822) q[0];
sx q[0];
rz(-1.21216) q[0];
sx q[0];
rz(0.50827311) q[0];
rz(-pi) q[1];
rz(-1.3959051) q[2];
sx q[2];
rz(-1.5064459) q[2];
sx q[2];
rz(-1.627587) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.7091211) q[1];
sx q[1];
rz(-2.2975031) q[1];
sx q[1];
rz(2.5120887) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9845731) q[3];
sx q[3];
rz(-2.5738705) q[3];
sx q[3];
rz(-2.1322676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.30194482) q[2];
sx q[2];
rz(-3.0366615) q[2];
sx q[2];
rz(-1.5860175) q[2];
rz(-2.1258866) q[3];
sx q[3];
rz(-2.000122) q[3];
sx q[3];
rz(-1.4122081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3509336) q[0];
sx q[0];
rz(-2.9582294) q[0];
sx q[0];
rz(-0.80107981) q[0];
rz(0.13310295) q[1];
sx q[1];
rz(-1.8201273) q[1];
sx q[1];
rz(-0.4020234) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19045705) q[0];
sx q[0];
rz(-1.2371089) q[0];
sx q[0];
rz(0.78898375) q[0];
rz(-1.2799995) q[2];
sx q[2];
rz(-1.5687402) q[2];
sx q[2];
rz(2.2711547) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.525156) q[1];
sx q[1];
rz(-1.0911988) q[1];
sx q[1];
rz(2.8308771) q[1];
rz(-pi) q[2];
rz(-1.4457213) q[3];
sx q[3];
rz(-2.2251724) q[3];
sx q[3];
rz(0.47548166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5037527) q[2];
sx q[2];
rz(-1.4053586) q[2];
sx q[2];
rz(-2.3657738) q[2];
rz(-1.4555629) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(1.0078526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(3.0966454) q[0];
sx q[0];
rz(-2.2438887) q[0];
sx q[0];
rz(-0.67895472) q[0];
rz(1.8567765) q[1];
sx q[1];
rz(-0.7130475) q[1];
sx q[1];
rz(1.7624034) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64039907) q[0];
sx q[0];
rz(-1.7441347) q[0];
sx q[0];
rz(-0.79357432) q[0];
rz(-pi) q[1];
rz(3.0581362) q[2];
sx q[2];
rz(-1.2922568) q[2];
sx q[2];
rz(-2.4065774) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0022565) q[1];
sx q[1];
rz(-0.69076194) q[1];
sx q[1];
rz(2.407807) q[1];
rz(-1.8399162) q[3];
sx q[3];
rz(-2.9342071) q[3];
sx q[3];
rz(-0.086750448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8046367) q[2];
sx q[2];
rz(-2.3865484) q[2];
sx q[2];
rz(-1.9026559) q[2];
rz(1.3321446) q[3];
sx q[3];
rz(-0.76056162) q[3];
sx q[3];
rz(0.24470394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181344) q[0];
sx q[0];
rz(-1.9754388) q[0];
sx q[0];
rz(0.13846692) q[0];
rz(-0.18374099) q[1];
sx q[1];
rz(-2.467149) q[1];
sx q[1];
rz(2.8750681) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040695) q[0];
sx q[0];
rz(-2.7752987) q[0];
sx q[0];
rz(1.1031325) q[0];
rz(-pi) q[1];
rz(-2.8710352) q[2];
sx q[2];
rz(-1.0204401) q[2];
sx q[2];
rz(1.5591043) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3764972) q[1];
sx q[1];
rz(-1.1516715) q[1];
sx q[1];
rz(3.0800099) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4146904) q[3];
sx q[3];
rz(-1.57734) q[3];
sx q[3];
rz(0.18101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0900241) q[2];
sx q[2];
rz(-1.8847621) q[2];
sx q[2];
rz(-2.9070692) q[2];
rz(1.6656434) q[3];
sx q[3];
rz(-1.6022976) q[3];
sx q[3];
rz(2.3852111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33673564) q[0];
sx q[0];
rz(-0.8661626) q[0];
sx q[0];
rz(-2.3418703) q[0];
rz(0.58894482) q[1];
sx q[1];
rz(-2.2942693) q[1];
sx q[1];
rz(-2.934093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8556535) q[0];
sx q[0];
rz(-2.7253599) q[0];
sx q[0];
rz(0.32815423) q[0];
rz(-pi) q[1];
rz(0.3090119) q[2];
sx q[2];
rz(-2.2591619) q[2];
sx q[2];
rz(2.7613044) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.55940699) q[1];
sx q[1];
rz(-1.9684125) q[1];
sx q[1];
rz(-2.29672) q[1];
rz(-pi) q[2];
rz(1.1010246) q[3];
sx q[3];
rz(-1.5964526) q[3];
sx q[3];
rz(-0.86650833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80692446) q[2];
sx q[2];
rz(-0.88403264) q[2];
sx q[2];
rz(0.36726382) q[2];
rz(-1.4372829) q[3];
sx q[3];
rz(-1.1742914) q[3];
sx q[3];
rz(1.1737163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87840286) q[0];
sx q[0];
rz(-1.0028361) q[0];
sx q[0];
rz(-0.81018418) q[0];
rz(-2.0267678) q[1];
sx q[1];
rz(-2.7572542) q[1];
sx q[1];
rz(-0.55327639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5779752) q[0];
sx q[0];
rz(-1.6262691) q[0];
sx q[0];
rz(-0.93240057) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79375879) q[2];
sx q[2];
rz(-1.1053893) q[2];
sx q[2];
rz(-0.88875729) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7825494) q[1];
sx q[1];
rz(-1.7573253) q[1];
sx q[1];
rz(-1.2839509) q[1];
rz(0.35240473) q[3];
sx q[3];
rz(-0.81768546) q[3];
sx q[3];
rz(-1.6116097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0193923) q[2];
sx q[2];
rz(-0.70851749) q[2];
sx q[2];
rz(-0.23966399) q[2];
rz(-1.8137118) q[3];
sx q[3];
rz(-2.0769104) q[3];
sx q[3];
rz(-2.6740668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0620621) q[0];
sx q[0];
rz(-1.7876328) q[0];
sx q[0];
rz(-2.7476516) q[0];
rz(2.486034) q[1];
sx q[1];
rz(-1.9656904) q[1];
sx q[1];
rz(0.94737731) q[1];
rz(-1.2389567) q[2];
sx q[2];
rz(-1.9658026) q[2];
sx q[2];
rz(2.6323242) q[2];
rz(-0.37674698) q[3];
sx q[3];
rz(-2.527311) q[3];
sx q[3];
rz(-0.75025076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
