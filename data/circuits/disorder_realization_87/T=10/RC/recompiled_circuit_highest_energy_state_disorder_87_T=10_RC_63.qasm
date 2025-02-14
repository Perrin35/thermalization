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
rz(1.4662161) q[0];
sx q[0];
rz(5.3386547) q[0];
sx q[0];
rz(9.6261779) q[0];
rz(-2.0693076) q[1];
sx q[1];
rz(-2.3787002) q[1];
sx q[1];
rz(1.720517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3823377) q[0];
sx q[0];
rz(-1.6256285) q[0];
sx q[0];
rz(-2.1817529) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1915053) q[2];
sx q[2];
rz(-0.7751152) q[2];
sx q[2];
rz(0.24573869) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0396467) q[1];
sx q[1];
rz(-1.1803487) q[1];
sx q[1];
rz(-2.4389308) q[1];
x q[2];
rz(-1.9626612) q[3];
sx q[3];
rz(-1.8800991) q[3];
sx q[3];
rz(0.40180692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.034915514) q[2];
sx q[2];
rz(-2.3372529) q[2];
sx q[2];
rz(0.2790645) q[2];
rz(1.2532824) q[3];
sx q[3];
rz(-2.8918355) q[3];
sx q[3];
rz(-0.15160027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6107553) q[0];
sx q[0];
rz(-2.294367) q[0];
sx q[0];
rz(2.8952059) q[0];
rz(-1.9465744) q[1];
sx q[1];
rz(-0.89391005) q[1];
sx q[1];
rz(-1.5066719) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55913299) q[0];
sx q[0];
rz(-1.8163067) q[0];
sx q[0];
rz(-2.0364291) q[0];
rz(-pi) q[1];
rz(0.042889281) q[2];
sx q[2];
rz(-1.873462) q[2];
sx q[2];
rz(1.8406701) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8283702) q[1];
sx q[1];
rz(-1.5130318) q[1];
sx q[1];
rz(1.934113) q[1];
rz(-pi) q[2];
rz(0.89119567) q[3];
sx q[3];
rz(-2.5053484) q[3];
sx q[3];
rz(-1.4435857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7404777) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(1.9096036) q[2];
rz(-1.1022107) q[3];
sx q[3];
rz(-0.004318459) q[3];
sx q[3];
rz(1.3365041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44322893) q[0];
sx q[0];
rz(-0.6969499) q[0];
sx q[0];
rz(-2.2337636) q[0];
rz(2.5746131) q[1];
sx q[1];
rz(-0.98273977) q[1];
sx q[1];
rz(-0.10279113) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11353569) q[0];
sx q[0];
rz(-2.5350219) q[0];
sx q[0];
rz(2.7675178) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0605761) q[2];
sx q[2];
rz(-1.4668494) q[2];
sx q[2];
rz(-2.8250717) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72702209) q[1];
sx q[1];
rz(-0.99491454) q[1];
sx q[1];
rz(-2.058279) q[1];
rz(-pi) q[2];
rz(2.1249173) q[3];
sx q[3];
rz(-1.301487) q[3];
sx q[3];
rz(-0.96162063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28321442) q[2];
sx q[2];
rz(-0.80867043) q[2];
sx q[2];
rz(1.7041448) q[2];
rz(-2.7490859) q[3];
sx q[3];
rz(-2.0065353) q[3];
sx q[3];
rz(-2.8431622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0482386) q[0];
sx q[0];
rz(-1.3627351) q[0];
sx q[0];
rz(1.9527973) q[0];
rz(-2.8554754) q[1];
sx q[1];
rz(-2.5240099) q[1];
sx q[1];
rz(1.5123222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4453569) q[0];
sx q[0];
rz(-2.0907479) q[0];
sx q[0];
rz(-0.32132863) q[0];
x q[1];
rz(0.76750038) q[2];
sx q[2];
rz(-1.8165117) q[2];
sx q[2];
rz(-1.9227288) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5315631) q[1];
sx q[1];
rz(-1.7563142) q[1];
sx q[1];
rz(2.1712135) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0779987) q[3];
sx q[3];
rz(-1.385687) q[3];
sx q[3];
rz(-2.0349658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35147038) q[2];
sx q[2];
rz(-2.0581547) q[2];
sx q[2];
rz(-0.67698395) q[2];
rz(1.8949159) q[3];
sx q[3];
rz(-1.8691984) q[3];
sx q[3];
rz(1.6251132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2077654) q[0];
sx q[0];
rz(-0.25660577) q[0];
sx q[0];
rz(-3.1166792) q[0];
rz(-2.086153) q[1];
sx q[1];
rz(-0.98506227) q[1];
sx q[1];
rz(2.8581462) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0830529) q[0];
sx q[0];
rz(-1.827936) q[0];
sx q[0];
rz(1.3015981) q[0];
x q[1];
rz(-1.7666498) q[2];
sx q[2];
rz(-1.699866) q[2];
sx q[2];
rz(-1.7650676) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39549669) q[1];
sx q[1];
rz(-2.0714134) q[1];
sx q[1];
rz(0.034828111) q[1];
rz(-2.2437566) q[3];
sx q[3];
rz(-1.3011609) q[3];
sx q[3];
rz(1.2854888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31263605) q[2];
sx q[2];
rz(-0.41746155) q[2];
sx q[2];
rz(-2.8157595) q[2];
rz(1.8587941) q[3];
sx q[3];
rz(-1.157434) q[3];
sx q[3];
rz(-1.5475984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40389898) q[0];
sx q[0];
rz(-1.6588545) q[0];
sx q[0];
rz(-2.7666336) q[0];
rz(-0.47863475) q[1];
sx q[1];
rz(-0.67239434) q[1];
sx q[1];
rz(-2.7714444) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3309114) q[0];
sx q[0];
rz(-0.025891455) q[0];
sx q[0];
rz(0.46790345) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0022718) q[2];
sx q[2];
rz(-1.754154) q[2];
sx q[2];
rz(-2.574711) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9795591) q[1];
sx q[1];
rz(-2.0311072) q[1];
sx q[1];
rz(-2.3737337) q[1];
rz(-2.7137854) q[3];
sx q[3];
rz(-1.5662346) q[3];
sx q[3];
rz(0.13739861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39773539) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(1.3810623) q[2];
rz(2.0928275) q[3];
sx q[3];
rz(-1.3449113) q[3];
sx q[3];
rz(0.63961187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8884856) q[0];
sx q[0];
rz(-2.3024004) q[0];
sx q[0];
rz(0.92415586) q[0];
rz(2.2249075) q[1];
sx q[1];
rz(-1.0662096) q[1];
sx q[1];
rz(-1.7402657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2271254) q[0];
sx q[0];
rz(-2.0642529) q[0];
sx q[0];
rz(1.6363793) q[0];
x q[1];
rz(1.7426874) q[2];
sx q[2];
rz(-2.2187833) q[2];
sx q[2];
rz(-1.6164219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1223334) q[1];
sx q[1];
rz(-0.75503765) q[1];
sx q[1];
rz(-2.9072059) q[1];
x q[2];
rz(-1.1772778) q[3];
sx q[3];
rz(-0.52214185) q[3];
sx q[3];
rz(2.4260184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9132729) q[2];
sx q[2];
rz(-1.6797804) q[2];
sx q[2];
rz(1.9165967) q[2];
rz(3.1332968) q[3];
sx q[3];
rz(-0.010992916) q[3];
sx q[3];
rz(-0.86709658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5834354) q[0];
sx q[0];
rz(-0.83034101) q[0];
sx q[0];
rz(-0.48026568) q[0];
rz(-0.14446124) q[1];
sx q[1];
rz(-2.2339349) q[1];
sx q[1];
rz(-2.2772148) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7378993) q[0];
sx q[0];
rz(-1.4072352) q[0];
sx q[0];
rz(3.0944194) q[0];
x q[1];
rz(1.6046674) q[2];
sx q[2];
rz(-1.665417) q[2];
sx q[2];
rz(2.692344) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56039366) q[1];
sx q[1];
rz(-2.542109) q[1];
sx q[1];
rz(-1.1349212) q[1];
rz(-1.2914574) q[3];
sx q[3];
rz(-1.8350826) q[3];
sx q[3];
rz(3.0077028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20624837) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(2.5065191) q[2];
rz(2.5687929) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(2.2662381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11005814) q[0];
sx q[0];
rz(-0.77339554) q[0];
sx q[0];
rz(0.12426678) q[0];
rz(0.36525137) q[1];
sx q[1];
rz(-1.4271913) q[1];
sx q[1];
rz(-1.431538) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5968558) q[0];
sx q[0];
rz(-1.7628306) q[0];
sx q[0];
rz(-1.3903244) q[0];
rz(-2.4489787) q[2];
sx q[2];
rz(-2.358105) q[2];
sx q[2];
rz(-0.88800511) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6811583) q[1];
sx q[1];
rz(-1.2660813) q[1];
sx q[1];
rz(2.4439993) q[1];
x q[2];
rz(-0.52732738) q[3];
sx q[3];
rz(-1.9550712) q[3];
sx q[3];
rz(-1.7155855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8753836) q[2];
sx q[2];
rz(-2.2201846) q[2];
sx q[2];
rz(1.9688152) q[2];
rz(1.4069936) q[3];
sx q[3];
rz(-1.809779) q[3];
sx q[3];
rz(1.9269358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239546) q[0];
sx q[0];
rz(-2.0068491) q[0];
sx q[0];
rz(-0.59463516) q[0];
rz(2.1356964) q[1];
sx q[1];
rz(-1.2024095) q[1];
sx q[1];
rz(2.2524021) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62458663) q[0];
sx q[0];
rz(-1.7910302) q[0];
sx q[0];
rz(3.0896679) q[0];
x q[1];
rz(2.0311622) q[2];
sx q[2];
rz(-1.5074919) q[2];
sx q[2];
rz(-0.6748131) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2584302) q[1];
sx q[1];
rz(-3.0297802) q[1];
sx q[1];
rz(-1.0055786) q[1];
x q[2];
rz(-0.066915705) q[3];
sx q[3];
rz(-0.51735462) q[3];
sx q[3];
rz(-0.87393119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1903926) q[2];
sx q[2];
rz(-1.9659646) q[2];
sx q[2];
rz(-0.54538837) q[2];
rz(-2.178318) q[3];
sx q[3];
rz(-1.9656209) q[3];
sx q[3];
rz(-1.6776599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49190285) q[0];
sx q[0];
rz(-2.2408673) q[0];
sx q[0];
rz(1.3229205) q[0];
rz(2.2979965) q[1];
sx q[1];
rz(-1.0782764) q[1];
sx q[1];
rz(-1.2596399) q[1];
rz(-1.1314992) q[2];
sx q[2];
rz(-1.7514624) q[2];
sx q[2];
rz(-3.024586) q[2];
rz(-2.4424408) q[3];
sx q[3];
rz(-2.4300162) q[3];
sx q[3];
rz(3.1292849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
