OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3040721) q[0];
sx q[0];
rz(3.8068258) q[0];
sx q[0];
rz(8.1991631) q[0];
rz(2.4576814) q[1];
sx q[1];
rz(-0.39165762) q[1];
sx q[1];
rz(-2.439523) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7312429) q[0];
sx q[0];
rz(-0.87603891) q[0];
sx q[0];
rz(-1.724051) q[0];
rz(-pi) q[1];
rz(0.021356301) q[2];
sx q[2];
rz(-1.735677) q[2];
sx q[2];
rz(-0.0065553824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42618034) q[1];
sx q[1];
rz(-2.3857834) q[1];
sx q[1];
rz(0.73985696) q[1];
x q[2];
rz(-2.9443378) q[3];
sx q[3];
rz(-0.88867696) q[3];
sx q[3];
rz(-1.6113987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2012653) q[2];
sx q[2];
rz(-1.6028812) q[2];
sx q[2];
rz(-2.8931457) q[2];
rz(-2.0072319) q[3];
sx q[3];
rz(-3.0076707) q[3];
sx q[3];
rz(1.1321446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092875384) q[0];
sx q[0];
rz(-2.5889914) q[0];
sx q[0];
rz(-1.5485113) q[0];
rz(0.50367194) q[1];
sx q[1];
rz(-1.9182938) q[1];
sx q[1];
rz(0.37685397) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5232613) q[0];
sx q[0];
rz(-1.0987765) q[0];
sx q[0];
rz(1.5923772) q[0];
rz(-0.81266788) q[2];
sx q[2];
rz(-2.118131) q[2];
sx q[2];
rz(-3.0603882) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.406817) q[1];
sx q[1];
rz(-1.8207014) q[1];
sx q[1];
rz(-2.6491449) q[1];
rz(-pi) q[2];
rz(-2.8089588) q[3];
sx q[3];
rz(-1.6432228) q[3];
sx q[3];
rz(2.4410332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5571931) q[2];
sx q[2];
rz(-2.445745) q[2];
sx q[2];
rz(2.2759571) q[2];
rz(1.668476) q[3];
sx q[3];
rz(-0.98554635) q[3];
sx q[3];
rz(-2.1540811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5675548) q[0];
sx q[0];
rz(-1.8542629) q[0];
sx q[0];
rz(-1.0512742) q[0];
rz(-1.6846664) q[1];
sx q[1];
rz(-1.3070062) q[1];
sx q[1];
rz(-1.6251224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74836377) q[0];
sx q[0];
rz(-0.81384515) q[0];
sx q[0];
rz(-2.7500344) q[0];
rz(-pi) q[1];
rz(-2.2411738) q[2];
sx q[2];
rz(-1.0068276) q[2];
sx q[2];
rz(2.0869067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1976002) q[1];
sx q[1];
rz(-1.9514582) q[1];
sx q[1];
rz(-2.073111) q[1];
x q[2];
rz(3.0805796) q[3];
sx q[3];
rz(-1.5170238) q[3];
sx q[3];
rz(0.53053571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.39525825) q[2];
sx q[2];
rz(-1.0552152) q[2];
sx q[2];
rz(-0.22221097) q[2];
rz(-0.5018417) q[3];
sx q[3];
rz(-1.7343438) q[3];
sx q[3];
rz(-0.80880222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2489081) q[0];
sx q[0];
rz(-0.48766708) q[0];
sx q[0];
rz(1.1647613) q[0];
rz(2.248863) q[1];
sx q[1];
rz(-1.7612532) q[1];
sx q[1];
rz(-2.8597615) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4087021) q[0];
sx q[0];
rz(-2.6441649) q[0];
sx q[0];
rz(-0.25411782) q[0];
rz(-pi) q[1];
rz(-2.149942) q[2];
sx q[2];
rz(-0.94303149) q[2];
sx q[2];
rz(0.53774688) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.033867) q[1];
sx q[1];
rz(-1.733779) q[1];
sx q[1];
rz(3.0708583) q[1];
x q[2];
rz(0.025679703) q[3];
sx q[3];
rz(-1.1246944) q[3];
sx q[3];
rz(-3.1334973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0622327) q[2];
sx q[2];
rz(-2.5793109) q[2];
sx q[2];
rz(2.6083561) q[2];
rz(-2.0594635) q[3];
sx q[3];
rz(-0.60492587) q[3];
sx q[3];
rz(2.3281039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.37394062) q[0];
sx q[0];
rz(-0.39893183) q[0];
sx q[0];
rz(1.8178513) q[0];
rz(-2.3228877) q[1];
sx q[1];
rz(-2.3981514) q[1];
sx q[1];
rz(-2.0735819) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5214106) q[0];
sx q[0];
rz(-0.17818923) q[0];
sx q[0];
rz(-2.6784269) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85535149) q[2];
sx q[2];
rz(-2.2673439) q[2];
sx q[2];
rz(-1.0698505) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1562742) q[1];
sx q[1];
rz(-0.66412369) q[1];
sx q[1];
rz(2.9742624) q[1];
rz(1.8858484) q[3];
sx q[3];
rz(-2.6482411) q[3];
sx q[3];
rz(0.82786614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1800804) q[2];
sx q[2];
rz(-2.027812) q[2];
sx q[2];
rz(-2.2799344) q[2];
rz(-2.1152451) q[3];
sx q[3];
rz(-1.9828826) q[3];
sx q[3];
rz(1.1177184) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3468129) q[0];
sx q[0];
rz(-0.81951278) q[0];
sx q[0];
rz(-2.2487707) q[0];
rz(2.9606441) q[1];
sx q[1];
rz(-2.4632958) q[1];
sx q[1];
rz(-0.13042626) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.82662) q[0];
sx q[0];
rz(-0.19222799) q[0];
sx q[0];
rz(-2.1703224) q[0];
x q[1];
rz(-0.67306913) q[2];
sx q[2];
rz(-2.4166738) q[2];
sx q[2];
rz(2.6672305) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1219668) q[1];
sx q[1];
rz(-1.5415915) q[1];
sx q[1];
rz(-1.8080416) q[1];
x q[2];
rz(-1.3731558) q[3];
sx q[3];
rz(-0.79953558) q[3];
sx q[3];
rz(-1.3621735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9170407) q[2];
sx q[2];
rz(-1.9150534) q[2];
sx q[2];
rz(-0.87635931) q[2];
rz(-1.8996436) q[3];
sx q[3];
rz(-0.94049898) q[3];
sx q[3];
rz(-2.131264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7414311) q[0];
sx q[0];
rz(-0.09859666) q[0];
sx q[0];
rz(0.36369351) q[0];
rz(-2.3997276) q[1];
sx q[1];
rz(-2.1153085) q[1];
sx q[1];
rz(2.738764) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14430732) q[0];
sx q[0];
rz(-1.8107521) q[0];
sx q[0];
rz(1.1474987) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9697519) q[2];
sx q[2];
rz(-1.2201933) q[2];
sx q[2];
rz(2.1828841) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6322569) q[1];
sx q[1];
rz(-1.4493353) q[1];
sx q[1];
rz(0.077510133) q[1];
rz(-pi) q[2];
rz(-0.63490156) q[3];
sx q[3];
rz(-1.3506119) q[3];
sx q[3];
rz(1.6600773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7824629) q[2];
sx q[2];
rz(-0.368258) q[2];
sx q[2];
rz(1.5472319) q[2];
rz(2.3063229) q[3];
sx q[3];
rz(-2.1850977) q[3];
sx q[3];
rz(2.6529151) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35992026) q[0];
sx q[0];
rz(-1.3363573) q[0];
sx q[0];
rz(0.51505995) q[0];
rz(1.7022279) q[1];
sx q[1];
rz(-1.8481588) q[1];
sx q[1];
rz(-2.7424367) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8989152) q[0];
sx q[0];
rz(-1.5689881) q[0];
sx q[0];
rz(3.010672) q[0];
rz(-pi) q[1];
rz(-1.3005343) q[2];
sx q[2];
rz(-2.3944602) q[2];
sx q[2];
rz(3.0556553) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6281575) q[1];
sx q[1];
rz(-1.6353459) q[1];
sx q[1];
rz(1.8316395) q[1];
x q[2];
rz(1.5720444) q[3];
sx q[3];
rz(-1.5051418) q[3];
sx q[3];
rz(1.1923238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1601552) q[2];
sx q[2];
rz(-1.8835386) q[2];
sx q[2];
rz(1.0464279) q[2];
rz(-2.4753172) q[3];
sx q[3];
rz(-2.8847238) q[3];
sx q[3];
rz(1.0506786) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1995354) q[0];
sx q[0];
rz(-0.10534795) q[0];
sx q[0];
rz(-0.60607213) q[0];
rz(2.3145158) q[1];
sx q[1];
rz(-1.3304293) q[1];
sx q[1];
rz(1.3333295) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5904986) q[0];
sx q[0];
rz(-2.7685389) q[0];
sx q[0];
rz(-2.8540552) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.016388) q[2];
sx q[2];
rz(-0.87832172) q[2];
sx q[2];
rz(0.25633474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.575653) q[1];
sx q[1];
rz(-0.90653949) q[1];
sx q[1];
rz(-2.9931426) q[1];
x q[2];
rz(-0.17154947) q[3];
sx q[3];
rz(-0.4643054) q[3];
sx q[3];
rz(-2.7363079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0102319) q[2];
sx q[2];
rz(-2.2183245) q[2];
sx q[2];
rz(-0.38491797) q[2];
rz(-0.63079232) q[3];
sx q[3];
rz(-1.8600978) q[3];
sx q[3];
rz(-1.3425286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25403062) q[0];
sx q[0];
rz(-1.7429202) q[0];
sx q[0];
rz(-1.4400462) q[0];
rz(0.74132672) q[1];
sx q[1];
rz(-1.4011551) q[1];
sx q[1];
rz(-2.8828566) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1848768) q[0];
sx q[0];
rz(-1.8167082) q[0];
sx q[0];
rz(-3.1374404) q[0];
rz(2.0912295) q[2];
sx q[2];
rz(-2.5602617) q[2];
sx q[2];
rz(0.54301013) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8347296) q[1];
sx q[1];
rz(-2.266699) q[1];
sx q[1];
rz(0.45180288) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7388938) q[3];
sx q[3];
rz(-0.81953632) q[3];
sx q[3];
rz(2.0961026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4197293) q[2];
sx q[2];
rz(-2.0294919) q[2];
sx q[2];
rz(-1.3191684) q[2];
rz(-2.3223274) q[3];
sx q[3];
rz(-2.0115439) q[3];
sx q[3];
rz(2.5949902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7405613) q[0];
sx q[0];
rz(-2.2673829) q[0];
sx q[0];
rz(-0.54947214) q[0];
rz(-1.7026547) q[1];
sx q[1];
rz(-0.66822744) q[1];
sx q[1];
rz(-2.6034036) q[1];
rz(-1.5835252) q[2];
sx q[2];
rz(-2.4314778) q[2];
sx q[2];
rz(1.1927803) q[2];
rz(-2.8945558) q[3];
sx q[3];
rz(-2.1185736) q[3];
sx q[3];
rz(2.0251956) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
