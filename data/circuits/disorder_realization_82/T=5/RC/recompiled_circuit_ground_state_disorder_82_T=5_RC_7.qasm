OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6717186) q[0];
sx q[0];
rz(-1.1550386) q[0];
sx q[0];
rz(-1.7116829) q[0];
rz(0.45541304) q[1];
sx q[1];
rz(5.6918511) q[1];
sx q[1];
rz(11.439352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10218638) q[0];
sx q[0];
rz(-1.8012905) q[0];
sx q[0];
rz(-3.1297157) q[0];
rz(-pi) q[1];
rz(-0.34268171) q[2];
sx q[2];
rz(-1.8215873) q[2];
sx q[2];
rz(-2.8042274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0309567) q[1];
sx q[1];
rz(-1.0070166) q[1];
sx q[1];
rz(-2.8085676) q[1];
rz(-pi) q[2];
rz(-1.5529049) q[3];
sx q[3];
rz(-1.3800637) q[3];
sx q[3];
rz(0.73303661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.11898018) q[2];
sx q[2];
rz(-0.70465124) q[2];
sx q[2];
rz(-1.653778) q[2];
rz(-0.37114272) q[3];
sx q[3];
rz(-1.1266339) q[3];
sx q[3];
rz(2.3625372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(0.30598518) q[0];
sx q[0];
rz(-1.7636517) q[0];
sx q[0];
rz(2.4775179) q[0];
rz(-0.56655073) q[1];
sx q[1];
rz(-1.605875) q[1];
sx q[1];
rz(2.9552592) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31748397) q[0];
sx q[0];
rz(-0.58360277) q[0];
sx q[0];
rz(1.9822796) q[0];
rz(0.26797023) q[2];
sx q[2];
rz(-1.9391962) q[2];
sx q[2];
rz(-3.1173988) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19263467) q[1];
sx q[1];
rz(-1.4073458) q[1];
sx q[1];
rz(-0.88731097) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6870895) q[3];
sx q[3];
rz(-2.6175559) q[3];
sx q[3];
rz(0.91479036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8191007) q[2];
sx q[2];
rz(-1.842061) q[2];
sx q[2];
rz(1.5619649) q[2];
rz(-1.9942079) q[3];
sx q[3];
rz(-2.8473144) q[3];
sx q[3];
rz(-1.3636205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086562432) q[0];
sx q[0];
rz(-1.6494305) q[0];
sx q[0];
rz(-0.30558875) q[0];
rz(-1.8606404) q[1];
sx q[1];
rz(-0.31550229) q[1];
sx q[1];
rz(-1.5234647) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.279583) q[0];
sx q[0];
rz(-1.209799) q[0];
sx q[0];
rz(-1.3230514) q[0];
x q[1];
rz(-2.8144437) q[2];
sx q[2];
rz(-0.047657813) q[2];
sx q[2];
rz(0.41027712) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0022920575) q[1];
sx q[1];
rz(-1.2131547) q[1];
sx q[1];
rz(0.73441084) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9039246) q[3];
sx q[3];
rz(-2.1690627) q[3];
sx q[3];
rz(-0.054758398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.65172255) q[2];
sx q[2];
rz(-1.8946596) q[2];
sx q[2];
rz(-0.22001246) q[2];
rz(-3.0618727) q[3];
sx q[3];
rz(-0.97697512) q[3];
sx q[3];
rz(-0.93130934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43059573) q[0];
sx q[0];
rz(-1.7544704) q[0];
sx q[0];
rz(0.0084477607) q[0];
rz(1.1620109) q[1];
sx q[1];
rz(-1.9879257) q[1];
sx q[1];
rz(-1.1147503) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82915562) q[0];
sx q[0];
rz(-0.36530802) q[0];
sx q[0];
rz(1.0201973) q[0];
x q[1];
rz(2.3643119) q[2];
sx q[2];
rz(-2.4156164) q[2];
sx q[2];
rz(2.7216146) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9544008) q[1];
sx q[1];
rz(-1.2174509) q[1];
sx q[1];
rz(-0.29275972) q[1];
rz(2.7211748) q[3];
sx q[3];
rz(-1.1913931) q[3];
sx q[3];
rz(2.8914352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4343425) q[2];
sx q[2];
rz(-1.0184526) q[2];
sx q[2];
rz(2.2960091) q[2];
rz(0.25104684) q[3];
sx q[3];
rz(-1.2716525) q[3];
sx q[3];
rz(-0.74560753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71395981) q[0];
sx q[0];
rz(-0.41249713) q[0];
sx q[0];
rz(1.0300256) q[0];
rz(2.5469942) q[1];
sx q[1];
rz(-1.4232891) q[1];
sx q[1];
rz(1.3535708) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6973482) q[0];
sx q[0];
rz(-1.5500907) q[0];
sx q[0];
rz(-0.022788825) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6313303) q[2];
sx q[2];
rz(-2.7760421) q[2];
sx q[2];
rz(-0.22399677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.36921484) q[1];
sx q[1];
rz(-2.2352043) q[1];
sx q[1];
rz(-2.496705) q[1];
rz(-pi) q[2];
rz(2.5491675) q[3];
sx q[3];
rz(-2.0090143) q[3];
sx q[3];
rz(-0.93417227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2305962) q[2];
sx q[2];
rz(-1.5065008) q[2];
sx q[2];
rz(-0.050749151) q[2];
rz(2.0283608) q[3];
sx q[3];
rz(-1.166393) q[3];
sx q[3];
rz(-0.39065233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0473061) q[0];
sx q[0];
rz(-2.081649) q[0];
sx q[0];
rz(-2.768709) q[0];
rz(-1.8376384) q[1];
sx q[1];
rz(-1.2645489) q[1];
sx q[1];
rz(-2.1162927) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2706755) q[0];
sx q[0];
rz(-1.1216721) q[0];
sx q[0];
rz(1.6555519) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2405459) q[2];
sx q[2];
rz(-1.7269968) q[2];
sx q[2];
rz(-1.3845342) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.655788) q[1];
sx q[1];
rz(-0.88690573) q[1];
sx q[1];
rz(0.25836103) q[1];
x q[2];
rz(-0.22475524) q[3];
sx q[3];
rz(-1.9850176) q[3];
sx q[3];
rz(-2.1638768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9103526) q[2];
sx q[2];
rz(-0.19008907) q[2];
sx q[2];
rz(-1.7115889) q[2];
rz(-1.6010239) q[3];
sx q[3];
rz(-0.68683306) q[3];
sx q[3];
rz(1.471126) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1533399) q[0];
sx q[0];
rz(-1.8415469) q[0];
sx q[0];
rz(-0.43149313) q[0];
rz(-1.2639812) q[1];
sx q[1];
rz(-1.5411721) q[1];
sx q[1];
rz(-0.65013179) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9857441) q[0];
sx q[0];
rz(-1.2835644) q[0];
sx q[0];
rz(1.4850856) q[0];
x q[1];
rz(-1.442914) q[2];
sx q[2];
rz(-0.68727109) q[2];
sx q[2];
rz(-2.5812347) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3172315) q[1];
sx q[1];
rz(-1.4556093) q[1];
sx q[1];
rz(0.86800602) q[1];
rz(-0.007644848) q[3];
sx q[3];
rz(-2.0286603) q[3];
sx q[3];
rz(2.1107822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.55107895) q[2];
sx q[2];
rz(-1.7503909) q[2];
sx q[2];
rz(-3.1371269) q[2];
rz(0.0095857754) q[3];
sx q[3];
rz(-1.8066581) q[3];
sx q[3];
rz(-0.34346223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15653217) q[0];
sx q[0];
rz(-2.8592906) q[0];
sx q[0];
rz(-2.9811133) q[0];
rz(1.7566682) q[1];
sx q[1];
rz(-2.3424708) q[1];
sx q[1];
rz(2.8327732) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64381448) q[0];
sx q[0];
rz(-2.9165932) q[0];
sx q[0];
rz(2.1589222) q[0];
x q[1];
rz(1.1868434) q[2];
sx q[2];
rz(-0.89659474) q[2];
sx q[2];
rz(0.37146682) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17559472) q[1];
sx q[1];
rz(-1.1958937) q[1];
sx q[1];
rz(-2.7213827) q[1];
x q[2];
rz(2.4229523) q[3];
sx q[3];
rz(-2.6860533) q[3];
sx q[3];
rz(0.18265858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.34362346) q[2];
sx q[2];
rz(-2.5324731) q[2];
sx q[2];
rz(-1.5999751) q[2];
rz(0.68495098) q[3];
sx q[3];
rz(-1.5545605) q[3];
sx q[3];
rz(-1.5233013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5110382) q[0];
sx q[0];
rz(-1.4383974) q[0];
sx q[0];
rz(-0.55139971) q[0];
rz(-2.664387) q[1];
sx q[1];
rz(-2.2255032) q[1];
sx q[1];
rz(2.3068857) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39838185) q[0];
sx q[0];
rz(-2.4373567) q[0];
sx q[0];
rz(2.387638) q[0];
rz(-1.536252) q[2];
sx q[2];
rz(-0.87022034) q[2];
sx q[2];
rz(-2.3242484) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7013423) q[1];
sx q[1];
rz(-1.6052488) q[1];
sx q[1];
rz(-0.95878307) q[1];
rz(-0.96179295) q[3];
sx q[3];
rz(-1.0474221) q[3];
sx q[3];
rz(0.23853049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4386374) q[2];
sx q[2];
rz(-2.1795887) q[2];
sx q[2];
rz(0.23923242) q[2];
rz(2.8171825) q[3];
sx q[3];
rz(-1.7041465) q[3];
sx q[3];
rz(-1.5391866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6259554) q[0];
sx q[0];
rz(-0.13787585) q[0];
sx q[0];
rz(0.71664083) q[0];
rz(-0.58249885) q[1];
sx q[1];
rz(-1.3526252) q[1];
sx q[1];
rz(2.8315721) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9649565) q[0];
sx q[0];
rz(-0.88439098) q[0];
sx q[0];
rz(2.6841037) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7298195) q[2];
sx q[2];
rz(-1.9343107) q[2];
sx q[2];
rz(2.1902254) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1056052) q[1];
sx q[1];
rz(-1.2140555) q[1];
sx q[1];
rz(-0.8143266) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4118092) q[3];
sx q[3];
rz(-2.457629) q[3];
sx q[3];
rz(-1.7548949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6325355) q[2];
sx q[2];
rz(-2.5639503) q[2];
sx q[2];
rz(1.7217815) q[2];
rz(1.4451197) q[3];
sx q[3];
rz(-2.4308379) q[3];
sx q[3];
rz(0.76326171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2390908) q[0];
sx q[0];
rz(-2.1949174) q[0];
sx q[0];
rz(1.3876023) q[0];
rz(1.4801964) q[1];
sx q[1];
rz(-1.4085242) q[1];
sx q[1];
rz(-2.9396802) q[1];
rz(0.26726503) q[2];
sx q[2];
rz(-1.6435192) q[2];
sx q[2];
rz(-2.7036713) q[2];
rz(0.46764163) q[3];
sx q[3];
rz(-2.1976039) q[3];
sx q[3];
rz(-0.80898367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
