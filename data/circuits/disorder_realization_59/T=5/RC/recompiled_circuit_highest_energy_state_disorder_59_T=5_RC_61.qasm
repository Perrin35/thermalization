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
rz(-2.2012329) q[0];
sx q[0];
rz(-1.1223963) q[0];
sx q[0];
rz(0.19413343) q[0];
rz(1.0951207) q[1];
sx q[1];
rz(-1.6835338) q[1];
sx q[1];
rz(-0.74498743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8391957) q[0];
sx q[0];
rz(-1.6212362) q[0];
sx q[0];
rz(0.10258672) q[0];
rz(-pi) q[1];
x q[1];
rz(2.216624) q[2];
sx q[2];
rz(-0.87069521) q[2];
sx q[2];
rz(0.89751228) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4575263) q[1];
sx q[1];
rz(-1.4722595) q[1];
sx q[1];
rz(-2.8981478) q[1];
rz(-pi) q[2];
rz(-1.1185094) q[3];
sx q[3];
rz(-1.2617526) q[3];
sx q[3];
rz(-0.13918258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8038883) q[2];
sx q[2];
rz(-1.8053728) q[2];
sx q[2];
rz(1.3275006) q[2];
rz(-1.368807) q[3];
sx q[3];
rz(-2.8076706) q[3];
sx q[3];
rz(-0.22148618) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9889481) q[0];
sx q[0];
rz(-1.4168318) q[0];
sx q[0];
rz(0.055135559) q[0];
rz(-2.4368743) q[1];
sx q[1];
rz(-1.5635419) q[1];
sx q[1];
rz(1.0447186) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44061384) q[0];
sx q[0];
rz(-1.199628) q[0];
sx q[0];
rz(1.2775902) q[0];
x q[1];
rz(-0.80514812) q[2];
sx q[2];
rz(-0.74294801) q[2];
sx q[2];
rz(2.2066903) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5767908) q[1];
sx q[1];
rz(-1.163244) q[1];
sx q[1];
rz(-0.75753404) q[1];
rz(1.509066) q[3];
sx q[3];
rz(-1.6654019) q[3];
sx q[3];
rz(-2.4199362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3220871) q[2];
sx q[2];
rz(-0.6531859) q[2];
sx q[2];
rz(1.8196003) q[2];
rz(-0.78754887) q[3];
sx q[3];
rz(-1.7371663) q[3];
sx q[3];
rz(-2.0866086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9287441) q[0];
sx q[0];
rz(-0.66900122) q[0];
sx q[0];
rz(-2.4231732) q[0];
rz(-0.33572117) q[1];
sx q[1];
rz(-1.6751143) q[1];
sx q[1];
rz(2.1885923) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18291423) q[0];
sx q[0];
rz(-2.1116858) q[0];
sx q[0];
rz(-3.0989585) q[0];
rz(-pi) q[1];
rz(1.8242055) q[2];
sx q[2];
rz(-2.6077787) q[2];
sx q[2];
rz(-1.6718503) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2359537) q[1];
sx q[1];
rz(-0.78879005) q[1];
sx q[1];
rz(-1.0324638) q[1];
x q[2];
rz(-0.6742874) q[3];
sx q[3];
rz(-2.4113048) q[3];
sx q[3];
rz(1.3004608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0042051729) q[2];
sx q[2];
rz(-2.0073828) q[2];
sx q[2];
rz(-2.4887264) q[2];
rz(0.79814923) q[3];
sx q[3];
rz(-1.3099202) q[3];
sx q[3];
rz(0.64364141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2586486) q[0];
sx q[0];
rz(-0.84027165) q[0];
sx q[0];
rz(2.3249481) q[0];
rz(2.4500997) q[1];
sx q[1];
rz(-1.34812) q[1];
sx q[1];
rz(-2.3930971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2484528) q[0];
sx q[0];
rz(-1.3723541) q[0];
sx q[0];
rz(-0.31427419) q[0];
rz(-1.9709936) q[2];
sx q[2];
rz(-0.23175254) q[2];
sx q[2];
rz(-1.0585001) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4481727) q[1];
sx q[1];
rz(-2.2058371) q[1];
sx q[1];
rz(-0.5902115) q[1];
rz(-pi) q[2];
rz(2.9660712) q[3];
sx q[3];
rz(-1.3427882) q[3];
sx q[3];
rz(1.0358178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3889435) q[2];
sx q[2];
rz(-0.77072531) q[2];
sx q[2];
rz(0.79461092) q[2];
rz(-1.7157582) q[3];
sx q[3];
rz(-2.1013997) q[3];
sx q[3];
rz(2.7690601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7727417) q[0];
sx q[0];
rz(-2.0617101) q[0];
sx q[0];
rz(-0.4976196) q[0];
rz(0.76958641) q[1];
sx q[1];
rz(-1.1520569) q[1];
sx q[1];
rz(-0.73046267) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.651938) q[0];
sx q[0];
rz(-1.7853072) q[0];
sx q[0];
rz(-0.69846054) q[0];
x q[1];
rz(0.25409296) q[2];
sx q[2];
rz(-2.1832018) q[2];
sx q[2];
rz(-1.8191561) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.11558696) q[1];
sx q[1];
rz(-2.3648944) q[1];
sx q[1];
rz(1.8724043) q[1];
rz(-pi) q[2];
rz(1.0845844) q[3];
sx q[3];
rz(-1.9745805) q[3];
sx q[3];
rz(2.9069834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89255014) q[2];
sx q[2];
rz(-1.6232792) q[2];
sx q[2];
rz(2.6882233) q[2];
rz(-1.8306277) q[3];
sx q[3];
rz(-2.9679208) q[3];
sx q[3];
rz(-3.0176676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740358) q[0];
sx q[0];
rz(-1.0834162) q[0];
sx q[0];
rz(1.7932844) q[0];
rz(-1.096161) q[1];
sx q[1];
rz(-1.6588914) q[1];
sx q[1];
rz(-2.5199264) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.877705) q[0];
sx q[0];
rz(-1.5866205) q[0];
sx q[0];
rz(-2.4272734) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7458199) q[2];
sx q[2];
rz(-1.4792551) q[2];
sx q[2];
rz(1.7965732) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.073879378) q[1];
sx q[1];
rz(-2.3150958) q[1];
sx q[1];
rz(1.2302876) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4765763) q[3];
sx q[3];
rz(-1.1248853) q[3];
sx q[3];
rz(0.2837821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2094476) q[2];
sx q[2];
rz(-2.3422082) q[2];
sx q[2];
rz(1.7032334) q[2];
rz(-0.98073331) q[3];
sx q[3];
rz(-1.8756198) q[3];
sx q[3];
rz(-1.2119306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.257523) q[0];
sx q[0];
rz(-1.7634089) q[0];
sx q[0];
rz(-0.32514611) q[0];
rz(2.6349321) q[1];
sx q[1];
rz(-2.1279361) q[1];
sx q[1];
rz(1.8258757) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0611826) q[0];
sx q[0];
rz(-1.4482486) q[0];
sx q[0];
rz(-2.2558801) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1998134) q[2];
sx q[2];
rz(-2.468839) q[2];
sx q[2];
rz(2.7655809) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81974492) q[1];
sx q[1];
rz(-0.79230601) q[1];
sx q[1];
rz(-0.56604947) q[1];
rz(1.6288888) q[3];
sx q[3];
rz(-0.63303715) q[3];
sx q[3];
rz(-0.35348383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.80960861) q[2];
sx q[2];
rz(-1.0386244) q[2];
sx q[2];
rz(1.6161551) q[2];
rz(1.7841313) q[3];
sx q[3];
rz(-1.4493891) q[3];
sx q[3];
rz(-1.6239032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1656533) q[0];
sx q[0];
rz(-1.140056) q[0];
sx q[0];
rz(0.64919382) q[0];
rz(0.80750418) q[1];
sx q[1];
rz(-1.1367926) q[1];
sx q[1];
rz(-1.3884707) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3102458) q[0];
sx q[0];
rz(-0.60605907) q[0];
sx q[0];
rz(1.5008657) q[0];
rz(-pi) q[1];
rz(0.138422) q[2];
sx q[2];
rz(-1.8571203) q[2];
sx q[2];
rz(0.22655205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5538881) q[1];
sx q[1];
rz(-2.9575051) q[1];
sx q[1];
rz(-1.192993) q[1];
x q[2];
rz(-1.3598594) q[3];
sx q[3];
rz(-1.1865133) q[3];
sx q[3];
rz(2.3944642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2041152) q[2];
sx q[2];
rz(-0.14675879) q[2];
sx q[2];
rz(0.80782962) q[2];
rz(2.4701123) q[3];
sx q[3];
rz(-1.5059794) q[3];
sx q[3];
rz(-1.5427264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9143739) q[0];
sx q[0];
rz(-0.53633538) q[0];
sx q[0];
rz(-0.45289034) q[0];
rz(-0.29531404) q[1];
sx q[1];
rz(-1.2295877) q[1];
sx q[1];
rz(1.7426851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2558289) q[0];
sx q[0];
rz(-1.654878) q[0];
sx q[0];
rz(-2.3039615) q[0];
x q[1];
rz(1.0571376) q[2];
sx q[2];
rz(-1.7073362) q[2];
sx q[2];
rz(1.650102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30946839) q[1];
sx q[1];
rz(-1.1946284) q[1];
sx q[1];
rz(2.6539567) q[1];
rz(-3.1196276) q[3];
sx q[3];
rz(-1.7679035) q[3];
sx q[3];
rz(-2.0052912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7704775) q[2];
sx q[2];
rz(-0.8728084) q[2];
sx q[2];
rz(2.5650043) q[2];
rz(-0.21814403) q[3];
sx q[3];
rz(-0.79272565) q[3];
sx q[3];
rz(-1.2208337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6169154) q[0];
sx q[0];
rz(-0.24743947) q[0];
sx q[0];
rz(-3.0531378) q[0];
rz(2.6103919) q[1];
sx q[1];
rz(-0.4041268) q[1];
sx q[1];
rz(-1.6204576) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5514497) q[0];
sx q[0];
rz(-1.1781173) q[0];
sx q[0];
rz(0.090242437) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4786104) q[2];
sx q[2];
rz(-0.47904821) q[2];
sx q[2];
rz(2.6625815) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38177165) q[1];
sx q[1];
rz(-2.9050937) q[1];
sx q[1];
rz(0.31707724) q[1];
rz(1.1559256) q[3];
sx q[3];
rz(-2.1262827) q[3];
sx q[3];
rz(-1.6270571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8497808) q[2];
sx q[2];
rz(-1.3498638) q[2];
sx q[2];
rz(1.6144217) q[2];
rz(1.269086) q[3];
sx q[3];
rz(-2.1554155) q[3];
sx q[3];
rz(2.0932978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.108718) q[0];
sx q[0];
rz(-1.1058818) q[0];
sx q[0];
rz(-0.82431128) q[0];
rz(-2.6678008) q[1];
sx q[1];
rz(-1.5693249) q[1];
sx q[1];
rz(1.5798461) q[1];
rz(-1.1444345) q[2];
sx q[2];
rz(-2.1254267) q[2];
sx q[2];
rz(2.8993901) q[2];
rz(-2.3403917) q[3];
sx q[3];
rz(-0.73794919) q[3];
sx q[3];
rz(-0.59384624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
