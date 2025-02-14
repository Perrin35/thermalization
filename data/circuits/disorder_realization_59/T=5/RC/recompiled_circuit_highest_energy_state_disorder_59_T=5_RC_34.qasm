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
rz(0.94035971) q[0];
sx q[0];
rz(-2.0191963) q[0];
sx q[0];
rz(2.9474592) q[0];
rz(-5.1880646) q[1];
sx q[1];
rz(1.4580589) q[1];
sx q[1];
rz(10.169765) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632089) q[0];
sx q[0];
rz(-1.4683405) q[0];
sx q[0];
rz(-1.5200903) q[0];
rz(-pi) q[1];
rz(2.3294747) q[2];
sx q[2];
rz(-2.0491164) q[2];
sx q[2];
rz(1.1253192) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2304334) q[1];
sx q[1];
rz(-1.3285561) q[1];
sx q[1];
rz(1.6723067) q[1];
rz(-2.2017415) q[3];
sx q[3];
rz(-0.54169023) q[3];
sx q[3];
rz(1.1506789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8038883) q[2];
sx q[2];
rz(-1.8053728) q[2];
sx q[2];
rz(-1.814092) q[2];
rz(1.368807) q[3];
sx q[3];
rz(-0.33392206) q[3];
sx q[3];
rz(-0.22148618) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9889481) q[0];
sx q[0];
rz(-1.7247609) q[0];
sx q[0];
rz(3.0864571) q[0];
rz(-0.70471835) q[1];
sx q[1];
rz(-1.5780508) q[1];
sx q[1];
rz(1.0447186) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44061384) q[0];
sx q[0];
rz(-1.199628) q[0];
sx q[0];
rz(-1.2775902) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56684871) q[2];
sx q[2];
rz(-1.0613614) q[2];
sx q[2];
rz(-1.2896565) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4028496) q[1];
sx q[1];
rz(-2.301005) q[1];
sx q[1];
rz(2.5806246) q[1];
x q[2];
rz(-3.0468076) q[3];
sx q[3];
rz(-1.6322502) q[3];
sx q[3];
rz(-0.8549785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.81950554) q[2];
sx q[2];
rz(-0.6531859) q[2];
sx q[2];
rz(-1.3219924) q[2];
rz(2.3540438) q[3];
sx q[3];
rz(-1.4044263) q[3];
sx q[3];
rz(-1.0549841) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9287441) q[0];
sx q[0];
rz(-0.66900122) q[0];
sx q[0];
rz(0.71841946) q[0];
rz(2.8058715) q[1];
sx q[1];
rz(-1.6751143) q[1];
sx q[1];
rz(2.1885923) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7756726) q[0];
sx q[0];
rz(-1.5342511) q[0];
sx q[0];
rz(2.1120872) q[0];
rz(1.0510873) q[2];
sx q[2];
rz(-1.4428836) q[2];
sx q[2];
rz(-3.0232883) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.905639) q[1];
sx q[1];
rz(-2.3528026) q[1];
sx q[1];
rz(-2.1091288) q[1];
x q[2];
rz(0.6742874) q[3];
sx q[3];
rz(-0.73028781) q[3];
sx q[3];
rz(-1.8411318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.0042051729) q[2];
sx q[2];
rz(-2.0073828) q[2];
sx q[2];
rz(0.65286621) q[2];
rz(0.79814923) q[3];
sx q[3];
rz(-1.3099202) q[3];
sx q[3];
rz(-2.4979512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2586486) q[0];
sx q[0];
rz(-2.301321) q[0];
sx q[0];
rz(0.81664455) q[0];
rz(2.4500997) q[1];
sx q[1];
rz(-1.7934727) q[1];
sx q[1];
rz(-0.74849558) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2740606) q[0];
sx q[0];
rz(-0.36991773) q[0];
sx q[0];
rz(2.5648613) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3567793) q[2];
sx q[2];
rz(-1.4811918) q[2];
sx q[2];
rz(-0.90286189) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15259057) q[1];
sx q[1];
rz(-2.3034597) q[1];
sx q[1];
rz(0.92392606) q[1];
rz(-pi) q[2];
rz(-2.2159202) q[3];
sx q[3];
rz(-0.28680772) q[3];
sx q[3];
rz(1.4405026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7526492) q[2];
sx q[2];
rz(-2.3708673) q[2];
sx q[2];
rz(-0.79461092) q[2];
rz(1.7157582) q[3];
sx q[3];
rz(-2.1013997) q[3];
sx q[3];
rz(-2.7690601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.36885095) q[0];
sx q[0];
rz(-2.0617101) q[0];
sx q[0];
rz(2.6439731) q[0];
rz(-2.3720062) q[1];
sx q[1];
rz(-1.1520569) q[1];
sx q[1];
rz(2.41113) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.651938) q[0];
sx q[0];
rz(-1.3562855) q[0];
sx q[0];
rz(-2.4431321) q[0];
rz(-1.914417) q[2];
sx q[2];
rz(-2.4848652) q[2];
sx q[2];
rz(1.7467787) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11558696) q[1];
sx q[1];
rz(-0.77669826) q[1];
sx q[1];
rz(1.8724043) q[1];
rz(-pi) q[2];
rz(0.45017193) q[3];
sx q[3];
rz(-1.1266337) q[3];
sx q[3];
rz(-1.1314363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2490425) q[2];
sx q[2];
rz(-1.6232792) q[2];
sx q[2];
rz(-2.6882233) q[2];
rz(1.8306277) q[3];
sx q[3];
rz(-0.1736719) q[3];
sx q[3];
rz(0.12392509) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3675568) q[0];
sx q[0];
rz(-1.0834162) q[0];
sx q[0];
rz(1.3483082) q[0];
rz(1.096161) q[1];
sx q[1];
rz(-1.6588914) q[1];
sx q[1];
rz(2.5199264) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2638877) q[0];
sx q[0];
rz(-1.5549721) q[0];
sx q[0];
rz(-0.71431922) q[0];
rz(2.9078324) q[2];
sx q[2];
rz(-2.7359161) q[2];
sx q[2];
rz(-2.7004227) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7324896) q[1];
sx q[1];
rz(-1.81899) q[1];
sx q[1];
rz(0.77381882) q[1];
rz(-1.0249025) q[3];
sx q[3];
rz(-2.1612447) q[3];
sx q[3];
rz(1.5284554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2094476) q[2];
sx q[2];
rz(-0.79938447) q[2];
sx q[2];
rz(-1.4383593) q[2];
rz(-0.98073331) q[3];
sx q[3];
rz(-1.8756198) q[3];
sx q[3];
rz(-1.2119306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88406968) q[0];
sx q[0];
rz(-1.3781837) q[0];
sx q[0];
rz(0.32514611) q[0];
rz(0.50666058) q[1];
sx q[1];
rz(-1.0136565) q[1];
sx q[1];
rz(-1.315717) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36111212) q[0];
sx q[0];
rz(-0.69421116) q[0];
sx q[0];
rz(-1.3785465) q[0];
x q[1];
rz(2.1431214) q[2];
sx q[2];
rz(-1.1954167) q[2];
sx q[2];
rz(1.4294415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.084242227) q[1];
sx q[1];
rz(-2.2154543) q[1];
sx q[1];
rz(2.0688384) q[1];
x q[2];
rz(-0.042577199) q[3];
sx q[3];
rz(-2.2025962) q[3];
sx q[3];
rz(-2.7160983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.80960861) q[2];
sx q[2];
rz(-1.0386244) q[2];
sx q[2];
rz(1.5254376) q[2];
rz(-1.7841313) q[3];
sx q[3];
rz(-1.4493891) q[3];
sx q[3];
rz(1.6239032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9759393) q[0];
sx q[0];
rz(-2.0015367) q[0];
sx q[0];
rz(2.4923988) q[0];
rz(2.3340885) q[1];
sx q[1];
rz(-2.0048001) q[1];
sx q[1];
rz(1.753122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8313468) q[0];
sx q[0];
rz(-2.5355336) q[0];
sx q[0];
rz(-1.640727) q[0];
x q[1];
rz(1.1325324) q[2];
sx q[2];
rz(-2.8243941) q[2];
sx q[2];
rz(-2.4567921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2040223) q[1];
sx q[1];
rz(-1.3998242) q[1];
sx q[1];
rz(3.0730166) q[1];
rz(-pi) q[2];
rz(1.3598594) q[3];
sx q[3];
rz(-1.9550793) q[3];
sx q[3];
rz(-0.74712844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2041152) q[2];
sx q[2];
rz(-0.14675879) q[2];
sx q[2];
rz(2.333763) q[2];
rz(2.4701123) q[3];
sx q[3];
rz(-1.6356133) q[3];
sx q[3];
rz(-1.5988662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2272187) q[0];
sx q[0];
rz(-0.53633538) q[0];
sx q[0];
rz(0.45289034) q[0];
rz(-2.8462786) q[1];
sx q[1];
rz(-1.9120049) q[1];
sx q[1];
rz(-1.3989075) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76052651) q[0];
sx q[0];
rz(-0.84080836) q[0];
sx q[0];
rz(-3.0286519) q[0];
rz(-pi) q[1];
rz(-2.0844551) q[2];
sx q[2];
rz(-1.7073362) q[2];
sx q[2];
rz(1.650102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4537482) q[1];
sx q[1];
rz(-2.0216989) q[1];
sx q[1];
rz(-1.1503673) q[1];
rz(-pi) q[2];
x q[2];
rz(0.021965071) q[3];
sx q[3];
rz(-1.3736892) q[3];
sx q[3];
rz(-1.1363014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37111515) q[2];
sx q[2];
rz(-0.8728084) q[2];
sx q[2];
rz(2.5650043) q[2];
rz(0.21814403) q[3];
sx q[3];
rz(-2.348867) q[3];
sx q[3];
rz(1.920759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6169154) q[0];
sx q[0];
rz(-0.24743947) q[0];
sx q[0];
rz(3.0531378) q[0];
rz(2.6103919) q[1];
sx q[1];
rz(-0.4041268) q[1];
sx q[1];
rz(-1.6204576) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1955521) q[0];
sx q[0];
rz(-1.6541535) q[0];
sx q[0];
rz(1.1766737) q[0];
x q[1];
rz(-2.0481061) q[2];
sx q[2];
rz(-1.6132406) q[2];
sx q[2];
rz(1.0099271) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38177165) q[1];
sx q[1];
rz(-2.9050937) q[1];
sx q[1];
rz(0.31707724) q[1];
rz(0.57595595) q[3];
sx q[3];
rz(-0.68000845) q[3];
sx q[3];
rz(-0.81871225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29181186) q[2];
sx q[2];
rz(-1.7917289) q[2];
sx q[2];
rz(-1.527171) q[2];
rz(1.8725066) q[3];
sx q[3];
rz(-0.98617712) q[3];
sx q[3];
rz(-1.0482949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0328746) q[0];
sx q[0];
rz(-2.0357108) q[0];
sx q[0];
rz(2.3172814) q[0];
rz(0.4737919) q[1];
sx q[1];
rz(-1.5693249) q[1];
sx q[1];
rz(1.5798461) q[1];
rz(1.9971581) q[2];
sx q[2];
rz(-2.1254267) q[2];
sx q[2];
rz(2.8993901) q[2];
rz(0.80120091) q[3];
sx q[3];
rz(-0.73794919) q[3];
sx q[3];
rz(-0.59384624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
