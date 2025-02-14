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
rz(2.6994045) q[0];
sx q[0];
rz(-1.9884041) q[0];
sx q[0];
rz(6.4118275) q[0];
rz(-4.6075912) q[1];
sx q[1];
rz(1.0990376) q[1];
sx q[1];
rz(8.3429835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2483467) q[0];
sx q[0];
rz(-1.5061646) q[0];
sx q[0];
rz(-1.5399702) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9456909) q[2];
sx q[2];
rz(-1.8842976) q[2];
sx q[2];
rz(1.9913265) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0073407081) q[1];
sx q[1];
rz(-1.6414343) q[1];
sx q[1];
rz(-1.411648) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2461189) q[3];
sx q[3];
rz(-2.8522308) q[3];
sx q[3];
rz(1.0288256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.094540207) q[2];
sx q[2];
rz(-0.69539842) q[2];
sx q[2];
rz(2.409234) q[2];
rz(-3.0426466) q[3];
sx q[3];
rz(-2.1061335) q[3];
sx q[3];
rz(-1.8100479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-3.052213) q[0];
sx q[0];
rz(-2.3591924) q[0];
sx q[0];
rz(-0.041444929) q[0];
rz(-2.4502358) q[1];
sx q[1];
rz(-1.3203011) q[1];
sx q[1];
rz(-0.88821214) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23032886) q[0];
sx q[0];
rz(-2.3155876) q[0];
sx q[0];
rz(2.2780096) q[0];
x q[1];
rz(2.4629243) q[2];
sx q[2];
rz(-1.1801085) q[2];
sx q[2];
rz(0.44579577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8687369) q[1];
sx q[1];
rz(-1.6140198) q[1];
sx q[1];
rz(0.24300139) q[1];
rz(-pi) q[2];
x q[2];
rz(2.668374) q[3];
sx q[3];
rz(-2.5788973) q[3];
sx q[3];
rz(-1.109888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9899675) q[2];
sx q[2];
rz(-1.3902367) q[2];
sx q[2];
rz(-1.3803231) q[2];
rz(-0.049792854) q[3];
sx q[3];
rz(-0.68792611) q[3];
sx q[3];
rz(0.5602347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3675073) q[0];
sx q[0];
rz(-0.15223509) q[0];
sx q[0];
rz(1.2138858) q[0];
rz(-1.3146575) q[1];
sx q[1];
rz(-2.1036802) q[1];
sx q[1];
rz(-1.5705869) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6187046) q[0];
sx q[0];
rz(-1.9490593) q[0];
sx q[0];
rz(1.8357651) q[0];
rz(-pi) q[1];
rz(-2.7309787) q[2];
sx q[2];
rz(-1.4272235) q[2];
sx q[2];
rz(-2.5559354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.27109001) q[1];
sx q[1];
rz(-1.6047641) q[1];
sx q[1];
rz(1.7351331) q[1];
rz(-pi) q[2];
rz(0.79292983) q[3];
sx q[3];
rz(-0.76628162) q[3];
sx q[3];
rz(-2.7310581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35459241) q[2];
sx q[2];
rz(-0.74942333) q[2];
sx q[2];
rz(2.836239) q[2];
rz(-0.8461771) q[3];
sx q[3];
rz(-1.927522) q[3];
sx q[3];
rz(-2.2727216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8531891) q[0];
sx q[0];
rz(-1.0193634) q[0];
sx q[0];
rz(2.1715721) q[0];
rz(0.24523973) q[1];
sx q[1];
rz(-1.0673362) q[1];
sx q[1];
rz(1.6397569) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062881447) q[0];
sx q[0];
rz(-0.20096603) q[0];
sx q[0];
rz(-0.68931343) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9059783) q[2];
sx q[2];
rz(-1.4019626) q[2];
sx q[2];
rz(1.1323358) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8322595) q[1];
sx q[1];
rz(-2.8974779) q[1];
sx q[1];
rz(-2.4629618) q[1];
rz(-pi) q[2];
rz(2.7795622) q[3];
sx q[3];
rz(-0.76043441) q[3];
sx q[3];
rz(2.8575051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5331427) q[2];
sx q[2];
rz(-1.5641944) q[2];
sx q[2];
rz(-1.2987785) q[2];
rz(-0.53457824) q[3];
sx q[3];
rz(-2.5375073) q[3];
sx q[3];
rz(2.4728313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088148549) q[0];
sx q[0];
rz(-1.7867333) q[0];
sx q[0];
rz(0.97766367) q[0];
rz(2.4560302) q[1];
sx q[1];
rz(-0.79427636) q[1];
sx q[1];
rz(0.96644863) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0880497) q[0];
sx q[0];
rz(-2.2849963) q[0];
sx q[0];
rz(-0.62936737) q[0];
x q[1];
rz(-0.56880237) q[2];
sx q[2];
rz(-1.9380271) q[2];
sx q[2];
rz(0.18062225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9271665) q[1];
sx q[1];
rz(-2.4220903) q[1];
sx q[1];
rz(1.7822766) q[1];
x q[2];
rz(0.54040106) q[3];
sx q[3];
rz(-0.62190074) q[3];
sx q[3];
rz(0.87582138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.55296772) q[2];
sx q[2];
rz(-1.3543465) q[2];
sx q[2];
rz(2.8589613) q[2];
rz(-2.1770554) q[3];
sx q[3];
rz(-1.5805565) q[3];
sx q[3];
rz(-2.7717822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10837567) q[0];
sx q[0];
rz(-2.7300457) q[0];
sx q[0];
rz(1.0785528) q[0];
rz(0.82303965) q[1];
sx q[1];
rz(-1.1230725) q[1];
sx q[1];
rz(-0.80026921) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9554515) q[0];
sx q[0];
rz(-1.5622883) q[0];
sx q[0];
rz(-1.5415989) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8326464) q[2];
sx q[2];
rz(-1.8712723) q[2];
sx q[2];
rz(2.7555675) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9335223) q[1];
sx q[1];
rz(-0.88729837) q[1];
sx q[1];
rz(-2.3166189) q[1];
rz(2.4752919) q[3];
sx q[3];
rz(-1.9283221) q[3];
sx q[3];
rz(0.038531664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0018953) q[2];
sx q[2];
rz(-1.8997833) q[2];
sx q[2];
rz(-3.042172) q[2];
rz(-2.5771778) q[3];
sx q[3];
rz(-1.5494917) q[3];
sx q[3];
rz(2.2216589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1426706) q[0];
sx q[0];
rz(-2.8746334) q[0];
sx q[0];
rz(0.024918407) q[0];
rz(-1.092356) q[1];
sx q[1];
rz(-1.9905636) q[1];
sx q[1];
rz(2.5999462) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6059581) q[0];
sx q[0];
rz(-0.97554699) q[0];
sx q[0];
rz(-2.7017977) q[0];
x q[1];
rz(-2.140803) q[2];
sx q[2];
rz(-1.8483682) q[2];
sx q[2];
rz(-2.9961207) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.147824) q[1];
sx q[1];
rz(-1.8986101) q[1];
sx q[1];
rz(-2.9891564) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8239519) q[3];
sx q[3];
rz(-2.5196155) q[3];
sx q[3];
rz(0.013381026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3148552) q[2];
sx q[2];
rz(-1.22437) q[2];
sx q[2];
rz(-2.0014191) q[2];
rz(-1.4402116) q[3];
sx q[3];
rz(-2.7482016) q[3];
sx q[3];
rz(-2.9668729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88406554) q[0];
sx q[0];
rz(-1.1648357) q[0];
sx q[0];
rz(-0.87673941) q[0];
rz(0.17403099) q[1];
sx q[1];
rz(-1.0935676) q[1];
sx q[1];
rz(-1.3678975) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0582188) q[0];
sx q[0];
rz(-1.030428) q[0];
sx q[0];
rz(-2.8182993) q[0];
rz(-pi) q[1];
rz(-2.1166685) q[2];
sx q[2];
rz(-2.721933) q[2];
sx q[2];
rz(2.1569606) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1845884) q[1];
sx q[1];
rz(-2.4589029) q[1];
sx q[1];
rz(2.1826988) q[1];
x q[2];
rz(-1.6927244) q[3];
sx q[3];
rz(-1.4953914) q[3];
sx q[3];
rz(0.46368515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6690663) q[2];
sx q[2];
rz(-1.3269227) q[2];
sx q[2];
rz(2.9533022) q[2];
rz(-1.7752198) q[3];
sx q[3];
rz(-1.3683189) q[3];
sx q[3];
rz(1.9395456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61757225) q[0];
sx q[0];
rz(-0.13372788) q[0];
sx q[0];
rz(0.66620859) q[0];
rz(-1.1353525) q[1];
sx q[1];
rz(-2.1609781) q[1];
sx q[1];
rz(1.1599783) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0157601) q[0];
sx q[0];
rz(-2.1938132) q[0];
sx q[0];
rz(0.71970223) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7889889) q[2];
sx q[2];
rz(-2.312783) q[2];
sx q[2];
rz(-0.11889501) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6659703) q[1];
sx q[1];
rz(-0.95487528) q[1];
sx q[1];
rz(3.1321976) q[1];
rz(-2.1362392) q[3];
sx q[3];
rz(-1.8400888) q[3];
sx q[3];
rz(-0.56391844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3975415) q[2];
sx q[2];
rz(-2.5752189) q[2];
sx q[2];
rz(-0.19691021) q[2];
rz(2.2714553) q[3];
sx q[3];
rz(-2.0700442) q[3];
sx q[3];
rz(1.8740694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63721913) q[0];
sx q[0];
rz(-1.0727896) q[0];
sx q[0];
rz(-0.037242591) q[0];
rz(1.0875018) q[1];
sx q[1];
rz(-1.0481513) q[1];
sx q[1];
rz(-2.9748999) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8213538) q[0];
sx q[0];
rz(-1.8379837) q[0];
sx q[0];
rz(2.4334879) q[0];
rz(2.6092763) q[2];
sx q[2];
rz(-0.54907832) q[2];
sx q[2];
rz(-0.86073574) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0089576978) q[1];
sx q[1];
rz(-1.8648498) q[1];
sx q[1];
rz(-1.4786665) q[1];
x q[2];
rz(-2.6981265) q[3];
sx q[3];
rz(-0.16032585) q[3];
sx q[3];
rz(0.043930862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7249666) q[2];
sx q[2];
rz(-2.5278957) q[2];
sx q[2];
rz(1.7566768) q[2];
rz(2.7409605) q[3];
sx q[3];
rz(-1.8568361) q[3];
sx q[3];
rz(-1.0111151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818125) q[0];
sx q[0];
rz(-2.0406944) q[0];
sx q[0];
rz(2.5115321) q[0];
rz(0.13314816) q[1];
sx q[1];
rz(-1.9825736) q[1];
sx q[1];
rz(2.0864743) q[1];
rz(0.038865969) q[2];
sx q[2];
rz(-1.7615223) q[2];
sx q[2];
rz(2.7675235) q[2];
rz(1.6937485) q[3];
sx q[3];
rz(-2.3364441) q[3];
sx q[3];
rz(0.91329109) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
