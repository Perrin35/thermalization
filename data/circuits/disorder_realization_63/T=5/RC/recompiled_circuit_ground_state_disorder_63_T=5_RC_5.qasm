OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7665793) q[0];
sx q[0];
rz(-1.2393247) q[0];
sx q[0];
rz(-2.0517321) q[0];
rz(1.0640979) q[1];
sx q[1];
rz(4.394726) q[1];
sx q[1];
rz(10.328007) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8860461) q[0];
sx q[0];
rz(-0.95661344) q[0];
sx q[0];
rz(-0.79844676) q[0];
x q[1];
rz(-0.8875928) q[2];
sx q[2];
rz(-1.5593865) q[2];
sx q[2];
rz(1.3225088) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9574388) q[1];
sx q[1];
rz(-1.2278964) q[1];
sx q[1];
rz(0.6253021) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12046195) q[3];
sx q[3];
rz(-1.0742004) q[3];
sx q[3];
rz(2.7169926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.696306) q[2];
sx q[2];
rz(-3.004965) q[2];
sx q[2];
rz(-2.7083) q[2];
rz(2.8516234) q[3];
sx q[3];
rz(-2.2765997) q[3];
sx q[3];
rz(-0.32052952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1746154) q[0];
sx q[0];
rz(-2.2791635) q[0];
sx q[0];
rz(-1.2906661) q[0];
rz(2.4621452) q[1];
sx q[1];
rz(-0.77630711) q[1];
sx q[1];
rz(-2.2490833) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7793625) q[0];
sx q[0];
rz(-0.17900544) q[0];
sx q[0];
rz(-1.711861) q[0];
rz(-pi) q[1];
rz(1.1721482) q[2];
sx q[2];
rz(-2.1855693) q[2];
sx q[2];
rz(3.1190256) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8645382) q[1];
sx q[1];
rz(-1.2965585) q[1];
sx q[1];
rz(-2.4692593) q[1];
rz(2.1492747) q[3];
sx q[3];
rz(-2.7035993) q[3];
sx q[3];
rz(2.7051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.14900011) q[2];
sx q[2];
rz(-2.0626103) q[2];
sx q[2];
rz(2.0089669) q[2];
rz(2.4752786) q[3];
sx q[3];
rz(-1.3601114) q[3];
sx q[3];
rz(-1.2247491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76729524) q[0];
sx q[0];
rz(-2.7485924) q[0];
sx q[0];
rz(2.1597916) q[0];
rz(-2.1994195) q[1];
sx q[1];
rz(-1.8323545) q[1];
sx q[1];
rz(0.91032496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2766586) q[0];
sx q[0];
rz(-2.4747643) q[0];
sx q[0];
rz(2.4163321) q[0];
rz(-pi) q[1];
rz(-1.6907755) q[2];
sx q[2];
rz(-1.8455122) q[2];
sx q[2];
rz(-0.071823013) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3094448) q[1];
sx q[1];
rz(-2.6845686) q[1];
sx q[1];
rz(0.8330658) q[1];
rz(-pi) q[2];
rz(2.0098088) q[3];
sx q[3];
rz(-1.3512423) q[3];
sx q[3];
rz(1.9569524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2008449) q[2];
sx q[2];
rz(-2.7767599) q[2];
sx q[2];
rz(2.9546837) q[2];
rz(2.9777891) q[3];
sx q[3];
rz(-2.2713594) q[3];
sx q[3];
rz(0.12266172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1214889) q[0];
sx q[0];
rz(-1.7262456) q[0];
sx q[0];
rz(-2.4439268) q[0];
rz(0.54667073) q[1];
sx q[1];
rz(-2.8488939) q[1];
sx q[1];
rz(1.9465416) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89681399) q[0];
sx q[0];
rz(-1.91511) q[0];
sx q[0];
rz(-0.94621511) q[0];
rz(-0.61364737) q[2];
sx q[2];
rz(-2.9093008) q[2];
sx q[2];
rz(0.42343806) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.66142144) q[1];
sx q[1];
rz(-1.464559) q[1];
sx q[1];
rz(-1.7376971) q[1];
rz(-pi) q[2];
rz(-3.0680858) q[3];
sx q[3];
rz(-1.8544844) q[3];
sx q[3];
rz(0.43108854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.67808548) q[2];
sx q[2];
rz(-1.4231851) q[2];
sx q[2];
rz(-0.19742337) q[2];
rz(3.0366963) q[3];
sx q[3];
rz(-1.1095122) q[3];
sx q[3];
rz(0.088002861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5941493) q[0];
sx q[0];
rz(-0.53590411) q[0];
sx q[0];
rz(2.1323668) q[0];
rz(3.1127473) q[1];
sx q[1];
rz(-1.4776769) q[1];
sx q[1];
rz(-0.65863329) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65153367) q[0];
sx q[0];
rz(-0.57350993) q[0];
sx q[0];
rz(-0.83322462) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20435996) q[2];
sx q[2];
rz(-1.0931226) q[2];
sx q[2];
rz(-0.43276873) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7004905) q[1];
sx q[1];
rz(-1.6630739) q[1];
sx q[1];
rz(-0.16903318) q[1];
rz(-pi) q[2];
rz(-2.4464408) q[3];
sx q[3];
rz(-1.3926695) q[3];
sx q[3];
rz(-0.069308829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.37904) q[2];
sx q[2];
rz(-0.54658824) q[2];
sx q[2];
rz(2.9075882) q[2];
rz(2.0512569) q[3];
sx q[3];
rz(-1.5749616) q[3];
sx q[3];
rz(2.0719297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91319084) q[0];
sx q[0];
rz(-0.15616067) q[0];
sx q[0];
rz(-1.2983904) q[0];
rz(-2.3550854) q[1];
sx q[1];
rz(-2.2857917) q[1];
sx q[1];
rz(1.7826805) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14058622) q[0];
sx q[0];
rz(-1.6504262) q[0];
sx q[0];
rz(-2.4796159) q[0];
rz(-pi) q[1];
rz(1.3471021) q[2];
sx q[2];
rz(-1.8321206) q[2];
sx q[2];
rz(1.4841207) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7715367) q[1];
sx q[1];
rz(-3.1395457) q[1];
sx q[1];
rz(0.63315804) q[1];
rz(-pi) q[2];
rz(-1.2511061) q[3];
sx q[3];
rz(-2.0123008) q[3];
sx q[3];
rz(2.6755345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7274373) q[2];
sx q[2];
rz(-1.6748019) q[2];
sx q[2];
rz(0.1725014) q[2];
rz(-2.6532459) q[3];
sx q[3];
rz(-1.1427053) q[3];
sx q[3];
rz(-1.1138227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1682424) q[0];
sx q[0];
rz(-0.872648) q[0];
sx q[0];
rz(2.8016222) q[0];
rz(-0.32644692) q[1];
sx q[1];
rz(-2.4531334) q[1];
sx q[1];
rz(-1.2350157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9654953) q[0];
sx q[0];
rz(-1.1998994) q[0];
sx q[0];
rz(2.4489016) q[0];
rz(-pi) q[1];
rz(-2.6831362) q[2];
sx q[2];
rz(-1.2453015) q[2];
sx q[2];
rz(0.95409648) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81081796) q[1];
sx q[1];
rz(-1.0913117) q[1];
sx q[1];
rz(-1.5768087) q[1];
x q[2];
rz(-0.49576393) q[3];
sx q[3];
rz(-2.7003482) q[3];
sx q[3];
rz(3.038681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1116703) q[2];
sx q[2];
rz(-0.80283529) q[2];
sx q[2];
rz(-0.55533448) q[2];
rz(2.9739144) q[3];
sx q[3];
rz(-1.6043112) q[3];
sx q[3];
rz(-0.47784561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4324206) q[0];
sx q[0];
rz(-0.92996159) q[0];
sx q[0];
rz(2.0322556) q[0];
rz(2.0093911) q[1];
sx q[1];
rz(-0.39674509) q[1];
sx q[1];
rz(-2.1999377) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50186673) q[0];
sx q[0];
rz(-1.4113771) q[0];
sx q[0];
rz(2.1306043) q[0];
rz(-2.0600165) q[2];
sx q[2];
rz(-0.3210381) q[2];
sx q[2];
rz(-2.5550317) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.68384993) q[1];
sx q[1];
rz(-0.93194973) q[1];
sx q[1];
rz(-2.9838377) q[1];
rz(-pi) q[2];
rz(1.3579943) q[3];
sx q[3];
rz(-1.9486685) q[3];
sx q[3];
rz(-1.6869643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38810101) q[2];
sx q[2];
rz(-1.0270303) q[2];
sx q[2];
rz(1.6027742) q[2];
rz(-1.7704891) q[3];
sx q[3];
rz(-1.1857827) q[3];
sx q[3];
rz(-0.051430844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082315363) q[0];
sx q[0];
rz(-2.5264854) q[0];
sx q[0];
rz(-2.1737461) q[0];
rz(-1.2713426) q[1];
sx q[1];
rz(-1.8495193) q[1];
sx q[1];
rz(-0.06180067) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74154749) q[0];
sx q[0];
rz(-1.9572659) q[0];
sx q[0];
rz(-0.58048141) q[0];
rz(-pi) q[1];
rz(0.65431739) q[2];
sx q[2];
rz(-0.32233244) q[2];
sx q[2];
rz(0.92609012) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.7605412) q[1];
sx q[1];
rz(-2.2534721) q[1];
sx q[1];
rz(-1.447601) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8093407) q[3];
sx q[3];
rz(-1.2945172) q[3];
sx q[3];
rz(-1.075923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9515848) q[2];
sx q[2];
rz(-2.9374359) q[2];
sx q[2];
rz(-0.07587138) q[2];
rz(-1.1672945) q[3];
sx q[3];
rz(-1.1018402) q[3];
sx q[3];
rz(-0.30437881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9645914) q[0];
sx q[0];
rz(-0.27934203) q[0];
sx q[0];
rz(-2.8569073) q[0];
rz(0.074706569) q[1];
sx q[1];
rz(-1.9120522) q[1];
sx q[1];
rz(-1.4178735) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79399949) q[0];
sx q[0];
rz(-2.3132304) q[0];
sx q[0];
rz(1.8338127) q[0];
rz(-pi) q[1];
rz(-1.8080797) q[2];
sx q[2];
rz(-0.56979376) q[2];
sx q[2];
rz(0.89536086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26894179) q[1];
sx q[1];
rz(-1.1564019) q[1];
sx q[1];
rz(0.68434836) q[1];
rz(-0.033662658) q[3];
sx q[3];
rz(-2.2164248) q[3];
sx q[3];
rz(1.0317623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8031926) q[2];
sx q[2];
rz(-2.000838) q[2];
sx q[2];
rz(1.8624064) q[2];
rz(-2.159436) q[3];
sx q[3];
rz(-1.4582062) q[3];
sx q[3];
rz(2.2568683) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2411156) q[0];
sx q[0];
rz(-1.8178839) q[0];
sx q[0];
rz(-1.5771014) q[0];
rz(2.1263532) q[1];
sx q[1];
rz(-1.1493586) q[1];
sx q[1];
rz(-1.1372067) q[1];
rz(1.8335291) q[2];
sx q[2];
rz(-1.4663525) q[2];
sx q[2];
rz(2.732085) q[2];
rz(-0.7393403) q[3];
sx q[3];
rz(-1.1237595) q[3];
sx q[3];
rz(2.9343284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
