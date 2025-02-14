OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7423695) q[0];
sx q[0];
rz(-0.41271451) q[0];
sx q[0];
rz(-2.3053919) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(6.3451938) q[1];
sx q[1];
rz(11.558029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22577408) q[0];
sx q[0];
rz(-1.4918707) q[0];
sx q[0];
rz(1.945334) q[0];
rz(-pi) q[1];
rz(2.7616829) q[2];
sx q[2];
rz(-1.7899982) q[2];
sx q[2];
rz(-0.70218147) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5656149) q[1];
sx q[1];
rz(-1.4139785) q[1];
sx q[1];
rz(-2.6601936) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7640231) q[3];
sx q[3];
rz(-1.591103) q[3];
sx q[3];
rz(1.8373674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8334373) q[2];
sx q[2];
rz(-0.94352949) q[2];
sx q[2];
rz(-0.65307871) q[2];
rz(3.0270882) q[3];
sx q[3];
rz(-1.3936035) q[3];
sx q[3];
rz(-2.0792927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36650518) q[0];
sx q[0];
rz(-2.1108284) q[0];
sx q[0];
rz(1.3436226) q[0];
rz(-1.1603181) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(2.5210099) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40784971) q[0];
sx q[0];
rz(-1.7622927) q[0];
sx q[0];
rz(-1.2194077) q[0];
x q[1];
rz(0.36850117) q[2];
sx q[2];
rz(-0.93612367) q[2];
sx q[2];
rz(2.2919185) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3239336) q[1];
sx q[1];
rz(-0.91485564) q[1];
sx q[1];
rz(1.1088015) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9602523) q[3];
sx q[3];
rz(-0.6529633) q[3];
sx q[3];
rz(-0.12784004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4429984) q[2];
sx q[2];
rz(-1.6220762) q[2];
sx q[2];
rz(-0.2300187) q[2];
rz(1.5551785) q[3];
sx q[3];
rz(-1.9929726) q[3];
sx q[3];
rz(0.27568278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8860633) q[0];
sx q[0];
rz(-0.38917381) q[0];
sx q[0];
rz(1.080876) q[0];
rz(2.4983662) q[1];
sx q[1];
rz(-2.3422362) q[1];
sx q[1];
rz(-0.08531514) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9832343) q[0];
sx q[0];
rz(-2.2081349) q[0];
sx q[0];
rz(2.4054224) q[0];
x q[1];
rz(2.8181452) q[2];
sx q[2];
rz(-0.5964884) q[2];
sx q[2];
rz(-1.8307387) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18135611) q[1];
sx q[1];
rz(-1.7367558) q[1];
sx q[1];
rz(1.7039429) q[1];
x q[2];
rz(-0.063488678) q[3];
sx q[3];
rz(-0.62219071) q[3];
sx q[3];
rz(0.49294642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2049415) q[2];
sx q[2];
rz(-0.0095657883) q[2];
sx q[2];
rz(0.19634518) q[2];
rz(0.28228545) q[3];
sx q[3];
rz(-2.0245656) q[3];
sx q[3];
rz(-1.0452247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1333756) q[0];
sx q[0];
rz(-0.5468002) q[0];
sx q[0];
rz(2.7561482) q[0];
rz(-1.4133981) q[1];
sx q[1];
rz(-1.9368659) q[1];
sx q[1];
rz(0.24994303) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9871368) q[0];
sx q[0];
rz(-2.3780795) q[0];
sx q[0];
rz(0.95434965) q[0];
x q[1];
rz(0.59993292) q[2];
sx q[2];
rz(-3.0351808) q[2];
sx q[2];
rz(-1.0628884) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7467667) q[1];
sx q[1];
rz(-1.7395081) q[1];
sx q[1];
rz(-0.49612237) q[1];
rz(-pi) q[2];
rz(1.3718666) q[3];
sx q[3];
rz(-2.3545111) q[3];
sx q[3];
rz(2.6347292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6803153) q[2];
sx q[2];
rz(-1.2948493) q[2];
sx q[2];
rz(0.2363905) q[2];
rz(-1.2605028) q[3];
sx q[3];
rz(-0.25006306) q[3];
sx q[3];
rz(-2.4642956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33609718) q[0];
sx q[0];
rz(-1.5376872) q[0];
sx q[0];
rz(2.6564823) q[0];
rz(-1.1803892) q[1];
sx q[1];
rz(-1.5472658) q[1];
sx q[1];
rz(-0.83650437) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.361995) q[0];
sx q[0];
rz(-1.9266202) q[0];
sx q[0];
rz(1.2339639) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9747462) q[2];
sx q[2];
rz(-1.8649668) q[2];
sx q[2];
rz(-2.6427513) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.066557601) q[1];
sx q[1];
rz(-1.5227942) q[1];
sx q[1];
rz(-0.90086277) q[1];
x q[2];
rz(0.14797609) q[3];
sx q[3];
rz(-2.736909) q[3];
sx q[3];
rz(-1.2860822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2430719) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(-0.76954976) q[2];
rz(-2.2122993) q[3];
sx q[3];
rz(-1.1309036) q[3];
sx q[3];
rz(-2.9088959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9418075) q[0];
sx q[0];
rz(-0.48518825) q[0];
sx q[0];
rz(2.3727681) q[0];
rz(-1.3517514) q[1];
sx q[1];
rz(-2.6685346) q[1];
sx q[1];
rz(-2.2519055) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8990435) q[0];
sx q[0];
rz(-2.7022323) q[0];
sx q[0];
rz(1.0602632) q[0];
rz(-pi) q[1];
rz(3.1402428) q[2];
sx q[2];
rz(-1.6151516) q[2];
sx q[2];
rz(2.7173017) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2592247) q[1];
sx q[1];
rz(-2.9881757) q[1];
sx q[1];
rz(1.2620366) q[1];
x q[2];
rz(1.6857288) q[3];
sx q[3];
rz(-1.1224261) q[3];
sx q[3];
rz(1.861972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6512904) q[2];
sx q[2];
rz(-1.7040665) q[2];
sx q[2];
rz(-1.5865405) q[2];
rz(-1.2706612) q[3];
sx q[3];
rz(-0.41897604) q[3];
sx q[3];
rz(-3.0095625) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1726058) q[0];
sx q[0];
rz(-1.1382599) q[0];
sx q[0];
rz(-1.1440811) q[0];
rz(0.83089685) q[1];
sx q[1];
rz(-1.3583207) q[1];
sx q[1];
rz(2.3544618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2427776) q[0];
sx q[0];
rz(-1.7022331) q[0];
sx q[0];
rz(-1.4351033) q[0];
rz(-pi) q[1];
rz(-2.2058555) q[2];
sx q[2];
rz(-0.79197403) q[2];
sx q[2];
rz(-2.9284649) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5980144) q[1];
sx q[1];
rz(-2.0435963) q[1];
sx q[1];
rz(0.531457) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5821286) q[3];
sx q[3];
rz(-1.2351523) q[3];
sx q[3];
rz(1.0192724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0234915) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(-0.033163158) q[2];
rz(-1.5432594) q[3];
sx q[3];
rz(-0.71591806) q[3];
sx q[3];
rz(-1.6395114) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.892266) q[0];
sx q[0];
rz(-1.5557657) q[0];
sx q[0];
rz(-1.3931042) q[0];
rz(-2.6847367) q[1];
sx q[1];
rz(-1.9970048) q[1];
sx q[1];
rz(1.1005864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8580061) q[0];
sx q[0];
rz(-1.6467983) q[0];
sx q[0];
rz(1.5690687) q[0];
rz(-pi) q[1];
rz(-0.92488168) q[2];
sx q[2];
rz(-2.1585228) q[2];
sx q[2];
rz(0.6494259) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7707899) q[1];
sx q[1];
rz(-0.90401559) q[1];
sx q[1];
rz(-0.043170269) q[1];
rz(1.4310775) q[3];
sx q[3];
rz(-1.9587687) q[3];
sx q[3];
rz(0.037484976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0988203) q[2];
sx q[2];
rz(-0.49751147) q[2];
sx q[2];
rz(-1.5607321) q[2];
rz(2.3780195) q[3];
sx q[3];
rz(-1.6288501) q[3];
sx q[3];
rz(1.0360576) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40912691) q[0];
sx q[0];
rz(-2.3865073) q[0];
sx q[0];
rz(-0.29148802) q[0];
rz(-1.2358707) q[1];
sx q[1];
rz(-1.9449077) q[1];
sx q[1];
rz(3.018697) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.204958) q[0];
sx q[0];
rz(-0.35047784) q[0];
sx q[0];
rz(0.15021439) q[0];
x q[1];
rz(1.8529296) q[2];
sx q[2];
rz(-1.3569433) q[2];
sx q[2];
rz(0.58889533) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1157411) q[1];
sx q[1];
rz(-0.23499712) q[1];
sx q[1];
rz(2.0629289) q[1];
x q[2];
rz(-1.4261774) q[3];
sx q[3];
rz(-1.5767617) q[3];
sx q[3];
rz(-0.98814135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62873944) q[2];
sx q[2];
rz(-1.3924007) q[2];
sx q[2];
rz(2.0780308) q[2];
rz(0.17744803) q[3];
sx q[3];
rz(-0.95824233) q[3];
sx q[3];
rz(-2.1232846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43750957) q[0];
sx q[0];
rz(-0.45314416) q[0];
sx q[0];
rz(-1.4310687) q[0];
rz(2.5408632) q[1];
sx q[1];
rz(-0.44980106) q[1];
sx q[1];
rz(-0.61753714) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73752943) q[0];
sx q[0];
rz(-1.4466219) q[0];
sx q[0];
rz(-1.7880718) q[0];
rz(-pi) q[1];
rz(-0.52458101) q[2];
sx q[2];
rz(-2.2426105) q[2];
sx q[2];
rz(1.8730522) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3869218) q[1];
sx q[1];
rz(-2.6082509) q[1];
sx q[1];
rz(-1.6120595) q[1];
x q[2];
rz(-1.2631868) q[3];
sx q[3];
rz(-2.121939) q[3];
sx q[3];
rz(1.2778417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2403229) q[2];
sx q[2];
rz(-2.1412886) q[2];
sx q[2];
rz(-2.1642302) q[2];
rz(2.0591002) q[3];
sx q[3];
rz(-2.2668224) q[3];
sx q[3];
rz(-0.055518363) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8025773) q[0];
sx q[0];
rz(-2.3727198) q[0];
sx q[0];
rz(2.4045237) q[0];
rz(-0.045724178) q[1];
sx q[1];
rz(-2.3241691) q[1];
sx q[1];
rz(-1.587422) q[1];
rz(-2.1004213) q[2];
sx q[2];
rz(-2.0176135) q[2];
sx q[2];
rz(1.8420646) q[2];
rz(-0.65683881) q[3];
sx q[3];
rz(-1.4359063) q[3];
sx q[3];
rz(-2.0834854) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
