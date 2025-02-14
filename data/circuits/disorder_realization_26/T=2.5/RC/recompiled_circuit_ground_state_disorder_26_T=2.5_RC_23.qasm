OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5717918) q[0];
sx q[0];
rz(-0.90919149) q[0];
sx q[0];
rz(1.5067014) q[0];
rz(-1.934496) q[1];
sx q[1];
rz(-0.49270448) q[1];
sx q[1];
rz(0.55941137) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052033612) q[0];
sx q[0];
rz(-1.8843643) q[0];
sx q[0];
rz(-1.3960209) q[0];
rz(-pi) q[1];
rz(-1.1790397) q[2];
sx q[2];
rz(-1.1572574) q[2];
sx q[2];
rz(-0.79038436) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2845498) q[1];
sx q[1];
rz(-2.8618097) q[1];
sx q[1];
rz(-1.5460874) q[1];
x q[2];
rz(0.97752882) q[3];
sx q[3];
rz(-1.9936322) q[3];
sx q[3];
rz(2.5893167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5514465) q[2];
sx q[2];
rz(-2.2821653) q[2];
sx q[2];
rz(-2.0476511) q[2];
rz(-1.9048196) q[3];
sx q[3];
rz(-1.9779132) q[3];
sx q[3];
rz(-0.26113233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7673489) q[0];
sx q[0];
rz(-1.2232895) q[0];
sx q[0];
rz(-0.71794024) q[0];
rz(-1.1942489) q[1];
sx q[1];
rz(-0.4078882) q[1];
sx q[1];
rz(-1.6494707) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2860752) q[0];
sx q[0];
rz(-1.4465032) q[0];
sx q[0];
rz(-3.0086424) q[0];
x q[1];
rz(2.7020383) q[2];
sx q[2];
rz(-1.9157404) q[2];
sx q[2];
rz(0.82122356) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1528984) q[1];
sx q[1];
rz(-2.0190658) q[1];
sx q[1];
rz(2.6397735) q[1];
rz(-pi) q[2];
rz(2.3158595) q[3];
sx q[3];
rz(-1.4108218) q[3];
sx q[3];
rz(-2.3478594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0933928) q[2];
sx q[2];
rz(-2.3651249) q[2];
sx q[2];
rz(-0.36273599) q[2];
rz(2.2705966) q[3];
sx q[3];
rz(-1.769915) q[3];
sx q[3];
rz(-1.8179651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0469623) q[0];
sx q[0];
rz(-1.8383263) q[0];
sx q[0];
rz(-2.5808425) q[0];
rz(-2.1669855) q[1];
sx q[1];
rz(-0.86169306) q[1];
sx q[1];
rz(-0.21751705) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.226885) q[0];
sx q[0];
rz(-0.98641005) q[0];
sx q[0];
rz(1.5938207) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0571805) q[2];
sx q[2];
rz(-1.660778) q[2];
sx q[2];
rz(-2.1291901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8629093) q[1];
sx q[1];
rz(-0.28537286) q[1];
sx q[1];
rz(-1.9888269) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60902304) q[3];
sx q[3];
rz(-1.6224738) q[3];
sx q[3];
rz(-0.9059815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64323032) q[2];
sx q[2];
rz(-0.78780323) q[2];
sx q[2];
rz(-2.1882449) q[2];
rz(-2.0638454) q[3];
sx q[3];
rz(-1.6809623) q[3];
sx q[3];
rz(-0.47422153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.137602) q[0];
sx q[0];
rz(-0.1739665) q[0];
sx q[0];
rz(-2.2250788) q[0];
rz(-2.7169531) q[1];
sx q[1];
rz(-1.3820796) q[1];
sx q[1];
rz(2.0133846) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5954428) q[0];
sx q[0];
rz(-2.2227123) q[0];
sx q[0];
rz(-0.94080399) q[0];
rz(-pi) q[1];
rz(-0.75707423) q[2];
sx q[2];
rz(-1.7221976) q[2];
sx q[2];
rz(-2.840207) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1188075) q[1];
sx q[1];
rz(-1.2712384) q[1];
sx q[1];
rz(-2.5980298) q[1];
rz(-pi) q[2];
rz(0.97336412) q[3];
sx q[3];
rz(-1.3829682) q[3];
sx q[3];
rz(-2.0300161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8757223) q[2];
sx q[2];
rz(-1.46572) q[2];
sx q[2];
rz(1.1992559) q[2];
rz(-1.9843598) q[3];
sx q[3];
rz(-2.3374989) q[3];
sx q[3];
rz(0.51586241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98738325) q[0];
sx q[0];
rz(-3.0976384) q[0];
sx q[0];
rz(-0.92535812) q[0];
rz(0.56272733) q[1];
sx q[1];
rz(-1.0147164) q[1];
sx q[1];
rz(-0.37511197) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0384509) q[0];
sx q[0];
rz(-1.4618357) q[0];
sx q[0];
rz(-0.54619638) q[0];
rz(-pi) q[1];
x q[1];
rz(1.135446) q[2];
sx q[2];
rz(-1.5223961) q[2];
sx q[2];
rz(-1.2788926) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7081333) q[1];
sx q[1];
rz(-2.9707751) q[1];
sx q[1];
rz(-0.14651919) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45467037) q[3];
sx q[3];
rz(-1.2068307) q[3];
sx q[3];
rz(-0.92586799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.90317) q[2];
sx q[2];
rz(-2.2245202) q[2];
sx q[2];
rz(-2.6440716) q[2];
rz(2.7072952) q[3];
sx q[3];
rz(-1.4562166) q[3];
sx q[3];
rz(0.98880497) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76181805) q[0];
sx q[0];
rz(-2.2803545) q[0];
sx q[0];
rz(-2.8616943) q[0];
rz(-2.1474536) q[1];
sx q[1];
rz(-0.86052624) q[1];
sx q[1];
rz(2.5631189) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46723235) q[0];
sx q[0];
rz(-2.0772837) q[0];
sx q[0];
rz(-1.1367528) q[0];
rz(-pi) q[1];
rz(-0.45409338) q[2];
sx q[2];
rz(-2.3857949) q[2];
sx q[2];
rz(-2.2735655) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0117409) q[1];
sx q[1];
rz(-1.7659582) q[1];
sx q[1];
rz(2.4052448) q[1];
x q[2];
rz(-0.26831823) q[3];
sx q[3];
rz(-1.473241) q[3];
sx q[3];
rz(0.25948157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0591187) q[2];
sx q[2];
rz(-2.5422577) q[2];
sx q[2];
rz(-2.2590051) q[2];
rz(2.5146218) q[3];
sx q[3];
rz(-2.4866703) q[3];
sx q[3];
rz(-2.2466808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6249228) q[0];
sx q[0];
rz(-2.1997917) q[0];
sx q[0];
rz(3.140977) q[0];
rz(2.2387538) q[1];
sx q[1];
rz(-0.87264624) q[1];
sx q[1];
rz(0.045086233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0233399) q[0];
sx q[0];
rz(-1.0824507) q[0];
sx q[0];
rz(-1.4834542) q[0];
x q[1];
rz(0.22901625) q[2];
sx q[2];
rz(-0.75521246) q[2];
sx q[2];
rz(1.3431741) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5645091) q[1];
sx q[1];
rz(-0.85751837) q[1];
sx q[1];
rz(-1.9803995) q[1];
x q[2];
rz(-0.99211971) q[3];
sx q[3];
rz(-1.0401772) q[3];
sx q[3];
rz(1.19095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2948598) q[2];
sx q[2];
rz(-0.555987) q[2];
sx q[2];
rz(2.2247458) q[2];
rz(-0.89546853) q[3];
sx q[3];
rz(-1.5119036) q[3];
sx q[3];
rz(-1.9303314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(3.004892) q[0];
sx q[0];
rz(-0.08575511) q[0];
sx q[0];
rz(-0.46491796) q[0];
rz(-1.1888986) q[1];
sx q[1];
rz(-1.7949972) q[1];
sx q[1];
rz(-2.7704923) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93286033) q[0];
sx q[0];
rz(-1.9813054) q[0];
sx q[0];
rz(-0.44637827) q[0];
rz(-3.1041652) q[2];
sx q[2];
rz(-2.3212395) q[2];
sx q[2];
rz(-1.4701173) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7830989) q[1];
sx q[1];
rz(-0.91121735) q[1];
sx q[1];
rz(0.5650608) q[1];
rz(1.281776) q[3];
sx q[3];
rz(-1.9403696) q[3];
sx q[3];
rz(0.7909382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4798639) q[2];
sx q[2];
rz(-2.3532545) q[2];
sx q[2];
rz(0.41024497) q[2];
rz(-1.5069626) q[3];
sx q[3];
rz(-0.40532902) q[3];
sx q[3];
rz(3.098367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0931382) q[0];
sx q[0];
rz(-2.6148836) q[0];
sx q[0];
rz(-0.17573892) q[0];
rz(-0.82548347) q[1];
sx q[1];
rz(-1.261542) q[1];
sx q[1];
rz(0.82297355) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1639481) q[0];
sx q[0];
rz(-2.2744482) q[0];
sx q[0];
rz(0.73935853) q[0];
rz(1.9972598) q[2];
sx q[2];
rz(-1.8606295) q[2];
sx q[2];
rz(2.5215931) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.041391782) q[1];
sx q[1];
rz(-0.90040246) q[1];
sx q[1];
rz(-1.1213699) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9814616) q[3];
sx q[3];
rz(-1.5370099) q[3];
sx q[3];
rz(0.50212348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2814111) q[2];
sx q[2];
rz(-1.9965636) q[2];
sx q[2];
rz(-0.54110503) q[2];
rz(-0.42630729) q[3];
sx q[3];
rz(-0.46190327) q[3];
sx q[3];
rz(-0.84686744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68778872) q[0];
sx q[0];
rz(-2.332088) q[0];
sx q[0];
rz(2.8059106) q[0];
rz(0.022739284) q[1];
sx q[1];
rz(-2.4191861) q[1];
sx q[1];
rz(2.4679599) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84647932) q[0];
sx q[0];
rz(-1.9178977) q[0];
sx q[0];
rz(-0.62008574) q[0];
x q[1];
rz(2.9433374) q[2];
sx q[2];
rz(-0.47706826) q[2];
sx q[2];
rz(-1.6828293) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9424098) q[1];
sx q[1];
rz(-1.7734151) q[1];
sx q[1];
rz(-0.97572787) q[1];
x q[2];
rz(-2.8073714) q[3];
sx q[3];
rz(-1.4386192) q[3];
sx q[3];
rz(-1.6652938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.444904) q[2];
sx q[2];
rz(-2.1477369) q[2];
sx q[2];
rz(2.135684) q[2];
rz(0.78399793) q[3];
sx q[3];
rz(-2.8204212) q[3];
sx q[3];
rz(1.706749) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0189331) q[0];
sx q[0];
rz(-1.4228595) q[0];
sx q[0];
rz(2.0056437) q[0];
rz(-0.38416531) q[1];
sx q[1];
rz(-1.9728248) q[1];
sx q[1];
rz(-1.493175) q[1];
rz(0.21239077) q[2];
sx q[2];
rz(-1.7707386) q[2];
sx q[2];
rz(2.7223827) q[2];
rz(0.85687153) q[3];
sx q[3];
rz(-1.8663434) q[3];
sx q[3];
rz(0.28874884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
