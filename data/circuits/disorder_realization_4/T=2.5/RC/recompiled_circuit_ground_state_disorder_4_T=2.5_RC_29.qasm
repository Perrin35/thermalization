OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57920116) q[0];
sx q[0];
rz(-0.76051036) q[0];
sx q[0];
rz(1.0150681) q[0];
rz(-1.1633582) q[1];
sx q[1];
rz(-2.7203163) q[1];
sx q[1];
rz(-1.2383229) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8684611) q[0];
sx q[0];
rz(-0.95306153) q[0];
sx q[0];
rz(2.4988089) q[0];
x q[1];
rz(-2.2767645) q[2];
sx q[2];
rz(-0.53542407) q[2];
sx q[2];
rz(-3.0147417) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0933669) q[1];
sx q[1];
rz(-1.2087355) q[1];
sx q[1];
rz(-2.0680548) q[1];
rz(-pi) q[2];
rz(-0.56391509) q[3];
sx q[3];
rz(-1.7132572) q[3];
sx q[3];
rz(-0.13232732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6282661) q[2];
sx q[2];
rz(-1.3506177) q[2];
sx q[2];
rz(-0.70293054) q[2];
rz(-2.2859196) q[3];
sx q[3];
rz(-0.84508768) q[3];
sx q[3];
rz(-1.2969016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6913476) q[0];
sx q[0];
rz(-2.0424728) q[0];
sx q[0];
rz(2.3968089) q[0];
rz(1.5738457) q[1];
sx q[1];
rz(-0.66190043) q[1];
sx q[1];
rz(0.94211284) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5013393) q[0];
sx q[0];
rz(-2.3177409) q[0];
sx q[0];
rz(1.2938611) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81824192) q[2];
sx q[2];
rz(-1.628792) q[2];
sx q[2];
rz(2.4596283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2750863) q[1];
sx q[1];
rz(-0.77436354) q[1];
sx q[1];
rz(1.7656209) q[1];
rz(-3.0511176) q[3];
sx q[3];
rz(-2.5535085) q[3];
sx q[3];
rz(0.17234853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88910237) q[2];
sx q[2];
rz(-0.95688755) q[2];
sx q[2];
rz(2.0460184) q[2];
rz(1.6628294) q[3];
sx q[3];
rz(-0.79083276) q[3];
sx q[3];
rz(1.7553294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11109322) q[0];
sx q[0];
rz(-1.4249304) q[0];
sx q[0];
rz(-1.7362562) q[0];
rz(1.7449215) q[1];
sx q[1];
rz(-1.7878572) q[1];
sx q[1];
rz(2.2415846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6557245) q[0];
sx q[0];
rz(-1.2085755) q[0];
sx q[0];
rz(2.3535504) q[0];
rz(-2.1271187) q[2];
sx q[2];
rz(-2.0790853) q[2];
sx q[2];
rz(-0.1137133) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0327889) q[1];
sx q[1];
rz(-2.3790763) q[1];
sx q[1];
rz(0.25188781) q[1];
rz(0.95696394) q[3];
sx q[3];
rz(-1.004289) q[3];
sx q[3];
rz(1.5885799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46093837) q[2];
sx q[2];
rz(-0.21549455) q[2];
sx q[2];
rz(2.1920965) q[2];
rz(0.62120581) q[3];
sx q[3];
rz(-2.216279) q[3];
sx q[3];
rz(-1.1533823) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7200274) q[0];
sx q[0];
rz(-1.8864487) q[0];
sx q[0];
rz(0.048728745) q[0];
rz(0.85982927) q[1];
sx q[1];
rz(-0.33165926) q[1];
sx q[1];
rz(-2.8299832) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4954715) q[0];
sx q[0];
rz(-0.55495431) q[0];
sx q[0];
rz(2.1053736) q[0];
rz(2.0987058) q[2];
sx q[2];
rz(-1.2706869) q[2];
sx q[2];
rz(-1.2712511) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3032461) q[1];
sx q[1];
rz(-0.7184808) q[1];
sx q[1];
rz(1.9796728) q[1];
x q[2];
rz(-2.8446557) q[3];
sx q[3];
rz(-1.1466168) q[3];
sx q[3];
rz(-2.3564828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.17778808) q[2];
sx q[2];
rz(-0.21169855) q[2];
sx q[2];
rz(-1.6507899) q[2];
rz(0.35495159) q[3];
sx q[3];
rz(-1.7345411) q[3];
sx q[3];
rz(-1.2995592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0519003) q[0];
sx q[0];
rz(-2.7841452) q[0];
sx q[0];
rz(1.8748913) q[0];
rz(3.025324) q[1];
sx q[1];
rz(-2.1513042) q[1];
sx q[1];
rz(1.5142534) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0564041) q[0];
sx q[0];
rz(-1.2582111) q[0];
sx q[0];
rz(2.2645135) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9562938) q[2];
sx q[2];
rz(-2.9045939) q[2];
sx q[2];
rz(2.1182107) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7092412) q[1];
sx q[1];
rz(-2.5951457) q[1];
sx q[1];
rz(1.4724456) q[1];
rz(-pi) q[2];
x q[2];
rz(2.299304) q[3];
sx q[3];
rz(-1.5429879) q[3];
sx q[3];
rz(-0.14735809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9194455) q[2];
sx q[2];
rz(-0.47480348) q[2];
sx q[2];
rz(1.9661281) q[2];
rz(2.2373824) q[3];
sx q[3];
rz(-1.892482) q[3];
sx q[3];
rz(-2.1632532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745875) q[0];
sx q[0];
rz(-2.8103204) q[0];
sx q[0];
rz(0.40147716) q[0];
rz(1.1270771) q[1];
sx q[1];
rz(-1.8311484) q[1];
sx q[1];
rz(0.81261596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3298523) q[0];
sx q[0];
rz(-2.5074258) q[0];
sx q[0];
rz(1.5777753) q[0];
rz(0.70131371) q[2];
sx q[2];
rz(-2.7879386) q[2];
sx q[2];
rz(2.3344699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5548812) q[1];
sx q[1];
rz(-2.5062392) q[1];
sx q[1];
rz(-2.79252) q[1];
rz(-0.74515588) q[3];
sx q[3];
rz(-2.0260915) q[3];
sx q[3];
rz(2.8139092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9408985) q[2];
sx q[2];
rz(-0.55018598) q[2];
sx q[2];
rz(1.5852488) q[2];
rz(-0.19041348) q[3];
sx q[3];
rz(-0.75473458) q[3];
sx q[3];
rz(1.2079027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.82454005) q[0];
sx q[0];
rz(-2.7523478) q[0];
sx q[0];
rz(0.12271605) q[0];
rz(1.9484733) q[1];
sx q[1];
rz(-0.90843186) q[1];
sx q[1];
rz(-0.76748031) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90119367) q[0];
sx q[0];
rz(-1.5195822) q[0];
sx q[0];
rz(-2.2177494) q[0];
rz(-pi) q[1];
rz(0.38184719) q[2];
sx q[2];
rz(-1.4918461) q[2];
sx q[2];
rz(0.71982924) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9373405) q[1];
sx q[1];
rz(-0.91751999) q[1];
sx q[1];
rz(-3.0863347) q[1];
x q[2];
rz(-0.40624491) q[3];
sx q[3];
rz(-1.5749388) q[3];
sx q[3];
rz(-0.7501445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7439338) q[2];
sx q[2];
rz(-3.1199516) q[2];
sx q[2];
rz(-1.5896612) q[2];
rz(1.4895561) q[3];
sx q[3];
rz(-1.7347615) q[3];
sx q[3];
rz(1.1154729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81812304) q[0];
sx q[0];
rz(-2.7796845) q[0];
sx q[0];
rz(0.37539151) q[0];
rz(0.88343945) q[1];
sx q[1];
rz(-1.5366303) q[1];
sx q[1];
rz(-0.72867957) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32266301) q[0];
sx q[0];
rz(-2.0897341) q[0];
sx q[0];
rz(-2.3170018) q[0];
rz(-pi) q[1];
rz(1.7394786) q[2];
sx q[2];
rz(-1.8034435) q[2];
sx q[2];
rz(-2.6017435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9364215) q[1];
sx q[1];
rz(-1.2146307) q[1];
sx q[1];
rz(0.89295279) q[1];
x q[2];
rz(-0.65547319) q[3];
sx q[3];
rz(-0.87067662) q[3];
sx q[3];
rz(-1.164013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6931307) q[2];
sx q[2];
rz(-1.5631661) q[2];
sx q[2];
rz(-2.3942088) q[2];
rz(-1.4061617) q[3];
sx q[3];
rz(-1.9236671) q[3];
sx q[3];
rz(1.5860484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2324227) q[0];
sx q[0];
rz(-2.2012043) q[0];
sx q[0];
rz(-0.7255834) q[0];
rz(-1.3369417) q[1];
sx q[1];
rz(-1.2533816) q[1];
sx q[1];
rz(0.57428378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0030356) q[0];
sx q[0];
rz(-1.4113562) q[0];
sx q[0];
rz(0.11388679) q[0];
rz(-pi) q[1];
rz(2.1663675) q[2];
sx q[2];
rz(-2.4831366) q[2];
sx q[2];
rz(2.0582046) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0813539) q[1];
sx q[1];
rz(-1.4574086) q[1];
sx q[1];
rz(-0.84073034) q[1];
rz(-2.3575555) q[3];
sx q[3];
rz(-2.0052387) q[3];
sx q[3];
rz(1.2617574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1672704) q[2];
sx q[2];
rz(-1.3784626) q[2];
sx q[2];
rz(2.4794225) q[2];
rz(-2.3769489) q[3];
sx q[3];
rz(-1.60631) q[3];
sx q[3];
rz(-1.8173328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.8793256) q[0];
sx q[0];
rz(-2.5536394) q[0];
sx q[0];
rz(1.7437438) q[0];
rz(1.0700048) q[1];
sx q[1];
rz(-2.013423) q[1];
sx q[1];
rz(-1.3319344) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71026826) q[0];
sx q[0];
rz(-1.5402176) q[0];
sx q[0];
rz(-1.7215183) q[0];
x q[1];
rz(-2.0831624) q[2];
sx q[2];
rz(-0.46736267) q[2];
sx q[2];
rz(0.99603727) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8269013) q[1];
sx q[1];
rz(-1.3413635) q[1];
sx q[1];
rz(-1.3342085) q[1];
rz(-pi) q[2];
rz(-0.58140786) q[3];
sx q[3];
rz(-0.96145844) q[3];
sx q[3];
rz(-1.7096303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8269044) q[2];
sx q[2];
rz(-2.0411699) q[2];
sx q[2];
rz(0.17987128) q[2];
rz(1.3966857) q[3];
sx q[3];
rz(-1.4254009) q[3];
sx q[3];
rz(0.21656187) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3708645) q[0];
sx q[0];
rz(-0.48038078) q[0];
sx q[0];
rz(0.90850716) q[0];
rz(0.52275672) q[1];
sx q[1];
rz(-2.0261384) q[1];
sx q[1];
rz(-1.1631858) q[1];
rz(-2.592587) q[2];
sx q[2];
rz(-0.71239757) q[2];
sx q[2];
rz(-3.0511643) q[2];
rz(2.497429) q[3];
sx q[3];
rz(-2.5732891) q[3];
sx q[3];
rz(2.7960232) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
