OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(-2.2990062) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(-1.0993212) q[1];
sx q[1];
rz(-1.6834747) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5320839) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(-0.30755933) q[0];
rz(-2.5277532) q[2];
sx q[2];
rz(-1.5868574) q[2];
sx q[2];
rz(1.2889372) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1229413) q[1];
sx q[1];
rz(-2.7416347) q[1];
sx q[1];
rz(-2.8040228) q[1];
rz(-pi) q[2];
rz(0.59431608) q[3];
sx q[3];
rz(-1.4694957) q[3];
sx q[3];
rz(0.017410226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6136916) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(2.9620985) q[2];
rz(1.2256631) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(-0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74801385) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(0.32546145) q[0];
rz(1.356396) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(1.9869841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5886473) q[0];
sx q[0];
rz(-1.5536904) q[0];
sx q[0];
rz(3.1228035) q[0];
rz(-pi) q[1];
rz(-0.94475586) q[2];
sx q[2];
rz(-1.895004) q[2];
sx q[2];
rz(2.7446483) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7696015) q[1];
sx q[1];
rz(-2.3725315) q[1];
sx q[1];
rz(0.12188697) q[1];
rz(-pi) q[2];
rz(1.3080018) q[3];
sx q[3];
rz(-1.7495219) q[3];
sx q[3];
rz(-2.5457515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(-2.2581805) q[2];
rz(-2.6702821) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(-0.78770351) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8283591) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(-1.6261684) q[0];
rz(2.5405163) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(-1.0916969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1611623) q[0];
sx q[0];
rz(-2.3083901) q[0];
sx q[0];
rz(-2.3479793) q[0];
rz(0.91471471) q[2];
sx q[2];
rz(-1.2351742) q[2];
sx q[2];
rz(0.4682954) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6669238) q[1];
sx q[1];
rz(-0.83819929) q[1];
sx q[1];
rz(-1.0679507) q[1];
x q[2];
rz(-3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(-1.2325866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.320257) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(-0.88095218) q[2];
rz(1.3736003) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3110733) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(2.7048892) q[0];
rz(-0.23315915) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(2.8312347) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6757641) q[0];
sx q[0];
rz(-2.7390263) q[0];
sx q[0];
rz(-2.7990544) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68508673) q[2];
sx q[2];
rz(-1.6685467) q[2];
sx q[2];
rz(-3.0435261) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.031361) q[1];
sx q[1];
rz(-2.7831315) q[1];
sx q[1];
rz(-1.5178174) q[1];
rz(1.849732) q[3];
sx q[3];
rz(-0.12860563) q[3];
sx q[3];
rz(-1.0517373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0115396) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(-1.0774353) q[2];
rz(0.056190101) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(-1.594054) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-2.8919343) q[0];
rz(-1.5646308) q[1];
sx q[1];
rz(-2.3639634) q[1];
sx q[1];
rz(2.2713984) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2244959) q[0];
sx q[0];
rz(-1.3670237) q[0];
sx q[0];
rz(-2.8208371) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7587897) q[2];
sx q[2];
rz(-2.2586939) q[2];
sx q[2];
rz(-0.57304136) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2493077) q[1];
sx q[1];
rz(-1.882949) q[1];
sx q[1];
rz(2.409163) q[1];
rz(-pi) q[2];
rz(0.16964511) q[3];
sx q[3];
rz(-1.0153474) q[3];
sx q[3];
rz(2.5559705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(0.67374054) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(-1.3195066) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34981397) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(0.2579903) q[0];
rz(0.42516431) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.4917096) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8404322) q[0];
sx q[0];
rz(-1.4727117) q[0];
sx q[0];
rz(-2.7105986) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6368124) q[2];
sx q[2];
rz(-0.74540388) q[2];
sx q[2];
rz(-0.2573075) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.7548435) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(2.2350603) q[1];
x q[2];
rz(0.10156472) q[3];
sx q[3];
rz(-1.2079117) q[3];
sx q[3];
rz(-0.18728072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.012718) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(-0.84364676) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(-0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3180852) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(-0.41123018) q[0];
rz(-0.86589083) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(3.1076028) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995178) q[0];
sx q[0];
rz(-2.3410428) q[0];
sx q[0];
rz(-2.9151025) q[0];
rz(0.88862822) q[2];
sx q[2];
rz(-2.3805328) q[2];
sx q[2];
rz(-1.6737446) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6778292) q[1];
sx q[1];
rz(-0.67968183) q[1];
sx q[1];
rz(-1.4903279) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4230698) q[3];
sx q[3];
rz(-0.96102321) q[3];
sx q[3];
rz(2.8019398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5380481) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(2.2650488) q[2];
rz(-2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2440764) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(-0.39392719) q[0];
rz(-0.36755964) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(-1.6961018) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1720393) q[0];
sx q[0];
rz(-1.5656099) q[0];
sx q[0];
rz(2.0041549) q[0];
x q[1];
rz(-2.5567899) q[2];
sx q[2];
rz(-1.0324761) q[2];
sx q[2];
rz(-2.1540097) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3280914) q[1];
sx q[1];
rz(-2.5655167) q[1];
sx q[1];
rz(0.2920132) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1962542) q[3];
sx q[3];
rz(-2.1350386) q[3];
sx q[3];
rz(-2.4953147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8470856) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(0.40714804) q[2];
rz(-1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(-2.8919162) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3354934) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.8898213) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(0.30977419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99684925) q[0];
sx q[0];
rz(-2.092917) q[0];
sx q[0];
rz(1.4247308) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6236213) q[2];
sx q[2];
rz(-0.17715684) q[2];
sx q[2];
rz(1.0747386) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27948353) q[1];
sx q[1];
rz(-0.21394193) q[1];
sx q[1];
rz(-1.2941542) q[1];
rz(2.5449949) q[3];
sx q[3];
rz(-1.784424) q[3];
sx q[3];
rz(1.240977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.70242515) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(1.9343728) q[2];
rz(-2.1045945) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(-0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050215125) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(-1.2058831) q[0];
rz(0.58569113) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(1.4996128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8382032) q[0];
sx q[0];
rz(-1.2801542) q[0];
sx q[0];
rz(0.10586664) q[0];
x q[1];
rz(-2.3775616) q[2];
sx q[2];
rz(-1.521763) q[2];
sx q[2];
rz(-1.3438091) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8780898) q[1];
sx q[1];
rz(-0.51551688) q[1];
sx q[1];
rz(1.131119) q[1];
rz(-pi) q[2];
rz(1.5564735) q[3];
sx q[3];
rz(-2.9687772) q[3];
sx q[3];
rz(0.9848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.88400921) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(0.94669) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(-1.0735738) q[2];
sx q[2];
rz(-1.9855269) q[2];
sx q[2];
rz(-1.2863458) q[2];
rz(-2.4214217) q[3];
sx q[3];
rz(-1.9577033) q[3];
sx q[3];
rz(-2.4693558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];