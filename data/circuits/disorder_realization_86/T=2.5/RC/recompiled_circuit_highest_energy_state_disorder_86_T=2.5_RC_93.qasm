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
rz(0.64492172) q[0];
sx q[0];
rz(2.6725197) q[0];
sx q[0];
rz(7.2735431) q[0];
rz(-1.1733836) q[1];
sx q[1];
rz(-0.85964179) q[1];
sx q[1];
rz(-2.2923861) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3135017) q[0];
sx q[0];
rz(-1.5679616) q[0];
sx q[0];
rz(-1.536973) q[0];
x q[1];
rz(1.9687378) q[2];
sx q[2];
rz(-0.74021562) q[2];
sx q[2];
rz(-2.5876683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9473331) q[1];
sx q[1];
rz(-0.68945706) q[1];
sx q[1];
rz(1.1442776) q[1];
rz(-pi) q[2];
rz(-1.0484656) q[3];
sx q[3];
rz(-1.4712787) q[3];
sx q[3];
rz(-0.24302788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2363756) q[2];
sx q[2];
rz(-1.507501) q[2];
sx q[2];
rz(0.53984731) q[2];
rz(-1.5763464) q[3];
sx q[3];
rz(-0.40894517) q[3];
sx q[3];
rz(1.4096155) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.09963116) q[0];
sx q[0];
rz(-0.67622447) q[0];
sx q[0];
rz(1.8145632) q[0];
rz(-2.3565893) q[1];
sx q[1];
rz(-1.7186586) q[1];
sx q[1];
rz(0.50055093) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47704044) q[0];
sx q[0];
rz(-2.121976) q[0];
sx q[0];
rz(1.6001979) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2039866) q[2];
sx q[2];
rz(-1.3932835) q[2];
sx q[2];
rz(2.8457763) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2296648) q[1];
sx q[1];
rz(-1.2357792) q[1];
sx q[1];
rz(2.204621) q[1];
x q[2];
rz(-2.331892) q[3];
sx q[3];
rz(-1.9040344) q[3];
sx q[3];
rz(-0.35077318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0850204) q[2];
sx q[2];
rz(-0.73740059) q[2];
sx q[2];
rz(0.48193398) q[2];
rz(-2.7815172) q[3];
sx q[3];
rz(-1.1432546) q[3];
sx q[3];
rz(2.9329407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8423186) q[0];
sx q[0];
rz(-1.1085008) q[0];
sx q[0];
rz(0.47766787) q[0];
rz(0.85992366) q[1];
sx q[1];
rz(-1.0483024) q[1];
sx q[1];
rz(2.8544676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6159218) q[0];
sx q[0];
rz(-1.6369372) q[0];
sx q[0];
rz(-1.69709) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3541917) q[2];
sx q[2];
rz(-1.9176716) q[2];
sx q[2];
rz(1.5225449) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3928007) q[1];
sx q[1];
rz(-0.8872036) q[1];
sx q[1];
rz(0.011286126) q[1];
rz(-pi) q[2];
rz(-1.1895387) q[3];
sx q[3];
rz(-1.4917456) q[3];
sx q[3];
rz(1.6565269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3890248) q[2];
sx q[2];
rz(-1.3501046) q[2];
sx q[2];
rz(-3.1217421) q[2];
rz(0.94414532) q[3];
sx q[3];
rz(-2.4343334) q[3];
sx q[3];
rz(2.7437362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860564) q[0];
sx q[0];
rz(-2.5113386) q[0];
sx q[0];
rz(1.3007042) q[0];
rz(-0.4862673) q[1];
sx q[1];
rz(-1.174289) q[1];
sx q[1];
rz(3.1062612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6636642) q[0];
sx q[0];
rz(-1.1424756) q[0];
sx q[0];
rz(3.1129063) q[0];
rz(-1.3372806) q[2];
sx q[2];
rz(-1.9583251) q[2];
sx q[2];
rz(3.1336409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17745397) q[1];
sx q[1];
rz(-1.576401) q[1];
sx q[1];
rz(0.023375794) q[1];
x q[2];
rz(-2.5132645) q[3];
sx q[3];
rz(-2.0355967) q[3];
sx q[3];
rz(-2.2006054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6639634) q[2];
sx q[2];
rz(-2.5783381) q[2];
sx q[2];
rz(-2.8832054) q[2];
rz(0.7102617) q[3];
sx q[3];
rz(-0.60451549) q[3];
sx q[3];
rz(1.5822423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6579987) q[0];
sx q[0];
rz(-0.76376629) q[0];
sx q[0];
rz(-0.64436954) q[0];
rz(-0.98980347) q[1];
sx q[1];
rz(-0.80392307) q[1];
sx q[1];
rz(-1.1092626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98609867) q[0];
sx q[0];
rz(-1.7664096) q[0];
sx q[0];
rz(-1.0367111) q[0];
rz(-3.101966) q[2];
sx q[2];
rz(-2.0840696) q[2];
sx q[2];
rz(0.64943343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7444744) q[1];
sx q[1];
rz(-2.3496685) q[1];
sx q[1];
rz(2.9379815) q[1];
x q[2];
rz(-0.30410366) q[3];
sx q[3];
rz(-2.3652206) q[3];
sx q[3];
rz(2.9714288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41205078) q[2];
sx q[2];
rz(-1.5645626) q[2];
sx q[2];
rz(2.5315419) q[2];
rz(0.080502056) q[3];
sx q[3];
rz(-2.9952315) q[3];
sx q[3];
rz(3.1350873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42008156) q[0];
sx q[0];
rz(-2.2048936) q[0];
sx q[0];
rz(-1.6625846) q[0];
rz(-1.1889907) q[1];
sx q[1];
rz(-0.34591302) q[1];
sx q[1];
rz(-1.4422653) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8540031) q[0];
sx q[0];
rz(-0.95335273) q[0];
sx q[0];
rz(2.94728) q[0];
x q[1];
rz(3.0411554) q[2];
sx q[2];
rz(-2.5578376) q[2];
sx q[2];
rz(0.69086087) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0008322) q[1];
sx q[1];
rz(-2.3545583) q[1];
sx q[1];
rz(2.0443022) q[1];
rz(0.53621323) q[3];
sx q[3];
rz(-1.271476) q[3];
sx q[3];
rz(1.0992536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6815765) q[2];
sx q[2];
rz(-0.67395335) q[2];
sx q[2];
rz(0.67503929) q[2];
rz(-0.33440822) q[3];
sx q[3];
rz(-1.6655917) q[3];
sx q[3];
rz(-0.36791754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6595031) q[0];
sx q[0];
rz(-0.98564321) q[0];
sx q[0];
rz(-0.12568812) q[0];
rz(0.19371678) q[1];
sx q[1];
rz(-2.2592762) q[1];
sx q[1];
rz(-1.151459) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2053368) q[0];
sx q[0];
rz(-2.057862) q[0];
sx q[0];
rz(-2.4508935) q[0];
rz(1.9415683) q[2];
sx q[2];
rz(-1.3091012) q[2];
sx q[2];
rz(0.38693025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9509778) q[1];
sx q[1];
rz(-1.1314609) q[1];
sx q[1];
rz(0.5885891) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56834055) q[3];
sx q[3];
rz(-2.3821444) q[3];
sx q[3];
rz(2.8041149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8573528) q[2];
sx q[2];
rz(-1.7809296) q[2];
sx q[2];
rz(-2.8947158) q[2];
rz(2.2829368) q[3];
sx q[3];
rz(-2.1571428) q[3];
sx q[3];
rz(2.8624559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7849279) q[0];
sx q[0];
rz(-0.44342884) q[0];
sx q[0];
rz(0.32522935) q[0];
rz(-2.6204956) q[1];
sx q[1];
rz(-2.5083713) q[1];
sx q[1];
rz(2.8281143) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26599182) q[0];
sx q[0];
rz(-0.82339215) q[0];
sx q[0];
rz(-0.70988795) q[0];
rz(-1.283139) q[2];
sx q[2];
rz(-0.23459841) q[2];
sx q[2];
rz(1.2755659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.65113803) q[1];
sx q[1];
rz(-0.50781194) q[1];
sx q[1];
rz(-0.74393028) q[1];
x q[2];
rz(1.0495124) q[3];
sx q[3];
rz(-1.3950328) q[3];
sx q[3];
rz(-2.5960032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43361214) q[2];
sx q[2];
rz(-1.2182451) q[2];
sx q[2];
rz(1.2797959) q[2];
rz(-0.72569877) q[3];
sx q[3];
rz(-1.9819219) q[3];
sx q[3];
rz(-2.3316135) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5028266) q[0];
sx q[0];
rz(-1.9483197) q[0];
sx q[0];
rz(-2.9916812) q[0];
rz(-1.8984849) q[1];
sx q[1];
rz(-1.5538235) q[1];
sx q[1];
rz(1.6601723) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7090764) q[0];
sx q[0];
rz(-0.066532739) q[0];
sx q[0];
rz(2.3155022) q[0];
rz(1.7779839) q[2];
sx q[2];
rz(-2.8958791) q[2];
sx q[2];
rz(0.65226511) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31563035) q[1];
sx q[1];
rz(-0.44188979) q[1];
sx q[1];
rz(2.619584) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0444938) q[3];
sx q[3];
rz(-2.8682531) q[3];
sx q[3];
rz(2.7784082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.82568613) q[2];
sx q[2];
rz(-1.3797727) q[2];
sx q[2];
rz(-2.1659577) q[2];
rz(-1.1286831) q[3];
sx q[3];
rz(-0.55207878) q[3];
sx q[3];
rz(-0.25477195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9752556) q[0];
sx q[0];
rz(-1.0003426) q[0];
sx q[0];
rz(-2.9794203) q[0];
rz(-1.8099248) q[1];
sx q[1];
rz(-2.5089896) q[1];
sx q[1];
rz(3.0955637) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23102681) q[0];
sx q[0];
rz(-1.5767158) q[0];
sx q[0];
rz(-1.5529412) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9580905) q[2];
sx q[2];
rz(-2.0047054) q[2];
sx q[2];
rz(2.6604685) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.522559) q[1];
sx q[1];
rz(-1.4231893) q[1];
sx q[1];
rz(-1.9760494) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8022553) q[3];
sx q[3];
rz(-1.7320314) q[3];
sx q[3];
rz(1.1649023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8056246) q[2];
sx q[2];
rz(-0.38463548) q[2];
sx q[2];
rz(-1.1495205) q[2];
rz(-0.51291054) q[3];
sx q[3];
rz(-1.3866813) q[3];
sx q[3];
rz(0.30470595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.649986) q[0];
sx q[0];
rz(-2.2254324) q[0];
sx q[0];
rz(-1.9944763) q[0];
rz(-0.06123771) q[1];
sx q[1];
rz(-1.1920659) q[1];
sx q[1];
rz(1.7658284) q[1];
rz(0.55228615) q[2];
sx q[2];
rz(-0.88197642) q[2];
sx q[2];
rz(-2.3563202) q[2];
rz(-1.9763532) q[3];
sx q[3];
rz(-0.83409393) q[3];
sx q[3];
rz(0.50611151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
