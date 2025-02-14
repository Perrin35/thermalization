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
rz(1.2827058) q[0];
sx q[0];
rz(-2.1252706) q[0];
sx q[0];
rz(-1.2555726) q[0];
rz(-1.3610871) q[1];
sx q[1];
rz(-0.95870107) q[1];
sx q[1];
rz(-2.5348742) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17598303) q[0];
sx q[0];
rz(-0.85411607) q[0];
sx q[0];
rz(0.59881439) q[0];
rz(3.1090281) q[2];
sx q[2];
rz(-1.3059907) q[2];
sx q[2];
rz(0.84003583) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.052362927) q[1];
sx q[1];
rz(-0.38859146) q[1];
sx q[1];
rz(-1.2227538) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64067556) q[3];
sx q[3];
rz(-2.9417613) q[3];
sx q[3];
rz(-2.8175648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9386193) q[2];
sx q[2];
rz(-1.9846658) q[2];
sx q[2];
rz(0.17091664) q[2];
rz(-2.0140698) q[3];
sx q[3];
rz(-0.5564965) q[3];
sx q[3];
rz(3.0881622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4776066) q[0];
sx q[0];
rz(-2.8150616) q[0];
sx q[0];
rz(0.68847454) q[0];
rz(1.3779878) q[1];
sx q[1];
rz(-2.1207899) q[1];
sx q[1];
rz(-3.0000906) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6571638) q[0];
sx q[0];
rz(-1.9701029) q[0];
sx q[0];
rz(-0.54175538) q[0];
rz(-2.7121041) q[2];
sx q[2];
rz(-1.4802684) q[2];
sx q[2];
rz(-1.8877826) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8162839) q[1];
sx q[1];
rz(-2.6806596) q[1];
sx q[1];
rz(-0.67616762) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9555952) q[3];
sx q[3];
rz(-1.4963004) q[3];
sx q[3];
rz(-0.39922478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7576393) q[2];
sx q[2];
rz(-0.44718224) q[2];
sx q[2];
rz(3.1128856) q[2];
rz(2.7590397) q[3];
sx q[3];
rz(-1.6304723) q[3];
sx q[3];
rz(0.20146519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460227) q[0];
sx q[0];
rz(-1.0423132) q[0];
sx q[0];
rz(2.7699455) q[0];
rz(0.21013513) q[1];
sx q[1];
rz(-2.7919283) q[1];
sx q[1];
rz(1.2079316) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2164566) q[0];
sx q[0];
rz(-2.0691923) q[0];
sx q[0];
rz(0.67277661) q[0];
x q[1];
rz(-1.5919331) q[2];
sx q[2];
rz(-1.653228) q[2];
sx q[2];
rz(1.3940514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4000798) q[1];
sx q[1];
rz(-0.99435213) q[1];
sx q[1];
rz(0.6735354) q[1];
rz(-0.34545107) q[3];
sx q[3];
rz(-1.6092704) q[3];
sx q[3];
rz(1.0055804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62260425) q[2];
sx q[2];
rz(-3.0486139) q[2];
sx q[2];
rz(-2.5034215) q[2];
rz(2.7291164) q[3];
sx q[3];
rz(-1.4700438) q[3];
sx q[3];
rz(-1.1023869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07688044) q[0];
sx q[0];
rz(-2.2314254) q[0];
sx q[0];
rz(1.8338715) q[0];
rz(-0.67963302) q[1];
sx q[1];
rz(-2.4145587) q[1];
sx q[1];
rz(-1.1839428) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1935496) q[0];
sx q[0];
rz(-1.8947161) q[0];
sx q[0];
rz(-0.69635038) q[0];
x q[1];
rz(-2.6908801) q[2];
sx q[2];
rz(-2.8808184) q[2];
sx q[2];
rz(2.5058172) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19540994) q[1];
sx q[1];
rz(-1.2210746) q[1];
sx q[1];
rz(-1.6433952) q[1];
rz(-pi) q[2];
x q[2];
rz(1.697942) q[3];
sx q[3];
rz(-1.4411488) q[3];
sx q[3];
rz(-1.6693019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1188941) q[2];
sx q[2];
rz(-2.9939632) q[2];
sx q[2];
rz(-1.0559319) q[2];
rz(0.28327709) q[3];
sx q[3];
rz(-1.2669468) q[3];
sx q[3];
rz(-1.7578846) q[3];
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
rz(-pi/2) q[0];
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
rz(-3.0530171) q[0];
sx q[0];
rz(-0.96240369) q[0];
sx q[0];
rz(-1.4085294) q[0];
rz(2.9119496) q[1];
sx q[1];
rz(-0.85669986) q[1];
sx q[1];
rz(-1.1913258) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2076805) q[0];
sx q[0];
rz(-1.7125074) q[0];
sx q[0];
rz(1.974154) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8764087) q[2];
sx q[2];
rz(-1.4508045) q[2];
sx q[2];
rz(0.12745276) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31684418) q[1];
sx q[1];
rz(-1.1543546) q[1];
sx q[1];
rz(1.29711) q[1];
rz(2.5727083) q[3];
sx q[3];
rz(-0.99939049) q[3];
sx q[3];
rz(-0.4898773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4952937) q[2];
sx q[2];
rz(-0.36094347) q[2];
sx q[2];
rz(1.4781282) q[2];
rz(0.16119257) q[3];
sx q[3];
rz(-1.5321923) q[3];
sx q[3];
rz(-0.78981367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6394871) q[0];
sx q[0];
rz(-2.7171071) q[0];
sx q[0];
rz(-2.2232527) q[0];
rz(-0.25686747) q[1];
sx q[1];
rz(-0.99595064) q[1];
sx q[1];
rz(-1.1161944) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74280069) q[0];
sx q[0];
rz(-2.0896119) q[0];
sx q[0];
rz(-2.1543845) q[0];
rz(2.8681197) q[2];
sx q[2];
rz(-1.6533706) q[2];
sx q[2];
rz(-1.7805748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0829385) q[1];
sx q[1];
rz(-1.3320597) q[1];
sx q[1];
rz(-1.7708562) q[1];
rz(2.5258216) q[3];
sx q[3];
rz(-0.29524657) q[3];
sx q[3];
rz(-2.8930882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9548698) q[2];
sx q[2];
rz(-0.8064417) q[2];
sx q[2];
rz(-2.7346129) q[2];
rz(2.5278029) q[3];
sx q[3];
rz(-1.611462) q[3];
sx q[3];
rz(0.49155244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75519049) q[0];
sx q[0];
rz(-0.84119216) q[0];
sx q[0];
rz(-2.2156773) q[0];
rz(2.6689802) q[1];
sx q[1];
rz(-2.0874529) q[1];
sx q[1];
rz(2.6714163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42379728) q[0];
sx q[0];
rz(-1.6542302) q[0];
sx q[0];
rz(-0.83441894) q[0];
rz(-pi) q[1];
rz(2.2165259) q[2];
sx q[2];
rz(-1.2534516) q[2];
sx q[2];
rz(2.6982174) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8308816) q[1];
sx q[1];
rz(-2.355) q[1];
sx q[1];
rz(1.6221694) q[1];
x q[2];
rz(-0.84104611) q[3];
sx q[3];
rz(-1.5435092) q[3];
sx q[3];
rz(-0.59852615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65972313) q[2];
sx q[2];
rz(-2.0707776) q[2];
sx q[2];
rz(2.3107963) q[2];
rz(0.080549084) q[3];
sx q[3];
rz(-2.060067) q[3];
sx q[3];
rz(-0.95675937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.047121) q[0];
sx q[0];
rz(-1.5249277) q[0];
sx q[0];
rz(-2.9834874) q[0];
rz(-2.415601) q[1];
sx q[1];
rz(-2.7114365) q[1];
sx q[1];
rz(-0.37011883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67664941) q[0];
sx q[0];
rz(-1.2386564) q[0];
sx q[0];
rz(1.667883) q[0];
rz(-pi) q[1];
rz(2.6441755) q[2];
sx q[2];
rz(-0.82225613) q[2];
sx q[2];
rz(-2.9679839) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.98277011) q[1];
sx q[1];
rz(-1.9590322) q[1];
sx q[1];
rz(-2.7169711) q[1];
x q[2];
rz(-2.1158989) q[3];
sx q[3];
rz(-1.926356) q[3];
sx q[3];
rz(-2.4073383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.98081723) q[2];
sx q[2];
rz(-1.2951916) q[2];
sx q[2];
rz(3.004461) q[2];
rz(1.49617) q[3];
sx q[3];
rz(-0.6936332) q[3];
sx q[3];
rz(0.48367286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5982323) q[0];
sx q[0];
rz(-2.1393263) q[0];
sx q[0];
rz(0.12635669) q[0];
rz(2.5670746) q[1];
sx q[1];
rz(-2.2210329) q[1];
sx q[1];
rz(-0.7472907) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80602801) q[0];
sx q[0];
rz(-1.3845516) q[0];
sx q[0];
rz(1.1769017) q[0];
x q[1];
rz(-0.70238446) q[2];
sx q[2];
rz(-1.1051854) q[2];
sx q[2];
rz(-1.3872176) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.14366985) q[1];
sx q[1];
rz(-1.2005245) q[1];
sx q[1];
rz(-2.0475494) q[1];
rz(-pi) q[2];
rz(-0.26227343) q[3];
sx q[3];
rz(-1.4424535) q[3];
sx q[3];
rz(2.195829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1025461) q[2];
sx q[2];
rz(-1.1209844) q[2];
sx q[2];
rz(-0.54753629) q[2];
rz(-1.3012137) q[3];
sx q[3];
rz(-1.1373212) q[3];
sx q[3];
rz(-0.12714061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8319156) q[0];
sx q[0];
rz(-1.7597821) q[0];
sx q[0];
rz(0.58050138) q[0];
rz(-2.8864587) q[1];
sx q[1];
rz(-2.3105123) q[1];
sx q[1];
rz(-1.1855804) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0190576) q[0];
sx q[0];
rz(-2.2154729) q[0];
sx q[0];
rz(0.77150795) q[0];
rz(-pi) q[1];
rz(-0.0015179356) q[2];
sx q[2];
rz(-1.9337144) q[2];
sx q[2];
rz(-1.5134821) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.285814) q[1];
sx q[1];
rz(-1.6596982) q[1];
sx q[1];
rz(-0.47089058) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85082741) q[3];
sx q[3];
rz(-1.4113857) q[3];
sx q[3];
rz(-0.042378332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.4127976) q[2];
sx q[2];
rz(-2.2788861) q[2];
sx q[2];
rz(1.0069138) q[2];
rz(-1.5236731) q[3];
sx q[3];
rz(-2.1606052) q[3];
sx q[3];
rz(-1.1716918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0278052) q[0];
sx q[0];
rz(-1.2600949) q[0];
sx q[0];
rz(-2.8366198) q[0];
rz(-1.2420568) q[1];
sx q[1];
rz(-1.8546974) q[1];
sx q[1];
rz(0.7484662) q[1];
rz(1.2075049) q[2];
sx q[2];
rz(-0.45798326) q[2];
sx q[2];
rz(1.9275673) q[2];
rz(0.92963582) q[3];
sx q[3];
rz(-2.2916678) q[3];
sx q[3];
rz(-2.5558932) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
