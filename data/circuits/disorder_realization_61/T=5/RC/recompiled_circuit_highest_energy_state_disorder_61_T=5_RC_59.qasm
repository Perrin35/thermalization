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
rz(0.8340303) q[0];
sx q[0];
rz(-2.7355255) q[0];
sx q[0];
rz(1.9598444) q[0];
rz(-1.4015247) q[1];
sx q[1];
rz(-2.6670246) q[1];
sx q[1];
rz(3.0085556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33677732) q[0];
sx q[0];
rz(-2.204737) q[0];
sx q[0];
rz(-2.8147349) q[0];
rz(-pi) q[1];
rz(-0.51271397) q[2];
sx q[2];
rz(-1.9641335) q[2];
sx q[2];
rz(-0.31189298) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.052357167) q[1];
sx q[1];
rz(-1.3468139) q[1];
sx q[1];
rz(1.2473945) q[1];
rz(1.0361133) q[3];
sx q[3];
rz(-2.4991192) q[3];
sx q[3];
rz(1.5758663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0066321) q[2];
sx q[2];
rz(-1.0198318) q[2];
sx q[2];
rz(-0.68438619) q[2];
rz(1.1413752) q[3];
sx q[3];
rz(-0.26641521) q[3];
sx q[3];
rz(-3.0548837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3619277) q[0];
sx q[0];
rz(-0.68988887) q[0];
sx q[0];
rz(-0.46098125) q[0];
rz(2.0224679) q[1];
sx q[1];
rz(-0.42367595) q[1];
sx q[1];
rz(2.8295595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49544612) q[0];
sx q[0];
rz(-1.5987885) q[0];
sx q[0];
rz(-1.5585446) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3316595) q[2];
sx q[2];
rz(-1.51568) q[2];
sx q[2];
rz(2.8365561) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2300517) q[1];
sx q[1];
rz(-1.5649867) q[1];
sx q[1];
rz(-1.7464386) q[1];
rz(-pi) q[2];
rz(-1.2106522) q[3];
sx q[3];
rz(-2.4575298) q[3];
sx q[3];
rz(1.6863914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95388874) q[2];
sx q[2];
rz(-1.1053332) q[2];
sx q[2];
rz(2.7375431) q[2];
rz(-0.36014253) q[3];
sx q[3];
rz(-1.6481684) q[3];
sx q[3];
rz(-2.4586316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3563511) q[0];
sx q[0];
rz(-1.0026362) q[0];
sx q[0];
rz(-0.32401618) q[0];
rz(1.0600545) q[1];
sx q[1];
rz(-2.8476604) q[1];
sx q[1];
rz(-2.4411328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4361079) q[0];
sx q[0];
rz(-1.7388845) q[0];
sx q[0];
rz(1.4072493) q[0];
rz(-pi) q[1];
rz(0.21493463) q[2];
sx q[2];
rz(-1.9086468) q[2];
sx q[2];
rz(-0.95906729) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6571275) q[1];
sx q[1];
rz(-1.336369) q[1];
sx q[1];
rz(-2.3695494) q[1];
rz(-pi) q[2];
rz(-1.992393) q[3];
sx q[3];
rz(-1.3918878) q[3];
sx q[3];
rz(0.94551555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5141653) q[2];
sx q[2];
rz(-0.24719396) q[2];
sx q[2];
rz(0.79666454) q[2];
rz(-1.1967777) q[3];
sx q[3];
rz(-1.5408885) q[3];
sx q[3];
rz(0.59320199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.45348039) q[0];
sx q[0];
rz(-2.1402833) q[0];
sx q[0];
rz(1.8185115) q[0];
rz(-2.6755013) q[1];
sx q[1];
rz(-2.2345462) q[1];
sx q[1];
rz(1.1493433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1150165) q[0];
sx q[0];
rz(-2.0043902) q[0];
sx q[0];
rz(-1.0267016) q[0];
rz(-0.39300275) q[2];
sx q[2];
rz(-0.98969029) q[2];
sx q[2];
rz(0.71497922) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9117078) q[1];
sx q[1];
rz(-2.9344892) q[1];
sx q[1];
rz(-2.9207499) q[1];
rz(1.9051311) q[3];
sx q[3];
rz(-1.4282246) q[3];
sx q[3];
rz(-2.4041374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.65679067) q[2];
sx q[2];
rz(-2.4718086) q[2];
sx q[2];
rz(0.82526866) q[2];
rz(-1.2800062) q[3];
sx q[3];
rz(-1.5214336) q[3];
sx q[3];
rz(0.72117225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.698302) q[0];
sx q[0];
rz(-1.8809141) q[0];
sx q[0];
rz(1.3662421) q[0];
rz(2.0514936) q[1];
sx q[1];
rz(-1.5049728) q[1];
sx q[1];
rz(-1.8025788) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60581029) q[0];
sx q[0];
rz(-0.92180862) q[0];
sx q[0];
rz(1.3507516) q[0];
rz(0.86947415) q[2];
sx q[2];
rz(-0.83842248) q[2];
sx q[2];
rz(-0.43763197) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80833188) q[1];
sx q[1];
rz(-1.3587316) q[1];
sx q[1];
rz(-1.6912837) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9129778) q[3];
sx q[3];
rz(-1.2552731) q[3];
sx q[3];
rz(-0.70759892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1947386) q[2];
sx q[2];
rz(-1.7297435) q[2];
sx q[2];
rz(2.8080688) q[2];
rz(0.58878318) q[3];
sx q[3];
rz(-2.6684561) q[3];
sx q[3];
rz(-2.3401882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3280585) q[0];
sx q[0];
rz(-1.7840339) q[0];
sx q[0];
rz(1.300977) q[0];
rz(2.1133555) q[1];
sx q[1];
rz(-2.6515549) q[1];
sx q[1];
rz(2.6992544) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9755632) q[0];
sx q[0];
rz(-2.0394169) q[0];
sx q[0];
rz(1.2965078) q[0];
rz(0.47719854) q[2];
sx q[2];
rz(-1.0893359) q[2];
sx q[2];
rz(1.9419326) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7117911) q[1];
sx q[1];
rz(-0.45908976) q[1];
sx q[1];
rz(-2.7613822) q[1];
rz(-pi) q[2];
rz(-0.30486561) q[3];
sx q[3];
rz(-2.8072522) q[3];
sx q[3];
rz(1.5078733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3698795) q[2];
sx q[2];
rz(-1.5481202) q[2];
sx q[2];
rz(0.4717353) q[2];
rz(-1.6996023) q[3];
sx q[3];
rz(-0.56157464) q[3];
sx q[3];
rz(-2.877318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1137375) q[0];
sx q[0];
rz(-1.4677445) q[0];
sx q[0];
rz(-2.9748528) q[0];
rz(2.3475504) q[1];
sx q[1];
rz(-1.200095) q[1];
sx q[1];
rz(-3.0418495) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8309495) q[0];
sx q[0];
rz(-1.8825873) q[0];
sx q[0];
rz(2.1556709) q[0];
rz(-pi) q[1];
rz(-1.0641618) q[2];
sx q[2];
rz(-0.68175137) q[2];
sx q[2];
rz(-0.1296986) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7078898) q[1];
sx q[1];
rz(-1.4447482) q[1];
sx q[1];
rz(-2.5688102) q[1];
rz(1.0519299) q[3];
sx q[3];
rz(-2.4104033) q[3];
sx q[3];
rz(-0.90381453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86228592) q[2];
sx q[2];
rz(-1.0655094) q[2];
sx q[2];
rz(0.12969895) q[2];
rz(-1.8779523) q[3];
sx q[3];
rz(-1.2010776) q[3];
sx q[3];
rz(-0.14550801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97067863) q[0];
sx q[0];
rz(-0.92249528) q[0];
sx q[0];
rz(2.7150735) q[0];
rz(2.049394) q[1];
sx q[1];
rz(-2.009232) q[1];
sx q[1];
rz(2.138864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4115017) q[0];
sx q[0];
rz(-0.97785946) q[0];
sx q[0];
rz(0.31296313) q[0];
rz(-pi) q[1];
rz(3.0363704) q[2];
sx q[2];
rz(-2.4159523) q[2];
sx q[2];
rz(1.2187472) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7510609) q[1];
sx q[1];
rz(-0.56274429) q[1];
sx q[1];
rz(0.90969245) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2360446) q[3];
sx q[3];
rz(-1.1268643) q[3];
sx q[3];
rz(0.027396552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2670474) q[2];
sx q[2];
rz(-0.55570498) q[2];
sx q[2];
rz(-3.0181273) q[2];
rz(-0.63001436) q[3];
sx q[3];
rz(-1.8450626) q[3];
sx q[3];
rz(2.4251078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89710871) q[0];
sx q[0];
rz(-0.69863313) q[0];
sx q[0];
rz(0.83936349) q[0];
rz(-0.15946236) q[1];
sx q[1];
rz(-2.1493252) q[1];
sx q[1];
rz(-0.92963591) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2004505) q[0];
sx q[0];
rz(-1.8243891) q[0];
sx q[0];
rz(-2.3932049) q[0];
rz(2.9603297) q[2];
sx q[2];
rz(-0.59872228) q[2];
sx q[2];
rz(-1.3186272) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8690455) q[1];
sx q[1];
rz(-0.89249498) q[1];
sx q[1];
rz(-0.1609623) q[1];
x q[2];
rz(1.2655018) q[3];
sx q[3];
rz(-2.4809894) q[3];
sx q[3];
rz(1.6317639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.9440426) q[2];
sx q[2];
rz(-1.9746747) q[2];
sx q[2];
rz(-2.3738532) q[2];
rz(-0.36457101) q[3];
sx q[3];
rz(-1.7578099) q[3];
sx q[3];
rz(0.069469623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0070873) q[0];
sx q[0];
rz(-0.61408478) q[0];
sx q[0];
rz(-2.2579204) q[0];
rz(-2.9332352) q[1];
sx q[1];
rz(-0.50269428) q[1];
sx q[1];
rz(1.3349894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55157298) q[0];
sx q[0];
rz(-3.0544222) q[0];
sx q[0];
rz(2.4159191) q[0];
rz(0.31317385) q[2];
sx q[2];
rz(-2.5921481) q[2];
sx q[2];
rz(-0.012481364) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.523004) q[1];
sx q[1];
rz(-2.4135313) q[1];
sx q[1];
rz(-2.2656206) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2741256) q[3];
sx q[3];
rz(-0.88109479) q[3];
sx q[3];
rz(2.3610601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.16365446) q[2];
sx q[2];
rz(-1.6341219) q[2];
sx q[2];
rz(0.046023544) q[2];
rz(0.16511354) q[3];
sx q[3];
rz(-0.29296994) q[3];
sx q[3];
rz(1.9273812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5788427) q[0];
sx q[0];
rz(-0.10012193) q[0];
sx q[0];
rz(-1.9520676) q[0];
rz(2.9764755) q[1];
sx q[1];
rz(-1.4585635) q[1];
sx q[1];
rz(-0.30452902) q[1];
rz(-1.2444185) q[2];
sx q[2];
rz(-1.9934629) q[2];
sx q[2];
rz(-1.6531682) q[2];
rz(2.0863927) q[3];
sx q[3];
rz(-1.1451117) q[3];
sx q[3];
rz(-2.3449863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
