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
rz(-1.6031185) q[0];
sx q[0];
rz(-1.7552019) q[0];
sx q[0];
rz(1.2586799) q[0];
rz(-1.7006682) q[1];
sx q[1];
rz(-0.52803334) q[1];
sx q[1];
rz(-0.33906403) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0345732) q[0];
sx q[0];
rz(-1.2385173) q[0];
sx q[0];
rz(0.1531202) q[0];
rz(2.6554675) q[2];
sx q[2];
rz(-0.74084254) q[2];
sx q[2];
rz(-1.6227826) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5524774) q[1];
sx q[1];
rz(-2.5259645) q[1];
sx q[1];
rz(1.6861123) q[1];
rz(-pi) q[2];
rz(1.6293832) q[3];
sx q[3];
rz(-2.226429) q[3];
sx q[3];
rz(2.7215304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6429834) q[2];
sx q[2];
rz(-2.6080841) q[2];
sx q[2];
rz(-1.3996997) q[2];
rz(-1.1534322) q[3];
sx q[3];
rz(-1.6250936) q[3];
sx q[3];
rz(-1.7216518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72700778) q[0];
sx q[0];
rz(-1.6731394) q[0];
sx q[0];
rz(2.6732095) q[0];
rz(0.98764694) q[1];
sx q[1];
rz(-1.3750319) q[1];
sx q[1];
rz(-2.1013026) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066369206) q[0];
sx q[0];
rz(-2.1589737) q[0];
sx q[0];
rz(2.4494314) q[0];
x q[1];
rz(-1.9557421) q[2];
sx q[2];
rz(-1.9296807) q[2];
sx q[2];
rz(2.5449616) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4978943) q[1];
sx q[1];
rz(-2.0308745) q[1];
sx q[1];
rz(-1.6804477) q[1];
x q[2];
rz(2.1352102) q[3];
sx q[3];
rz(-1.4335535) q[3];
sx q[3];
rz(-2.8999527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4686244) q[2];
sx q[2];
rz(-1.1138223) q[2];
sx q[2];
rz(1.683453) q[2];
rz(-0.80312076) q[3];
sx q[3];
rz(-1.8773269) q[3];
sx q[3];
rz(-1.7709812) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9010381) q[0];
sx q[0];
rz(-0.63837186) q[0];
sx q[0];
rz(-1.4580131) q[0];
rz(1.1395678) q[1];
sx q[1];
rz(-1.5007867) q[1];
sx q[1];
rz(2.4511888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6891606) q[0];
sx q[0];
rz(-2.3827973) q[0];
sx q[0];
rz(-0.62907416) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1369526) q[2];
sx q[2];
rz(-1.829471) q[2];
sx q[2];
rz(-2.0438922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1377167) q[1];
sx q[1];
rz(-2.5270577) q[1];
sx q[1];
rz(2.4412066) q[1];
rz(-2.2182057) q[3];
sx q[3];
rz(-1.5316803) q[3];
sx q[3];
rz(-0.95992328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0475426) q[2];
sx q[2];
rz(-1.2529145) q[2];
sx q[2];
rz(2.3599153) q[2];
rz(2.9983669) q[3];
sx q[3];
rz(-2.4476624) q[3];
sx q[3];
rz(0.028698746) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52611065) q[0];
sx q[0];
rz(-0.019793864) q[0];
sx q[0];
rz(0.97446781) q[0];
rz(0.027007667) q[1];
sx q[1];
rz(-1.6504495) q[1];
sx q[1];
rz(2.9496121) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0483069) q[0];
sx q[0];
rz(-1.9294717) q[0];
sx q[0];
rz(-1.8762588) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5938283) q[2];
sx q[2];
rz(-0.27850391) q[2];
sx q[2];
rz(0.26879569) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0897191) q[1];
sx q[1];
rz(-2.8579144) q[1];
sx q[1];
rz(-1.5121626) q[1];
rz(0.2918029) q[3];
sx q[3];
rz(-2.553396) q[3];
sx q[3];
rz(3.0384529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.30569046) q[2];
sx q[2];
rz(-1.4985936) q[2];
sx q[2];
rz(2.4198325) q[2];
rz(2.6137958) q[3];
sx q[3];
rz(-0.58776394) q[3];
sx q[3];
rz(1.2036948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0498306) q[0];
sx q[0];
rz(-1.993655) q[0];
sx q[0];
rz(-2.9312396) q[0];
rz(2.6166088) q[1];
sx q[1];
rz(-2.4139082) q[1];
sx q[1];
rz(-2.4997589) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1686629) q[0];
sx q[0];
rz(-2.7073858) q[0];
sx q[0];
rz(-0.28247873) q[0];
rz(-2.9027356) q[2];
sx q[2];
rz(-1.3625267) q[2];
sx q[2];
rz(1.3778401) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.32165471) q[1];
sx q[1];
rz(-1.7366855) q[1];
sx q[1];
rz(-0.031601918) q[1];
rz(-1.8689496) q[3];
sx q[3];
rz(-1.5437727) q[3];
sx q[3];
rz(-2.2397796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.019023808) q[2];
sx q[2];
rz(-1.6876612) q[2];
sx q[2];
rz(-0.34073487) q[2];
rz(2.6797898) q[3];
sx q[3];
rz(-1.267649) q[3];
sx q[3];
rz(0.99892282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13071624) q[0];
sx q[0];
rz(-1.0557446) q[0];
sx q[0];
rz(-0.51489949) q[0];
rz(1.8574832) q[1];
sx q[1];
rz(-1.037037) q[1];
sx q[1];
rz(2.5426755) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8209957) q[0];
sx q[0];
rz(-0.55984523) q[0];
sx q[0];
rz(1.2843156) q[0];
rz(1.7755505) q[2];
sx q[2];
rz(-1.6806112) q[2];
sx q[2];
rz(1.5653953) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1753483) q[1];
sx q[1];
rz(-1.6754901) q[1];
sx q[1];
rz(2.1123516) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1647366) q[3];
sx q[3];
rz(-0.58771509) q[3];
sx q[3];
rz(2.5669861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9942223) q[2];
sx q[2];
rz(-1.3753336) q[2];
sx q[2];
rz(-0.71225172) q[2];
rz(1.1143403) q[3];
sx q[3];
rz(-2.2388191) q[3];
sx q[3];
rz(-2.4207992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095605843) q[0];
sx q[0];
rz(-1.1271891) q[0];
sx q[0];
rz(-1.7186694) q[0];
rz(-1.0320041) q[1];
sx q[1];
rz(-1.2675266) q[1];
sx q[1];
rz(1.9292018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40875567) q[0];
sx q[0];
rz(-1.2877595) q[0];
sx q[0];
rz(0.7247337) q[0];
x q[1];
rz(2.5676377) q[2];
sx q[2];
rz(-0.35517737) q[2];
sx q[2];
rz(-1.4188302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29323762) q[1];
sx q[1];
rz(-2.2841286) q[1];
sx q[1];
rz(2.7687702) q[1];
rz(-0.29799975) q[3];
sx q[3];
rz(-1.2825507) q[3];
sx q[3];
rz(0.73747613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.57264868) q[2];
sx q[2];
rz(-2.0062916) q[2];
sx q[2];
rz(0.14986077) q[2];
rz(3.0585238) q[3];
sx q[3];
rz(-0.84851256) q[3];
sx q[3];
rz(1.8547295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.5158373) q[0];
sx q[0];
rz(-2.5560162) q[0];
sx q[0];
rz(2.9199912) q[0];
rz(-1.8009456) q[1];
sx q[1];
rz(-2.4602175) q[1];
sx q[1];
rz(-0.63366205) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1213729) q[0];
sx q[0];
rz(-0.069761666) q[0];
sx q[0];
rz(-2.1345977) q[0];
rz(-pi) q[1];
rz(-2.7886798) q[2];
sx q[2];
rz(-2.196759) q[2];
sx q[2];
rz(-2.7865648) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3089655) q[1];
sx q[1];
rz(-1.6436952) q[1];
sx q[1];
rz(-2.4935007) q[1];
rz(-pi) q[2];
rz(2.996534) q[3];
sx q[3];
rz(-0.71509711) q[3];
sx q[3];
rz(1.3655658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0076912045) q[2];
sx q[2];
rz(-2.2654221) q[2];
sx q[2];
rz(0.42789856) q[2];
rz(2.2513572) q[3];
sx q[3];
rz(-2.4211703) q[3];
sx q[3];
rz(0.015457411) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7240922) q[0];
sx q[0];
rz(-0.98560792) q[0];
sx q[0];
rz(-0.86522657) q[0];
rz(2.4100307) q[1];
sx q[1];
rz(-2.1831436) q[1];
sx q[1];
rz(2.1581214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.959592) q[0];
sx q[0];
rz(-2.1361217) q[0];
sx q[0];
rz(1.8046787) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8268565) q[2];
sx q[2];
rz(-1.6552002) q[2];
sx q[2];
rz(-0.39324444) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76170761) q[1];
sx q[1];
rz(-1.7049676) q[1];
sx q[1];
rz(0.11685808) q[1];
x q[2];
rz(2.8494469) q[3];
sx q[3];
rz(-2.2315624) q[3];
sx q[3];
rz(-3.0432971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5359155) q[2];
sx q[2];
rz(-1.8226049) q[2];
sx q[2];
rz(0.21161045) q[2];
rz(-2.0896046) q[3];
sx q[3];
rz(-0.34422031) q[3];
sx q[3];
rz(-2.2229164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6720523) q[0];
sx q[0];
rz(-0.34228995) q[0];
sx q[0];
rz(2.4218609) q[0];
rz(-1.4620048) q[1];
sx q[1];
rz(-2.3639677) q[1];
sx q[1];
rz(3.0696312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0746243) q[0];
sx q[0];
rz(-0.82199854) q[0];
sx q[0];
rz(-2.7496265) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9182704) q[2];
sx q[2];
rz(-0.90969786) q[2];
sx q[2];
rz(-1.2375792) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0256361) q[1];
sx q[1];
rz(-1.9731338) q[1];
sx q[1];
rz(1.3570804) q[1];
rz(-pi) q[2];
rz(-2.3945599) q[3];
sx q[3];
rz(-0.45901925) q[3];
sx q[3];
rz(-1.6667509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0422684) q[2];
sx q[2];
rz(-0.86809731) q[2];
sx q[2];
rz(-0.68359366) q[2];
rz(-1.4156703) q[3];
sx q[3];
rz(-2.1779124) q[3];
sx q[3];
rz(-1.8583309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0409446) q[0];
sx q[0];
rz(-1.9293979) q[0];
sx q[0];
rz(1.9579493) q[0];
rz(0.64507858) q[1];
sx q[1];
rz(-1.8408142) q[1];
sx q[1];
rz(-2.849978) q[1];
rz(1.245247) q[2];
sx q[2];
rz(-1.799188) q[2];
sx q[2];
rz(-1.3273018) q[2];
rz(0.24786077) q[3];
sx q[3];
rz(-1.9209721) q[3];
sx q[3];
rz(2.1192931) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
