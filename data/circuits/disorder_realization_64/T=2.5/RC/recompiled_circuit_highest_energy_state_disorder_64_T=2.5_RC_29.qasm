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
rz(-0.023802726) q[0];
sx q[0];
rz(4.2476141) q[0];
sx q[0];
rz(10.201693) q[0];
rz(1.39224) q[1];
sx q[1];
rz(-1.3148146) q[1];
sx q[1];
rz(2.1652752) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2433246) q[0];
sx q[0];
rz(-1.7043566) q[0];
sx q[0];
rz(-1.1569174) q[0];
rz(1.637865) q[2];
sx q[2];
rz(-2.133495) q[2];
sx q[2];
rz(0.85516155) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7853422) q[1];
sx q[1];
rz(-1.1467198) q[1];
sx q[1];
rz(-0.29845684) q[1];
rz(-pi) q[2];
rz(-1.5209274) q[3];
sx q[3];
rz(-0.5736151) q[3];
sx q[3];
rz(0.22465868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6355847) q[2];
sx q[2];
rz(-3.0878461) q[2];
sx q[2];
rz(-2.7347943) q[2];
rz(-0.16945101) q[3];
sx q[3];
rz(-0.5295161) q[3];
sx q[3];
rz(2.0690401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88103831) q[0];
sx q[0];
rz(-0.23242234) q[0];
sx q[0];
rz(-0.0037923092) q[0];
rz(-3.0637528) q[1];
sx q[1];
rz(-0.65996116) q[1];
sx q[1];
rz(2.8357764) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77340375) q[0];
sx q[0];
rz(-2.0819252) q[0];
sx q[0];
rz(-0.21171932) q[0];
x q[1];
rz(-2.3897116) q[2];
sx q[2];
rz(-1.6691748) q[2];
sx q[2];
rz(1.4954612) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.393049) q[1];
sx q[1];
rz(-2.3094588) q[1];
sx q[1];
rz(1.0219177) q[1];
rz(-pi) q[2];
rz(-1.634785) q[3];
sx q[3];
rz(-1.5100237) q[3];
sx q[3];
rz(-1.2343386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1514312) q[2];
sx q[2];
rz(-1.3913245) q[2];
sx q[2];
rz(-0.64275098) q[2];
rz(-0.07240545) q[3];
sx q[3];
rz(-1.0591155) q[3];
sx q[3];
rz(1.7806627) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0533503) q[0];
sx q[0];
rz(-0.14277661) q[0];
sx q[0];
rz(2.9852168) q[0];
rz(-0.0414255) q[1];
sx q[1];
rz(-2.5138469) q[1];
sx q[1];
rz(1.5511537) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32937059) q[0];
sx q[0];
rz(-0.58766627) q[0];
sx q[0];
rz(2.7633322) q[0];
x q[1];
rz(3.0535544) q[2];
sx q[2];
rz(-2.4836342) q[2];
sx q[2];
rz(-0.41706271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98387466) q[1];
sx q[1];
rz(-1.0441171) q[1];
sx q[1];
rz(-1.4406564) q[1];
rz(-0.82156397) q[3];
sx q[3];
rz(-0.85431495) q[3];
sx q[3];
rz(-2.2301073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7330043) q[2];
sx q[2];
rz(-2.6027347) q[2];
sx q[2];
rz(-3.1206701) q[2];
rz(0.18713348) q[3];
sx q[3];
rz(-0.20550607) q[3];
sx q[3];
rz(3.0278897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1784096) q[0];
sx q[0];
rz(-0.33247501) q[0];
sx q[0];
rz(2.7030429) q[0];
rz(-1.6167538) q[1];
sx q[1];
rz(-0.33477819) q[1];
sx q[1];
rz(2.8964892) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97814752) q[0];
sx q[0];
rz(-2.4171099) q[0];
sx q[0];
rz(-1.6969157) q[0];
rz(1.4133873) q[2];
sx q[2];
rz(-1.9345967) q[2];
sx q[2];
rz(-0.25870332) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0742798) q[1];
sx q[1];
rz(-1.7099705) q[1];
sx q[1];
rz(1.0275296) q[1];
rz(-pi) q[2];
rz(2.80527) q[3];
sx q[3];
rz(-2.3592279) q[3];
sx q[3];
rz(-0.7366283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.601292) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(-0.33176804) q[2];
rz(2.6541384) q[3];
sx q[3];
rz(-2.09477) q[3];
sx q[3];
rz(-2.148518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8496534) q[0];
sx q[0];
rz(-1.4537469) q[0];
sx q[0];
rz(0.77350235) q[0];
rz(-1.1812814) q[1];
sx q[1];
rz(-0.14207323) q[1];
sx q[1];
rz(1.389651) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7252325) q[0];
sx q[0];
rz(-1.9342039) q[0];
sx q[0];
rz(0.4298707) q[0];
rz(-pi) q[1];
rz(0.90221407) q[2];
sx q[2];
rz(-2.3093975) q[2];
sx q[2];
rz(-2.8813643) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90926814) q[1];
sx q[1];
rz(-2.0438015) q[1];
sx q[1];
rz(-1.6632516) q[1];
rz(-0.64719154) q[3];
sx q[3];
rz(-0.29509896) q[3];
sx q[3];
rz(1.5403252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9224077) q[2];
sx q[2];
rz(-1.9318523) q[2];
sx q[2];
rz(-0.47214559) q[2];
rz(1.2989429) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(2.3310272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2209114) q[0];
sx q[0];
rz(-2.8784316) q[0];
sx q[0];
rz(2.8780908) q[0];
rz(-1.1031411) q[1];
sx q[1];
rz(-1.8270854) q[1];
sx q[1];
rz(0.37364328) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140681) q[0];
sx q[0];
rz(-1.6434945) q[0];
sx q[0];
rz(-3.0811653) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0194026) q[2];
sx q[2];
rz(-0.8475248) q[2];
sx q[2];
rz(-3.0360589) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8220351) q[1];
sx q[1];
rz(-0.6614092) q[1];
sx q[1];
rz(-1.0378077) q[1];
rz(-pi) q[2];
rz(-0.5881891) q[3];
sx q[3];
rz(-2.0547227) q[3];
sx q[3];
rz(1.5744792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11334795) q[2];
sx q[2];
rz(-2.9714606) q[2];
sx q[2];
rz(2.587758) q[2];
rz(-1.7438186) q[3];
sx q[3];
rz(-2.5388986) q[3];
sx q[3];
rz(0.3127313) q[3];
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
rz(-2.5572307) q[0];
sx q[0];
rz(-2.1350242) q[0];
sx q[0];
rz(1.9118017) q[0];
rz(0.23904414) q[1];
sx q[1];
rz(-1.5115279) q[1];
sx q[1];
rz(-2.8412433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2588033) q[0];
sx q[0];
rz(-0.20016328) q[0];
sx q[0];
rz(-0.61997719) q[0];
rz(-pi) q[1];
rz(-2.0652169) q[2];
sx q[2];
rz(-2.2110614) q[2];
sx q[2];
rz(2.6312021) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3546037) q[1];
sx q[1];
rz(-1.5704078) q[1];
sx q[1];
rz(1.5860737) q[1];
x q[2];
rz(2.2836779) q[3];
sx q[3];
rz(-1.5522209) q[3];
sx q[3];
rz(1.9481079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.131669) q[2];
sx q[2];
rz(-1.4574304) q[2];
sx q[2];
rz(-2.8966676) q[2];
rz(-2.6217672) q[3];
sx q[3];
rz(-0.86331415) q[3];
sx q[3];
rz(-2.4533217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7898665) q[0];
sx q[0];
rz(-0.34857294) q[0];
sx q[0];
rz(1.1630195) q[0];
rz(0.066896833) q[1];
sx q[1];
rz(-1.6480548) q[1];
sx q[1];
rz(2.1287207) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5236008) q[0];
sx q[0];
rz(-2.3396684) q[0];
sx q[0];
rz(0.84419925) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9628721) q[2];
sx q[2];
rz(-1.0105437) q[2];
sx q[2];
rz(0.58748875) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3607531) q[1];
sx q[1];
rz(-1.8855125) q[1];
sx q[1];
rz(2.7937522) q[1];
rz(0.68655218) q[3];
sx q[3];
rz(-1.2408537) q[3];
sx q[3];
rz(-0.8984962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0248727) q[2];
sx q[2];
rz(-0.97815424) q[2];
sx q[2];
rz(-0.32279521) q[2];
rz(-2.5358477) q[3];
sx q[3];
rz(-0.79234684) q[3];
sx q[3];
rz(-2.7927223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.087273) q[0];
sx q[0];
rz(-0.1467341) q[0];
sx q[0];
rz(3.1291381) q[0];
rz(-2.3948578) q[1];
sx q[1];
rz(-0.92339271) q[1];
sx q[1];
rz(-2.8616203) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0736948) q[0];
sx q[0];
rz(-1.6215186) q[0];
sx q[0];
rz(-1.6832666) q[0];
x q[1];
rz(2.0048281) q[2];
sx q[2];
rz(-1.2723337) q[2];
sx q[2];
rz(-0.091574319) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88290434) q[1];
sx q[1];
rz(-1.8778442) q[1];
sx q[1];
rz(-0.76489246) q[1];
x q[2];
rz(2.8948726) q[3];
sx q[3];
rz(-1.2182052) q[3];
sx q[3];
rz(2.9143434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3296457) q[2];
sx q[2];
rz(-0.75355607) q[2];
sx q[2];
rz(-0.21198708) q[2];
rz(0.82344615) q[3];
sx q[3];
rz(-1.6971089) q[3];
sx q[3];
rz(0.25920355) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14281808) q[0];
sx q[0];
rz(-3.0840254) q[0];
sx q[0];
rz(-0.69277358) q[0];
rz(-2.5686) q[1];
sx q[1];
rz(-1.7968105) q[1];
sx q[1];
rz(0.43100345) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8739024) q[0];
sx q[0];
rz(-2.5866716) q[0];
sx q[0];
rz(0.11124994) q[0];
rz(-2.0153322) q[2];
sx q[2];
rz(-1.017414) q[2];
sx q[2];
rz(2.0972507) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5332408) q[1];
sx q[1];
rz(-0.87366381) q[1];
sx q[1];
rz(1.240154) q[1];
x q[2];
rz(-1.0209084) q[3];
sx q[3];
rz(-1.8326933) q[3];
sx q[3];
rz(2.6481215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4219605) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(-2.5813622) q[2];
rz(-2.6719921) q[3];
sx q[3];
rz(-0.40265366) q[3];
sx q[3];
rz(-0.71389055) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8568759) q[0];
sx q[0];
rz(-1.4061883) q[0];
sx q[0];
rz(-1.0549369) q[0];
rz(0.84125413) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(1.2011436) q[2];
sx q[2];
rz(-1.2751725) q[2];
sx q[2];
rz(-1.0775492) q[2];
rz(1.5359405) q[3];
sx q[3];
rz(-2.4823454) q[3];
sx q[3];
rz(-2.1128826) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
