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
rz(-2.9450077) q[0];
sx q[0];
rz(-0.44214806) q[0];
sx q[0];
rz(-2.2801939) q[0];
rz(1.5098894) q[1];
sx q[1];
rz(-1.8169401) q[1];
sx q[1];
rz(3.1060001) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067757879) q[0];
sx q[0];
rz(-1.375421) q[0];
sx q[0];
rz(-0.34123502) q[0];
rz(-pi) q[1];
rz(-2.8350949) q[2];
sx q[2];
rz(-1.674429) q[2];
sx q[2];
rz(0.95970861) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1110927) q[1];
sx q[1];
rz(-0.72734088) q[1];
sx q[1];
rz(0.05383454) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89081141) q[3];
sx q[3];
rz(-1.5847795) q[3];
sx q[3];
rz(-2.3774862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.55277905) q[2];
sx q[2];
rz(-1.3014883) q[2];
sx q[2];
rz(2.0471052) q[2];
rz(-1.6257446) q[3];
sx q[3];
rz(-2.2690014) q[3];
sx q[3];
rz(-2.998013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0482408) q[0];
sx q[0];
rz(-2.1513262) q[0];
sx q[0];
rz(0.37407237) q[0];
rz(-2.2834942) q[1];
sx q[1];
rz(-2.2007807) q[1];
sx q[1];
rz(2.0657952) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6158524) q[0];
sx q[0];
rz(-1.6739556) q[0];
sx q[0];
rz(0.10248689) q[0];
rz(-pi) q[1];
rz(-2.7611742) q[2];
sx q[2];
rz(-2.5606692) q[2];
sx q[2];
rz(0.38512938) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.43861345) q[1];
sx q[1];
rz(-1.9697048) q[1];
sx q[1];
rz(1.6970509) q[1];
x q[2];
rz(-0.91173633) q[3];
sx q[3];
rz(-2.1359332) q[3];
sx q[3];
rz(2.679058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6263803) q[2];
sx q[2];
rz(-0.43312803) q[2];
sx q[2];
rz(2.058378) q[2];
rz(1.8488688) q[3];
sx q[3];
rz(-1.1336528) q[3];
sx q[3];
rz(-0.92866549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.2556297) q[0];
sx q[0];
rz(-0.77115458) q[0];
sx q[0];
rz(-2.6166925) q[0];
rz(1.6820172) q[1];
sx q[1];
rz(-2.4039905) q[1];
sx q[1];
rz(-0.14579138) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4065259) q[0];
sx q[0];
rz(-0.13938381) q[0];
sx q[0];
rz(2.4295619) q[0];
rz(-pi) q[1];
rz(-2.6078659) q[2];
sx q[2];
rz(-1.0366129) q[2];
sx q[2];
rz(1.6426517) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.98602146) q[1];
sx q[1];
rz(-1.1386445) q[1];
sx q[1];
rz(1.4714824) q[1];
rz(-0.84050982) q[3];
sx q[3];
rz(-2.0183992) q[3];
sx q[3];
rz(2.2882035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8968481) q[2];
sx q[2];
rz(-1.9795828) q[2];
sx q[2];
rz(0.16540089) q[2];
rz(-0.8616972) q[3];
sx q[3];
rz(-0.69178897) q[3];
sx q[3];
rz(-2.944788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2794613) q[0];
sx q[0];
rz(-2.617351) q[0];
sx q[0];
rz(-2.4133546) q[0];
rz(-0.32314745) q[1];
sx q[1];
rz(-2.1073982) q[1];
sx q[1];
rz(-2.1023777) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0511) q[0];
sx q[0];
rz(-1.9673825) q[0];
sx q[0];
rz(-1.1086247) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8604554) q[2];
sx q[2];
rz(-0.95401154) q[2];
sx q[2];
rz(0.58831085) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7936662) q[1];
sx q[1];
rz(-1.6632342) q[1];
sx q[1];
rz(-0.4853918) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81890727) q[3];
sx q[3];
rz(-1.0207796) q[3];
sx q[3];
rz(-0.14980042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49179658) q[2];
sx q[2];
rz(-2.0911262) q[2];
sx q[2];
rz(2.6226079) q[2];
rz(2.0670048) q[3];
sx q[3];
rz(-1.855775) q[3];
sx q[3];
rz(-2.6463267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.0786667) q[0];
sx q[0];
rz(-0.92200297) q[0];
sx q[0];
rz(-0.17157383) q[0];
rz(1.0942787) q[1];
sx q[1];
rz(-1.8405874) q[1];
sx q[1];
rz(2.9069854) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07148347) q[0];
sx q[0];
rz(-0.30710328) q[0];
sx q[0];
rz(-1.135541) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0472946) q[2];
sx q[2];
rz(-0.50912947) q[2];
sx q[2];
rz(-1.5056339) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.68069862) q[1];
sx q[1];
rz(-2.1354224) q[1];
sx q[1];
rz(2.2959397) q[1];
rz(-pi) q[2];
rz(1.9453641) q[3];
sx q[3];
rz(-0.8959594) q[3];
sx q[3];
rz(2.5684211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2285063) q[2];
sx q[2];
rz(-1.191491) q[2];
sx q[2];
rz(1.7388434) q[2];
rz(2.5168822) q[3];
sx q[3];
rz(-2.1731845) q[3];
sx q[3];
rz(1.3431842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45873555) q[0];
sx q[0];
rz(-0.0077489297) q[0];
sx q[0];
rz(-0.32753456) q[0];
rz(-3.1134743) q[1];
sx q[1];
rz(-1.1437806) q[1];
sx q[1];
rz(2.1498674) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1234489) q[0];
sx q[0];
rz(-2.2010692) q[0];
sx q[0];
rz(-2.2176377) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6985833) q[2];
sx q[2];
rz(-0.24925079) q[2];
sx q[2];
rz(-0.64068782) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8918482) q[1];
sx q[1];
rz(-1.5288903) q[1];
sx q[1];
rz(-0.81467198) q[1];
rz(-pi) q[2];
rz(-1.5434274) q[3];
sx q[3];
rz(-2.0413766) q[3];
sx q[3];
rz(-1.1740299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.907054) q[2];
sx q[2];
rz(-1.3776642) q[2];
sx q[2];
rz(-1.482359) q[2];
rz(-2.0756857) q[3];
sx q[3];
rz(-1.1803455) q[3];
sx q[3];
rz(-1.0970241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7557573) q[0];
sx q[0];
rz(-0.11180728) q[0];
sx q[0];
rz(0.29543153) q[0];
rz(-2.9040728) q[1];
sx q[1];
rz(-2.2078881) q[1];
sx q[1];
rz(0.74737731) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.462446) q[0];
sx q[0];
rz(-2.0751795) q[0];
sx q[0];
rz(0.64179582) q[0];
x q[1];
rz(0.68589034) q[2];
sx q[2];
rz(-0.7524366) q[2];
sx q[2];
rz(2.8519832) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.566639) q[1];
sx q[1];
rz(-0.63640928) q[1];
sx q[1];
rz(-2.2242935) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1959869) q[3];
sx q[3];
rz(-2.0829933) q[3];
sx q[3];
rz(-2.1846445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4733009) q[2];
sx q[2];
rz(-1.1857727) q[2];
sx q[2];
rz(0.2549003) q[2];
rz(0.21236803) q[3];
sx q[3];
rz(-2.2771211) q[3];
sx q[3];
rz(0.00096850639) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2905228) q[0];
sx q[0];
rz(-1.2238598) q[0];
sx q[0];
rz(0.12390027) q[0];
rz(2.2663785) q[1];
sx q[1];
rz(-0.83299914) q[1];
sx q[1];
rz(0.55878729) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8494069) q[0];
sx q[0];
rz(-1.8208139) q[0];
sx q[0];
rz(1.0305189) q[0];
rz(-pi) q[1];
rz(-0.28382878) q[2];
sx q[2];
rz(-2.0528194) q[2];
sx q[2];
rz(1.986077) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.24148252) q[1];
sx q[1];
rz(-2.0962976) q[1];
sx q[1];
rz(1.8837758) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93620091) q[3];
sx q[3];
rz(-1.9125328) q[3];
sx q[3];
rz(0.5638916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.045804068) q[2];
sx q[2];
rz(-0.49171058) q[2];
sx q[2];
rz(-2.4208505) q[2];
rz(-2.7785684) q[3];
sx q[3];
rz(-1.9178773) q[3];
sx q[3];
rz(-1.5117517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95661288) q[0];
sx q[0];
rz(-0.19421254) q[0];
sx q[0];
rz(0.089056253) q[0];
rz(2.7748499) q[1];
sx q[1];
rz(-2.2373503) q[1];
sx q[1];
rz(-1.6820224) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90208611) q[0];
sx q[0];
rz(-1.7988822) q[0];
sx q[0];
rz(2.2091764) q[0];
rz(-pi) q[1];
rz(0.681293) q[2];
sx q[2];
rz(-0.80901399) q[2];
sx q[2];
rz(1.9495277) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4521317) q[1];
sx q[1];
rz(-2.2575827) q[1];
sx q[1];
rz(0.48312123) q[1];
rz(-pi) q[2];
rz(2.5404471) q[3];
sx q[3];
rz(-0.69128762) q[3];
sx q[3];
rz(-0.14562452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9350932) q[2];
sx q[2];
rz(-1.7743899) q[2];
sx q[2];
rz(1.2500786) q[2];
rz(1.9153204) q[3];
sx q[3];
rz(-0.63392249) q[3];
sx q[3];
rz(-2.5026076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5117689) q[0];
sx q[0];
rz(-2.0084232) q[0];
sx q[0];
rz(2.0015707) q[0];
rz(-2.5584768) q[1];
sx q[1];
rz(-1.4864328) q[1];
sx q[1];
rz(-2.4737632) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9662922) q[0];
sx q[0];
rz(-0.81476962) q[0];
sx q[0];
rz(-2.6867784) q[0];
rz(0.23867757) q[2];
sx q[2];
rz(-0.88310034) q[2];
sx q[2];
rz(2.4872045) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3558423) q[1];
sx q[1];
rz(-1.7652581) q[1];
sx q[1];
rz(2.0114312) q[1];
x q[2];
rz(1.1125426) q[3];
sx q[3];
rz(-0.83247165) q[3];
sx q[3];
rz(2.3979681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2976611) q[2];
sx q[2];
rz(-2.6466978) q[2];
sx q[2];
rz(0.75221357) q[2];
rz(-0.080605896) q[3];
sx q[3];
rz(-1.7919431) q[3];
sx q[3];
rz(-2.5940671) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3362296) q[0];
sx q[0];
rz(-1.6148051) q[0];
sx q[0];
rz(-0.057407277) q[0];
rz(0.84857955) q[1];
sx q[1];
rz(-0.72863693) q[1];
sx q[1];
rz(-1.4562664) q[1];
rz(2.3586629) q[2];
sx q[2];
rz(-2.30884) q[2];
sx q[2];
rz(-2.5306551) q[2];
rz(2.4918588) q[3];
sx q[3];
rz(-0.55745535) q[3];
sx q[3];
rz(-1.0563323) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
