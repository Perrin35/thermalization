OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2744098) q[0];
sx q[0];
rz(-0.46719587) q[0];
sx q[0];
rz(2.245477) q[0];
rz(2.3236302) q[1];
sx q[1];
rz(-1.5939413) q[1];
sx q[1];
rz(0.98734468) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86292279) q[0];
sx q[0];
rz(-0.76919829) q[0];
sx q[0];
rz(-1.0095566) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57968037) q[2];
sx q[2];
rz(-2.2319921) q[2];
sx q[2];
rz(1.5722317) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.27459135) q[1];
sx q[1];
rz(-1.4030292) q[1];
sx q[1];
rz(0.28047362) q[1];
rz(-pi) q[2];
rz(2.7475376) q[3];
sx q[3];
rz(-2.7480304) q[3];
sx q[3];
rz(1.2979899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92153111) q[2];
sx q[2];
rz(-1.0342197) q[2];
sx q[2];
rz(0.49911487) q[2];
rz(-2.3257997) q[3];
sx q[3];
rz(-0.59320265) q[3];
sx q[3];
rz(-0.35713404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.65207425) q[0];
sx q[0];
rz(-0.3638142) q[0];
sx q[0];
rz(0.2336842) q[0];
rz(-0.09985996) q[1];
sx q[1];
rz(-1.4870817) q[1];
sx q[1];
rz(-1.5335836) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8692538) q[0];
sx q[0];
rz(-1.4492479) q[0];
sx q[0];
rz(-0.097313332) q[0];
x q[1];
rz(-2.7614715) q[2];
sx q[2];
rz(-2.673966) q[2];
sx q[2];
rz(-1.1884226) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.93615426) q[1];
sx q[1];
rz(-2.1742651) q[1];
sx q[1];
rz(-0.11118576) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0879806) q[3];
sx q[3];
rz(-1.2948961) q[3];
sx q[3];
rz(1.538572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.097215501) q[2];
sx q[2];
rz(-0.88572398) q[2];
sx q[2];
rz(-3.0866747) q[2];
rz(2.4954259) q[3];
sx q[3];
rz(-1.5271527) q[3];
sx q[3];
rz(3.0687148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1845448) q[0];
sx q[0];
rz(-0.2453198) q[0];
sx q[0];
rz(2.3967337) q[0];
rz(-1.2767731) q[1];
sx q[1];
rz(-1.0294015) q[1];
sx q[1];
rz(0.863711) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.540032) q[0];
sx q[0];
rz(-1.1313725) q[0];
sx q[0];
rz(-1.5063398) q[0];
rz(2.6361739) q[2];
sx q[2];
rz(-0.8552455) q[2];
sx q[2];
rz(-0.16929132) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59434592) q[1];
sx q[1];
rz(-2.1859096) q[1];
sx q[1];
rz(-2.1662461) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66074547) q[3];
sx q[3];
rz(-1.163394) q[3];
sx q[3];
rz(-1.7254207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.471571) q[2];
sx q[2];
rz(-1.6100641) q[2];
sx q[2];
rz(-0.395533) q[2];
rz(2.1192571) q[3];
sx q[3];
rz(-2.6774355) q[3];
sx q[3];
rz(-2.6497604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96524298) q[0];
sx q[0];
rz(-1.9597766) q[0];
sx q[0];
rz(-3.0856207) q[0];
rz(-0.61198992) q[1];
sx q[1];
rz(-1.5925708) q[1];
sx q[1];
rz(2.2956119) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6406785) q[0];
sx q[0];
rz(-1.0553103) q[0];
sx q[0];
rz(-2.8436996) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95161423) q[2];
sx q[2];
rz(-1.2834833) q[2];
sx q[2];
rz(3.0653022) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8236134) q[1];
sx q[1];
rz(-2.7732964) q[1];
sx q[1];
rz(2.5846534) q[1];
rz(-pi) q[2];
rz(-2.0457129) q[3];
sx q[3];
rz(-1.2327415) q[3];
sx q[3];
rz(0.89382225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0116288) q[2];
sx q[2];
rz(-1.6946946) q[2];
sx q[2];
rz(-0.72032991) q[2];
rz(2.7885041) q[3];
sx q[3];
rz(-1.947764) q[3];
sx q[3];
rz(-0.24875719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40097749) q[0];
sx q[0];
rz(-1.688513) q[0];
sx q[0];
rz(0.98689669) q[0];
rz(-1.997021) q[1];
sx q[1];
rz(-0.83895504) q[1];
sx q[1];
rz(-2.1867337) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3564044) q[0];
sx q[0];
rz(-2.4739728) q[0];
sx q[0];
rz(-1.6130527) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0667162) q[2];
sx q[2];
rz(-0.62037797) q[2];
sx q[2];
rz(-1.7091144) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2209619) q[1];
sx q[1];
rz(-2.536533) q[1];
sx q[1];
rz(-0.69843881) q[1];
rz(-2.0935861) q[3];
sx q[3];
rz(-0.29055933) q[3];
sx q[3];
rz(2.235242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5130875) q[2];
sx q[2];
rz(-2.2272765) q[2];
sx q[2];
rz(1.0465735) q[2];
rz(2.4322677) q[3];
sx q[3];
rz(-1.9966634) q[3];
sx q[3];
rz(0.63372248) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3777622) q[0];
sx q[0];
rz(-1.7338294) q[0];
sx q[0];
rz(0.27107006) q[0];
rz(3.0625878) q[1];
sx q[1];
rz(-2.7075691) q[1];
sx q[1];
rz(-1.431042) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5749211) q[0];
sx q[0];
rz(-1.780176) q[0];
sx q[0];
rz(1.7185332) q[0];
x q[1];
rz(-0.76381229) q[2];
sx q[2];
rz(-0.71893636) q[2];
sx q[2];
rz(2.1670408) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9985314) q[1];
sx q[1];
rz(-2.240671) q[1];
sx q[1];
rz(1.8865442) q[1];
rz(1.1523139) q[3];
sx q[3];
rz(-1.7177204) q[3];
sx q[3];
rz(0.72673405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.048910353) q[2];
sx q[2];
rz(-1.1834669) q[2];
sx q[2];
rz(-3.1177706) q[2];
rz(-1.228099) q[3];
sx q[3];
rz(-2.7619599) q[3];
sx q[3];
rz(-2.9983799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.989885) q[0];
sx q[0];
rz(-2.8360974) q[0];
sx q[0];
rz(-0.09224961) q[0];
rz(0.30091885) q[1];
sx q[1];
rz(-1.9124799) q[1];
sx q[1];
rz(-2.0613964) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0940822) q[0];
sx q[0];
rz(-1.9157529) q[0];
sx q[0];
rz(0.10070078) q[0];
rz(0.87026922) q[2];
sx q[2];
rz(-0.53127161) q[2];
sx q[2];
rz(-0.16933717) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.791886) q[1];
sx q[1];
rz(-0.14794803) q[1];
sx q[1];
rz(0.011758864) q[1];
rz(-pi) q[2];
rz(1.2545414) q[3];
sx q[3];
rz(-0.90083921) q[3];
sx q[3];
rz(3.0681821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13176189) q[2];
sx q[2];
rz(-0.35085446) q[2];
sx q[2];
rz(2.5133613) q[2];
rz(-0.96588165) q[3];
sx q[3];
rz(-1.5768257) q[3];
sx q[3];
rz(-0.32111827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1107165) q[0];
sx q[0];
rz(-0.70699152) q[0];
sx q[0];
rz(-1.4068756) q[0];
rz(3.0137317) q[1];
sx q[1];
rz(-1.7117932) q[1];
sx q[1];
rz(2.0157287) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0681852) q[0];
sx q[0];
rz(-0.74755423) q[0];
sx q[0];
rz(-0.27184591) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1190542) q[2];
sx q[2];
rz(-1.5875419) q[2];
sx q[2];
rz(2.708205) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8496823) q[1];
sx q[1];
rz(-2.9314802) q[1];
sx q[1];
rz(0.30847524) q[1];
rz(-pi) q[2];
rz(-2.5223658) q[3];
sx q[3];
rz(-1.9398749) q[3];
sx q[3];
rz(-2.823373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0109978) q[2];
sx q[2];
rz(-2.0086292) q[2];
sx q[2];
rz(0.093718378) q[2];
rz(-1.4334009) q[3];
sx q[3];
rz(-2.8580557) q[3];
sx q[3];
rz(-1.7482429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-0.691064) q[0];
sx q[0];
rz(-2.0557623) q[0];
sx q[0];
rz(-1.0040671) q[0];
rz(1.1035236) q[1];
sx q[1];
rz(-2.6595778) q[1];
sx q[1];
rz(1.4080661) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226949) q[0];
sx q[0];
rz(-2.0859045) q[0];
sx q[0];
rz(1.5638419) q[0];
rz(0.31230052) q[2];
sx q[2];
rz(-1.6767297) q[2];
sx q[2];
rz(0.16367975) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6390761) q[1];
sx q[1];
rz(-2.6682819) q[1];
sx q[1];
rz(1.972354) q[1];
x q[2];
rz(-3.0549307) q[3];
sx q[3];
rz(-2.6209957) q[3];
sx q[3];
rz(2.749032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7496926) q[2];
sx q[2];
rz(-1.5211279) q[2];
sx q[2];
rz(-2.5938972) q[2];
rz(0.31217602) q[3];
sx q[3];
rz(-1.7269937) q[3];
sx q[3];
rz(-1.7267905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3373229) q[0];
sx q[0];
rz(-1.4708568) q[0];
sx q[0];
rz(-2.3626589) q[0];
rz(2.7670822) q[1];
sx q[1];
rz(-1.51314) q[1];
sx q[1];
rz(2.3096854) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4171281) q[0];
sx q[0];
rz(-2.3848371) q[0];
sx q[0];
rz(2.7555076) q[0];
rz(0.7124165) q[2];
sx q[2];
rz(-2.2210026) q[2];
sx q[2];
rz(-0.67931108) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.9495957) q[1];
sx q[1];
rz(-1.0284541) q[1];
sx q[1];
rz(2.3344912) q[1];
x q[2];
rz(0.79672265) q[3];
sx q[3];
rz(-1.5925538) q[3];
sx q[3];
rz(0.12783229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1067918) q[2];
sx q[2];
rz(-2.782244) q[2];
sx q[2];
rz(0.0084776004) q[2];
rz(2.7360385) q[3];
sx q[3];
rz(-1.6885992) q[3];
sx q[3];
rz(1.4439553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1494898) q[0];
sx q[0];
rz(-1.533951) q[0];
sx q[0];
rz(0.78716192) q[0];
rz(1.0943195) q[1];
sx q[1];
rz(-2.7631187) q[1];
sx q[1];
rz(-0.43597058) q[1];
rz(1.3588253) q[2];
sx q[2];
rz(-1.4478085) q[2];
sx q[2];
rz(0.87926452) q[2];
rz(0.094219128) q[3];
sx q[3];
rz(-1.5658656) q[3];
sx q[3];
rz(1.6938536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
