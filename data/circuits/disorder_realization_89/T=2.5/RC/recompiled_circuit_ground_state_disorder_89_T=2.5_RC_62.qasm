OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9417579) q[0];
sx q[0];
rz(-2.4500442) q[0];
sx q[0];
rz(-2.3681695) q[0];
rz(-0.091104658) q[1];
sx q[1];
rz(-2.3478822) q[1];
sx q[1];
rz(-1.8082126) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0036732) q[0];
sx q[0];
rz(-0.9124476) q[0];
sx q[0];
rz(-1.5594894) q[0];
rz(1.0716419) q[2];
sx q[2];
rz(-1.7737651) q[2];
sx q[2];
rz(2.2405193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7385834) q[1];
sx q[1];
rz(-2.3883551) q[1];
sx q[1];
rz(-2.1639216) q[1];
x q[2];
rz(-1.638014) q[3];
sx q[3];
rz(-1.5537694) q[3];
sx q[3];
rz(0.57651455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8495463) q[2];
sx q[2];
rz(-1.3214194) q[2];
sx q[2];
rz(-0.41479659) q[2];
rz(1.241812) q[3];
sx q[3];
rz(-1.1153699) q[3];
sx q[3];
rz(2.950086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6280129) q[0];
sx q[0];
rz(-3.0256671) q[0];
sx q[0];
rz(0.43981788) q[0];
rz(3.0385333) q[1];
sx q[1];
rz(-2.8524046) q[1];
sx q[1];
rz(1.9724847) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77020184) q[0];
sx q[0];
rz(-1.7097211) q[0];
sx q[0];
rz(-2.2004805) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98054072) q[2];
sx q[2];
rz(-1.613918) q[2];
sx q[2];
rz(1.0999964) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7767765) q[1];
sx q[1];
rz(-2.4668982) q[1];
sx q[1];
rz(-1.2801835) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7984932) q[3];
sx q[3];
rz(-2.4191354) q[3];
sx q[3];
rz(1.3355153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.64572) q[2];
sx q[2];
rz(-1.3912667) q[2];
sx q[2];
rz(-2.8953841) q[2];
rz(1.5496893) q[3];
sx q[3];
rz(-0.61384765) q[3];
sx q[3];
rz(1.7483819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98591268) q[0];
sx q[0];
rz(-2.0014626) q[0];
sx q[0];
rz(-0.21173665) q[0];
rz(-0.011292975) q[1];
sx q[1];
rz(-0.79835049) q[1];
sx q[1];
rz(-1.4439772) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2560476) q[0];
sx q[0];
rz(-3.1064337) q[0];
sx q[0];
rz(3.1395182) q[0];
x q[1];
rz(-1.774774) q[2];
sx q[2];
rz(-0.24572554) q[2];
sx q[2];
rz(1.3947006) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1423751) q[1];
sx q[1];
rz(-0.96227598) q[1];
sx q[1];
rz(2.0318248) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6908247) q[3];
sx q[3];
rz(-1.6524466) q[3];
sx q[3];
rz(-1.9729561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4766562) q[2];
sx q[2];
rz(-1.0756476) q[2];
sx q[2];
rz(-0.54534379) q[2];
rz(2.2319345) q[3];
sx q[3];
rz(-0.46234149) q[3];
sx q[3];
rz(1.9285412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.42105168) q[0];
sx q[0];
rz(-0.67890972) q[0];
sx q[0];
rz(-0.82756591) q[0];
rz(-1.958581) q[1];
sx q[1];
rz(-1.8667826) q[1];
sx q[1];
rz(-1.8756728) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2823845) q[0];
sx q[0];
rz(-3.0050238) q[0];
sx q[0];
rz(-0.84117667) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6235731) q[2];
sx q[2];
rz(-0.20773331) q[2];
sx q[2];
rz(0.29333255) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3193533) q[1];
sx q[1];
rz(-0.6262928) q[1];
sx q[1];
rz(-1.6250454) q[1];
rz(-0.16014853) q[3];
sx q[3];
rz(-1.3831098) q[3];
sx q[3];
rz(1.2040862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2563235) q[2];
sx q[2];
rz(-1.6689391) q[2];
sx q[2];
rz(3.0598158) q[2];
rz(0.30803382) q[3];
sx q[3];
rz(-2.6576198) q[3];
sx q[3];
rz(-1.6035732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1643739) q[0];
sx q[0];
rz(-2.866221) q[0];
sx q[0];
rz(-3.0710869) q[0];
rz(-2.8071857) q[1];
sx q[1];
rz(-0.53135482) q[1];
sx q[1];
rz(0.45323429) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95246661) q[0];
sx q[0];
rz(-1.0116018) q[0];
sx q[0];
rz(0.97049148) q[0];
rz(-2.9137827) q[2];
sx q[2];
rz(-1.1880298) q[2];
sx q[2];
rz(-1.5441976) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6649225) q[1];
sx q[1];
rz(-1.8174662) q[1];
sx q[1];
rz(-1.4566521) q[1];
rz(-pi) q[2];
rz(2.5284834) q[3];
sx q[3];
rz(-2.1804348) q[3];
sx q[3];
rz(0.10559374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8411023) q[2];
sx q[2];
rz(-1.6359685) q[2];
sx q[2];
rz(2.9312768) q[2];
rz(-2.7503843) q[3];
sx q[3];
rz(-2.936383) q[3];
sx q[3];
rz(-2.6989663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6415569) q[0];
sx q[0];
rz(-2.0099202) q[0];
sx q[0];
rz(0.17182194) q[0];
rz(-1.7956644) q[1];
sx q[1];
rz(-2.185952) q[1];
sx q[1];
rz(1.1489493) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35691914) q[0];
sx q[0];
rz(-1.5789512) q[0];
sx q[0];
rz(0.11155931) q[0];
rz(-pi) q[1];
rz(-0.8900155) q[2];
sx q[2];
rz(-0.35951722) q[2];
sx q[2];
rz(2.7744164) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.90075452) q[1];
sx q[1];
rz(-0.60610702) q[1];
sx q[1];
rz(-1.89487) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9257718) q[3];
sx q[3];
rz(-0.53484155) q[3];
sx q[3];
rz(0.25391423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.020784) q[2];
sx q[2];
rz(-0.9641996) q[2];
sx q[2];
rz(-0.42050335) q[2];
rz(-1.6546107) q[3];
sx q[3];
rz(-1.9872811) q[3];
sx q[3];
rz(1.1544863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7789975) q[0];
sx q[0];
rz(-1.2263466) q[0];
sx q[0];
rz(-0.37780365) q[0];
rz(1.7516349) q[1];
sx q[1];
rz(-1.708834) q[1];
sx q[1];
rz(-2.7094944) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4946232) q[0];
sx q[0];
rz(-1.7996117) q[0];
sx q[0];
rz(-1.7499313) q[0];
rz(-0.092707002) q[2];
sx q[2];
rz(-1.5376629) q[2];
sx q[2];
rz(0.2011782) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3687146) q[1];
sx q[1];
rz(-2.0324824) q[1];
sx q[1];
rz(-1.6108455) q[1];
rz(1.0749519) q[3];
sx q[3];
rz(-2.4044839) q[3];
sx q[3];
rz(0.51966156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3720588) q[2];
sx q[2];
rz(-0.62166119) q[2];
sx q[2];
rz(2.9080234) q[2];
rz(1.6893049) q[3];
sx q[3];
rz(-1.7100311) q[3];
sx q[3];
rz(-1.5805894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772783) q[0];
sx q[0];
rz(-2.1121341) q[0];
sx q[0];
rz(2.8218063) q[0];
rz(-2.2794967) q[1];
sx q[1];
rz(-1.2513221) q[1];
sx q[1];
rz(-1.7822942) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0523543) q[0];
sx q[0];
rz(-1.0066832) q[0];
sx q[0];
rz(-2.4648239) q[0];
rz(-pi) q[1];
rz(-1.3699652) q[2];
sx q[2];
rz(-1.4322831) q[2];
sx q[2];
rz(0.85182933) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8492125) q[1];
sx q[1];
rz(-1.74382) q[1];
sx q[1];
rz(0.72489691) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2044851) q[3];
sx q[3];
rz(-1.7965297) q[3];
sx q[3];
rz(-0.87701517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84997815) q[2];
sx q[2];
rz(-2.7225967) q[2];
sx q[2];
rz(2.6739547) q[2];
rz(-2.6575798) q[3];
sx q[3];
rz(-1.7734807) q[3];
sx q[3];
rz(1.8638994) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17396946) q[0];
sx q[0];
rz(-1.4894217) q[0];
sx q[0];
rz(-1.1349348) q[0];
rz(-0.060215503) q[1];
sx q[1];
rz(-0.73679149) q[1];
sx q[1];
rz(-1.5312451) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0031207) q[0];
sx q[0];
rz(-2.0299737) q[0];
sx q[0];
rz(-0.91158406) q[0];
rz(-1.0410454) q[2];
sx q[2];
rz(-0.95472958) q[2];
sx q[2];
rz(0.043827961) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.169226) q[1];
sx q[1];
rz(-1.833583) q[1];
sx q[1];
rz(1.6714833) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1020622) q[3];
sx q[3];
rz(-2.3590204) q[3];
sx q[3];
rz(1.4631997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.327329) q[2];
sx q[2];
rz(-1.9244104) q[2];
sx q[2];
rz(2.4129756) q[2];
rz(-2.9260351) q[3];
sx q[3];
rz(-1.7451127) q[3];
sx q[3];
rz(-2.047915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88290596) q[0];
sx q[0];
rz(-0.031143324) q[0];
sx q[0];
rz(2.8236142) q[0];
rz(1.3975573) q[1];
sx q[1];
rz(-1.3293068) q[1];
sx q[1];
rz(2.8584282) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508789) q[0];
sx q[0];
rz(-0.15199272) q[0];
sx q[0];
rz(-2.1426523) q[0];
rz(2.3874902) q[2];
sx q[2];
rz(-1.5148544) q[2];
sx q[2];
rz(-2.2365582) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.52022987) q[1];
sx q[1];
rz(-0.89530357) q[1];
sx q[1];
rz(-0.9654926) q[1];
rz(1.3628325) q[3];
sx q[3];
rz(-1.4803518) q[3];
sx q[3];
rz(-1.2048723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2240923) q[2];
sx q[2];
rz(-1.516569) q[2];
sx q[2];
rz(3.0523114) q[2];
rz(-1.8381522) q[3];
sx q[3];
rz(-0.83804122) q[3];
sx q[3];
rz(0.97344056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78171171) q[0];
sx q[0];
rz(-0.64890535) q[0];
sx q[0];
rz(-2.1606408) q[0];
rz(-2.5557062) q[1];
sx q[1];
rz(-1.2955019) q[1];
sx q[1];
rz(-1.6265709) q[1];
rz(-0.70102767) q[2];
sx q[2];
rz(-0.27576896) q[2];
sx q[2];
rz(-2.452843) q[2];
rz(-1.1896776) q[3];
sx q[3];
rz(-1.8101235) q[3];
sx q[3];
rz(2.9265612) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
