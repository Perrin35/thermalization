OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(2.7913845) q[0];
sx q[0];
rz(12.933001) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(1.4555414) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35858425) q[0];
sx q[0];
rz(-1.9134247) q[0];
sx q[0];
rz(2.0185828) q[0];
rz(-pi) q[1];
rz(-0.86906616) q[2];
sx q[2];
rz(-1.367374) q[2];
sx q[2];
rz(0.8806526) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9526706) q[1];
sx q[1];
rz(-1.023472) q[1];
sx q[1];
rz(0.31618936) q[1];
x q[2];
rz(3.0936712) q[3];
sx q[3];
rz(-1.5600292) q[3];
sx q[3];
rz(2.3445232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.1047487) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(-2.74995) q[2];
rz(-0.13970217) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7249107) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(1.8649944) q[0];
rz(0.87031594) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(-1.2044027) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8613974) q[0];
sx q[0];
rz(-1.9635927) q[0];
sx q[0];
rz(1.1108206) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9559829) q[2];
sx q[2];
rz(-2.2679272) q[2];
sx q[2];
rz(2.6999161) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.345563) q[1];
sx q[1];
rz(-1.9351164) q[1];
sx q[1];
rz(1.4644535) q[1];
x q[2];
rz(0.28537206) q[3];
sx q[3];
rz(-2.223569) q[3];
sx q[3];
rz(-0.6245581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.63699547) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(-1.2191999) q[2];
rz(-2.7870264) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(-1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5511659) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(2.5235126) q[0];
rz(-0.016050054) q[1];
sx q[1];
rz(-2.2647104) q[1];
sx q[1];
rz(-1.9504257) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5844849) q[0];
sx q[0];
rz(-0.35042052) q[0];
sx q[0];
rz(-2.0273897) q[0];
rz(-pi) q[1];
rz(-2.8580655) q[2];
sx q[2];
rz(-1.8054188) q[2];
sx q[2];
rz(-1.4153751) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1515854) q[1];
sx q[1];
rz(-1.3602435) q[1];
sx q[1];
rz(-1.7819808) q[1];
rz(-pi) q[2];
rz(-1.972584) q[3];
sx q[3];
rz(-2.0766633) q[3];
sx q[3];
rz(-1.5617621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.96737343) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(2.1276786) q[2];
rz(-1.3570471) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(-2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8106666) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(2.9372835) q[0];
rz(1.3775685) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(-0.92299443) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2134593) q[0];
sx q[0];
rz(-1.8753578) q[0];
sx q[0];
rz(0.34780985) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5590018) q[2];
sx q[2];
rz(-1.6147385) q[2];
sx q[2];
rz(-0.61461385) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35034414) q[1];
sx q[1];
rz(-1.6748168) q[1];
sx q[1];
rz(1.9700325) q[1];
rz(-pi) q[2];
rz(0.8664341) q[3];
sx q[3];
rz(-1.5878864) q[3];
sx q[3];
rz(-0.41364241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6874281) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(0.50951177) q[2];
rz(-2.4984958) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(-2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.61280695) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(-0.91941419) q[0];
rz(-0.98584229) q[1];
sx q[1];
rz(-1.3580094) q[1];
sx q[1];
rz(-1.3607508) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082598535) q[0];
sx q[0];
rz(-1.7599802) q[0];
sx q[0];
rz(0.22342213) q[0];
x q[1];
rz(0.947457) q[2];
sx q[2];
rz(-0.91295487) q[2];
sx q[2];
rz(2.5156977) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9779224) q[1];
sx q[1];
rz(-1.1174035) q[1];
sx q[1];
rz(-1.8015253) q[1];
x q[2];
rz(2.1702483) q[3];
sx q[3];
rz(-1.6347307) q[3];
sx q[3];
rz(1.5980699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.78648606) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(-0.11486593) q[2];
rz(0.45587513) q[3];
sx q[3];
rz(-1.8084278) q[3];
sx q[3];
rz(-1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.24546656) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(1.8967569) q[0];
rz(-0.70760977) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(-0.27522603) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21101418) q[0];
sx q[0];
rz(-1.0761346) q[0];
sx q[0];
rz(-0.37822322) q[0];
x q[1];
rz(-2.2116824) q[2];
sx q[2];
rz(-1.9800131) q[2];
sx q[2];
rz(0.90604679) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.59906193) q[1];
sx q[1];
rz(-1.9179357) q[1];
sx q[1];
rz(0.51862006) q[1];
rz(-1.6861378) q[3];
sx q[3];
rz(-2.5177258) q[3];
sx q[3];
rz(1.3184402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.747921) q[2];
sx q[2];
rz(-1.5966406) q[2];
sx q[2];
rz(-0.45864027) q[2];
rz(0.7115055) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(0.59035629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8608619) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(1.3635427) q[0];
rz(-1.7215615) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(0.71969676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8450711) q[0];
sx q[0];
rz(-2.4356027) q[0];
sx q[0];
rz(0.4476053) q[0];
rz(0.72549707) q[2];
sx q[2];
rz(-0.91425397) q[2];
sx q[2];
rz(0.072349116) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.67932018) q[1];
sx q[1];
rz(-1.5264891) q[1];
sx q[1];
rz(-2.5965152) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68088716) q[3];
sx q[3];
rz(-1.7828935) q[3];
sx q[3];
rz(-1.6915481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.482243) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(0.68816319) q[2];
rz(-0.041953772) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(-0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(2.9550609) q[0];
rz(2.5371011) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(1.6465181) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.625531) q[0];
sx q[0];
rz(-2.7357258) q[0];
sx q[0];
rz(0.56674515) q[0];
rz(-pi) q[1];
rz(-0.69627701) q[2];
sx q[2];
rz(-1.7211203) q[2];
sx q[2];
rz(0.50762774) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.013156803) q[1];
sx q[1];
rz(-1.6746579) q[1];
sx q[1];
rz(0.26058148) q[1];
x q[2];
rz(-2.2349615) q[3];
sx q[3];
rz(-2.2615221) q[3];
sx q[3];
rz(0.074113473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8994393) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(0.014766679) q[2];
rz(1.8642558) q[3];
sx q[3];
rz(-1.8321962) q[3];
sx q[3];
rz(1.7285255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16167851) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(2.263608) q[0];
rz(0.19628482) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(-1.3605798) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36688766) q[0];
sx q[0];
rz(-0.57292143) q[0];
sx q[0];
rz(1.6402871) q[0];
rz(-pi) q[1];
rz(-2.1979245) q[2];
sx q[2];
rz(-1.5329648) q[2];
sx q[2];
rz(-0.19247069) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66879241) q[1];
sx q[1];
rz(-1.9023412) q[1];
sx q[1];
rz(-1.214295) q[1];
rz(1.2336041) q[3];
sx q[3];
rz(-2.7003151) q[3];
sx q[3];
rz(-2.6233167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62375162) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(-2.2488135) q[2];
rz(3.1023074) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(0.75004309) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84353012) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(-3.0950586) q[0];
rz(-2.2325366) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(-2.4216901) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4601446) q[0];
sx q[0];
rz(-1.3272734) q[0];
sx q[0];
rz(3.0597276) q[0];
rz(-pi) q[1];
rz(-1.5653531) q[2];
sx q[2];
rz(-1.2619702) q[2];
sx q[2];
rz(2.960161) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7457909) q[1];
sx q[1];
rz(-1.5862641) q[1];
sx q[1];
rz(1.323911) q[1];
x q[2];
rz(0.87749691) q[3];
sx q[3];
rz(-2.174456) q[3];
sx q[3];
rz(-0.391215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4204734) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(3.1151248) q[2];
rz(-1.3039533) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(1.8855689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532886) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(2.9251255) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(-0.70710612) q[2];
sx q[2];
rz(-1.8137992) q[2];
sx q[2];
rz(0.9546311) q[2];
rz(-2.1574216) q[3];
sx q[3];
rz(-1.5979206) q[3];
sx q[3];
rz(0.4978705) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];