OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(3.7319558) q[0];
sx q[0];
rz(9.0537602) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(1.376027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5698619) q[0];
sx q[0];
rz(-2.1751582) q[0];
sx q[0];
rz(0.5423003) q[0];
x q[1];
rz(2.5901428) q[2];
sx q[2];
rz(-0.79556888) q[2];
sx q[2];
rz(1.7099107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0378117) q[1];
sx q[1];
rz(-1.4191287) q[1];
sx q[1];
rz(-1.2396461) q[1];
x q[2];
rz(0.052470603) q[3];
sx q[3];
rz(-0.82004181) q[3];
sx q[3];
rz(2.6657871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.084289) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(-2.1526745) q[2];
rz(-0.75254285) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9343524) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(-1.9616615) q[0];
rz(-0.99769366) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(2.4172799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1441919) q[0];
sx q[0];
rz(-1.8840944) q[0];
sx q[0];
rz(0.86667592) q[0];
rz(-pi) q[1];
rz(2.9421147) q[2];
sx q[2];
rz(-1.6530767) q[2];
sx q[2];
rz(0.85180887) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.38072941) q[1];
sx q[1];
rz(-2.7298096) q[1];
sx q[1];
rz(-1.0357344) q[1];
x q[2];
rz(1.3734666) q[3];
sx q[3];
rz(-1.8036246) q[3];
sx q[3];
rz(0.80430921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2362242) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(0.36188564) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1419462) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(-1.746159) q[0];
rz(-2.6793001) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(1.9225072) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31419471) q[0];
sx q[0];
rz(-2.6959246) q[0];
sx q[0];
rz(-2.218194) q[0];
rz(0.8692603) q[2];
sx q[2];
rz(-1.3420891) q[2];
sx q[2];
rz(1.9321835) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3238941) q[1];
sx q[1];
rz(-0.59191275) q[1];
sx q[1];
rz(-1.1281668) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3027906) q[3];
sx q[3];
rz(-1.5558814) q[3];
sx q[3];
rz(0.55004317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42157713) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(1.6960309) q[2];
rz(2.5727663) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(-3.0310757) q[3];
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
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-0.51825994) q[0];
rz(0.7154243) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(0.82675654) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402156) q[0];
sx q[0];
rz(-2.4670521) q[0];
sx q[0];
rz(-2.4504689) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0024198) q[2];
sx q[2];
rz(-0.77251245) q[2];
sx q[2];
rz(-0.13538361) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7552232) q[1];
sx q[1];
rz(-1.9519023) q[1];
sx q[1];
rz(1.2632881) q[1];
rz(1.8825674) q[3];
sx q[3];
rz(-0.87083737) q[3];
sx q[3];
rz(1.242897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.15239079) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(1.3789122) q[2];
rz(3.0692696) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(-1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(2.0157053) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(-0.28516969) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0641159) q[0];
sx q[0];
rz(-0.47668326) q[0];
sx q[0];
rz(-1.6633196) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4816149) q[2];
sx q[2];
rz(-0.72313213) q[2];
sx q[2];
rz(-2.8358104) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3045029) q[1];
sx q[1];
rz(-1.2577406) q[1];
sx q[1];
rz(0.00043991107) q[1];
rz(-0.40163715) q[3];
sx q[3];
rz(-1.4816227) q[3];
sx q[3];
rz(3.0107486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8482762) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(2.6021393) q[2];
rz(-0.30682492) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834615) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(-0.15701292) q[0];
rz(-0.69333386) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(-1.3670115) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088928662) q[0];
sx q[0];
rz(-3.0658709) q[0];
sx q[0];
rz(1.1195539) q[0];
x q[1];
rz(-0.96732803) q[2];
sx q[2];
rz(-0.8562932) q[2];
sx q[2];
rz(0.29010233) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4945592) q[1];
sx q[1];
rz(-1.4899947) q[1];
sx q[1];
rz(1.2809491) q[1];
rz(-pi) q[2];
rz(-0.27512392) q[3];
sx q[3];
rz(-0.36558662) q[3];
sx q[3];
rz(0.58055731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29629016) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(0.40346754) q[2];
rz(-2.6599595) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(0.51923716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-0.8738628) q[0];
rz(2.6938687) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.9708995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0543594) q[0];
sx q[0];
rz(-1.5938252) q[0];
sx q[0];
rz(-1.5769632) q[0];
rz(-pi) q[1];
rz(0.41202338) q[2];
sx q[2];
rz(-3.014866) q[2];
sx q[2];
rz(2.3840981) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9782941) q[1];
sx q[1];
rz(-2.1566609) q[1];
sx q[1];
rz(-0.75978029) q[1];
rz(-pi) q[2];
rz(0.4865173) q[3];
sx q[3];
rz(-0.58972893) q[3];
sx q[3];
rz(-0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0447023) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(-0.61075413) q[2];
rz(-2.6664873) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(-2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(-0.29712594) q[0];
rz(1.3946474) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(0.64613211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0028249) q[0];
sx q[0];
rz(-1.7768304) q[0];
sx q[0];
rz(3.0384008) q[0];
x q[1];
rz(-2.6693194) q[2];
sx q[2];
rz(-1.4204645) q[2];
sx q[2];
rz(-0.78778247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7568126) q[1];
sx q[1];
rz(-1.9128748) q[1];
sx q[1];
rz(0.53201075) q[1];
x q[2];
rz(2.331278) q[3];
sx q[3];
rz(-1.1018361) q[3];
sx q[3];
rz(1.1427493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(2.6867552) q[2];
rz(-2.440195) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69944537) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(2.3440857) q[0];
rz(2.6240255) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(-0.10841766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3603044) q[0];
sx q[0];
rz(-2.4753248) q[0];
sx q[0];
rz(1.6643307) q[0];
rz(-pi) q[1];
rz(-1.8234532) q[2];
sx q[2];
rz(-0.18441072) q[2];
sx q[2];
rz(-0.84469634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.11549599) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(-2.2309169) q[1];
x q[2];
rz(-2.2750862) q[3];
sx q[3];
rz(-1.7139072) q[3];
sx q[3];
rz(-1.0673616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.3468578) q[2];
sx q[2];
rz(-0.33995315) q[2];
rz(-2.7231976) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.6091992) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(2.4627731) q[0];
rz(-0.36418307) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(3.0864339) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63294166) q[0];
sx q[0];
rz(-1.5061437) q[0];
sx q[0];
rz(0.37737329) q[0];
rz(-pi) q[1];
rz(1.7807547) q[2];
sx q[2];
rz(-1.0210438) q[2];
sx q[2];
rz(1.2325665) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.060505796) q[1];
sx q[1];
rz(-2.2844237) q[1];
sx q[1];
rz(0.75018261) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1147898) q[3];
sx q[3];
rz(-1.6037233) q[3];
sx q[3];
rz(-1.959521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98383343) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(0.49017635) q[2];
rz(-0.13752078) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(-0.2086808) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(-1.6677042) q[2];
sx q[2];
rz(-1.2599535) q[2];
sx q[2];
rz(-2.8355666) q[2];
rz(0.18110885) q[3];
sx q[3];
rz(-0.24387471) q[3];
sx q[3];
rz(-2.7004164) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
