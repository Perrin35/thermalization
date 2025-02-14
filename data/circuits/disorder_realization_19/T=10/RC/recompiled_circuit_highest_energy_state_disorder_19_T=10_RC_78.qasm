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
rz(0.28111464) q[0];
sx q[0];
rz(-1.3507564) q[0];
sx q[0];
rz(-0.95066345) q[0];
rz(-2.9188393) q[1];
sx q[1];
rz(-0.25382257) q[1];
sx q[1];
rz(-2.4737127) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.107936) q[0];
sx q[0];
rz(-2.78068) q[0];
sx q[0];
rz(-0.080789195) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4761042) q[2];
sx q[2];
rz(-0.166278) q[2];
sx q[2];
rz(-0.68782297) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4227071) q[1];
sx q[1];
rz(-1.4816007) q[1];
sx q[1];
rz(-1.6929106) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5667772) q[3];
sx q[3];
rz(-1.7815456) q[3];
sx q[3];
rz(-0.8523418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.604598) q[2];
sx q[2];
rz(-0.94702417) q[2];
sx q[2];
rz(-0.14744082) q[2];
rz(-1.3747181) q[3];
sx q[3];
rz(-1.7539897) q[3];
sx q[3];
rz(1.714777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073567063) q[0];
sx q[0];
rz(-1.0140714) q[0];
sx q[0];
rz(0.2036988) q[0];
rz(0.112946) q[1];
sx q[1];
rz(-0.68165439) q[1];
sx q[1];
rz(3.122094) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18960139) q[0];
sx q[0];
rz(-0.87965779) q[0];
sx q[0];
rz(-2.7777872) q[0];
x q[1];
rz(-1.1910559) q[2];
sx q[2];
rz(-1.4342208) q[2];
sx q[2];
rz(-1.4771763) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9128117) q[1];
sx q[1];
rz(-1.2626909) q[1];
sx q[1];
rz(-2.5768215) q[1];
rz(-pi) q[2];
rz(-0.56657378) q[3];
sx q[3];
rz(-2.1687897) q[3];
sx q[3];
rz(-1.8600885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.514275) q[2];
sx q[2];
rz(-1.5760199) q[2];
sx q[2];
rz(0.19228284) q[2];
rz(-2.9228041) q[3];
sx q[3];
rz(-1.8407121) q[3];
sx q[3];
rz(-1.3505664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.669765) q[0];
sx q[0];
rz(-2.8596523) q[0];
sx q[0];
rz(2.6523253) q[0];
rz(1.4504112) q[1];
sx q[1];
rz(-1.1510808) q[1];
sx q[1];
rz(-0.53389186) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0301186) q[0];
sx q[0];
rz(-1.5970802) q[0];
sx q[0];
rz(2.8199955) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13907532) q[2];
sx q[2];
rz(-1.5267045) q[2];
sx q[2];
rz(-1.9382375) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3855672) q[1];
sx q[1];
rz(-1.4293164) q[1];
sx q[1];
rz(-1.1668851) q[1];
rz(-pi) q[2];
rz(1.3891641) q[3];
sx q[3];
rz(-2.2881621) q[3];
sx q[3];
rz(-2.7434012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78407946) q[2];
sx q[2];
rz(-1.0369022) q[2];
sx q[2];
rz(0.64884031) q[2];
rz(-0.31323788) q[3];
sx q[3];
rz(-2.1514838) q[3];
sx q[3];
rz(-1.311897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637777) q[0];
sx q[0];
rz(-2.2721993) q[0];
sx q[0];
rz(0.75611269) q[0];
rz(2.6847878) q[1];
sx q[1];
rz(-2.017338) q[1];
sx q[1];
rz(0.66064984) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74213492) q[0];
sx q[0];
rz(-0.44299865) q[0];
sx q[0];
rz(0.42794064) q[0];
rz(-pi) q[1];
rz(-1.2006423) q[2];
sx q[2];
rz(-1.3602019) q[2];
sx q[2];
rz(2.7538607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7380437) q[1];
sx q[1];
rz(-1.3894086) q[1];
sx q[1];
rz(2.2970143) q[1];
rz(2.8023671) q[3];
sx q[3];
rz(-2.0008828) q[3];
sx q[3];
rz(0.11232377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.75055355) q[2];
sx q[2];
rz(-1.4918574) q[2];
sx q[2];
rz(2.1569815) q[2];
rz(1.5084958) q[3];
sx q[3];
rz(-1.0909785) q[3];
sx q[3];
rz(2.7896037) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9462747) q[0];
sx q[0];
rz(-1.0349422) q[0];
sx q[0];
rz(-2.4786733) q[0];
rz(-1.4074869) q[1];
sx q[1];
rz(-2.2588142) q[1];
sx q[1];
rz(2.1972806) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51296455) q[0];
sx q[0];
rz(-1.7288393) q[0];
sx q[0];
rz(2.4751644) q[0];
rz(-pi) q[1];
rz(-0.64651368) q[2];
sx q[2];
rz(-2.1120089) q[2];
sx q[2];
rz(0.806964) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8710526) q[1];
sx q[1];
rz(-0.72977282) q[1];
sx q[1];
rz(-2.5403028) q[1];
rz(-pi) q[2];
rz(1.8098894) q[3];
sx q[3];
rz(-2.0951562) q[3];
sx q[3];
rz(0.6839377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0974836) q[2];
sx q[2];
rz(-2.5554843) q[2];
sx q[2];
rz(3.0666472) q[2];
rz(0.99503851) q[3];
sx q[3];
rz(-2.0375662) q[3];
sx q[3];
rz(1.9869355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5950045) q[0];
sx q[0];
rz(-1.1052479) q[0];
sx q[0];
rz(2.8357491) q[0];
rz(2.2587237) q[1];
sx q[1];
rz(-0.63778937) q[1];
sx q[1];
rz(-0.84552228) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99793154) q[0];
sx q[0];
rz(-0.73135644) q[0];
sx q[0];
rz(1.2242009) q[0];
rz(3.1404308) q[2];
sx q[2];
rz(-1.8550623) q[2];
sx q[2];
rz(-2.6068316) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9915139) q[1];
sx q[1];
rz(-2.2631915) q[1];
sx q[1];
rz(-2.451339) q[1];
rz(-1.7843397) q[3];
sx q[3];
rz(-2.41666) q[3];
sx q[3];
rz(-1.2675347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4800097) q[2];
sx q[2];
rz(-2.2955387) q[2];
sx q[2];
rz(2.0704849) q[2];
rz(2.8969911) q[3];
sx q[3];
rz(-2.1961803) q[3];
sx q[3];
rz(-1.5379813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59413183) q[0];
sx q[0];
rz(-2.5762711) q[0];
sx q[0];
rz(0.93455899) q[0];
rz(2.4681828) q[1];
sx q[1];
rz(-2.5118561) q[1];
sx q[1];
rz(-3.0455132) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5058511) q[0];
sx q[0];
rz(-1.4560501) q[0];
sx q[0];
rz(0.026057292) q[0];
rz(-1.3911909) q[2];
sx q[2];
rz(-1.7471004) q[2];
sx q[2];
rz(1.4063032) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1818349) q[1];
sx q[1];
rz(-1.3980306) q[1];
sx q[1];
rz(-2.8448425) q[1];
rz(-pi) q[2];
rz(0.096287829) q[3];
sx q[3];
rz(-0.91595338) q[3];
sx q[3];
rz(1.8055467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38976321) q[2];
sx q[2];
rz(-1.4233754) q[2];
sx q[2];
rz(-0.54330379) q[2];
rz(3.1144888) q[3];
sx q[3];
rz(-2.6684941) q[3];
sx q[3];
rz(1.0075587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84697023) q[0];
sx q[0];
rz(-2.4610418) q[0];
sx q[0];
rz(2.3934225) q[0];
rz(1.6698042) q[1];
sx q[1];
rz(-1.4133778) q[1];
sx q[1];
rz(-2.4459623) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78139898) q[0];
sx q[0];
rz(-1.8273786) q[0];
sx q[0];
rz(0.11959038) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3514481) q[2];
sx q[2];
rz(-2.2168014) q[2];
sx q[2];
rz(1.9771119) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4804617) q[1];
sx q[1];
rz(-0.5782402) q[1];
sx q[1];
rz(-2.4938512) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4616667) q[3];
sx q[3];
rz(-1.6798476) q[3];
sx q[3];
rz(-2.8449279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8814322) q[2];
sx q[2];
rz(-2.3023534) q[2];
sx q[2];
rz(-3.0879171) q[2];
rz(-1.7565049) q[3];
sx q[3];
rz(-1.799492) q[3];
sx q[3];
rz(-2.9746941) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78228918) q[0];
sx q[0];
rz(-1.259869) q[0];
sx q[0];
rz(-2.2480929) q[0];
rz(2.9529849) q[1];
sx q[1];
rz(-0.74917561) q[1];
sx q[1];
rz(0.20208727) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2938439) q[0];
sx q[0];
rz(-1.5408775) q[0];
sx q[0];
rz(-2.317753) q[0];
rz(-pi) q[1];
rz(-1.484732) q[2];
sx q[2];
rz(-2.616639) q[2];
sx q[2];
rz(-0.88366547) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.34101453) q[1];
sx q[1];
rz(-2.1034001) q[1];
sx q[1];
rz(2.8982603) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5361183) q[3];
sx q[3];
rz(-2.6993963) q[3];
sx q[3];
rz(-2.1911606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.14822745) q[2];
sx q[2];
rz(-1.8348285) q[2];
sx q[2];
rz(0.6692878) q[2];
rz(2.3547442) q[3];
sx q[3];
rz(-1.9653178) q[3];
sx q[3];
rz(2.440786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65449077) q[0];
sx q[0];
rz(-2.632532) q[0];
sx q[0];
rz(1.85602) q[0];
rz(-0.14699832) q[1];
sx q[1];
rz(-1.1451984) q[1];
sx q[1];
rz(-0.95050341) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72785801) q[0];
sx q[0];
rz(-2.0986027) q[0];
sx q[0];
rz(2.1377273) q[0];
rz(-pi) q[1];
rz(-3.0790607) q[2];
sx q[2];
rz(-1.034063) q[2];
sx q[2];
rz(1.5471396) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9626395) q[1];
sx q[1];
rz(-0.8118642) q[1];
sx q[1];
rz(-2.35747) q[1];
rz(-3.0026423) q[3];
sx q[3];
rz(-1.6873056) q[3];
sx q[3];
rz(-0.22288915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46108437) q[2];
sx q[2];
rz(-0.66404873) q[2];
sx q[2];
rz(0.39592478) q[2];
rz(0.33216533) q[3];
sx q[3];
rz(-1.1484523) q[3];
sx q[3];
rz(-2.240326) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0300776) q[0];
sx q[0];
rz(-0.60965309) q[0];
sx q[0];
rz(-1.6994221) q[0];
rz(-3.1051927) q[1];
sx q[1];
rz(-1.8489238) q[1];
sx q[1];
rz(1.363149) q[1];
rz(-1.2716952) q[2];
sx q[2];
rz(-0.80812412) q[2];
sx q[2];
rz(-3.0199188) q[2];
rz(-0.66988173) q[3];
sx q[3];
rz(-1.6139779) q[3];
sx q[3];
rz(-1.1077751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
