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
rz(1.2443378) q[0];
sx q[0];
rz(-0.79255784) q[0];
sx q[0];
rz(2.8886524) q[0];
rz(5.5090299) q[1];
sx q[1];
rz(2.7634662) q[1];
sx q[1];
rz(9.0378349) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88322631) q[0];
sx q[0];
rz(-0.41825003) q[0];
sx q[0];
rz(-2.2453488) q[0];
rz(-pi) q[1];
rz(3.0091506) q[2];
sx q[2];
rz(-0.58183607) q[2];
sx q[2];
rz(-1.0755607) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9305715) q[1];
sx q[1];
rz(-1.4912349) q[1];
sx q[1];
rz(-3.0573175) q[1];
rz(1.0798595) q[3];
sx q[3];
rz(-2.1055974) q[3];
sx q[3];
rz(-0.25568572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0953377) q[2];
sx q[2];
rz(-0.75901186) q[2];
sx q[2];
rz(1.7389899) q[2];
rz(1.7272353) q[3];
sx q[3];
rz(-0.71310133) q[3];
sx q[3];
rz(-2.6426017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6801179) q[0];
sx q[0];
rz(-2.6533227) q[0];
sx q[0];
rz(-0.71037355) q[0];
rz(-1.0164227) q[1];
sx q[1];
rz(-1.8417532) q[1];
sx q[1];
rz(-2.3764835) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23952627) q[0];
sx q[0];
rz(-1.8886811) q[0];
sx q[0];
rz(2.5796298) q[0];
rz(-2.2369628) q[2];
sx q[2];
rz(-2.1470765) q[2];
sx q[2];
rz(1.2699256) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.73280947) q[1];
sx q[1];
rz(-1.2460695) q[1];
sx q[1];
rz(0.88692437) q[1];
rz(-0.14987544) q[3];
sx q[3];
rz(-1.9016163) q[3];
sx q[3];
rz(-0.56043032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0244828) q[2];
sx q[2];
rz(-2.4958002) q[2];
sx q[2];
rz(-0.70466858) q[2];
rz(-0.17413983) q[3];
sx q[3];
rz(-0.87248674) q[3];
sx q[3];
rz(2.7599938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0524549) q[0];
sx q[0];
rz(-1.8280886) q[0];
sx q[0];
rz(-0.88743368) q[0];
rz(1.8916091) q[1];
sx q[1];
rz(-1.9307815) q[1];
sx q[1];
rz(1.3608305) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75653532) q[0];
sx q[0];
rz(-0.65246118) q[0];
sx q[0];
rz(0.65076179) q[0];
rz(-pi) q[1];
rz(2.0949006) q[2];
sx q[2];
rz(-2.9827318) q[2];
sx q[2];
rz(-2.9422974) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2584784) q[1];
sx q[1];
rz(-0.7683903) q[1];
sx q[1];
rz(1.7151296) q[1];
x q[2];
rz(1.4668457) q[3];
sx q[3];
rz(-2.7505005) q[3];
sx q[3];
rz(1.1268738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36797324) q[2];
sx q[2];
rz(-2.7757288) q[2];
sx q[2];
rz(-1.492929) q[2];
rz(-0.44935539) q[3];
sx q[3];
rz(-1.8466693) q[3];
sx q[3];
rz(-1.9213283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0538977) q[0];
sx q[0];
rz(-0.59006417) q[0];
sx q[0];
rz(-1.4087403) q[0];
rz(2.159481) q[1];
sx q[1];
rz(-1.5487919) q[1];
sx q[1];
rz(-1.7353479) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3228673) q[0];
sx q[0];
rz(-1.7850842) q[0];
sx q[0];
rz(-3.0279798) q[0];
rz(-0.92806973) q[2];
sx q[2];
rz(-3.0756223) q[2];
sx q[2];
rz(-0.82974354) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50735578) q[1];
sx q[1];
rz(-2.0869531) q[1];
sx q[1];
rz(2.8576351) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.222417) q[3];
sx q[3];
rz(-0.67871415) q[3];
sx q[3];
rz(-1.7308337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9881607) q[2];
sx q[2];
rz(-2.1250696) q[2];
sx q[2];
rz(-2.9856258) q[2];
rz(-0.65034741) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(-0.11421886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0735737) q[0];
sx q[0];
rz(-1.2696215) q[0];
sx q[0];
rz(-0.20075783) q[0];
rz(-1.9643895) q[1];
sx q[1];
rz(-0.8526082) q[1];
sx q[1];
rz(-1.1600561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2748799) q[0];
sx q[0];
rz(-1.794941) q[0];
sx q[0];
rz(0.17205162) q[0];
rz(-pi) q[1];
rz(1.6408553) q[2];
sx q[2];
rz(-1.4558917) q[2];
sx q[2];
rz(-0.15945486) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1694035) q[1];
sx q[1];
rz(-0.47220017) q[1];
sx q[1];
rz(-0.71581033) q[1];
rz(-1.2358642) q[3];
sx q[3];
rz(-2.5865002) q[3];
sx q[3];
rz(2.5344283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7572299) q[2];
sx q[2];
rz(-0.94567662) q[2];
sx q[2];
rz(0.45822701) q[2];
rz(0.630817) q[3];
sx q[3];
rz(-1.2354555) q[3];
sx q[3];
rz(-0.49921504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45796564) q[0];
sx q[0];
rz(-0.85245913) q[0];
sx q[0];
rz(3.1347347) q[0];
rz(0.231617) q[1];
sx q[1];
rz(-1.3791142) q[1];
sx q[1];
rz(-2.0793656) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6303321) q[0];
sx q[0];
rz(-2.8046135) q[0];
sx q[0];
rz(-2.760051) q[0];
x q[1];
rz(1.3590603) q[2];
sx q[2];
rz(-1.0811624) q[2];
sx q[2];
rz(2.8608866) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5353229) q[1];
sx q[1];
rz(-2.4237666) q[1];
sx q[1];
rz(2.7503783) q[1];
rz(-2.9763664) q[3];
sx q[3];
rz(-2.2731785) q[3];
sx q[3];
rz(-0.14152292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.14083938) q[2];
sx q[2];
rz(-1.8893628) q[2];
sx q[2];
rz(-1.8737277) q[2];
rz(2.2800692) q[3];
sx q[3];
rz(-1.0637161) q[3];
sx q[3];
rz(1.4388194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93290257) q[0];
sx q[0];
rz(-0.33354315) q[0];
sx q[0];
rz(2.9260337) q[0];
rz(-2.6453099) q[1];
sx q[1];
rz(-1.9535306) q[1];
sx q[1];
rz(0.15942474) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8474583) q[0];
sx q[0];
rz(-2.205664) q[0];
sx q[0];
rz(-3.1383187) q[0];
rz(2.6706598) q[2];
sx q[2];
rz(-1.3976946) q[2];
sx q[2];
rz(-1.7488232) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9262808) q[1];
sx q[1];
rz(-0.23880733) q[1];
sx q[1];
rz(1.396149) q[1];
rz(-pi) q[2];
rz(0.88626363) q[3];
sx q[3];
rz(-0.66826398) q[3];
sx q[3];
rz(-1.22235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7364007) q[2];
sx q[2];
rz(-2.0173732) q[2];
sx q[2];
rz(-2.6250725) q[2];
rz(-1.9308331) q[3];
sx q[3];
rz(-1.4529994) q[3];
sx q[3];
rz(1.3794544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6571758) q[0];
sx q[0];
rz(-0.076198904) q[0];
sx q[0];
rz(-2.6013689) q[0];
rz(1.6172488) q[1];
sx q[1];
rz(-1.0018307) q[1];
sx q[1];
rz(1.8744972) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1138685) q[0];
sx q[0];
rz(-2.0696215) q[0];
sx q[0];
rz(-1.8201102) q[0];
rz(-pi) q[1];
rz(2.3643963) q[2];
sx q[2];
rz(-2.1146449) q[2];
sx q[2];
rz(1.0842619) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6798576) q[1];
sx q[1];
rz(-1.2038132) q[1];
sx q[1];
rz(-1.7253897) q[1];
rz(2.0319875) q[3];
sx q[3];
rz(-2.4465585) q[3];
sx q[3];
rz(2.024533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8871062) q[2];
sx q[2];
rz(-1.1595414) q[2];
sx q[2];
rz(-0.17733388) q[2];
rz(-1.0698498) q[3];
sx q[3];
rz(-2.0900574) q[3];
sx q[3];
rz(0.18226084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7808481) q[0];
sx q[0];
rz(-1.1540664) q[0];
sx q[0];
rz(3.0408707) q[0];
rz(-1.992647) q[1];
sx q[1];
rz(-1.9457685) q[1];
sx q[1];
rz(-1.9675868) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35648221) q[0];
sx q[0];
rz(-2.504341) q[0];
sx q[0];
rz(-2.8204945) q[0];
x q[1];
rz(0.10020013) q[2];
sx q[2];
rz(-0.76876193) q[2];
sx q[2];
rz(0.59610808) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8619949) q[1];
sx q[1];
rz(-1.5513541) q[1];
sx q[1];
rz(1.5330381) q[1];
rz(-pi) q[2];
rz(2.2618746) q[3];
sx q[3];
rz(-1.6501325) q[3];
sx q[3];
rz(-0.21896958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9775057) q[2];
sx q[2];
rz(-1.6567433) q[2];
sx q[2];
rz(-0.40317765) q[2];
rz(-0.66344231) q[3];
sx q[3];
rz(-0.57938975) q[3];
sx q[3];
rz(0.0042393953) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73096257) q[0];
sx q[0];
rz(-2.0616489) q[0];
sx q[0];
rz(-0.078911111) q[0];
rz(2.2321189) q[1];
sx q[1];
rz(-2.0957004) q[1];
sx q[1];
rz(0.39631072) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27061227) q[0];
sx q[0];
rz(-1.5798693) q[0];
sx q[0];
rz(-1.5551644) q[0];
x q[1];
rz(0.0070856241) q[2];
sx q[2];
rz(-0.79372294) q[2];
sx q[2];
rz(-2.7929579) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2036613) q[1];
sx q[1];
rz(-1.0962558) q[1];
sx q[1];
rz(1.8349667) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5509866) q[3];
sx q[3];
rz(-2.1175623) q[3];
sx q[3];
rz(1.3133698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4497455) q[2];
sx q[2];
rz(-0.88374603) q[2];
sx q[2];
rz(0.79908243) q[2];
rz(0.07829047) q[3];
sx q[3];
rz(-1.8872063) q[3];
sx q[3];
rz(-2.4962795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4770724) q[0];
sx q[0];
rz(-1.399853) q[0];
sx q[0];
rz(0.54367263) q[0];
rz(1.1935344) q[1];
sx q[1];
rz(-1.8652893) q[1];
sx q[1];
rz(0.51020772) q[1];
rz(1.9075248) q[2];
sx q[2];
rz(-0.48067817) q[2];
sx q[2];
rz(-0.59412063) q[2];
rz(-2.067749) q[3];
sx q[3];
rz(-0.90917773) q[3];
sx q[3];
rz(-2.3346693) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
