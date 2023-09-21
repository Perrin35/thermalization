OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7678087) q[0];
sx q[0];
rz(5.8435506) q[0];
sx q[0];
rz(6.2018659) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(-2.3587956) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37047526) q[0];
sx q[0];
rz(-1.3210216) q[0];
sx q[0];
rz(0.063637861) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9847203) q[2];
sx q[2];
rz(-1.1638767) q[2];
sx q[2];
rz(-3.0770609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0040972) q[1];
sx q[1];
rz(-1.4427408) q[1];
sx q[1];
rz(2.0069684) q[1];
rz(-pi) q[2];
rz(0.26773914) q[3];
sx q[3];
rz(-0.32666884) q[3];
sx q[3];
rz(0.79145811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89447442) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(0.11581126) q[2];
rz(-1.5995021) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(-2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88749921) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(2.9462573) q[0];
rz(-2.7665566) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(-2.9017752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1447906) q[0];
sx q[0];
rz(-1.7273434) q[0];
sx q[0];
rz(-2.8926204) q[0];
x q[1];
rz(0.96599483) q[2];
sx q[2];
rz(-2.7808393) q[2];
sx q[2];
rz(-0.15168562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.61566831) q[1];
sx q[1];
rz(-0.79155542) q[1];
sx q[1];
rz(3.0676003) q[1];
rz(-pi) q[2];
rz(0.014157045) q[3];
sx q[3];
rz(-1.5544958) q[3];
sx q[3];
rz(2.1099159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5043162) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(1.3298539) q[2];
rz(1.3416393) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(2.4734316) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(-1.0659165) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0695755) q[0];
sx q[0];
rz(-1.5298651) q[0];
sx q[0];
rz(1.8191562) q[0];
rz(-pi) q[1];
x q[1];
rz(1.714528) q[2];
sx q[2];
rz(-1.020069) q[2];
sx q[2];
rz(1.3155754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4855811) q[1];
sx q[1];
rz(-2.1335568) q[1];
sx q[1];
rz(2.807711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12001868) q[3];
sx q[3];
rz(-0.91306049) q[3];
sx q[3];
rz(2.0002055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.187591) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(-0.032547396) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2739094) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(2.4348863) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(-1.6548086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.810881) q[0];
sx q[0];
rz(-1.8510305) q[0];
sx q[0];
rz(0.34513721) q[0];
x q[1];
rz(-1.2232259) q[2];
sx q[2];
rz(-0.98993694) q[2];
sx q[2];
rz(-0.52085224) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1176016) q[1];
sx q[1];
rz(-1.7204493) q[1];
sx q[1];
rz(2.8531122) q[1];
rz(-0.82641043) q[3];
sx q[3];
rz(-1.5757757) q[3];
sx q[3];
rz(1.6299562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5965745) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(-2.13307) q[2];
rz(1.0962983) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(3.0116459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0356045) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(1.6812356) q[0];
rz(1.5885072) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(-0.016074093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0468633) q[0];
sx q[0];
rz(-2.097058) q[0];
sx q[0];
rz(-0.88386436) q[0];
x q[1];
rz(0.0083382567) q[2];
sx q[2];
rz(-1.1401046) q[2];
sx q[2];
rz(-2.8991933) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.15522038) q[1];
sx q[1];
rz(-1.9934137) q[1];
sx q[1];
rz(-0.95814725) q[1];
rz(1.5662976) q[3];
sx q[3];
rz(-1.9330977) q[3];
sx q[3];
rz(3.0294861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0323223) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(2.5197022) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(-0.2203075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9428403) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(-2.2348256) q[0];
rz(-2.3268907) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(1.2247359) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82299267) q[0];
sx q[0];
rz(-0.66245125) q[0];
sx q[0];
rz(1.9255161) q[0];
rz(-pi) q[1];
rz(-2.0133063) q[2];
sx q[2];
rz(-1.759521) q[2];
sx q[2];
rz(-0.21274266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.743222) q[1];
sx q[1];
rz(-2.1757158) q[1];
sx q[1];
rz(-0.066992316) q[1];
rz(-pi) q[2];
rz(1.46035) q[3];
sx q[3];
rz(-1.9697646) q[3];
sx q[3];
rz(-0.16050592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.99888745) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(-0.93969807) q[2];
rz(0.21329221) q[3];
sx q[3];
rz(-2.8010938) q[3];
sx q[3];
rz(1.3903769) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1619103) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(-0.58037037) q[0];
rz(-1.0549818) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(-2.4408128) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0133007) q[0];
sx q[0];
rz(-1.4585146) q[0];
sx q[0];
rz(-1.1178455) q[0];
rz(-2.8220196) q[2];
sx q[2];
rz(-1.459889) q[2];
sx q[2];
rz(-3.0802397) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5138445) q[1];
sx q[1];
rz(-0.57250896) q[1];
sx q[1];
rz(1.4036914) q[1];
rz(0.84007646) q[3];
sx q[3];
rz(-1.5080161) q[3];
sx q[3];
rz(1.9907794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3646399) q[2];
sx q[2];
rz(-0.31034714) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(0.74470216) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3547524) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(1.4165075) q[0];
rz(1.7658866) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(-2.2498806) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4448924) q[0];
sx q[0];
rz(-0.80388821) q[0];
sx q[0];
rz(-0.21960396) q[0];
rz(-pi) q[1];
rz(2.0217998) q[2];
sx q[2];
rz(-2.1037256) q[2];
sx q[2];
rz(-1.230513) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85372323) q[1];
sx q[1];
rz(-2.3239377) q[1];
sx q[1];
rz(-0.57662782) q[1];
rz(-pi) q[2];
rz(-1.9705087) q[3];
sx q[3];
rz(-2.4418695) q[3];
sx q[3];
rz(1.2801998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4884168) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(2.3573504) q[2];
rz(0.50576058) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(-0.23323664) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6655675) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(-0.72203565) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(-1.3649712) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9997172) q[0];
sx q[0];
rz(-1.5293855) q[0];
sx q[0];
rz(2.3977604) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7609673) q[2];
sx q[2];
rz(-2.724078) q[2];
sx q[2];
rz(1.4830358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0627928) q[1];
sx q[1];
rz(-1.9140869) q[1];
sx q[1];
rz(2.417056) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9421259) q[3];
sx q[3];
rz(-1.8623127) q[3];
sx q[3];
rz(-0.22649543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0329131) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(-1.2314679) q[2];
rz(0.03406295) q[3];
sx q[3];
rz(-1.27682) q[3];
sx q[3];
rz(-2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.05474) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(-1.6636794) q[0];
rz(-1.0832896) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(-2.1733984) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6598845) q[0];
sx q[0];
rz(-1.1904926) q[0];
sx q[0];
rz(2.7754521) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7231611) q[2];
sx q[2];
rz(-1.1088587) q[2];
sx q[2];
rz(-0.07428169) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8755175) q[1];
sx q[1];
rz(-2.7873758) q[1];
sx q[1];
rz(-2.7060899) q[1];
rz(-1.5700941) q[3];
sx q[3];
rz(-0.32214221) q[3];
sx q[3];
rz(2.567167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0782464) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(3.1402804) q[2];
rz(-1.1408268) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(-1.7989981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.39682) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(-0.36623476) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(-2.3920849) q[2];
sx q[2];
rz(-1.8029677) q[2];
sx q[2];
rz(0.36063902) q[2];
rz(2.1288539) q[3];
sx q[3];
rz(-1.2663208) q[3];
sx q[3];
rz(-2.3940621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];