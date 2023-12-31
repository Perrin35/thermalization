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
rz(-0.43963471) q[0];
sx q[0];
rz(3.0602732) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(-2.3587956) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37047526) q[0];
sx q[0];
rz(-1.820571) q[0];
sx q[0];
rz(3.0779548) q[0];
rz(1.159367) q[2];
sx q[2];
rz(-1.4268268) q[2];
sx q[2];
rz(-1.6978482) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.3738457) q[1];
sx q[1];
rz(-1.1384374) q[1];
sx q[1];
rz(3.0004764) q[1];
rz(-pi) q[2];
rz(-0.31580117) q[3];
sx q[3];
rz(-1.6557906) q[3];
sx q[3];
rz(-1.0335361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2471182) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(3.0257814) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(-2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2540934) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(-2.9462573) q[0];
rz(-2.7665566) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(0.23981747) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46562425) q[0];
sx q[0];
rz(-1.8166607) q[0];
sx q[0];
rz(1.7322391) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9302836) q[2];
sx q[2];
rz(-1.8654056) q[2];
sx q[2];
rz(-0.78795563) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5259243) q[1];
sx q[1];
rz(-0.79155542) q[1];
sx q[1];
rz(3.0676003) q[1];
rz(-pi) q[2];
rz(-1.5544942) q[3];
sx q[3];
rz(-1.5566412) q[3];
sx q[3];
rz(-0.53888884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6372765) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(1.8117388) q[2];
rz(1.7999533) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2770237) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(0.66816107) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(-1.0659165) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50915584) q[0];
sx q[0];
rz(-1.3226489) q[0];
sx q[0];
rz(0.042225348) q[0];
rz(-pi) q[1];
rz(0.22914825) q[2];
sx q[2];
rz(-2.5742968) q[2];
sx q[2];
rz(1.0457525) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26820688) q[1];
sx q[1];
rz(-1.2899439) q[1];
sx q[1];
rz(-0.98209776) q[1];
rz(-pi) q[2];
rz(-2.2320281) q[3];
sx q[3];
rz(-1.665691) q[3];
sx q[3];
rz(0.50300099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(-2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86768326) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(-0.70670635) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(-1.6548086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.810881) q[0];
sx q[0];
rz(-1.2905621) q[0];
sx q[0];
rz(2.7964554) q[0];
rz(-pi) q[1];
rz(2.6629278) q[2];
sx q[2];
rz(-0.66648167) q[2];
sx q[2];
rz(2.037231) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6390037) q[1];
sx q[1];
rz(-1.8559615) q[1];
sx q[1];
rz(-1.7267978) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1348226) q[3];
sx q[3];
rz(-2.3151708) q[3];
sx q[3];
rz(3.0778459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5965745) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-1.0085227) q[2];
rz(-1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(-0.1299468) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(1.460357) q[0];
rz(-1.5885072) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(0.016074093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8668629) q[0];
sx q[0];
rz(-2.1511973) q[0];
sx q[0];
rz(-0.6443364) q[0];
rz(-pi) q[1];
rz(-1.1400914) q[2];
sx q[2];
rz(-1.5632196) q[2];
sx q[2];
rz(1.3318782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.88735089) q[1];
sx q[1];
rz(-0.72853959) q[1];
sx q[1];
rz(2.2345047) q[1];
rz(-1.5752951) q[3];
sx q[3];
rz(-1.2084949) q[3];
sx q[3];
rz(0.11210657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1092704) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(-0.62189046) q[2];
rz(-1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.19875232) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(-2.2348256) q[0];
rz(0.81470195) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(1.2247359) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46366102) q[0];
sx q[0];
rz(-1.3555129) q[0];
sx q[0];
rz(-0.93925516) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0133063) q[2];
sx q[2];
rz(-1.759521) q[2];
sx q[2];
rz(-0.21274266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8606557) q[1];
sx q[1];
rz(-0.60815629) q[1];
sx q[1];
rz(-1.4742736) q[1];
x q[2];
rz(-0.40116061) q[3];
sx q[3];
rz(-1.6725371) q[3];
sx q[3];
rz(-1.3672369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99888745) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(-0.93969807) q[2];
rz(2.9283004) q[3];
sx q[3];
rz(-2.8010938) q[3];
sx q[3];
rz(1.7512158) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9796824) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(2.5612223) q[0];
rz(1.0549818) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(-0.70077983) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0133007) q[0];
sx q[0];
rz(-1.4585146) q[0];
sx q[0];
rz(-2.0237472) q[0];
x q[1];
rz(1.4540265) q[2];
sx q[2];
rz(-1.8883369) q[2];
sx q[2];
rz(1.5460528) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.339387) q[1];
sx q[1];
rz(-1.6610258) q[1];
sx q[1];
rz(-2.1369364) q[1];
rz(0.084214597) q[3];
sx q[3];
rz(-0.84184064) q[3];
sx q[3];
rz(0.47615151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3646399) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(-2.3968905) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72043334) q[0];
sx q[0];
rz(-1.4132858) q[0];
sx q[0];
rz(-0.79173761) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1197929) q[2];
sx q[2];
rz(-1.0378671) q[2];
sx q[2];
rz(-1.9110796) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6139639) q[1];
sx q[1];
rz(-0.91273897) q[1];
sx q[1];
rz(1.0440473) q[1];
rz(2.2303921) q[3];
sx q[3];
rz(-1.824114) q[3];
sx q[3];
rz(0.60318702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4884168) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(2.3573504) q[2];
rz(-0.50576058) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(0.23323664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4760251) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(-0.72203565) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(-1.7766215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7576335) q[0];
sx q[0];
rz(-0.74476349) q[0];
sx q[0];
rz(0.061116771) q[0];
x q[1];
rz(-0.083655595) q[2];
sx q[2];
rz(-1.9803279) q[2];
sx q[2];
rz(1.6905897) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.939146) q[1];
sx q[1];
rz(-2.2447963) q[1];
sx q[1];
rz(1.1254315) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19946675) q[3];
sx q[3];
rz(-1.27928) q[3];
sx q[3];
rz(0.22649543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1086796) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(-1.2314679) q[2];
rz(3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(-2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.05474) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(-1.6636794) q[0];
rz(-2.058303) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(-0.96819425) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4614842) q[0];
sx q[0];
rz(-0.52163863) q[0];
sx q[0];
rz(0.84043829) q[0];
rz(2.2553026) q[2];
sx q[2];
rz(-2.5286525) q[2];
sx q[2];
rz(2.431589) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4149975) q[1];
sx q[1];
rz(-1.2508878) q[1];
sx q[1];
rz(-1.725561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5714985) q[3];
sx q[3];
rz(-0.32214221) q[3];
sx q[3];
rz(2.567167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0633462) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(-3.1402804) q[2];
rz(-1.1408268) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(1.7447727) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(0.36623476) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(0.33404074) q[2];
sx q[2];
rz(-2.3636849) q[2];
sx q[2];
rz(-1.452527) q[2];
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
