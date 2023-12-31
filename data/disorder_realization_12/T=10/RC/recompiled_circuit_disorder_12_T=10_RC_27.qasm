OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(-2.4175329) q[0];
sx q[0];
rz(1.6568503) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(0.88820052) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6523525) q[0];
sx q[0];
rz(-1.3999192) q[0];
sx q[0];
rz(1.6160374) q[0];
rz(1.8445831) q[2];
sx q[2];
rz(-0.55371504) q[2];
sx q[2];
rz(-0.48534976) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.29477316) q[1];
sx q[1];
rz(-1.5702489) q[1];
sx q[1];
rz(-2.9746387) q[1];
rz(-pi) q[2];
rz(-2.4429951) q[3];
sx q[3];
rz(-1.5454096) q[3];
sx q[3];
rz(-2.4217525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28329864) q[2];
sx q[2];
rz(-2.7259939) q[2];
sx q[2];
rz(1.1179914) q[2];
rz(-0.14532146) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(-0.045923559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5730826) q[0];
sx q[0];
rz(-1.2643603) q[0];
sx q[0];
rz(-2.4867687) q[0];
rz(-1.2163935) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(-0.12589802) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4100285) q[0];
sx q[0];
rz(-2.0429987) q[0];
sx q[0];
rz(-0.41759755) q[0];
rz(-pi) q[1];
rz(0.6217896) q[2];
sx q[2];
rz(-0.99025531) q[2];
sx q[2];
rz(0.34678005) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.241908) q[1];
sx q[1];
rz(-3.0113314) q[1];
sx q[1];
rz(-1.3379407) q[1];
x q[2];
rz(0.84729362) q[3];
sx q[3];
rz(-1.5627268) q[3];
sx q[3];
rz(1.1524259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0931603) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(2.753567) q[2];
rz(1.7175425) q[3];
sx q[3];
rz(-0.63801304) q[3];
sx q[3];
rz(2.5542636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5557142) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(3.0622603) q[0];
rz(-3.0575867) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(-1.9817339) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9281611) q[0];
sx q[0];
rz(-1.0614938) q[0];
sx q[0];
rz(-0.094497735) q[0];
x q[1];
rz(0.50767501) q[2];
sx q[2];
rz(-2.1360364) q[2];
sx q[2];
rz(-0.6914247) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0283454) q[1];
sx q[1];
rz(-2.5302437) q[1];
sx q[1];
rz(-2.4436185) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4781811) q[3];
sx q[3];
rz(-2.5941656) q[3];
sx q[3];
rz(0.80143354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.40538654) q[2];
sx q[2];
rz(-1.8182886) q[2];
sx q[2];
rz(0.1082871) q[2];
rz(0.64374271) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.165034) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(-1.4820341) q[0];
rz(-0.7011134) q[1];
sx q[1];
rz(-1.2524403) q[1];
sx q[1];
rz(-2.8569417) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50243044) q[0];
sx q[0];
rz(-2.6579755) q[0];
sx q[0];
rz(1.731427) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95732032) q[2];
sx q[2];
rz(-2.7942604) q[2];
sx q[2];
rz(-0.50328244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.22420158) q[1];
sx q[1];
rz(-2.9926139) q[1];
sx q[1];
rz(-1.0148744) q[1];
x q[2];
rz(0.90162006) q[3];
sx q[3];
rz(-1.1812783) q[3];
sx q[3];
rz(-0.31182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.231679) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(-0.56751928) q[2];
rz(-2.7275758) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.652997) q[0];
sx q[0];
rz(-0.96019205) q[0];
sx q[0];
rz(1.3265142) q[0];
rz(-1.9891706) q[1];
sx q[1];
rz(-1.3782586) q[1];
sx q[1];
rz(-0.93793905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91468231) q[0];
sx q[0];
rz(-1.4786647) q[0];
sx q[0];
rz(3.1196306) q[0];
rz(-pi) q[1];
rz(-2.1300975) q[2];
sx q[2];
rz(-2.0213631) q[2];
sx q[2];
rz(-0.089103854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93463072) q[1];
sx q[1];
rz(-3.0983371) q[1];
sx q[1];
rz(0.91911493) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7004847) q[3];
sx q[3];
rz(-1.9228336) q[3];
sx q[3];
rz(2.6854533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9662629) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(1.0413292) q[2];
rz(-1.8390309) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(-1.8235824) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(2.5842216) q[0];
rz(0.56466651) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(-0.55647892) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0628478) q[0];
sx q[0];
rz(-1.9468465) q[0];
sx q[0];
rz(1.3486805) q[0];
rz(-1.3715903) q[2];
sx q[2];
rz(-2.4161985) q[2];
sx q[2];
rz(2.9273916) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4887052) q[1];
sx q[1];
rz(-2.8294551) q[1];
sx q[1];
rz(2.1991792) q[1];
rz(-pi) q[2];
rz(-2.3658386) q[3];
sx q[3];
rz(-1.5343101) q[3];
sx q[3];
rz(-3.0690103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6094728) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(1.2247941) q[2];
rz(-1.5103643) q[3];
sx q[3];
rz(-1.3204201) q[3];
sx q[3];
rz(-1.640655) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0790134) q[0];
sx q[0];
rz(-1.2591079) q[0];
sx q[0];
rz(3.0986837) q[0];
rz(-0.91730109) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(0.50061217) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6635054) q[0];
sx q[0];
rz(-2.9604719) q[0];
sx q[0];
rz(1.7630793) q[0];
rz(2.1130354) q[2];
sx q[2];
rz(-1.2985897) q[2];
sx q[2];
rz(-0.82424639) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73270479) q[1];
sx q[1];
rz(-1.9703456) q[1];
sx q[1];
rz(1.7016181) q[1];
rz(-2.385419) q[3];
sx q[3];
rz(-1.2740967) q[3];
sx q[3];
rz(0.77615661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6033972) q[2];
sx q[2];
rz(-2.0624702) q[2];
sx q[2];
rz(-2.4659757) q[2];
rz(-2.7006941) q[3];
sx q[3];
rz(-1.7023804) q[3];
sx q[3];
rz(2.6202257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065141) q[0];
sx q[0];
rz(-0.32442176) q[0];
sx q[0];
rz(2.0741529) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.503711) q[1];
sx q[1];
rz(-1.6759466) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7536403) q[0];
sx q[0];
rz(-2.6590829) q[0];
sx q[0];
rz(1.2778736) q[0];
x q[1];
rz(1.6797811) q[2];
sx q[2];
rz(-1.3853067) q[2];
sx q[2];
rz(0.89599228) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.88361909) q[1];
sx q[1];
rz(-2.5273364) q[1];
sx q[1];
rz(-2.1780464) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.500324) q[3];
sx q[3];
rz(-2.0740168) q[3];
sx q[3];
rz(-0.36676952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33891588) q[2];
sx q[2];
rz(-2.2854476) q[2];
sx q[2];
rz(-1.476293) q[2];
rz(-2.2680797) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(2.6749558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8840238) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(2.4043758) q[0];
rz(-3.1228512) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(2.2163056) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46426526) q[0];
sx q[0];
rz(-0.87809169) q[0];
sx q[0];
rz(-2.4753184) q[0];
rz(1.4383573) q[2];
sx q[2];
rz(-0.91184154) q[2];
sx q[2];
rz(3.1140285) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9463897) q[1];
sx q[1];
rz(-2.9159947) q[1];
sx q[1];
rz(0.619508) q[1];
rz(-2.0426644) q[3];
sx q[3];
rz(-0.77583757) q[3];
sx q[3];
rz(-2.8692506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3725738) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(-1.0037237) q[2];
rz(-3.051493) q[3];
sx q[3];
rz(-0.026244791) q[3];
sx q[3];
rz(1.1130921) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039624778) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(1.2596624) q[0];
rz(0.26578495) q[1];
sx q[1];
rz(-2.5140285) q[1];
sx q[1];
rz(-0.75751799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3525317) q[0];
sx q[0];
rz(-1.4333945) q[0];
sx q[0];
rz(-1.3760516) q[0];
rz(2.9847758) q[2];
sx q[2];
rz(-2.1553401) q[2];
sx q[2];
rz(2.2371694) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1265035) q[1];
sx q[1];
rz(-1.9941829) q[1];
sx q[1];
rz(-0.24366118) q[1];
x q[2];
rz(-3.1086139) q[3];
sx q[3];
rz(-2.4313201) q[3];
sx q[3];
rz(0.23746333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1511128) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(2.5027067) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(2.1209774) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4929852) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(1.6408625) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(2.5095148) q[2];
sx q[2];
rz(-2.6161604) q[2];
sx q[2];
rz(-2.5642774) q[2];
rz(-1.100148) q[3];
sx q[3];
rz(-2.5082519) q[3];
sx q[3];
rz(-2.8520907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
