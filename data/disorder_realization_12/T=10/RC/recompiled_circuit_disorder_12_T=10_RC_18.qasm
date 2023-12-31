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
rz(3.8656524) q[0];
sx q[0];
rz(11.081628) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(-2.2533921) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0523383) q[0];
sx q[0];
rz(-1.6153781) q[0];
sx q[0];
rz(0.17104878) q[0];
rz(-pi) q[1];
rz(-1.0339123) q[2];
sx q[2];
rz(-1.7134588) q[2];
sx q[2];
rz(2.2906274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2792714) q[1];
sx q[1];
rz(-0.16695484) q[1];
sx q[1];
rz(3.1382986) q[1];
x q[2];
rz(-1.5376484) q[3];
sx q[3];
rz(-2.2691233) q[3];
sx q[3];
rz(-0.87227277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.858294) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(-1.1179914) q[2];
rz(2.9962712) q[3];
sx q[3];
rz(-1.5829007) q[3];
sx q[3];
rz(0.045923559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5685101) q[0];
sx q[0];
rz(-1.2643603) q[0];
sx q[0];
rz(-0.65482393) q[0];
rz(1.2163935) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(0.12589802) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1032216) q[0];
sx q[0];
rz(-1.9403337) q[0];
sx q[0];
rz(-1.0612556) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89181487) q[2];
sx q[2];
rz(-1.0620772) q[2];
sx q[2];
rz(-1.5985135) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0071348) q[1];
sx q[1];
rz(-1.6975228) q[1];
sx q[1];
rz(-0.030220672) q[1];
x q[2];
rz(-1.5829854) q[3];
sx q[3];
rz(-0.72353957) q[3];
sx q[3];
rz(-0.42750588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0931603) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(-2.753567) q[2];
rz(1.7175425) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(-2.5542636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5858784) q[0];
sx q[0];
rz(-0.782574) q[0];
sx q[0];
rz(-0.079332381) q[0];
rz(-0.084005984) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(1.1598587) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2134316) q[0];
sx q[0];
rz(-2.0800989) q[0];
sx q[0];
rz(-3.0470949) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94295393) q[2];
sx q[2];
rz(-1.9938333) q[2];
sx q[2];
rz(1.1689651) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0594271) q[1];
sx q[1];
rz(-1.9485928) q[1];
sx q[1];
rz(-2.6487745) q[1];
rz(-1.0252762) q[3];
sx q[3];
rz(-1.6189515) q[3];
sx q[3];
rz(0.84850509) q[3];
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
rz(-0.64374271) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(-2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.165034) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(-1.6595586) q[0];
rz(-0.7011134) q[1];
sx q[1];
rz(-1.2524403) q[1];
sx q[1];
rz(-2.8569417) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6391622) q[0];
sx q[0];
rz(-2.6579755) q[0];
sx q[0];
rz(1.4101656) q[0];
rz(-pi) q[1];
rz(1.2830164) q[2];
sx q[2];
rz(-1.7680401) q[2];
sx q[2];
rz(-1.6523199) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9173911) q[1];
sx q[1];
rz(-0.14897878) q[1];
sx q[1];
rz(1.0148744) q[1];
rz(-pi) q[2];
rz(-2.1553667) q[3];
sx q[3];
rz(-2.3826736) q[3];
sx q[3];
rz(-0.81134568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9099137) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(0.56751928) q[2];
rz(-0.41401687) q[3];
sx q[3];
rz(-2.0040138) q[3];
sx q[3];
rz(-2.0295985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48859566) q[0];
sx q[0];
rz(-0.96019205) q[0];
sx q[0];
rz(1.3265142) q[0];
rz(-1.1524221) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(2.2036536) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2269103) q[0];
sx q[0];
rz(-1.662928) q[0];
sx q[0];
rz(-3.1196306) q[0];
x q[1];
rz(-2.3100501) q[2];
sx q[2];
rz(-0.70280308) q[2];
sx q[2];
rz(-2.0895095) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2873958) q[1];
sx q[1];
rz(-1.5445659) q[1];
sx q[1];
rz(-1.5363974) q[1];
rz(-1.7004847) q[3];
sx q[3];
rz(-1.9228336) q[3];
sx q[3];
rz(-0.45613939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1753297) q[2];
sx q[2];
rz(-2.9523409) q[2];
sx q[2];
rz(-1.0413292) q[2];
rz(1.3025618) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(-1.8235824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43907169) q[0];
sx q[0];
rz(-1.9938001) q[0];
sx q[0];
rz(0.55737108) q[0];
rz(0.56466651) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(2.5851137) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0628478) q[0];
sx q[0];
rz(-1.9468465) q[0];
sx q[0];
rz(1.7929121) q[0];
x q[1];
rz(1.7700023) q[2];
sx q[2];
rz(-0.7253941) q[2];
sx q[2];
rz(0.2142011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6528875) q[1];
sx q[1];
rz(-0.31213752) q[1];
sx q[1];
rz(-2.1991792) q[1];
x q[2];
rz(0.052080215) q[3];
sx q[3];
rz(-0.77643231) q[3];
sx q[3];
rz(1.6805502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6094728) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(1.2247941) q[2];
rz(-1.5103643) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(-1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0625793) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(0.042908948) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(-0.50061217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0450889) q[0];
sx q[0];
rz(-1.6052264) q[0];
sx q[0];
rz(-1.3929429) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0285573) q[2];
sx q[2];
rz(-1.8430029) q[2];
sx q[2];
rz(-0.82424639) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73270479) q[1];
sx q[1];
rz(-1.171247) q[1];
sx q[1];
rz(-1.7016181) q[1];
x q[2];
rz(1.9686437) q[3];
sx q[3];
rz(-0.85507353) q[3];
sx q[3];
rz(-2.0778823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5381955) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-0.67561692) q[2];
rz(2.7006941) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(-0.52136695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0764517) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(2.0741529) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(1.6759466) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38795234) q[0];
sx q[0];
rz(-0.48250972) q[0];
sx q[0];
rz(-1.8637191) q[0];
rz(-pi) q[1];
rz(2.6162716) q[2];
sx q[2];
rz(-2.9267731) q[2];
sx q[2];
rz(-0.36052442) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2579736) q[1];
sx q[1];
rz(-2.5273364) q[1];
sx q[1];
rz(-0.96354624) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74387868) q[3];
sx q[3];
rz(-0.79259593) q[3];
sx q[3];
rz(0.63046968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8026768) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(-1.476293) q[2];
rz(0.87351292) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(2.6749558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2575689) q[0];
sx q[0];
rz(-0.57618657) q[0];
sx q[0];
rz(-2.4043758) q[0];
rz(3.1228512) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(-2.2163056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6773274) q[0];
sx q[0];
rz(-2.263501) q[0];
sx q[0];
rz(2.4753184) q[0];
x q[1];
rz(0.66321744) q[2];
sx q[2];
rz(-1.466201) q[2];
sx q[2];
rz(1.6246206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.23183322) q[1];
sx q[1];
rz(-1.4405466) q[1];
sx q[1];
rz(2.9568683) q[1];
rz(-2.2889745) q[3];
sx q[3];
rz(-1.8947621) q[3];
sx q[3];
rz(0.94911239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3725738) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(-1.0037237) q[2];
rz(-0.090099661) q[3];
sx q[3];
rz(-0.026244791) q[3];
sx q[3];
rz(2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039624778) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(-1.2596624) q[0];
rz(-0.26578495) q[1];
sx q[1];
rz(-2.5140285) q[1];
sx q[1];
rz(0.75751799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.789061) q[0];
sx q[0];
rz(-1.4333945) q[0];
sx q[0];
rz(-1.3760516) q[0];
x q[1];
rz(-2.1610356) q[2];
sx q[2];
rz(-1.440181) q[2];
sx q[2];
rz(-2.5622501) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47093686) q[1];
sx q[1];
rz(-2.6568012) q[1];
sx q[1];
rz(-1.0792586) q[1];
x q[2];
rz(0.032978756) q[3];
sx q[3];
rz(-0.71027256) q[3];
sx q[3];
rz(2.9041293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1511128) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(-2.5027067) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(-2.1209774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4929852) q[0];
sx q[0];
rz(-0.47825559) q[0];
sx q[0];
rz(-0.91260845) q[0];
rz(1.5007301) q[1];
sx q[1];
rz(-0.91703569) q[1];
sx q[1];
rz(-1.348319) q[1];
rz(2.5095148) q[2];
sx q[2];
rz(-2.6161604) q[2];
sx q[2];
rz(-2.5642774) q[2];
rz(0.99132514) q[3];
sx q[3];
rz(-1.2990868) q[3];
sx q[3];
rz(2.2494863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
