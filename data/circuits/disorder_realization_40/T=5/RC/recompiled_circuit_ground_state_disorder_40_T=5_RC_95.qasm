OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6177144) q[0];
sx q[0];
rz(-0.59433794) q[0];
sx q[0];
rz(-2.8327827) q[0];
rz(0.29769695) q[1];
sx q[1];
rz(4.3990064) q[1];
sx q[1];
rz(10.036751) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2997433) q[0];
sx q[0];
rz(-2.1349944) q[0];
sx q[0];
rz(-0.10625221) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9813479) q[2];
sx q[2];
rz(-1.8547684) q[2];
sx q[2];
rz(-0.71061963) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4472297) q[1];
sx q[1];
rz(-1.2389394) q[1];
sx q[1];
rz(-2.993078) q[1];
rz(1.6211007) q[3];
sx q[3];
rz(-0.87891776) q[3];
sx q[3];
rz(1.8722201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1674126) q[2];
sx q[2];
rz(-1.1777425) q[2];
sx q[2];
rz(-0.56498945) q[2];
rz(-0.51088339) q[3];
sx q[3];
rz(-2.9027945) q[3];
sx q[3];
rz(1.6142982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086455) q[0];
sx q[0];
rz(-2.8605509) q[0];
sx q[0];
rz(2.1458022) q[0];
rz(0.58798724) q[1];
sx q[1];
rz(-2.7016787) q[1];
sx q[1];
rz(-1.8221375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4026227) q[0];
sx q[0];
rz(-1.8085294) q[0];
sx q[0];
rz(-2.0031702) q[0];
x q[1];
rz(0.075602268) q[2];
sx q[2];
rz(-1.8972862) q[2];
sx q[2];
rz(-1.9734427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4012137) q[1];
sx q[1];
rz(-1.3644364) q[1];
sx q[1];
rz(2.3310082) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4493345) q[3];
sx q[3];
rz(-1.0912885) q[3];
sx q[3];
rz(-0.012727199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1002645) q[2];
sx q[2];
rz(-1.0039971) q[2];
sx q[2];
rz(3.1190994) q[2];
rz(-3.0371173) q[3];
sx q[3];
rz(-1.6040809) q[3];
sx q[3];
rz(-0.87048602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3290038) q[0];
sx q[0];
rz(-2.1109695) q[0];
sx q[0];
rz(1.6382244) q[0];
rz(-3.0606048) q[1];
sx q[1];
rz(-2.4826725) q[1];
sx q[1];
rz(1.170105) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2572944) q[0];
sx q[0];
rz(-2.5405875) q[0];
sx q[0];
rz(-2.8544507) q[0];
rz(-pi) q[1];
rz(1.3721714) q[2];
sx q[2];
rz(-0.93588698) q[2];
sx q[2];
rz(2.9609307) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31331949) q[1];
sx q[1];
rz(-2.0707284) q[1];
sx q[1];
rz(-3.1089422) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25024869) q[3];
sx q[3];
rz(-2.0648533) q[3];
sx q[3];
rz(1.2376461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4795503) q[2];
sx q[2];
rz(-0.68632555) q[2];
sx q[2];
rz(0.043206841) q[2];
rz(-2.9442545) q[3];
sx q[3];
rz(-2.1102326) q[3];
sx q[3];
rz(-0.432338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2546805) q[0];
sx q[0];
rz(-2.0557025) q[0];
sx q[0];
rz(0.32549724) q[0];
rz(0.85572851) q[1];
sx q[1];
rz(-0.41494644) q[1];
sx q[1];
rz(2.0554481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0551712) q[0];
sx q[0];
rz(-2.060411) q[0];
sx q[0];
rz(-1.7167164) q[0];
x q[1];
rz(0.83502533) q[2];
sx q[2];
rz(-2.1673598) q[2];
sx q[2];
rz(1.2460097) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.73472527) q[1];
sx q[1];
rz(-1.4063837) q[1];
sx q[1];
rz(1.0477935) q[1];
x q[2];
rz(0.086768199) q[3];
sx q[3];
rz(-1.7132856) q[3];
sx q[3];
rz(-2.2693107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2073652) q[2];
sx q[2];
rz(-0.97565979) q[2];
sx q[2];
rz(-0.18386851) q[2];
rz(1.81987) q[3];
sx q[3];
rz(-2.4403641) q[3];
sx q[3];
rz(0.13300657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70165271) q[0];
sx q[0];
rz(-0.30783215) q[0];
sx q[0];
rz(0.93691784) q[0];
rz(2.2266455) q[1];
sx q[1];
rz(-0.85275424) q[1];
sx q[1];
rz(1.459704) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57758439) q[0];
sx q[0];
rz(-1.47808) q[0];
sx q[0];
rz(-2.2683539) q[0];
x q[1];
rz(2.1906846) q[2];
sx q[2];
rz(-1.0703329) q[2];
sx q[2];
rz(-1.9096979) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.062511584) q[1];
sx q[1];
rz(-2.3214503) q[1];
sx q[1];
rz(-0.13343498) q[1];
x q[2];
rz(1.3077626) q[3];
sx q[3];
rz(-1.1888148) q[3];
sx q[3];
rz(0.39417496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0437643) q[2];
sx q[2];
rz(-1.8057258) q[2];
sx q[2];
rz(0.14399993) q[2];
rz(-2.0740267) q[3];
sx q[3];
rz(-0.34393603) q[3];
sx q[3];
rz(2.4501154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3095793) q[0];
sx q[0];
rz(-0.80736512) q[0];
sx q[0];
rz(0.76835865) q[0];
rz(-2.2840624) q[1];
sx q[1];
rz(-2.0950967) q[1];
sx q[1];
rz(-2.5879477) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93282774) q[0];
sx q[0];
rz(-1.1195445) q[0];
sx q[0];
rz(-2.8594467) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3411936) q[2];
sx q[2];
rz(-2.3586949) q[2];
sx q[2];
rz(1.4318493) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.92451292) q[1];
sx q[1];
rz(-1.1486774) q[1];
sx q[1];
rz(-0.98084992) q[1];
rz(-0.97691128) q[3];
sx q[3];
rz(-2.2008228) q[3];
sx q[3];
rz(-1.6225369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7553317) q[2];
sx q[2];
rz(-2.7109881) q[2];
sx q[2];
rz(0.44580305) q[2];
rz(-1.8949932) q[3];
sx q[3];
rz(-1.7720902) q[3];
sx q[3];
rz(0.8500475) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069139473) q[0];
sx q[0];
rz(-0.9599762) q[0];
sx q[0];
rz(3.0015216) q[0];
rz(2.2135997) q[1];
sx q[1];
rz(-1.8047787) q[1];
sx q[1];
rz(1.450052) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6312184) q[0];
sx q[0];
rz(-1.4112177) q[0];
sx q[0];
rz(0.11103156) q[0];
x q[1];
rz(1.9533402) q[2];
sx q[2];
rz(-1.7859283) q[2];
sx q[2];
rz(2.6184788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.99388945) q[1];
sx q[1];
rz(-2.1793587) q[1];
sx q[1];
rz(-2.1046361) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18168707) q[3];
sx q[3];
rz(-1.4188926) q[3];
sx q[3];
rz(0.10703281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.037584) q[2];
sx q[2];
rz(-0.17825492) q[2];
sx q[2];
rz(-3.0577216) q[2];
rz(-2.1508079) q[3];
sx q[3];
rz(-1.1889941) q[3];
sx q[3];
rz(-1.8535463) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24638076) q[0];
sx q[0];
rz(-1.3061433) q[0];
sx q[0];
rz(2.7847248) q[0];
rz(1.6995947) q[1];
sx q[1];
rz(-2.5394963) q[1];
sx q[1];
rz(-3.1325565) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9401902) q[0];
sx q[0];
rz(-1.3497258) q[0];
sx q[0];
rz(1.8044492) q[0];
rz(-pi) q[1];
rz(-2.6918654) q[2];
sx q[2];
rz(-2.3894261) q[2];
sx q[2];
rz(-0.2470242) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3956093) q[1];
sx q[1];
rz(-1.8425178) q[1];
sx q[1];
rz(-0.094825788) q[1];
rz(0.44329109) q[3];
sx q[3];
rz(-1.6066243) q[3];
sx q[3];
rz(1.2092502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8845727) q[2];
sx q[2];
rz(-1.7947861) q[2];
sx q[2];
rz(-0.70551562) q[2];
rz(2.8722615) q[3];
sx q[3];
rz(-2.0651385) q[3];
sx q[3];
rz(-2.1721325) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2483098) q[0];
sx q[0];
rz(-0.80535424) q[0];
sx q[0];
rz(-0.70948187) q[0];
rz(-2.4995038) q[1];
sx q[1];
rz(-2.700192) q[1];
sx q[1];
rz(-0.4253687) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72634436) q[0];
sx q[0];
rz(-1.1443024) q[0];
sx q[0];
rz(0.70856673) q[0];
rz(-pi) q[1];
rz(0.9832731) q[2];
sx q[2];
rz(-0.66907489) q[2];
sx q[2];
rz(-2.4226505) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8249564) q[1];
sx q[1];
rz(-1.4805838) q[1];
sx q[1];
rz(-2.6949203) q[1];
x q[2];
rz(-2.6900979) q[3];
sx q[3];
rz(-0.62314829) q[3];
sx q[3];
rz(2.5939306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2702668) q[2];
sx q[2];
rz(-2.9623803) q[2];
sx q[2];
rz(-2.6628009) q[2];
rz(-2.791413) q[3];
sx q[3];
rz(-1.9637354) q[3];
sx q[3];
rz(-2.71463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4737074) q[0];
sx q[0];
rz(-0.29692867) q[0];
sx q[0];
rz(-0.83734751) q[0];
rz(-0.26652023) q[1];
sx q[1];
rz(-1.8217249) q[1];
sx q[1];
rz(-1.3574319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6476319) q[0];
sx q[0];
rz(-1.5742366) q[0];
sx q[0];
rz(1.6572862) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0827882) q[2];
sx q[2];
rz(-3.0072525) q[2];
sx q[2];
rz(-2.5299046) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7818004) q[1];
sx q[1];
rz(-0.6524274) q[1];
sx q[1];
rz(-0.020571938) q[1];
rz(-pi) q[2];
rz(-0.32517631) q[3];
sx q[3];
rz(-2.0132228) q[3];
sx q[3];
rz(-2.5669702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9799161) q[2];
sx q[2];
rz(-0.53407532) q[2];
sx q[2];
rz(2.7034289) q[2];
rz(-2.2257889) q[3];
sx q[3];
rz(-1.5235498) q[3];
sx q[3];
rz(0.80317909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5220779) q[0];
sx q[0];
rz(-1.8997471) q[0];
sx q[0];
rz(2.1258623) q[0];
rz(3.0370514) q[1];
sx q[1];
rz(-1.3019982) q[1];
sx q[1];
rz(-1.7560538) q[1];
rz(1.9011433) q[2];
sx q[2];
rz(-1.5120244) q[2];
sx q[2];
rz(2.8373847) q[2];
rz(-2.2833179) q[3];
sx q[3];
rz(-2.0701874) q[3];
sx q[3];
rz(1.7449995) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
