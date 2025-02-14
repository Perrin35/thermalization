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
rz(2.9313791) q[0];
sx q[0];
rz(5.9271521) q[0];
sx q[0];
rz(7.9071101) q[0];
rz(1.8310945) q[1];
sx q[1];
rz(3.9929462) q[1];
sx q[1];
rz(8.5265018) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0097281) q[0];
sx q[0];
rz(-2.4180331) q[0];
sx q[0];
rz(1.0148763) q[0];
rz(-2.8509628) q[2];
sx q[2];
rz(-2.0465474) q[2];
sx q[2];
rz(1.0889183) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15816244) q[1];
sx q[1];
rz(-1.3884516) q[1];
sx q[1];
rz(1.4901864) q[1];
rz(-pi) q[2];
rz(0.0086011767) q[3];
sx q[3];
rz(-1.7349173) q[3];
sx q[3];
rz(-0.65879956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5400759) q[2];
sx q[2];
rz(-2.8105152) q[2];
sx q[2];
rz(0.69851056) q[2];
rz(1.6554333) q[3];
sx q[3];
rz(-2.555002) q[3];
sx q[3];
rz(1.1521888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0450714) q[0];
sx q[0];
rz(-2.4424545) q[0];
sx q[0];
rz(-2.844098) q[0];
rz(-2.7104764) q[1];
sx q[1];
rz(-2.2747206) q[1];
sx q[1];
rz(-1.8284304) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8968412) q[0];
sx q[0];
rz(-0.8936106) q[0];
sx q[0];
rz(1.1245825) q[0];
rz(-pi) q[1];
rz(1.976491) q[2];
sx q[2];
rz(-1.5303474) q[2];
sx q[2];
rz(1.8581529) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0370223) q[1];
sx q[1];
rz(-2.3442019) q[1];
sx q[1];
rz(3.131011) q[1];
rz(-2.6752279) q[3];
sx q[3];
rz(-0.27606872) q[3];
sx q[3];
rz(1.6735145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3678652) q[2];
sx q[2];
rz(-0.44168681) q[2];
sx q[2];
rz(2.3750677) q[2];
rz(0.68486989) q[3];
sx q[3];
rz(-1.6133512) q[3];
sx q[3];
rz(-0.49055704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8399452) q[0];
sx q[0];
rz(-1.3293581) q[0];
sx q[0];
rz(-2.8369821) q[0];
rz(-2.1412762) q[1];
sx q[1];
rz(-1.1023003) q[1];
sx q[1];
rz(-0.89835483) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72422623) q[0];
sx q[0];
rz(-0.53970102) q[0];
sx q[0];
rz(-1.4090572) q[0];
x q[1];
rz(-1.1834621) q[2];
sx q[2];
rz(-0.74046052) q[2];
sx q[2];
rz(-1.0473521) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.5872941) q[1];
sx q[1];
rz(-2.0600494) q[1];
sx q[1];
rz(-2.4913408) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4130745) q[3];
sx q[3];
rz(-2.6135859) q[3];
sx q[3];
rz(-0.64380336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0072713) q[2];
sx q[2];
rz(-1.4192702) q[2];
sx q[2];
rz(2.814494) q[2];
rz(1.6359811) q[3];
sx q[3];
rz(-1.6749195) q[3];
sx q[3];
rz(-1.8065642) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026767749) q[0];
sx q[0];
rz(-2.6752495) q[0];
sx q[0];
rz(2.603671) q[0];
rz(0.7750569) q[1];
sx q[1];
rz(-1.8102976) q[1];
sx q[1];
rz(-0.34781003) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7737601) q[0];
sx q[0];
rz(-1.1080386) q[0];
sx q[0];
rz(-2.5823043) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6762907) q[2];
sx q[2];
rz(-1.4660133) q[2];
sx q[2];
rz(-1.8035165) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8533142) q[1];
sx q[1];
rz(-0.90196672) q[1];
sx q[1];
rz(2.1200256) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9640173) q[3];
sx q[3];
rz(-1.7361509) q[3];
sx q[3];
rz(-2.4158784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9378918) q[2];
sx q[2];
rz(-2.4163279) q[2];
sx q[2];
rz(2.4791278) q[2];
rz(-2.0857816) q[3];
sx q[3];
rz(-2.8231088) q[3];
sx q[3];
rz(1.8752347) q[3];
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
rz(2.1556959) q[0];
sx q[0];
rz(-0.53845423) q[0];
sx q[0];
rz(2.1180617) q[0];
rz(-2.0899978) q[1];
sx q[1];
rz(-1.6513499) q[1];
sx q[1];
rz(-0.016062707) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4996289) q[0];
sx q[0];
rz(-2.1132541) q[0];
sx q[0];
rz(-0.15799518) q[0];
rz(0.93463411) q[2];
sx q[2];
rz(-0.98630465) q[2];
sx q[2];
rz(-2.5713657) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50243584) q[1];
sx q[1];
rz(-2.0896488) q[1];
sx q[1];
rz(-1.3154037) q[1];
x q[2];
rz(1.8422801) q[3];
sx q[3];
rz(-1.5010271) q[3];
sx q[3];
rz(1.2742364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84231275) q[2];
sx q[2];
rz(-2.2613342) q[2];
sx q[2];
rz(-0.24097815) q[2];
rz(-2.0321417) q[3];
sx q[3];
rz(-1.8532608) q[3];
sx q[3];
rz(-2.747005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1269162) q[0];
sx q[0];
rz(-0.83634818) q[0];
sx q[0];
rz(-2.8447004) q[0];
rz(-1.1709921) q[1];
sx q[1];
rz(-1.820727) q[1];
sx q[1];
rz(2.6062633) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3936345) q[0];
sx q[0];
rz(-1.7553333) q[0];
sx q[0];
rz(-0.18963295) q[0];
rz(-pi) q[1];
rz(-0.29285999) q[2];
sx q[2];
rz(-2.3032078) q[2];
sx q[2];
rz(2.5590961) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6397093) q[1];
sx q[1];
rz(-1.090126) q[1];
sx q[1];
rz(-2.9021846) q[1];
rz(-pi) q[2];
rz(-3.0020079) q[3];
sx q[3];
rz(-1.4142766) q[3];
sx q[3];
rz(-2.8536882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31763306) q[2];
sx q[2];
rz(-2.3667658) q[2];
sx q[2];
rz(0.63014692) q[2];
rz(1.0050425) q[3];
sx q[3];
rz(-0.21868394) q[3];
sx q[3];
rz(2.4439243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9610577) q[0];
sx q[0];
rz(-0.068004161) q[0];
sx q[0];
rz(1.8649752) q[0];
rz(-2.0969773) q[1];
sx q[1];
rz(-1.4199384) q[1];
sx q[1];
rz(2.9081664) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2058207) q[0];
sx q[0];
rz(-2.290831) q[0];
sx q[0];
rz(-0.056179742) q[0];
rz(-pi) q[1];
rz(0.40005513) q[2];
sx q[2];
rz(-0.68184417) q[2];
sx q[2];
rz(-1.8703731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6957558) q[1];
sx q[1];
rz(-1.5323116) q[1];
sx q[1];
rz(0.74798788) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0393321) q[3];
sx q[3];
rz(-0.58057154) q[3];
sx q[3];
rz(-1.5619265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6650247) q[2];
sx q[2];
rz(-1.2580405) q[2];
sx q[2];
rz(-0.27102077) q[2];
rz(-0.34936163) q[3];
sx q[3];
rz(-0.91785279) q[3];
sx q[3];
rz(2.5640986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9409598) q[0];
sx q[0];
rz(-2.2043493) q[0];
sx q[0];
rz(1.3295133) q[0];
rz(0.26893523) q[1];
sx q[1];
rz(-2.4766998) q[1];
sx q[1];
rz(1.7609133) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5912) q[0];
sx q[0];
rz(-1.7730646) q[0];
sx q[0];
rz(-1.3161385) q[0];
rz(3.1166638) q[2];
sx q[2];
rz(-2.5453794) q[2];
sx q[2];
rz(-2.2298857) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8963651) q[1];
sx q[1];
rz(-0.62313634) q[1];
sx q[1];
rz(-1.289403) q[1];
rz(-pi) q[2];
rz(0.56048067) q[3];
sx q[3];
rz(-1.6874325) q[3];
sx q[3];
rz(2.9612535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.55714503) q[2];
sx q[2];
rz(-0.86981213) q[2];
sx q[2];
rz(2.4073041) q[2];
rz(-1.0555438) q[3];
sx q[3];
rz(-1.5493834) q[3];
sx q[3];
rz(-0.9744823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.4771117) q[0];
sx q[0];
rz(-0.25150126) q[0];
sx q[0];
rz(0.86790458) q[0];
rz(-1.8382629) q[1];
sx q[1];
rz(-1.5918599) q[1];
sx q[1];
rz(2.910639) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69567861) q[0];
sx q[0];
rz(-1.1584846) q[0];
sx q[0];
rz(-2.5492909) q[0];
rz(0.43416331) q[2];
sx q[2];
rz(-1.4462398) q[2];
sx q[2];
rz(-1.9366154) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24915299) q[1];
sx q[1];
rz(-2.0241535) q[1];
sx q[1];
rz(1.813375) q[1];
x q[2];
rz(-0.6720613) q[3];
sx q[3];
rz(-1.5717744) q[3];
sx q[3];
rz(0.031571183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93297282) q[2];
sx q[2];
rz(-1.4459556) q[2];
sx q[2];
rz(0.43819532) q[2];
rz(-2.4402601) q[3];
sx q[3];
rz(-1.067679) q[3];
sx q[3];
rz(2.7856538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1040795) q[0];
sx q[0];
rz(-2.6689745) q[0];
sx q[0];
rz(-2.0613101) q[0];
rz(2.089962) q[1];
sx q[1];
rz(-2.0104505) q[1];
sx q[1];
rz(1.0952605) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8488377) q[0];
sx q[0];
rz(-2.0803323) q[0];
sx q[0];
rz(-1.900702) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11028408) q[2];
sx q[2];
rz(-0.37026893) q[2];
sx q[2];
rz(2.8831589) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.238823) q[1];
sx q[1];
rz(-2.2939957) q[1];
sx q[1];
rz(0.1478424) q[1];
x q[2];
rz(2.0729823) q[3];
sx q[3];
rz(-2.1081703) q[3];
sx q[3];
rz(-2.2110232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0206535) q[2];
sx q[2];
rz(-1.8755269) q[2];
sx q[2];
rz(2.2361163) q[2];
rz(2.4261273) q[3];
sx q[3];
rz(-2.4162636) q[3];
sx q[3];
rz(0.65413094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.41659551) q[0];
sx q[0];
rz(-1.7323957) q[0];
sx q[0];
rz(0.61687627) q[0];
rz(1.1299409) q[1];
sx q[1];
rz(-1.2810974) q[1];
sx q[1];
rz(1.7249736) q[1];
rz(0.89492013) q[2];
sx q[2];
rz(-2.2883359) q[2];
sx q[2];
rz(-2.2872912) q[2];
rz(0.15177095) q[3];
sx q[3];
rz(-1.3625154) q[3];
sx q[3];
rz(-1.2073928) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
