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
rz(-2.5858606) q[0];
sx q[0];
rz(-1.2795804) q[0];
sx q[0];
rz(0.32787856) q[0];
rz(0.15283395) q[1];
sx q[1];
rz(-0.489355) q[1];
sx q[1];
rz(-2.1305003) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7578825) q[0];
sx q[0];
rz(-2.5284736) q[0];
sx q[0];
rz(-2.8611373) q[0];
x q[1];
rz(1.3921803) q[2];
sx q[2];
rz(-1.1682604) q[2];
sx q[2];
rz(2.2312763) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3950306) q[1];
sx q[1];
rz(-0.85828188) q[1];
sx q[1];
rz(1.0814352) q[1];
rz(-pi) q[2];
rz(1.1527083) q[3];
sx q[3];
rz(-1.7888513) q[3];
sx q[3];
rz(0.59948987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8622387) q[2];
sx q[2];
rz(-2.3591159) q[2];
sx q[2];
rz(-1.5581101) q[2];
rz(0.33485788) q[3];
sx q[3];
rz(-1.0840651) q[3];
sx q[3];
rz(0.28086942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8849477) q[0];
sx q[0];
rz(-0.56459752) q[0];
sx q[0];
rz(-2.8096492) q[0];
rz(-0.360082) q[1];
sx q[1];
rz(-1.8372476) q[1];
sx q[1];
rz(0.28775451) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1126039) q[0];
sx q[0];
rz(-1.8191914) q[0];
sx q[0];
rz(1.7775787) q[0];
x q[1];
rz(-0.73189484) q[2];
sx q[2];
rz(-1.8127155) q[2];
sx q[2];
rz(0.11920028) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33720484) q[1];
sx q[1];
rz(-1.241695) q[1];
sx q[1];
rz(-1.2044524) q[1];
rz(-pi) q[2];
rz(-1.4075564) q[3];
sx q[3];
rz(-1.7989915) q[3];
sx q[3];
rz(1.850224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45431367) q[2];
sx q[2];
rz(-1.6087029) q[2];
sx q[2];
rz(-1.7506556) q[2];
rz(2.2823997) q[3];
sx q[3];
rz(-1.2733368) q[3];
sx q[3];
rz(0.19101645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39621064) q[0];
sx q[0];
rz(-1.5353545) q[0];
sx q[0];
rz(2.9111653) q[0];
rz(-0.51741171) q[1];
sx q[1];
rz(-2.0473862) q[1];
sx q[1];
rz(0.28883019) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5342425) q[0];
sx q[0];
rz(-1.6301883) q[0];
sx q[0];
rz(-1.1618105) q[0];
rz(-pi) q[1];
rz(-0.14634653) q[2];
sx q[2];
rz(-0.47292559) q[2];
sx q[2];
rz(2.2951916) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6284184) q[1];
sx q[1];
rz(-2.7283154) q[1];
sx q[1];
rz(-1.7032313) q[1];
rz(-pi) q[2];
rz(0.90694859) q[3];
sx q[3];
rz(-0.28093279) q[3];
sx q[3];
rz(0.30981608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0383703) q[2];
sx q[2];
rz(-2.4028845) q[2];
sx q[2];
rz(0.040741097) q[2];
rz(-3.0401958) q[3];
sx q[3];
rz(-1.8056168) q[3];
sx q[3];
rz(-1.066677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7539702) q[0];
sx q[0];
rz(-3.1356223) q[0];
sx q[0];
rz(-0.23400865) q[0];
rz(-0.19800828) q[1];
sx q[1];
rz(-2.104069) q[1];
sx q[1];
rz(2.1580946) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6925544) q[0];
sx q[0];
rz(-2.6561167) q[0];
sx q[0];
rz(-1.4794769) q[0];
x q[1];
rz(2.6279672) q[2];
sx q[2];
rz(-1.6503929) q[2];
sx q[2];
rz(-0.74515504) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2605839) q[1];
sx q[1];
rz(-2.1420292) q[1];
sx q[1];
rz(0.84239475) q[1];
x q[2];
rz(0.03201501) q[3];
sx q[3];
rz(-1.181662) q[3];
sx q[3];
rz(1.8097046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46821758) q[2];
sx q[2];
rz(-2.8523291) q[2];
sx q[2];
rz(-2.3667228) q[2];
rz(-1.8153927) q[3];
sx q[3];
rz(-1.1893136) q[3];
sx q[3];
rz(-0.51138043) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73959094) q[0];
sx q[0];
rz(-0.87257659) q[0];
sx q[0];
rz(-3.1091029) q[0];
rz(1.2292817) q[1];
sx q[1];
rz(-0.928343) q[1];
sx q[1];
rz(-0.083018735) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036839824) q[0];
sx q[0];
rz(-2.6304465) q[0];
sx q[0];
rz(-2.6698861) q[0];
rz(-pi) q[1];
rz(1.4013702) q[2];
sx q[2];
rz(-1.3063141) q[2];
sx q[2];
rz(1.5096017) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28369812) q[1];
sx q[1];
rz(-2.2304728) q[1];
sx q[1];
rz(-0.96925756) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.015542726) q[3];
sx q[3];
rz(-2.2076026) q[3];
sx q[3];
rz(-1.4525082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3518389) q[2];
sx q[2];
rz(-2.3408076) q[2];
sx q[2];
rz(-2.9808673) q[2];
rz(1.1927346) q[3];
sx q[3];
rz(-0.21218097) q[3];
sx q[3];
rz(-2.3737657) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47650325) q[0];
sx q[0];
rz(-1.1192717) q[0];
sx q[0];
rz(0.29092586) q[0];
rz(1.7581958) q[1];
sx q[1];
rz(-1.8427126) q[1];
sx q[1];
rz(-0.49984041) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2384278) q[0];
sx q[0];
rz(-0.060783371) q[0];
sx q[0];
rz(-0.77068909) q[0];
rz(-pi) q[1];
rz(0.55576365) q[2];
sx q[2];
rz(-2.2190208) q[2];
sx q[2];
rz(2.6157635) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.11092) q[1];
sx q[1];
rz(-0.77092147) q[1];
sx q[1];
rz(-2.4708807) q[1];
rz(-pi) q[2];
rz(0.860943) q[3];
sx q[3];
rz(-1.1490886) q[3];
sx q[3];
rz(0.56906869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9342039) q[2];
sx q[2];
rz(-2.5298205) q[2];
sx q[2];
rz(-0.70154166) q[2];
rz(-2.1675341) q[3];
sx q[3];
rz(-0.8166703) q[3];
sx q[3];
rz(2.2402703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3153673) q[0];
sx q[0];
rz(-2.3235445) q[0];
sx q[0];
rz(0.65959626) q[0];
rz(-1.9024128) q[1];
sx q[1];
rz(-2.0678803) q[1];
sx q[1];
rz(2.0674131) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2240018) q[0];
sx q[0];
rz(-1.5842284) q[0];
sx q[0];
rz(1.55634) q[0];
rz(-pi) q[1];
rz(-0.72788179) q[2];
sx q[2];
rz(-0.54833503) q[2];
sx q[2];
rz(0.20073433) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2905137) q[1];
sx q[1];
rz(-2.004262) q[1];
sx q[1];
rz(2.9032533) q[1];
x q[2];
rz(0.76403615) q[3];
sx q[3];
rz(-1.8978528) q[3];
sx q[3];
rz(2.5491109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2913975) q[2];
sx q[2];
rz(-0.75866282) q[2];
sx q[2];
rz(-2.5568753) q[2];
rz(1.0662474) q[3];
sx q[3];
rz(-1.9138347) q[3];
sx q[3];
rz(-2.7339981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9503815) q[0];
sx q[0];
rz(-1.9712912) q[0];
sx q[0];
rz(-0.48072746) q[0];
rz(-1.2113384) q[1];
sx q[1];
rz(-1.6119266) q[1];
sx q[1];
rz(0.617625) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6473501) q[0];
sx q[0];
rz(-1.2456919) q[0];
sx q[0];
rz(3.1371497) q[0];
x q[1];
rz(-2.5343393) q[2];
sx q[2];
rz(-2.0130081) q[2];
sx q[2];
rz(-2.2970125) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1943975) q[1];
sx q[1];
rz(-0.36809599) q[1];
sx q[1];
rz(-0.93666623) q[1];
rz(1.5354352) q[3];
sx q[3];
rz(-1.3786982) q[3];
sx q[3];
rz(2.0336322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9628613) q[2];
sx q[2];
rz(-1.3946673) q[2];
sx q[2];
rz(2.2507131) q[2];
rz(2.1122872) q[3];
sx q[3];
rz(-1.3903214) q[3];
sx q[3];
rz(-1.5911969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35128281) q[0];
sx q[0];
rz(-2.2064378) q[0];
sx q[0];
rz(1.9679605) q[0];
rz(-1.1890746) q[1];
sx q[1];
rz(-2.3937841) q[1];
sx q[1];
rz(0.48114166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8136776) q[0];
sx q[0];
rz(-2.5338533) q[0];
sx q[0];
rz(2.743676) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6716154) q[2];
sx q[2];
rz(-1.013213) q[2];
sx q[2];
rz(0.65396152) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9607842) q[1];
sx q[1];
rz(-2.167302) q[1];
sx q[1];
rz(2.9636613) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88925006) q[3];
sx q[3];
rz(-2.4947531) q[3];
sx q[3];
rz(1.7689887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22566191) q[2];
sx q[2];
rz(-1.2292726) q[2];
sx q[2];
rz(-1.4813598) q[2];
rz(0.046317421) q[3];
sx q[3];
rz(-1.9527083) q[3];
sx q[3];
rz(-1.2272629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34761053) q[0];
sx q[0];
rz(-1.0337669) q[0];
sx q[0];
rz(-0.069742918) q[0];
rz(-2.8522885) q[1];
sx q[1];
rz(-1.8569088) q[1];
sx q[1];
rz(2.0893673) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69221321) q[0];
sx q[0];
rz(-1.2468296) q[0];
sx q[0];
rz(2.762336) q[0];
rz(1.3585111) q[2];
sx q[2];
rz(-2.4683778) q[2];
sx q[2];
rz(-0.5897943) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1919162) q[1];
sx q[1];
rz(-1.9468309) q[1];
sx q[1];
rz(0.38825775) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9931562) q[3];
sx q[3];
rz(-1.9655272) q[3];
sx q[3];
rz(-1.0843866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59794402) q[2];
sx q[2];
rz(-1.9660549) q[2];
sx q[2];
rz(3.0808466) q[2];
rz(2.8525823) q[3];
sx q[3];
rz(-1.776639) q[3];
sx q[3];
rz(1.2750767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-2.65171) q[0];
sx q[0];
rz(-1.4986421) q[0];
sx q[0];
rz(-1.0810252) q[0];
rz(2.091485) q[1];
sx q[1];
rz(-2.0790015) q[1];
sx q[1];
rz(1.9008295) q[1];
rz(2.0403258) q[2];
sx q[2];
rz(-1.7055837) q[2];
sx q[2];
rz(-2.3619426) q[2];
rz(-0.28239863) q[3];
sx q[3];
rz(-1.1544607) q[3];
sx q[3];
rz(-0.37731597) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
