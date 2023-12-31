OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(-1.7320002) q[0];
sx q[0];
rz(1.4341266) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(6.8847818) q[1];
sx q[1];
rz(9.8431982) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9040065) q[0];
sx q[0];
rz(-1.307784) q[0];
sx q[0];
rz(-1.0884652) q[0];
x q[1];
rz(1.7069874) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(1.9905123) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99416713) q[1];
sx q[1];
rz(-2.6174183) q[1];
sx q[1];
rz(-1.0679354) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1316142) q[3];
sx q[3];
rz(-1.7712799) q[3];
sx q[3];
rz(0.85103121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.819954) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(2.3036172) q[2];
rz(-0.49301246) q[3];
sx q[3];
rz(-0.27291441) q[3];
sx q[3];
rz(3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74719602) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(-1.2778506) q[0];
rz(-0.17678075) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(0.4321672) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1457739) q[0];
sx q[0];
rz(-2.5403025) q[0];
sx q[0];
rz(-2.6390618) q[0];
rz(-pi) q[1];
rz(-0.400153) q[2];
sx q[2];
rz(-1.9674941) q[2];
sx q[2];
rz(2.951159) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9464311) q[1];
sx q[1];
rz(-0.45743194) q[1];
sx q[1];
rz(-1.3034348) q[1];
rz(-pi) q[2];
rz(-2.1357972) q[3];
sx q[3];
rz(-1.5497991) q[3];
sx q[3];
rz(-0.82783031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5923578) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(0.48669997) q[2];
rz(-1.7633847) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(2.6087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.10087) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(0.41734636) q[0];
rz(1.6529282) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(0.506385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7931472) q[0];
sx q[0];
rz(-1.2198824) q[0];
sx q[0];
rz(0.1883513) q[0];
x q[1];
rz(-0.73451368) q[2];
sx q[2];
rz(-1.3909512) q[2];
sx q[2];
rz(0.70914662) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3370034) q[1];
sx q[1];
rz(-2.5924006) q[1];
sx q[1];
rz(-2.5227491) q[1];
rz(-1.9403463) q[3];
sx q[3];
rz(-1.4003716) q[3];
sx q[3];
rz(0.91526645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.69904077) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(-2.55012) q[2];
rz(2.5555723) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.1699003) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(2.3024094) q[0];
rz(0.025578586) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(-1.5930088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65028134) q[0];
sx q[0];
rz(-2.6822753) q[0];
sx q[0];
rz(-1.7772872) q[0];
rz(-pi) q[1];
rz(-2.3392802) q[2];
sx q[2];
rz(-1.017184) q[2];
sx q[2];
rz(-0.72788903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6099445) q[1];
sx q[1];
rz(-0.87024401) q[1];
sx q[1];
rz(2.080426) q[1];
rz(-pi) q[2];
rz(2.7110093) q[3];
sx q[3];
rz(-1.8292556) q[3];
sx q[3];
rz(0.57627288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30535355) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(-0.099686064) q[2];
rz(0.95885197) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(-1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56617671) q[0];
sx q[0];
rz(-1.7594936) q[0];
sx q[0];
rz(2.8856522) q[0];
rz(2.6804965) q[1];
sx q[1];
rz(-2.0979116) q[1];
sx q[1];
rz(-0.76006132) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3292023) q[0];
sx q[0];
rz(-0.15604067) q[0];
sx q[0];
rz(2.4183256) q[0];
x q[1];
rz(1.9854457) q[2];
sx q[2];
rz(-1.6290602) q[2];
sx q[2];
rz(2.0410048) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.45436817) q[1];
sx q[1];
rz(-1.6974653) q[1];
sx q[1];
rz(-1.6997937) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7759336) q[3];
sx q[3];
rz(-0.95698157) q[3];
sx q[3];
rz(-0.97096503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.57050675) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(2.467353) q[2];
rz(2.9267866) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(-3.1242483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53133416) q[0];
sx q[0];
rz(-1.4704309) q[0];
sx q[0];
rz(-1.9624788) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(1.0669605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93766312) q[0];
sx q[0];
rz(-1.5852889) q[0];
sx q[0];
rz(-0.020676215) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7166769) q[2];
sx q[2];
rz(-1.9746466) q[2];
sx q[2];
rz(-2.0770819) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8300007) q[1];
sx q[1];
rz(-2.316906) q[1];
sx q[1];
rz(-0.064706116) q[1];
rz(1.4076036) q[3];
sx q[3];
rz(-0.7810775) q[3];
sx q[3];
rz(1.2490602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.71010464) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(-2.5816494) q[2];
rz(2.4152749) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(-2.8360951) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2840435) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(-1.42111) q[0];
rz(-0.20206085) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(-0.85817671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7875123) q[0];
sx q[0];
rz(-1.4424099) q[0];
sx q[0];
rz(1.7187198) q[0];
rz(-1.467534) q[2];
sx q[2];
rz(-0.41245663) q[2];
sx q[2];
rz(-2.7445284) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8691751) q[1];
sx q[1];
rz(-3.1088447) q[1];
sx q[1];
rz(0.48519022) q[1];
x q[2];
rz(-2.5704727) q[3];
sx q[3];
rz(-1.8652328) q[3];
sx q[3];
rz(1.7366228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.28785607) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(-1.8161592) q[2];
rz(0.89007968) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(0.98852283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.4637852) q[0];
sx q[0];
rz(-0.57254922) q[0];
sx q[0];
rz(-0.37471399) q[0];
rz(-0.97887865) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(1.3495061) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6357248) q[0];
sx q[0];
rz(-2.3946107) q[0];
sx q[0];
rz(-2.5542459) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5478412) q[2];
sx q[2];
rz(-2.288753) q[2];
sx q[2];
rz(-2.9771476) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3186036) q[1];
sx q[1];
rz(-2.9530596) q[1];
sx q[1];
rz(1.3033426) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54578636) q[3];
sx q[3];
rz(-1.7958399) q[3];
sx q[3];
rz(-0.033566098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65138856) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(-2.6728969) q[2];
rz(-1.1941341) q[3];
sx q[3];
rz(-1.8374551) q[3];
sx q[3];
rz(1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63012183) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(-1.2619031) q[0];
rz(0.17503861) q[1];
sx q[1];
rz(-1.1418399) q[1];
sx q[1];
rz(-1.6040241) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81048548) q[0];
sx q[0];
rz(-1.8750422) q[0];
sx q[0];
rz(-2.9097793) q[0];
x q[1];
rz(1.6749009) q[2];
sx q[2];
rz(-1.1541919) q[2];
sx q[2];
rz(0.57550752) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9242212) q[1];
sx q[1];
rz(-1.6842168) q[1];
sx q[1];
rz(2.4356615) q[1];
rz(-pi) q[2];
rz(-2.2213307) q[3];
sx q[3];
rz(-2.5839621) q[3];
sx q[3];
rz(-0.44616163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.039915446) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(-2.3802479) q[2];
rz(-0.90041655) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(-0.049023978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3392357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(-0.21690579) q[0];
rz(2.5096109) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(0.95473081) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8666329) q[0];
sx q[0];
rz(-1.8909591) q[0];
sx q[0];
rz(-0.76036705) q[0];
x q[1];
rz(-2.5108699) q[2];
sx q[2];
rz(-0.64881697) q[2];
sx q[2];
rz(2.0231252) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.54143822) q[1];
sx q[1];
rz(-0.46240515) q[1];
sx q[1];
rz(-0.17141436) q[1];
rz(-2.6984152) q[3];
sx q[3];
rz(-2.157353) q[3];
sx q[3];
rz(-0.12936684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1251936) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(2.1968502) q[2];
rz(-2.7567806) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(0.95705664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-0.64086296) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(2.2568933) q[1];
sx q[1];
rz(-2.2365166) q[1];
sx q[1];
rz(2.8832163) q[1];
rz(-0.62467081) q[2];
sx q[2];
rz(-0.93745898) q[2];
sx q[2];
rz(-2.6425101) q[2];
rz(-2.0142043) q[3];
sx q[3];
rz(-1.756712) q[3];
sx q[3];
rz(-2.0851019) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
