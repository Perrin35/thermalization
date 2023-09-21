OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(-1.6337448) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(2.642282) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7802785) q[0];
sx q[0];
rz(-1.8457992) q[0];
sx q[0];
rz(-1.6367903) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7668031) q[2];
sx q[2];
rz(-2.4110846) q[2];
sx q[2];
rz(-2.1240049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5159113) q[1];
sx q[1];
rz(-1.7209007) q[1];
sx q[1];
rz(-2.4019269) q[1];
x q[2];
rz(-0.026704549) q[3];
sx q[3];
rz(-1.4674125) q[3];
sx q[3];
rz(-2.0863233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.22380655) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(-1.0144368) q[2];
rz(-0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6681799) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(-2.1372674) q[0];
rz(-1.6218119) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(-1.0027592) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6111345) q[0];
sx q[0];
rz(-0.76871745) q[0];
sx q[0];
rz(-2.4918633) q[0];
rz(-0.82680306) q[2];
sx q[2];
rz(-2.8892592) q[2];
sx q[2];
rz(-1.2834872) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36675378) q[1];
sx q[1];
rz(-0.60791053) q[1];
sx q[1];
rz(1.8444091) q[1];
rz(-pi) q[2];
rz(-2.5428883) q[3];
sx q[3];
rz(-1.8667392) q[3];
sx q[3];
rz(-2.242089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8721547) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(2.9906452) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(-3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8931483) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(2.7422854) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(2.2580106) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040366216) q[0];
sx q[0];
rz(-1.4063615) q[0];
sx q[0];
rz(-0.40584392) q[0];
x q[1];
rz(-1.1728889) q[2];
sx q[2];
rz(-1.8430317) q[2];
sx q[2];
rz(0.18323252) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98595847) q[1];
sx q[1];
rz(-1.2859207) q[1];
sx q[1];
rz(3.1180711) q[1];
rz(-pi) q[2];
rz(-2.482445) q[3];
sx q[3];
rz(-1.4775606) q[3];
sx q[3];
rz(-2.6574709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.8847382) q[2];
sx q[2];
rz(-0.42923129) q[2];
rz(-2.1515576) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(0.86301962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(0.4317959) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(-0.61169949) q[0];
rz(-1.1071831) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(-2.591419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69458354) q[0];
sx q[0];
rz(-0.83568474) q[0];
sx q[0];
rz(-1.6741333) q[0];
x q[1];
rz(1.6646531) q[2];
sx q[2];
rz(-1.2774733) q[2];
sx q[2];
rz(0.23767995) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.53766996) q[1];
sx q[1];
rz(-0.99891716) q[1];
sx q[1];
rz(-0.13725431) q[1];
x q[2];
rz(0.22927852) q[3];
sx q[3];
rz(-1.5203272) q[3];
sx q[3];
rz(-2.5073187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6198373) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(-0.27553976) q[2];
rz(-0.11166212) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(-0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4438641) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(2.2648947) q[0];
rz(-2.450401) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(-2.2263288) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6633776) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(2.1992654) q[0];
rz(-3.0643164) q[2];
sx q[2];
rz(-0.52646598) q[2];
sx q[2];
rz(0.72409814) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1716869) q[1];
sx q[1];
rz(-1.7261506) q[1];
sx q[1];
rz(-0.032226493) q[1];
rz(-pi) q[2];
rz(0.74113412) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(-0.65628624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(2.6780224) q[2];
rz(-2.5772337) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(1.0790496) q[0];
rz(2.5462529) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(0.39658305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81178938) q[0];
sx q[0];
rz(-1.8237231) q[0];
sx q[0];
rz(1.1008218) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55690342) q[2];
sx q[2];
rz(-0.54585005) q[2];
sx q[2];
rz(2.5715373) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.035316) q[1];
sx q[1];
rz(-2.4748487) q[1];
sx q[1];
rz(-1.4186526) q[1];
x q[2];
rz(1.8144238) q[3];
sx q[3];
rz(-1.0001839) q[3];
sx q[3];
rz(2.7057196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98809272) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(-1.68613) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-2.3144408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39847386) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(2.356785) q[0];
rz(-1.2706884) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(-1.126359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92857498) q[0];
sx q[0];
rz(-1.2733766) q[0];
sx q[0];
rz(1.6638882) q[0];
x q[1];
rz(1.3190086) q[2];
sx q[2];
rz(-2.6195824) q[2];
sx q[2];
rz(0.28282794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.78325242) q[1];
sx q[1];
rz(-1.1249152) q[1];
sx q[1];
rz(-1.7794442) q[1];
rz(1.3133615) q[3];
sx q[3];
rz(-2.7830887) q[3];
sx q[3];
rz(0.076971019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9324947) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(-1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37041935) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(-0.11288189) q[0];
rz(2.1408634) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(0.032756068) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2643124) q[0];
sx q[0];
rz(-2.7240629) q[0];
sx q[0];
rz(-2.2420922) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44185408) q[2];
sx q[2];
rz(-1.0519069) q[2];
sx q[2];
rz(1.9372802) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9084839) q[1];
sx q[1];
rz(-1.9875803) q[1];
sx q[1];
rz(0.064518708) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0524391) q[3];
sx q[3];
rz(-1.1925863) q[3];
sx q[3];
rz(-1.13812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.96413606) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(2.5194871) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(1.5266248) q[0];
rz(-2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(0.48535767) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94333121) q[0];
sx q[0];
rz(-1.7051464) q[0];
sx q[0];
rz(1.6249379) q[0];
x q[1];
rz(-1.4114686) q[2];
sx q[2];
rz(-0.97852856) q[2];
sx q[2];
rz(-1.2194022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8977412) q[1];
sx q[1];
rz(-0.85695367) q[1];
sx q[1];
rz(0.41390151) q[1];
rz(-0.58166196) q[3];
sx q[3];
rz(-1.186944) q[3];
sx q[3];
rz(1.4130842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0316524) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(-0.1594485) q[2];
rz(-1.4032646) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(-3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6195246) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(-1.2592978) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(2.5433345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78865096) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.3593332) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9293289) q[2];
sx q[2];
rz(-1.4258254) q[2];
sx q[2];
rz(-0.14130172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.64893374) q[1];
sx q[1];
rz(-1.4760541) q[1];
sx q[1];
rz(1.9365063) q[1];
rz(-0.15828295) q[3];
sx q[3];
rz(-1.7332819) q[3];
sx q[3];
rz(-0.76009258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0593807) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(2.5478798) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3020637) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(-1.9051753) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(2.9700206) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(-0.39246172) q[3];
sx q[3];
rz(-0.93467181) q[3];
sx q[3];
rz(-1.296464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
