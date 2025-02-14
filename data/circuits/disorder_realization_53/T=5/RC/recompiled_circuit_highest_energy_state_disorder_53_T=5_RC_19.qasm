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
rz(-0.1102912) q[0];
sx q[0];
rz(-1.2854486) q[0];
sx q[0];
rz(-0.31874803) q[0];
rz(1.1822074) q[1];
sx q[1];
rz(-1.6637586) q[1];
sx q[1];
rz(1.9786662) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6642374) q[0];
sx q[0];
rz(-1.2384982) q[0];
sx q[0];
rz(2.8927781) q[0];
rz(3.0672795) q[2];
sx q[2];
rz(-0.19467672) q[2];
sx q[2];
rz(-3.1271324) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9388814) q[1];
sx q[1];
rz(-1.8716646) q[1];
sx q[1];
rz(1.2510514) q[1];
x q[2];
rz(1.5908904) q[3];
sx q[3];
rz(-0.99599518) q[3];
sx q[3];
rz(1.2736831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.17243871) q[2];
sx q[2];
rz(-0.46010751) q[2];
sx q[2];
rz(3.0058506) q[2];
rz(-2.8067449) q[3];
sx q[3];
rz(-1.2232774) q[3];
sx q[3];
rz(2.7443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927476) q[0];
sx q[0];
rz(-2.1362342) q[0];
sx q[0];
rz(-2.1951065) q[0];
rz(1.954156) q[1];
sx q[1];
rz(-0.58381909) q[1];
sx q[1];
rz(-0.28396398) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0683354) q[0];
sx q[0];
rz(-0.26480246) q[0];
sx q[0];
rz(2.3786484) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2912441) q[2];
sx q[2];
rz(-0.7182622) q[2];
sx q[2];
rz(0.58092434) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.52580035) q[1];
sx q[1];
rz(-2.6944571) q[1];
sx q[1];
rz(-2.8195803) q[1];
rz(-0.3881298) q[3];
sx q[3];
rz(-1.8947269) q[3];
sx q[3];
rz(0.38979724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2076608) q[2];
sx q[2];
rz(-1.8457103) q[2];
sx q[2];
rz(-0.80061039) q[2];
rz(-1.8396395) q[3];
sx q[3];
rz(-1.1336361) q[3];
sx q[3];
rz(3.083526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12450739) q[0];
sx q[0];
rz(-0.43211102) q[0];
sx q[0];
rz(0.79297638) q[0];
rz(-2.4555581) q[1];
sx q[1];
rz(-1.7162836) q[1];
sx q[1];
rz(2.3763903) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8350462) q[0];
sx q[0];
rz(-1.5126499) q[0];
sx q[0];
rz(-1.3594199) q[0];
rz(-pi) q[1];
rz(0.55770959) q[2];
sx q[2];
rz(-0.06076014) q[2];
sx q[2];
rz(2.9094158) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0211693) q[1];
sx q[1];
rz(-0.83008728) q[1];
sx q[1];
rz(-0.3786901) q[1];
rz(-0.63777615) q[3];
sx q[3];
rz(-0.90472639) q[3];
sx q[3];
rz(0.037925743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5746295) q[2];
sx q[2];
rz(-0.69090635) q[2];
sx q[2];
rz(1.8951269) q[2];
rz(1.9555107) q[3];
sx q[3];
rz(-1.7639152) q[3];
sx q[3];
rz(0.20868364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935811) q[0];
sx q[0];
rz(-1.8210541) q[0];
sx q[0];
rz(-3.0082974) q[0];
rz(-2.8721299) q[1];
sx q[1];
rz(-0.91454426) q[1];
sx q[1];
rz(2.6752313) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3379221) q[0];
sx q[0];
rz(-2.161177) q[0];
sx q[0];
rz(-1.7147816) q[0];
rz(-pi) q[1];
rz(1.3106636) q[2];
sx q[2];
rz(-0.64412737) q[2];
sx q[2];
rz(0.82714236) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34750858) q[1];
sx q[1];
rz(-1.7219436) q[1];
sx q[1];
rz(0.61194365) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6367988) q[3];
sx q[3];
rz(-1.0579234) q[3];
sx q[3];
rz(2.6309225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.93126297) q[2];
sx q[2];
rz(-2.7208021) q[2];
sx q[2];
rz(-2.8974864) q[2];
rz(-1.3611475) q[3];
sx q[3];
rz(-0.94435349) q[3];
sx q[3];
rz(-1.2066427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72531438) q[0];
sx q[0];
rz(-2.2940574) q[0];
sx q[0];
rz(-3.0552979) q[0];
rz(1.7591954) q[1];
sx q[1];
rz(-2.081213) q[1];
sx q[1];
rz(-1.28654) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3495958) q[0];
sx q[0];
rz(-1.1832058) q[0];
sx q[0];
rz(2.9583065) q[0];
rz(-pi) q[1];
rz(2.9164739) q[2];
sx q[2];
rz(-2.3500366) q[2];
sx q[2];
rz(2.5560372) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4205192) q[1];
sx q[1];
rz(-2.2700078) q[1];
sx q[1];
rz(3.0341604) q[1];
rz(-pi) q[2];
rz(2.6161353) q[3];
sx q[3];
rz(-2.1920929) q[3];
sx q[3];
rz(-2.3182825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4061819) q[2];
sx q[2];
rz(-1.4124796) q[2];
sx q[2];
rz(-2.0162876) q[2];
rz(-2.3451037) q[3];
sx q[3];
rz(-0.73553604) q[3];
sx q[3];
rz(0.19851941) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7859802) q[0];
sx q[0];
rz(-0.2934083) q[0];
sx q[0];
rz(-0.32817131) q[0];
rz(2.7525821) q[1];
sx q[1];
rz(-1.5711454) q[1];
sx q[1];
rz(1.4555812) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9756115) q[0];
sx q[0];
rz(-1.7811547) q[0];
sx q[0];
rz(-1.4882404) q[0];
rz(-1.7057538) q[2];
sx q[2];
rz(-1.696387) q[2];
sx q[2];
rz(-1.7474945) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.38992369) q[1];
sx q[1];
rz(-1.8018541) q[1];
sx q[1];
rz(2.2242145) q[1];
rz(-pi) q[2];
rz(1.2120655) q[3];
sx q[3];
rz(-1.609612) q[3];
sx q[3];
rz(-1.9129378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97633156) q[2];
sx q[2];
rz(-0.92504048) q[2];
sx q[2];
rz(0.46249214) q[2];
rz(-2.1180604) q[3];
sx q[3];
rz(-2.0868802) q[3];
sx q[3];
rz(0.83824497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8031215) q[0];
sx q[0];
rz(-1.4305038) q[0];
sx q[0];
rz(-2.9679003) q[0];
rz(-2.9202785) q[1];
sx q[1];
rz(-0.68562713) q[1];
sx q[1];
rz(1.3632704) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6741329) q[0];
sx q[0];
rz(-1.1674351) q[0];
sx q[0];
rz(-2.9441071) q[0];
rz(-pi) q[1];
rz(-2.3214746) q[2];
sx q[2];
rz(-3.0228399) q[2];
sx q[2];
rz(1.135716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9611836) q[1];
sx q[1];
rz(-2.3272657) q[1];
sx q[1];
rz(-0.66988887) q[1];
rz(-0.84529384) q[3];
sx q[3];
rz(-0.35652439) q[3];
sx q[3];
rz(2.561352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26019105) q[2];
sx q[2];
rz(-1.1611725) q[2];
sx q[2];
rz(-1.4455522) q[2];
rz(-2.3675303) q[3];
sx q[3];
rz(-0.71662199) q[3];
sx q[3];
rz(1.4708446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.562302) q[0];
sx q[0];
rz(-2.2023872) q[0];
sx q[0];
rz(-0.24765177) q[0];
rz(0.66894764) q[1];
sx q[1];
rz(-1.1890143) q[1];
sx q[1];
rz(2.155969) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0730608) q[0];
sx q[0];
rz(-1.685623) q[0];
sx q[0];
rz(-2.5825744) q[0];
x q[1];
rz(1.2560166) q[2];
sx q[2];
rz(-1.0295964) q[2];
sx q[2];
rz(1.7371295) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0528763) q[1];
sx q[1];
rz(-0.50584882) q[1];
sx q[1];
rz(-0.90152503) q[1];
rz(2.5466315) q[3];
sx q[3];
rz(-1.56968) q[3];
sx q[3];
rz(0.11754712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0011657) q[2];
sx q[2];
rz(-2.3900034) q[2];
sx q[2];
rz(-1.2933732) q[2];
rz(-1.2398531) q[3];
sx q[3];
rz(-0.69609061) q[3];
sx q[3];
rz(-0.70824879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86466113) q[0];
sx q[0];
rz(-1.7129352) q[0];
sx q[0];
rz(0.96555936) q[0];
rz(2.7652265) q[1];
sx q[1];
rz(-1.8195567) q[1];
sx q[1];
rz(1.9115062) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9133643) q[0];
sx q[0];
rz(-1.3051864) q[0];
sx q[0];
rz(-1.6227386) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9006017) q[2];
sx q[2];
rz(-1.1993186) q[2];
sx q[2];
rz(-2.7895751) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4415746) q[1];
sx q[1];
rz(-0.61005615) q[1];
sx q[1];
rz(3.0185591) q[1];
x q[2];
rz(-3.0885124) q[3];
sx q[3];
rz(-0.81392589) q[3];
sx q[3];
rz(0.9313213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85912117) q[2];
sx q[2];
rz(-2.4716447) q[2];
sx q[2];
rz(-3.0100789) q[2];
rz(-2.6604743) q[3];
sx q[3];
rz(-0.93004623) q[3];
sx q[3];
rz(2.6432945) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1199101) q[0];
sx q[0];
rz(-0.57761884) q[0];
sx q[0];
rz(-0.48900327) q[0];
rz(-1.8148212) q[1];
sx q[1];
rz(-1.860447) q[1];
sx q[1];
rz(-0.16924032) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76194644) q[0];
sx q[0];
rz(-2.6874306) q[0];
sx q[0];
rz(1.2962925) q[0];
rz(-pi) q[1];
rz(1.660295) q[2];
sx q[2];
rz(-1.3907593) q[2];
sx q[2];
rz(1.5370739) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0161017) q[1];
sx q[1];
rz(-0.91052848) q[1];
sx q[1];
rz(1.0339526) q[1];
rz(-0.27948252) q[3];
sx q[3];
rz(-1.8516083) q[3];
sx q[3];
rz(-1.8008417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5883611) q[2];
sx q[2];
rz(-1.0855805) q[2];
sx q[2];
rz(2.9300743) q[2];
rz(-1.616098) q[3];
sx q[3];
rz(-2.1414521) q[3];
sx q[3];
rz(-0.26743993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78032988) q[0];
sx q[0];
rz(-1.7050671) q[0];
sx q[0];
rz(0.29062511) q[0];
rz(0.81193874) q[1];
sx q[1];
rz(-0.62807905) q[1];
sx q[1];
rz(-1.5608578) q[1];
rz(-0.042150368) q[2];
sx q[2];
rz(-0.43402815) q[2];
sx q[2];
rz(2.8708906) q[2];
rz(-0.62150443) q[3];
sx q[3];
rz(-0.92172289) q[3];
sx q[3];
rz(-0.48529101) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
