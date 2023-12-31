OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5470619) q[0];
sx q[0];
rz(4.2630258) q[0];
sx q[0];
rz(6.1017258) q[0];
rz(-1.0815066) q[1];
sx q[1];
rz(-2.4681611) q[1];
sx q[1];
rz(1.0531309) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7979413) q[0];
sx q[0];
rz(-2.6290253) q[0];
sx q[0];
rz(2.0128065) q[0];
rz(-pi) q[1];
rz(1.2744781) q[2];
sx q[2];
rz(-2.1726492) q[2];
sx q[2];
rz(-1.5822496) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2145572) q[1];
sx q[1];
rz(-1.7662449) q[1];
sx q[1];
rz(2.3656225) q[1];
x q[2];
rz(1.1374723) q[3];
sx q[3];
rz(-0.37131272) q[3];
sx q[3];
rz(2.5759047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9782605) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(1.367761) q[2];
rz(-2.1286428) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(-3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0435836) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(2.5464771) q[0];
rz(-2.0960506) q[1];
sx q[1];
rz(-1.7273993) q[1];
sx q[1];
rz(-1.5140623) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1986952) q[0];
sx q[0];
rz(-1.1637582) q[0];
sx q[0];
rz(-2.1172949) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83346955) q[2];
sx q[2];
rz(-1.5916628) q[2];
sx q[2];
rz(-1.2534864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.086386911) q[1];
sx q[1];
rz(-1.3761763) q[1];
sx q[1];
rz(-2.9662532) q[1];
x q[2];
rz(-1.4684832) q[3];
sx q[3];
rz(-1.6204837) q[3];
sx q[3];
rz(1.4440086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4864768) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(-0.79616037) q[2];
rz(-2.1697309) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(-2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650836) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(0.064095108) q[0];
rz(0.31072101) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(-1.4583189) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70143914) q[0];
sx q[0];
rz(-0.15471409) q[0];
sx q[0];
rz(-1.2604777) q[0];
x q[1];
rz(-0.73297357) q[2];
sx q[2];
rz(-2.3419215) q[2];
sx q[2];
rz(3.1250931) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4255193) q[1];
sx q[1];
rz(-0.41285535) q[1];
sx q[1];
rz(-0.13456657) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0890373) q[3];
sx q[3];
rz(-0.31523809) q[3];
sx q[3];
rz(1.0759575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.98214275) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(1.3726161) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(0.94846559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963592) q[0];
sx q[0];
rz(-2.5722752) q[0];
sx q[0];
rz(1.1244208) q[0];
rz(1.9212978) q[1];
sx q[1];
rz(-1.9492457) q[1];
sx q[1];
rz(1.4845928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35614466) q[0];
sx q[0];
rz(-1.6011642) q[0];
sx q[0];
rz(1.5201475) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9494483) q[2];
sx q[2];
rz(-2.0803183) q[2];
sx q[2];
rz(-2.0083049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.209219) q[1];
sx q[1];
rz(-2.0062431) q[1];
sx q[1];
rz(2.2030764) q[1];
x q[2];
rz(1.901239) q[3];
sx q[3];
rz(-2.2664865) q[3];
sx q[3];
rz(0.89360039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1241887) q[2];
sx q[2];
rz(-1.2839395) q[2];
sx q[2];
rz(-1.0162639) q[2];
rz(1.7381564) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(-2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50773412) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(2.741709) q[0];
rz(1.9790861) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(2.9679325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6337316) q[0];
sx q[0];
rz(-1.5426239) q[0];
sx q[0];
rz(0.072576056) q[0];
rz(-2.6661751) q[2];
sx q[2];
rz(-0.43736514) q[2];
sx q[2];
rz(1.4471444) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.225519) q[1];
sx q[1];
rz(-1.5991296) q[1];
sx q[1];
rz(-1.4654935) q[1];
x q[2];
rz(-2.5912656) q[3];
sx q[3];
rz(-0.384207) q[3];
sx q[3];
rz(-1.3889544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6614723) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(2.373467) q[2];
rz(-0.85401946) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(-2.1693726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0356692) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(1.4703898) q[0];
rz(0.51180965) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(-1.1434198) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15018806) q[0];
sx q[0];
rz(-0.17230573) q[0];
sx q[0];
rz(-2.6323787) q[0];
rz(-pi) q[1];
rz(-0.46164718) q[2];
sx q[2];
rz(-0.6558154) q[2];
sx q[2];
rz(-2.0222424) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1089576) q[1];
sx q[1];
rz(-0.40954486) q[1];
sx q[1];
rz(-1.1213379) q[1];
rz(-0.30541909) q[3];
sx q[3];
rz(-2.0913887) q[3];
sx q[3];
rz(0.49314317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64289552) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(1.5303622) q[2];
rz(-1.6879843) q[3];
sx q[3];
rz(-1.0691103) q[3];
sx q[3];
rz(-2.8924275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.11944184) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(-1.4402333) q[0];
rz(0.72921905) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(1.1332606) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92111174) q[0];
sx q[0];
rz(-2.8163914) q[0];
sx q[0];
rz(0.13334206) q[0];
rz(-pi) q[1];
rz(-1.5845756) q[2];
sx q[2];
rz(-0.85859495) q[2];
sx q[2];
rz(-0.53879246) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82517351) q[1];
sx q[1];
rz(-0.9269886) q[1];
sx q[1];
rz(0.34367798) q[1];
rz(-pi) q[2];
rz(0.53448581) q[3];
sx q[3];
rz(-2.839698) q[3];
sx q[3];
rz(1.2045977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40522727) q[2];
sx q[2];
rz(-1.1232802) q[2];
sx q[2];
rz(-2.6531632) q[2];
rz(-1.3119665) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(2.5002938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.064780386) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(-3.0704165) q[0];
rz(3.1094303) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(-1.9326928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8457444) q[0];
sx q[0];
rz(-2.309572) q[0];
sx q[0];
rz(2.8315298) q[0];
rz(2.2279943) q[2];
sx q[2];
rz(-0.29619869) q[2];
sx q[2];
rz(-1.2072472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68122411) q[1];
sx q[1];
rz(-1.7011233) q[1];
sx q[1];
rz(3.1254915) q[1];
rz(0.62187059) q[3];
sx q[3];
rz(-1.8727881) q[3];
sx q[3];
rz(-2.5627476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.20748392) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(1.5709546) q[2];
rz(2.2682244) q[3];
sx q[3];
rz(-1.7375172) q[3];
sx q[3];
rz(0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97380012) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(0.31162509) q[0];
rz(-2.3198126) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(1.6315546) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25699297) q[0];
sx q[0];
rz(-1.8879226) q[0];
sx q[0];
rz(1.7811799) q[0];
rz(-pi) q[1];
rz(2.9456003) q[2];
sx q[2];
rz(-1.2336944) q[2];
sx q[2];
rz(0.54346426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3703287) q[1];
sx q[1];
rz(-2.1762098) q[1];
sx q[1];
rz(1.3243115) q[1];
rz(-pi) q[2];
rz(-2.4771677) q[3];
sx q[3];
rz(-2.8492152) q[3];
sx q[3];
rz(-1.2053306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(-0.81595016) q[2];
rz(2.6319035) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(-1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2508535) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(0.2302641) q[0];
rz(-2.5157805) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(2.4831916) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9068245) q[0];
sx q[0];
rz(-1.0848197) q[0];
sx q[0];
rz(-2.4776239) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67450972) q[2];
sx q[2];
rz(-1.3042547) q[2];
sx q[2];
rz(-0.36905497) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8734332) q[1];
sx q[1];
rz(-2.6323316) q[1];
sx q[1];
rz(1.2251653) q[1];
rz(0.55024054) q[3];
sx q[3];
rz(-0.60866683) q[3];
sx q[3];
rz(-0.88702162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4663503) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(2.8038483) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.7815536) q[3];
sx q[3];
rz(1.9256928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52453775) q[0];
sx q[0];
rz(-2.239997) q[0];
sx q[0];
rz(-0.064185113) q[0];
rz(-1.2383923) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
rz(0.94454371) q[2];
sx q[2];
rz(-1.4251475) q[2];
sx q[2];
rz(0.057694358) q[2];
rz(1.6347377) q[3];
sx q[3];
rz(-2.1422374) q[3];
sx q[3];
rz(-0.16850785) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
