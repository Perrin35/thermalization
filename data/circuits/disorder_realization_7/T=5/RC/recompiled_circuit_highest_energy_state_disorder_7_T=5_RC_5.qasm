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
rz(-3.0149674) q[0];
sx q[0];
rz(-1.5746483) q[0];
sx q[0];
rz(-2.6259165) q[0];
rz(0.22817837) q[1];
sx q[1];
rz(5.5362267) q[1];
sx q[1];
rz(9.8510392) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58881271) q[0];
sx q[0];
rz(-0.79040397) q[0];
sx q[0];
rz(-1.8145723) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2834163) q[2];
sx q[2];
rz(-2.0567908) q[2];
sx q[2];
rz(-0.11746841) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.961024) q[1];
sx q[1];
rz(-0.60446179) q[1];
sx q[1];
rz(0.1015517) q[1];
rz(0.98661749) q[3];
sx q[3];
rz(-1.296935) q[3];
sx q[3];
rz(3.1308208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2446186) q[2];
sx q[2];
rz(-1.3099058) q[2];
sx q[2];
rz(-0.2429602) q[2];
rz(-2.7729559) q[3];
sx q[3];
rz(-0.60550767) q[3];
sx q[3];
rz(-1.2322371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826137) q[0];
sx q[0];
rz(-0.22268) q[0];
sx q[0];
rz(-2.7742703) q[0];
rz(-0.79633725) q[1];
sx q[1];
rz(-2.0834736) q[1];
sx q[1];
rz(-0.64250362) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8272907) q[0];
sx q[0];
rz(-1.1513984) q[0];
sx q[0];
rz(0.51879306) q[0];
rz(-pi) q[1];
rz(-2.5359277) q[2];
sx q[2];
rz(-2.4993651) q[2];
sx q[2];
rz(-2.5544744) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31256235) q[1];
sx q[1];
rz(-2.1305269) q[1];
sx q[1];
rz(-1.0360495) q[1];
rz(2.0525371) q[3];
sx q[3];
rz(-1.6828487) q[3];
sx q[3];
rz(2.5555536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.502304) q[2];
sx q[2];
rz(-0.93561155) q[2];
sx q[2];
rz(1.3703692) q[2];
rz(2.9141407) q[3];
sx q[3];
rz(-1.8919614) q[3];
sx q[3];
rz(1.9500218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.3962536) q[0];
sx q[0];
rz(-0.17530137) q[0];
sx q[0];
rz(1.9434209) q[0];
rz(2.0612969) q[1];
sx q[1];
rz(-2.9273169) q[1];
sx q[1];
rz(-0.094873039) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8137275) q[0];
sx q[0];
rz(-1.3133606) q[0];
sx q[0];
rz(1.0125005) q[0];
x q[1];
rz(-2.7327161) q[2];
sx q[2];
rz(-0.95252242) q[2];
sx q[2];
rz(-0.28291075) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.97040365) q[1];
sx q[1];
rz(-0.53817777) q[1];
sx q[1];
rz(1.6506877) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8646144) q[3];
sx q[3];
rz(-1.2001032) q[3];
sx q[3];
rz(-1.0007559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0723116) q[2];
sx q[2];
rz(-0.60439622) q[2];
sx q[2];
rz(-0.4064694) q[2];
rz(1.1786002) q[3];
sx q[3];
rz(-1.4402163) q[3];
sx q[3];
rz(-1.1627722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56226319) q[0];
sx q[0];
rz(-0.13450204) q[0];
sx q[0];
rz(2.80559) q[0];
rz(2.621189) q[1];
sx q[1];
rz(-2.2774179) q[1];
sx q[1];
rz(-0.10428183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0020427) q[0];
sx q[0];
rz(-2.6536396) q[0];
sx q[0];
rz(-0.81172184) q[0];
rz(-pi) q[1];
rz(-2.9851172) q[2];
sx q[2];
rz(-1.556852) q[2];
sx q[2];
rz(-2.5283611) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6809685) q[1];
sx q[1];
rz(-1.9379341) q[1];
sx q[1];
rz(-2.2904094) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6881677) q[3];
sx q[3];
rz(-1.1201829) q[3];
sx q[3];
rz(-2.4291458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5229554) q[2];
sx q[2];
rz(-1.0982265) q[2];
sx q[2];
rz(-1.1445507) q[2];
rz(2.5942904) q[3];
sx q[3];
rz(-1.2283044) q[3];
sx q[3];
rz(1.087629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8892141) q[0];
sx q[0];
rz(-2.9504898) q[0];
sx q[0];
rz(0.55437535) q[0];
rz(0.91122183) q[1];
sx q[1];
rz(-1.2634042) q[1];
sx q[1];
rz(-0.30141452) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2269991) q[0];
sx q[0];
rz(-0.98557878) q[0];
sx q[0];
rz(0.27497681) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6658556) q[2];
sx q[2];
rz(-0.43125501) q[2];
sx q[2];
rz(0.83275821) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.44228783) q[1];
sx q[1];
rz(-1.7831452) q[1];
sx q[1];
rz(-1.9599171) q[1];
x q[2];
rz(2.2603792) q[3];
sx q[3];
rz(-2.6699319) q[3];
sx q[3];
rz(-2.415433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1201799) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(-0.0069228355) q[2];
rz(-2.210468) q[3];
sx q[3];
rz(-1.1313063) q[3];
sx q[3];
rz(2.0324223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329426) q[0];
sx q[0];
rz(-2.3045492) q[0];
sx q[0];
rz(1.6078) q[0];
rz(-0.7026698) q[1];
sx q[1];
rz(-0.53121316) q[1];
sx q[1];
rz(2.839397) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.142573) q[0];
sx q[0];
rz(-2.4732865) q[0];
sx q[0];
rz(-1.7354986) q[0];
rz(2.7415761) q[2];
sx q[2];
rz(-0.75655327) q[2];
sx q[2];
rz(0.053002593) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6249891) q[1];
sx q[1];
rz(-1.7440197) q[1];
sx q[1];
rz(0.94329639) q[1];
x q[2];
rz(-1.8994529) q[3];
sx q[3];
rz(-0.16020003) q[3];
sx q[3];
rz(0.70806187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2493784) q[2];
sx q[2];
rz(-2.0912781) q[2];
sx q[2];
rz(-0.92602473) q[2];
rz(1.2146436) q[3];
sx q[3];
rz(-0.49321431) q[3];
sx q[3];
rz(-2.8652969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4183384) q[0];
sx q[0];
rz(-2.3746018) q[0];
sx q[0];
rz(-2.0330698) q[0];
rz(-2.8077937) q[1];
sx q[1];
rz(-1.8148986) q[1];
sx q[1];
rz(1.0858067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1340116) q[0];
sx q[0];
rz(-1.3659048) q[0];
sx q[0];
rz(-1.4896859) q[0];
x q[1];
rz(-0.18098197) q[2];
sx q[2];
rz(-1.9662734) q[2];
sx q[2];
rz(-0.68812319) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1841504) q[1];
sx q[1];
rz(-2.085902) q[1];
sx q[1];
rz(0.95179291) q[1];
rz(0.9341888) q[3];
sx q[3];
rz(-1.459834) q[3];
sx q[3];
rz(-2.8257089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51599017) q[2];
sx q[2];
rz(-0.36449271) q[2];
sx q[2];
rz(-1.1642574) q[2];
rz(2.8734251) q[3];
sx q[3];
rz(-1.8987013) q[3];
sx q[3];
rz(0.75781649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5371573) q[0];
sx q[0];
rz(-0.51979655) q[0];
sx q[0];
rz(-1.1489768) q[0];
rz(2.8474498) q[1];
sx q[1];
rz(-1.5831169) q[1];
sx q[1];
rz(-2.099096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1746848) q[0];
sx q[0];
rz(-2.3174441) q[0];
sx q[0];
rz(-1.1134201) q[0];
rz(-pi) q[1];
rz(-0.22886054) q[2];
sx q[2];
rz(-2.3386152) q[2];
sx q[2];
rz(-2.2265018) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.52400541) q[1];
sx q[1];
rz(-1.397314) q[1];
sx q[1];
rz(-2.7702727) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5339666) q[3];
sx q[3];
rz(-2.4772518) q[3];
sx q[3];
rz(0.0078594154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59757549) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(1.7768804) q[2];
rz(-0.65258604) q[3];
sx q[3];
rz(-2.3502974) q[3];
sx q[3];
rz(0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6259916) q[0];
sx q[0];
rz(-2.7643272) q[0];
sx q[0];
rz(3.1242477) q[0];
rz(0.01677244) q[1];
sx q[1];
rz(-0.78712946) q[1];
sx q[1];
rz(0.15388547) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090899793) q[0];
sx q[0];
rz(-2.8150301) q[0];
sx q[0];
rz(2.9039277) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4520313) q[2];
sx q[2];
rz(-0.84367263) q[2];
sx q[2];
rz(0.52983701) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9964136) q[1];
sx q[1];
rz(-2.6657678) q[1];
sx q[1];
rz(0.68527542) q[1];
rz(-pi) q[2];
rz(-3.1176223) q[3];
sx q[3];
rz(-1.6466738) q[3];
sx q[3];
rz(-0.45437231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.953557) q[2];
sx q[2];
rz(-2.3253658) q[2];
sx q[2];
rz(0.44450644) q[2];
rz(-0.19017531) q[3];
sx q[3];
rz(-1.5580274) q[3];
sx q[3];
rz(-1.7688513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.001215) q[0];
sx q[0];
rz(-1.8920521) q[0];
sx q[0];
rz(-1.0706527) q[0];
rz(-3.095678) q[1];
sx q[1];
rz(-1.4713947) q[1];
sx q[1];
rz(-0.79107034) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96595461) q[0];
sx q[0];
rz(-1.6583867) q[0];
sx q[0];
rz(-1.9071294) q[0];
rz(-0.39761646) q[2];
sx q[2];
rz(-1.3221581) q[2];
sx q[2];
rz(1.4225117) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2998878) q[1];
sx q[1];
rz(-1.9652525) q[1];
sx q[1];
rz(1.9576555) q[1];
x q[2];
rz(-1.3338575) q[3];
sx q[3];
rz(-0.51722368) q[3];
sx q[3];
rz(0.57843971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.644824) q[2];
sx q[2];
rz(-0.18459979) q[2];
sx q[2];
rz(-1.6688639) q[2];
rz(1.5276927) q[3];
sx q[3];
rz(-2.0852641) q[3];
sx q[3];
rz(-1.9014026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.938217) q[0];
sx q[0];
rz(-1.4980409) q[0];
sx q[0];
rz(-0.71312755) q[0];
rz(-0.51364246) q[1];
sx q[1];
rz(-1.4331663) q[1];
sx q[1];
rz(-1.6930361) q[1];
rz(2.7265139) q[2];
sx q[2];
rz(-0.54176258) q[2];
sx q[2];
rz(-1.1272346) q[2];
rz(-2.33251) q[3];
sx q[3];
rz(-1.2590903) q[3];
sx q[3];
rz(2.1303582) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
