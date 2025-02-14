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
rz(0.02778223) q[0];
sx q[0];
rz(-1.2454998) q[0];
sx q[0];
rz(1.0863289) q[0];
rz(-1.0983941) q[1];
sx q[1];
rz(-1.2135222) q[1];
sx q[1];
rz(-1.9953802) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6667696) q[0];
sx q[0];
rz(-2.7753029) q[0];
sx q[0];
rz(-1.1046871) q[0];
x q[1];
rz(0.5733344) q[2];
sx q[2];
rz(-2.4681915) q[2];
sx q[2];
rz(1.2902989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.024373) q[1];
sx q[1];
rz(-2.0659323) q[1];
sx q[1];
rz(-0.58039033) q[1];
rz(-pi) q[2];
rz(2.9153053) q[3];
sx q[3];
rz(-2.0620966) q[3];
sx q[3];
rz(1.1102344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.22831297) q[2];
sx q[2];
rz(-1.5659119) q[2];
sx q[2];
rz(2.7794465) q[2];
rz(-0.95237887) q[3];
sx q[3];
rz(-0.55914545) q[3];
sx q[3];
rz(2.4366116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.1700965) q[0];
sx q[0];
rz(-2.9269452) q[0];
sx q[0];
rz(-1.5317408) q[0];
rz(1.2486628) q[1];
sx q[1];
rz(-1.5208533) q[1];
sx q[1];
rz(-0.85266399) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.296761) q[0];
sx q[0];
rz(-1.0687901) q[0];
sx q[0];
rz(0.39433483) q[0];
rz(1.1617848) q[2];
sx q[2];
rz(-0.30050052) q[2];
sx q[2];
rz(-0.40846339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.98140796) q[1];
sx q[1];
rz(-1.7716737) q[1];
sx q[1];
rz(2.1698879) q[1];
x q[2];
rz(0.38532704) q[3];
sx q[3];
rz(-1.0443976) q[3];
sx q[3];
rz(-0.86643636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2109005) q[2];
sx q[2];
rz(-0.96052581) q[2];
sx q[2];
rz(-1.1962013) q[2];
rz(3.1148552) q[3];
sx q[3];
rz(-0.49632159) q[3];
sx q[3];
rz(-1.668821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6579599) q[0];
sx q[0];
rz(-2.8909029) q[0];
sx q[0];
rz(-2.2657917) q[0];
rz(-0.35975131) q[1];
sx q[1];
rz(-1.3641554) q[1];
sx q[1];
rz(-1.0035427) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26311603) q[0];
sx q[0];
rz(-2.5777393) q[0];
sx q[0];
rz(0.13335689) q[0];
x q[1];
rz(2.7458756) q[2];
sx q[2];
rz(-2.1747602) q[2];
sx q[2];
rz(-0.25091083) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4565125) q[1];
sx q[1];
rz(-0.69275038) q[1];
sx q[1];
rz(0.22489203) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67387132) q[3];
sx q[3];
rz(-2.0099075) q[3];
sx q[3];
rz(1.278745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2063736) q[2];
sx q[2];
rz(-0.71722764) q[2];
sx q[2];
rz(-0.75500542) q[2];
rz(-0.098091789) q[3];
sx q[3];
rz(-2.5710929) q[3];
sx q[3];
rz(-0.3977972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0841242) q[0];
sx q[0];
rz(-1.6300772) q[0];
sx q[0];
rz(0.64495069) q[0];
rz(-2.0461252) q[1];
sx q[1];
rz(-2.7332833) q[1];
sx q[1];
rz(-1.1114978) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31331691) q[0];
sx q[0];
rz(-3.0544825) q[0];
sx q[0];
rz(0.75096186) q[0];
rz(-1.6531282) q[2];
sx q[2];
rz(-1.4042743) q[2];
sx q[2];
rz(0.67556591) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.20852986) q[1];
sx q[1];
rz(-0.25756076) q[1];
sx q[1];
rz(-2.527586) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74097775) q[3];
sx q[3];
rz(-1.7347226) q[3];
sx q[3];
rz(1.0865097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5081386) q[2];
sx q[2];
rz(-2.3544669) q[2];
sx q[2];
rz(0.48466361) q[2];
rz(0.44114068) q[3];
sx q[3];
rz(-1.7287247) q[3];
sx q[3];
rz(-1.2539366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50705528) q[0];
sx q[0];
rz(-0.84369722) q[0];
sx q[0];
rz(1.4265358) q[0];
rz(-0.85975319) q[1];
sx q[1];
rz(-0.25329241) q[1];
sx q[1];
rz(-0.8304798) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7380421) q[0];
sx q[0];
rz(-2.134424) q[0];
sx q[0];
rz(0.72569816) q[0];
rz(-1.6489588) q[2];
sx q[2];
rz(-2.507308) q[2];
sx q[2];
rz(1.2595578) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1089835) q[1];
sx q[1];
rz(-0.84798893) q[1];
sx q[1];
rz(-1.7574786) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36019616) q[3];
sx q[3];
rz(-1.7203698) q[3];
sx q[3];
rz(-2.5392578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44437235) q[2];
sx q[2];
rz(-1.0151007) q[2];
sx q[2];
rz(2.3946136) q[2];
rz(2.8367481) q[3];
sx q[3];
rz(-1.3957142) q[3];
sx q[3];
rz(-1.9029118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9656669) q[0];
sx q[0];
rz(-1.3015863) q[0];
sx q[0];
rz(2.9081705) q[0];
rz(2.6790791) q[1];
sx q[1];
rz(-1.9513444) q[1];
sx q[1];
rz(-2.7375284) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3041026) q[0];
sx q[0];
rz(-1.7902052) q[0];
sx q[0];
rz(-0.93774937) q[0];
rz(2.7163804) q[2];
sx q[2];
rz(-1.325144) q[2];
sx q[2];
rz(2.6157891) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7884454) q[1];
sx q[1];
rz(-0.57203469) q[1];
sx q[1];
rz(2.9414909) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3498106) q[3];
sx q[3];
rz(-0.95181247) q[3];
sx q[3];
rz(-0.19863752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.81249753) q[2];
sx q[2];
rz(-1.7915244) q[2];
sx q[2];
rz(-0.65328944) q[2];
rz(1.6434044) q[3];
sx q[3];
rz(-2.523962) q[3];
sx q[3];
rz(2.5352855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46635982) q[0];
sx q[0];
rz(-1.9283858) q[0];
sx q[0];
rz(2.8208222) q[0];
rz(-0.65387154) q[1];
sx q[1];
rz(-2.8120698) q[1];
sx q[1];
rz(3.0810862) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0573567) q[0];
sx q[0];
rz(-2.8222146) q[0];
sx q[0];
rz(2.1288298) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67530151) q[2];
sx q[2];
rz(-2.2847698) q[2];
sx q[2];
rz(1.7973822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.49566073) q[1];
sx q[1];
rz(-2.0510608) q[1];
sx q[1];
rz(-0.86841327) q[1];
x q[2];
rz(1.1641422) q[3];
sx q[3];
rz(-1.2429626) q[3];
sx q[3];
rz(-1.145112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74762809) q[2];
sx q[2];
rz(-1.4614033) q[2];
sx q[2];
rz(2.4192269) q[2];
rz(-2.8584621) q[3];
sx q[3];
rz(-2.2524998) q[3];
sx q[3];
rz(-0.50962555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79001456) q[0];
sx q[0];
rz(-0.55792648) q[0];
sx q[0];
rz(0.18490069) q[0];
rz(0.50390538) q[1];
sx q[1];
rz(-1.3495219) q[1];
sx q[1];
rz(-0.46461836) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5123925) q[0];
sx q[0];
rz(-0.2357395) q[0];
sx q[0];
rz(0.78099536) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9275168) q[2];
sx q[2];
rz(-0.42595562) q[2];
sx q[2];
rz(1.730012) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1211627) q[1];
sx q[1];
rz(-1.1926225) q[1];
sx q[1];
rz(-0.44831033) q[1];
rz(1.7691602) q[3];
sx q[3];
rz(-0.86097211) q[3];
sx q[3];
rz(-1.3876624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9087312) q[2];
sx q[2];
rz(-1.6249012) q[2];
sx q[2];
rz(1.0535343) q[2];
rz(-2.4145224) q[3];
sx q[3];
rz(-1.9126242) q[3];
sx q[3];
rz(0.54522902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19186774) q[0];
sx q[0];
rz(-1.1911012) q[0];
sx q[0];
rz(-0.18549347) q[0];
rz(2.5718451) q[1];
sx q[1];
rz(-0.92420095) q[1];
sx q[1];
rz(-1.9843019) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7737024) q[0];
sx q[0];
rz(-1.4285402) q[0];
sx q[0];
rz(-0.41928798) q[0];
rz(0.050008372) q[2];
sx q[2];
rz(-0.93863025) q[2];
sx q[2];
rz(-0.061246733) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4450106) q[1];
sx q[1];
rz(-2.5392476) q[1];
sx q[1];
rz(2.7578043) q[1];
x q[2];
rz(-2.3367711) q[3];
sx q[3];
rz(-1.8880883) q[3];
sx q[3];
rz(-2.1498093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40413228) q[2];
sx q[2];
rz(-0.22455939) q[2];
sx q[2];
rz(-1.4004716) q[2];
rz(0.36462668) q[3];
sx q[3];
rz(-0.76321634) q[3];
sx q[3];
rz(2.3188685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6309361) q[0];
sx q[0];
rz(-0.83844227) q[0];
sx q[0];
rz(-0.034570845) q[0];
rz(2.8055387) q[1];
sx q[1];
rz(-2.2269109) q[1];
sx q[1];
rz(2.5539982) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66389342) q[0];
sx q[0];
rz(-1.4283533) q[0];
sx q[0];
rz(-0.42147984) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4108836) q[2];
sx q[2];
rz(-1.8160422) q[2];
sx q[2];
rz(2.2718549) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61240367) q[1];
sx q[1];
rz(-1.8029787) q[1];
sx q[1];
rz(-0.043223765) q[1];
rz(-2.5925007) q[3];
sx q[3];
rz(-0.7585511) q[3];
sx q[3];
rz(2.9545934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33362886) q[2];
sx q[2];
rz(-0.89666349) q[2];
sx q[2];
rz(2.4133852) q[2];
rz(2.7306469) q[3];
sx q[3];
rz(-0.22308895) q[3];
sx q[3];
rz(-1.8691011) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67644607) q[0];
sx q[0];
rz(-1.5855753) q[0];
sx q[0];
rz(1.3517071) q[0];
rz(-2.7017055) q[1];
sx q[1];
rz(-0.37275795) q[1];
sx q[1];
rz(2.8511924) q[1];
rz(1.2705304) q[2];
sx q[2];
rz(-1.9004442) q[2];
sx q[2];
rz(1.6729542) q[2];
rz(1.5241809) q[3];
sx q[3];
rz(-1.6773994) q[3];
sx q[3];
rz(1.027153) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
