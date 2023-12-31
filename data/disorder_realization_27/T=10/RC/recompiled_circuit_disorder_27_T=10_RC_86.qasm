OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.0032089631) q[0];
sx q[0];
rz(-0.15455833) q[0];
sx q[0];
rz(0.69252339) q[0];
rz(5.0737557) q[1];
sx q[1];
rz(4.3901246) q[1];
sx q[1];
rz(7.6683383) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6801075) q[0];
sx q[0];
rz(-2.2930817) q[0];
sx q[0];
rz(1.1397584) q[0];
x q[1];
rz(-1.295747) q[2];
sx q[2];
rz(-1.0359456) q[2];
sx q[2];
rz(1.5616852) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79917158) q[1];
sx q[1];
rz(-2.85596) q[1];
sx q[1];
rz(0.61113961) q[1];
rz(-0.65269835) q[3];
sx q[3];
rz(-1.1678809) q[3];
sx q[3];
rz(2.9626915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8866855) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(-2.936426) q[2];
rz(0.77130476) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(2.0390959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40760621) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(-0.45390391) q[0];
rz(-1.0247963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(-1.227238) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1821546) q[0];
sx q[0];
rz(-1.4417138) q[0];
sx q[0];
rz(3.0653619) q[0];
rz(-pi) q[1];
x q[1];
rz(1.111755) q[2];
sx q[2];
rz(-2.0437721) q[2];
sx q[2];
rz(2.3943261) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.40185803) q[1];
sx q[1];
rz(-1.2793853) q[1];
sx q[1];
rz(1.3009562) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4144054) q[3];
sx q[3];
rz(-1.7692493) q[3];
sx q[3];
rz(-0.95201492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1002905) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(2.5741637) q[2];
rz(-0.36519095) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.658618) q[0];
sx q[0];
rz(-2.5768319) q[0];
sx q[0];
rz(2.2429402) q[0];
rz(2.1458416) q[1];
sx q[1];
rz(-1.5834705) q[1];
sx q[1];
rz(-2.8083037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.916327) q[0];
sx q[0];
rz(-0.81575459) q[0];
sx q[0];
rz(2.4832721) q[0];
x q[1];
rz(2.0731508) q[2];
sx q[2];
rz(-0.2012673) q[2];
sx q[2];
rz(-1.7114491) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4132727) q[1];
sx q[1];
rz(-1.4188758) q[1];
sx q[1];
rz(-0.9539714) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7379052) q[3];
sx q[3];
rz(-1.0567769) q[3];
sx q[3];
rz(-0.46883632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4553392) q[2];
sx q[2];
rz(-1.3505961) q[2];
sx q[2];
rz(-1.9906445) q[2];
rz(2.3006556) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(1.9870728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.44160098) q[0];
sx q[0];
rz(-1.561152) q[0];
sx q[0];
rz(0.68471318) q[0];
rz(-1.0355863) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(1.9365786) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9158463) q[0];
sx q[0];
rz(-1.5374743) q[0];
sx q[0];
rz(-0.87945625) q[0];
x q[1];
rz(2.2150546) q[2];
sx q[2];
rz(-1.4069923) q[2];
sx q[2];
rz(2.0870199) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64441427) q[1];
sx q[1];
rz(-0.89343151) q[1];
sx q[1];
rz(-2.7888984) q[1];
rz(1.4865925) q[3];
sx q[3];
rz(-2.4349672) q[3];
sx q[3];
rz(0.018761793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8923607) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(-0.37115804) q[2];
rz(-1.7403729) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(-2.0223117) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.086833) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(3.0084685) q[0];
rz(-2.1482824) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(0.55508074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18300444) q[0];
sx q[0];
rz(-2.8259813) q[0];
sx q[0];
rz(-1.6118227) q[0];
rz(-0.06818469) q[2];
sx q[2];
rz(-1.9848595) q[2];
sx q[2];
rz(-2.806864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3913369) q[1];
sx q[1];
rz(-1.5899982) q[1];
sx q[1];
rz(-2.0160497) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0606907) q[3];
sx q[3];
rz(-1.6228075) q[3];
sx q[3];
rz(-2.4018283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30620265) q[2];
sx q[2];
rz(-2.1357048) q[2];
sx q[2];
rz(-3.0026657) q[2];
rz(0.94240087) q[3];
sx q[3];
rz(-1.5709632) q[3];
sx q[3];
rz(-2.5901103) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5979364) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(-0.57975769) q[0];
rz(-0.12750164) q[1];
sx q[1];
rz(-1.189905) q[1];
sx q[1];
rz(-1.5396083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67028763) q[0];
sx q[0];
rz(-1.3677214) q[0];
sx q[0];
rz(0.85104403) q[0];
x q[1];
rz(2.8069205) q[2];
sx q[2];
rz(-0.13813189) q[2];
sx q[2];
rz(-1.7931256) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3547937) q[1];
sx q[1];
rz(-1.8537632) q[1];
sx q[1];
rz(-1.1250886) q[1];
rz(-pi) q[2];
rz(1.2195915) q[3];
sx q[3];
rz(-0.71577365) q[3];
sx q[3];
rz(-2.2367665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5876028) q[2];
sx q[2];
rz(-2.8911399) q[2];
sx q[2];
rz(-0.26947752) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(3.0814734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6948029) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(-2.457298) q[0];
rz(3.0220095) q[1];
sx q[1];
rz(-1.8493098) q[1];
sx q[1];
rz(0.51876846) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25012384) q[0];
sx q[0];
rz(-2.2995673) q[0];
sx q[0];
rz(2.9617589) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3709626) q[2];
sx q[2];
rz(-1.7297041) q[2];
sx q[2];
rz(-2.4385914) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8956633) q[1];
sx q[1];
rz(-2.3024125) q[1];
sx q[1];
rz(-2.500446) q[1];
x q[2];
rz(-1.5090844) q[3];
sx q[3];
rz(-2.7402174) q[3];
sx q[3];
rz(-2.9369831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8873022) q[2];
sx q[2];
rz(-1.3672978) q[2];
sx q[2];
rz(0.56345144) q[2];
rz(-3.0900132) q[3];
sx q[3];
rz(-1.1392461) q[3];
sx q[3];
rz(-2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96034399) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(-1.7425591) q[0];
rz(2.3545806) q[1];
sx q[1];
rz(-2.0128638) q[1];
sx q[1];
rz(-2.3972437) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58957419) q[0];
sx q[0];
rz(-2.9840143) q[0];
sx q[0];
rz(-2.7872859) q[0];
rz(-1.3049576) q[2];
sx q[2];
rz(-2.611534) q[2];
sx q[2];
rz(-2.838138) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5163755) q[1];
sx q[1];
rz(-2.5853734) q[1];
sx q[1];
rz(2.7079627) q[1];
x q[2];
rz(2.6870071) q[3];
sx q[3];
rz(-1.5919519) q[3];
sx q[3];
rz(1.3190312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4259592) q[2];
sx q[2];
rz(-1.3954433) q[2];
sx q[2];
rz(-2.5320833) q[2];
rz(-0.65731796) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(-2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1881926) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(1.6375861) q[0];
rz(1.2414744) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(-2.3666568) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06668815) q[0];
sx q[0];
rz(-0.79830352) q[0];
sx q[0];
rz(-1.0134646) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8911792) q[2];
sx q[2];
rz(-0.46386007) q[2];
sx q[2];
rz(-1.8360209) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1800268) q[1];
sx q[1];
rz(-2.0225836) q[1];
sx q[1];
rz(-1.3614484) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9990066) q[3];
sx q[3];
rz(-1.3538059) q[3];
sx q[3];
rz(-0.91689527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9541786) q[2];
sx q[2];
rz(-0.22970197) q[2];
sx q[2];
rz(-2.9837218) q[2];
rz(1.9291417) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(-1.8635748) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.071844) q[0];
sx q[0];
rz(-0.97244111) q[0];
sx q[0];
rz(-0.21433314) q[0];
rz(-2.4841323) q[1];
sx q[1];
rz(-0.22413707) q[1];
sx q[1];
rz(-1.0459895) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8911274) q[0];
sx q[0];
rz(-1.9292826) q[0];
sx q[0];
rz(1.6842708) q[0];
rz(1.5418391) q[2];
sx q[2];
rz(-2.1615897) q[2];
sx q[2];
rz(2.7383885) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5221094) q[1];
sx q[1];
rz(-1.2088747) q[1];
sx q[1];
rz(-0.98720179) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5973741) q[3];
sx q[3];
rz(-1.1655032) q[3];
sx q[3];
rz(-1.6534896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41632286) q[2];
sx q[2];
rz(-0.060083397) q[2];
sx q[2];
rz(2.0521169) q[2];
rz(1.5754835) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(0.65264788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8626704) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(-1.6090341) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(-1.8794466) q[2];
sx q[2];
rz(-1.1790922) q[2];
sx q[2];
rz(-2.14711) q[2];
rz(-3.0388721) q[3];
sx q[3];
rz(-0.552388) q[3];
sx q[3];
rz(-1.296476) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
