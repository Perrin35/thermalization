OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(-0.78200114) q[0];
sx q[0];
rz(1.8703823) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9936375) q[0];
sx q[0];
rz(-1.7339098) q[0];
sx q[0];
rz(-1.9440584) q[0];
x q[1];
rz(-1.5288562) q[2];
sx q[2];
rz(-1.1239927) q[2];
sx q[2];
rz(2.1697793) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.672294) q[1];
sx q[1];
rz(-1.0743595) q[1];
sx q[1];
rz(-1.4908355) q[1];
rz(-pi) q[2];
rz(1.9418342) q[3];
sx q[3];
rz(-1.8098117) q[3];
sx q[3];
rz(-2.1109076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15443054) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(0.74938613) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(-0.40482503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4352903) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(0.97066561) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(2.326139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1708508) q[0];
sx q[0];
rz(-2.4903957) q[0];
sx q[0];
rz(-0.09911508) q[0];
x q[1];
rz(2.8380727) q[2];
sx q[2];
rz(-2.3887861) q[2];
sx q[2];
rz(-2.7941861) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.26317877) q[1];
sx q[1];
rz(-2.8385332) q[1];
sx q[1];
rz(-0.11942272) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29173298) q[3];
sx q[3];
rz(-0.66666224) q[3];
sx q[3];
rz(-0.093689703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4619535) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(-1.1535545) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(-0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84045029) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.7069838) q[1];
sx q[1];
rz(-2.1420746) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32854983) q[0];
sx q[0];
rz(-1.8694436) q[0];
sx q[0];
rz(2.9850328) q[0];
x q[1];
rz(-2.524316) q[2];
sx q[2];
rz(-2.1356138) q[2];
sx q[2];
rz(1.9895983) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42400186) q[1];
sx q[1];
rz(-1.0832936) q[1];
sx q[1];
rz(-1.2815777) q[1];
rz(-pi) q[2];
rz(-1.7129094) q[3];
sx q[3];
rz(-0.51321533) q[3];
sx q[3];
rz(1.3877447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6040566) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(-0.58004722) q[2];
rz(-0.81702685) q[3];
sx q[3];
rz(-1.3823119) q[3];
sx q[3];
rz(1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.375305) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(2.2312009) q[0];
rz(-2.6903649) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(-0.27483637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3954454) q[0];
sx q[0];
rz(-1.6593885) q[0];
sx q[0];
rz(-3.0625312) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.574013) q[2];
sx q[2];
rz(-0.46041691) q[2];
sx q[2];
rz(-2.761063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.30157629) q[1];
sx q[1];
rz(-0.77743545) q[1];
sx q[1];
rz(0.90228723) q[1];
rz(1.4196017) q[3];
sx q[3];
rz(-2.4393743) q[3];
sx q[3];
rz(0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3466907) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(1.5083195) q[2];
rz(-1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(-2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4836924) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(2.143798) q[0];
rz(2.9580341) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(1.6246187) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51380101) q[0];
sx q[0];
rz(-0.9538981) q[0];
sx q[0];
rz(-0.4517201) q[0];
x q[1];
rz(0.60913182) q[2];
sx q[2];
rz(-2.5172148) q[2];
sx q[2];
rz(-0.49027157) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6301873) q[1];
sx q[1];
rz(-1.1668219) q[1];
sx q[1];
rz(-1.1955111) q[1];
rz(-pi) q[2];
x q[2];
rz(2.153271) q[3];
sx q[3];
rz(-2.0867996) q[3];
sx q[3];
rz(-0.9048681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8020442) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(2.8175763) q[2];
rz(-1.8185395) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(1.2639686) q[0];
rz(-0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(2.8009159) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37939385) q[0];
sx q[0];
rz(-1.6025935) q[0];
sx q[0];
rz(-1.6367153) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8315115) q[2];
sx q[2];
rz(-0.78044621) q[2];
sx q[2];
rz(2.9030637) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1250455) q[1];
sx q[1];
rz(-0.53495896) q[1];
sx q[1];
rz(-0.46334456) q[1];
x q[2];
rz(-0.86990279) q[3];
sx q[3];
rz(-2.5330336) q[3];
sx q[3];
rz(-0.22939798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.95057758) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(2.3699956) q[2];
rz(-2.5937882) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724378) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(0.34564885) q[0];
rz(-0.06282839) q[1];
sx q[1];
rz(-0.47880104) q[1];
sx q[1];
rz(-0.46494928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3195254) q[0];
sx q[0];
rz(-1.7626581) q[0];
sx q[0];
rz(-1.096154) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99636997) q[2];
sx q[2];
rz(-1.1466221) q[2];
sx q[2];
rz(3.0976354) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17864922) q[1];
sx q[1];
rz(-1.6766251) q[1];
sx q[1];
rz(0.21957285) q[1];
rz(1.5246478) q[3];
sx q[3];
rz(-1.5449636) q[3];
sx q[3];
rz(-1.8125364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1376301) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(-2.8239992) q[2];
rz(-2.5701304) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6222318) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(-2.8572594) q[0];
rz(-0.55150664) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(-0.078358738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6358444) q[0];
sx q[0];
rz(-2.3730179) q[0];
sx q[0];
rz(1.4710674) q[0];
rz(0.80348357) q[2];
sx q[2];
rz(-2.8065971) q[2];
sx q[2];
rz(1.4575046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2625418) q[1];
sx q[1];
rz(-1.3971551) q[1];
sx q[1];
rz(0.91211984) q[1];
rz(-0.034268495) q[3];
sx q[3];
rz(-1.4141603) q[3];
sx q[3];
rz(-0.53965118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4006965) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(-2.890214) q[2];
rz(2.5583983) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3437929) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(0.051368512) q[0];
rz(2.2180166) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(2.267568) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6459991) q[0];
sx q[0];
rz(-0.80276239) q[0];
sx q[0];
rz(-1.5587224) q[0];
rz(-pi) q[1];
rz(-1.1633515) q[2];
sx q[2];
rz(-1.2375087) q[2];
sx q[2];
rz(-2.8826706) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8928788) q[1];
sx q[1];
rz(-1.6422179) q[1];
sx q[1];
rz(2.0715269) q[1];
x q[2];
rz(0.60871082) q[3];
sx q[3];
rz(-1.8551644) q[3];
sx q[3];
rz(-1.1876719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.727227) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(2.5218463) q[2];
rz(-1.9571346) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(-1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0062362) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(0.7243048) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.9627409) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26319474) q[0];
sx q[0];
rz(-0.34391719) q[0];
sx q[0];
rz(-2.2180936) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12702282) q[2];
sx q[2];
rz(-1.6633031) q[2];
sx q[2];
rz(-0.18735838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2420173) q[1];
sx q[1];
rz(-1.3841108) q[1];
sx q[1];
rz(-0.94355299) q[1];
rz(-2.6188649) q[3];
sx q[3];
rz(-1.7676815) q[3];
sx q[3];
rz(-2.6678391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.05802352) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(-0.89938346) q[2];
rz(-2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.6806867) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(-0.75795603) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(-1.1307217) q[2];
sx q[2];
rz(-1.6025087) q[2];
sx q[2];
rz(-0.58379731) q[2];
rz(1.422613) q[3];
sx q[3];
rz(-1.2978745) q[3];
sx q[3];
rz(0.10980448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
