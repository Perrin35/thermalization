OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2744098) q[0];
sx q[0];
rz(-0.46719587) q[0];
sx q[0];
rz(2.245477) q[0];
rz(-0.81796247) q[1];
sx q[1];
rz(-1.5476513) q[1];
sx q[1];
rz(-0.98734468) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2786699) q[0];
sx q[0];
rz(-0.76919829) q[0];
sx q[0];
rz(1.0095566) q[0];
x q[1];
rz(0.9573663) q[2];
sx q[2];
rz(-0.84954689) q[2];
sx q[2];
rz(-2.3867875) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8670013) q[1];
sx q[1];
rz(-1.7385635) q[1];
sx q[1];
rz(-2.861119) q[1];
rz(0.36611661) q[3];
sx q[3];
rz(-1.7185655) q[3];
sx q[3];
rz(0.63945668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2200615) q[2];
sx q[2];
rz(-2.107373) q[2];
sx q[2];
rz(-2.6424778) q[2];
rz(-0.81579298) q[3];
sx q[3];
rz(-0.59320265) q[3];
sx q[3];
rz(0.35713404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4895184) q[0];
sx q[0];
rz(-2.7777785) q[0];
sx q[0];
rz(2.9079085) q[0];
rz(0.09985996) q[1];
sx q[1];
rz(-1.4870817) q[1];
sx q[1];
rz(1.5335836) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8692538) q[0];
sx q[0];
rz(-1.6923447) q[0];
sx q[0];
rz(0.097313332) q[0];
rz(0.38012114) q[2];
sx q[2];
rz(-0.46762662) q[2];
sx q[2];
rz(1.1884226) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2054384) q[1];
sx q[1];
rz(-2.1742651) q[1];
sx q[1];
rz(0.11118576) q[1];
rz(2.8321974) q[3];
sx q[3];
rz(-1.1077048) q[3];
sx q[3];
rz(2.9675296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0443772) q[2];
sx q[2];
rz(-0.88572398) q[2];
sx q[2];
rz(0.054917939) q[2];
rz(2.4954259) q[3];
sx q[3];
rz(-1.5271527) q[3];
sx q[3];
rz(3.0687148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1845448) q[0];
sx q[0];
rz(-0.2453198) q[0];
sx q[0];
rz(0.74485892) q[0];
rz(1.2767731) q[1];
sx q[1];
rz(-1.0294015) q[1];
sx q[1];
rz(-0.863711) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0833763) q[0];
sx q[0];
rz(-1.6291219) q[0];
sx q[0];
rz(-2.7013679) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3529007) q[2];
sx q[2];
rz(-1.1967107) q[2];
sx q[2];
rz(-1.7497334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7924713) q[1];
sx q[1];
rz(-2.0466261) q[1];
sx q[1];
rz(2.4351067) q[1];
rz(1.0705804) q[3];
sx q[3];
rz(-0.97220927) q[3];
sx q[3];
rz(-0.14414302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.471571) q[2];
sx q[2];
rz(-1.5315285) q[2];
sx q[2];
rz(-0.395533) q[2];
rz(-2.1192571) q[3];
sx q[3];
rz(-2.6774355) q[3];
sx q[3];
rz(-0.49183229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96524298) q[0];
sx q[0];
rz(-1.1818161) q[0];
sx q[0];
rz(3.0856207) q[0];
rz(-0.61198992) q[1];
sx q[1];
rz(-1.5490218) q[1];
sx q[1];
rz(0.8459808) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056123646) q[0];
sx q[0];
rz(-0.58857354) q[0];
sx q[0];
rz(1.0928297) q[0];
rz(-2.0417287) q[2];
sx q[2];
rz(-0.67457565) q[2];
sx q[2];
rz(1.2686307) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27056672) q[1];
sx q[1];
rz(-1.2602087) q[1];
sx q[1];
rz(1.7720211) q[1];
x q[2];
rz(-0.9153644) q[3];
sx q[3];
rz(-2.5662321) q[3];
sx q[3];
rz(-1.2498472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0116288) q[2];
sx q[2];
rz(-1.6946946) q[2];
sx q[2];
rz(-0.72032991) q[2];
rz(2.7885041) q[3];
sx q[3];
rz(-1.1938286) q[3];
sx q[3];
rz(0.24875719) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40097749) q[0];
sx q[0];
rz(-1.688513) q[0];
sx q[0];
rz(0.98689669) q[0];
rz(1.1445716) q[1];
sx q[1];
rz(-2.3026376) q[1];
sx q[1];
rz(2.1867337) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83897774) q[0];
sx q[0];
rz(-0.90388008) q[0];
sx q[0];
rz(0.033292183) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32771941) q[2];
sx q[2];
rz(-2.1075038) q[2];
sx q[2];
rz(-0.84568727) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.1249117) q[1];
sx q[1];
rz(-2.0215258) q[1];
sx q[1];
rz(-1.1523516) q[1];
rz(-pi) q[2];
rz(-1.0480065) q[3];
sx q[3];
rz(-2.8510333) q[3];
sx q[3];
rz(2.235242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6285051) q[2];
sx q[2];
rz(-0.91431618) q[2];
sx q[2];
rz(-2.0950192) q[2];
rz(-0.70932499) q[3];
sx q[3];
rz(-1.1449292) q[3];
sx q[3];
rz(2.5078702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3777622) q[0];
sx q[0];
rz(-1.4077633) q[0];
sx q[0];
rz(-2.8705226) q[0];
rz(3.0625878) q[1];
sx q[1];
rz(-0.43402356) q[1];
sx q[1];
rz(-1.7105506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9535429) q[0];
sx q[0];
rz(-2.8859647) q[0];
sx q[0];
rz(0.60582692) q[0];
rz(-2.1151401) q[2];
sx q[2];
rz(-1.0751131) q[2];
sx q[2];
rz(-3.0718671) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.62785) q[1];
sx q[1];
rz(-1.3249036) q[1];
sx q[1];
rz(-0.69468433) q[1];
x q[2];
rz(-1.1523139) q[3];
sx q[3];
rz(-1.7177204) q[3];
sx q[3];
rz(-0.72673405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0926823) q[2];
sx q[2];
rz(-1.9581257) q[2];
sx q[2];
rz(-3.1177706) q[2];
rz(-1.9134936) q[3];
sx q[3];
rz(-2.7619599) q[3];
sx q[3];
rz(2.9983799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.989885) q[0];
sx q[0];
rz(-2.8360974) q[0];
sx q[0];
rz(-3.049343) q[0];
rz(-0.30091885) q[1];
sx q[1];
rz(-1.9124799) q[1];
sx q[1];
rz(2.0613964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5841519) q[0];
sx q[0];
rz(-1.6655465) q[0];
sx q[0];
rz(1.2242203) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2713234) q[2];
sx q[2];
rz(-2.610321) q[2];
sx q[2];
rz(-0.16933717) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34970666) q[1];
sx q[1];
rz(-2.9936446) q[1];
sx q[1];
rz(-0.011758864) q[1];
x q[2];
rz(1.2545414) q[3];
sx q[3];
rz(-2.2407534) q[3];
sx q[3];
rz(0.073410598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0098308) q[2];
sx q[2];
rz(-2.7907382) q[2];
sx q[2];
rz(0.62823137) q[2];
rz(-2.175711) q[3];
sx q[3];
rz(-1.5647669) q[3];
sx q[3];
rz(-0.32111827) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0308762) q[0];
sx q[0];
rz(-0.70699152) q[0];
sx q[0];
rz(1.734717) q[0];
rz(0.12786099) q[1];
sx q[1];
rz(-1.7117932) q[1];
sx q[1];
rz(1.1258639) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28982535) q[0];
sx q[0];
rz(-0.85677108) q[0];
sx q[0];
rz(-1.8147574) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0225384) q[2];
sx q[2];
rz(-1.5875419) q[2];
sx q[2];
rz(-0.43338768) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.606877) q[1];
sx q[1];
rz(-1.770853) q[1];
sx q[1];
rz(1.6354531) q[1];
x q[2];
rz(-0.61922686) q[3];
sx q[3];
rz(-1.2017177) q[3];
sx q[3];
rz(-2.823373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.13059482) q[2];
sx q[2];
rz(-1.1329634) q[2];
sx q[2];
rz(-3.0478743) q[2];
rz(-1.4334009) q[3];
sx q[3];
rz(-2.8580557) q[3];
sx q[3];
rz(-1.7482429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4505287) q[0];
sx q[0];
rz(-2.0557623) q[0];
sx q[0];
rz(2.1375256) q[0];
rz(1.1035236) q[1];
sx q[1];
rz(-2.6595778) q[1];
sx q[1];
rz(-1.7335266) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226949) q[0];
sx q[0];
rz(-1.0556882) q[0];
sx q[0];
rz(1.5777508) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33317487) q[2];
sx q[2];
rz(-0.3292203) q[2];
sx q[2];
rz(1.7235989) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6390761) q[1];
sx q[1];
rz(-0.47331077) q[1];
sx q[1];
rz(-1.972354) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6203815) q[3];
sx q[3];
rz(-2.0892429) q[3];
sx q[3];
rz(-2.848846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7496926) q[2];
sx q[2];
rz(-1.5211279) q[2];
sx q[2];
rz(0.54769546) q[2];
rz(-0.31217602) q[3];
sx q[3];
rz(-1.4145989) q[3];
sx q[3];
rz(1.4148022) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80426973) q[0];
sx q[0];
rz(-1.4708568) q[0];
sx q[0];
rz(-0.7789337) q[0];
rz(0.3745105) q[1];
sx q[1];
rz(-1.51314) q[1];
sx q[1];
rz(-2.3096854) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4171281) q[0];
sx q[0];
rz(-2.3848371) q[0];
sx q[0];
rz(-0.3860851) q[0];
rz(-pi) q[1];
rz(0.7124165) q[2];
sx q[2];
rz(-2.2210026) q[2];
sx q[2];
rz(-0.67931108) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9495957) q[1];
sx q[1];
rz(-2.1131385) q[1];
sx q[1];
rz(-0.80710141) q[1];
x q[2];
rz(1.5396773) q[3];
sx q[3];
rz(-0.77431576) q[3];
sx q[3];
rz(1.4207135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1067918) q[2];
sx q[2];
rz(-0.35934862) q[2];
sx q[2];
rz(0.0084776004) q[2];
rz(0.40555412) q[3];
sx q[3];
rz(-1.4529934) q[3];
sx q[3];
rz(1.4439553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99210284) q[0];
sx q[0];
rz(-1.6076417) q[0];
sx q[0];
rz(-2.3544307) q[0];
rz(-2.0472732) q[1];
sx q[1];
rz(-2.7631187) q[1];
sx q[1];
rz(-0.43597058) q[1];
rz(2.1020066) q[2];
sx q[2];
rz(-2.8969904) q[2];
sx q[2];
rz(1.9319509) q[2];
rz(0.094219128) q[3];
sx q[3];
rz(-1.5658656) q[3];
sx q[3];
rz(1.6938536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
