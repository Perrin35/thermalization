OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7423695) q[0];
sx q[0];
rz(2.7288781) q[0];
sx q[0];
rz(5.4469845) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(6.3451938) q[1];
sx q[1];
rz(11.558029) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5429687) q[0];
sx q[0];
rz(-2.7592139) q[0];
sx q[0];
rz(-1.3578869) q[0];
x q[1];
rz(1.8062334) q[2];
sx q[2];
rz(-1.2004235) q[2];
sx q[2];
rz(-0.78200227) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5656149) q[1];
sx q[1];
rz(-1.7276141) q[1];
sx q[1];
rz(-0.48139907) q[1];
x q[2];
rz(1.3775695) q[3];
sx q[3];
rz(-1.5504897) q[3];
sx q[3];
rz(1.3042252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8334373) q[2];
sx q[2];
rz(-2.1980632) q[2];
sx q[2];
rz(0.65307871) q[2];
rz(-3.0270882) q[3];
sx q[3];
rz(-1.7479892) q[3];
sx q[3];
rz(-2.0792927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7750875) q[0];
sx q[0];
rz(-2.1108284) q[0];
sx q[0];
rz(1.3436226) q[0];
rz(1.1603181) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(0.62058273) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7337429) q[0];
sx q[0];
rz(-1.7622927) q[0];
sx q[0];
rz(-1.922185) q[0];
x q[1];
rz(2.0257904) q[2];
sx q[2];
rz(-0.7209076) q[2];
sx q[2];
rz(1.4269331) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3239336) q[1];
sx q[1];
rz(-0.91485564) q[1];
sx q[1];
rz(-2.0327912) q[1];
rz(1.9602523) q[3];
sx q[3];
rz(-0.6529633) q[3];
sx q[3];
rz(0.12784004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6985942) q[2];
sx q[2];
rz(-1.6220762) q[2];
sx q[2];
rz(0.2300187) q[2];
rz(-1.5551785) q[3];
sx q[3];
rz(-1.9929726) q[3];
sx q[3];
rz(2.8659099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8860633) q[0];
sx q[0];
rz(-0.38917381) q[0];
sx q[0];
rz(2.0607167) q[0];
rz(-2.4983662) q[1];
sx q[1];
rz(-0.79935646) q[1];
sx q[1];
rz(-0.08531514) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15835831) q[0];
sx q[0];
rz(-2.2081349) q[0];
sx q[0];
rz(-2.4054224) q[0];
x q[1];
rz(1.3582466) q[2];
sx q[2];
rz(-1.0091127) q[2];
sx q[2];
rz(1.4457653) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9602365) q[1];
sx q[1];
rz(-1.4048368) q[1];
sx q[1];
rz(1.7039429) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.078104) q[3];
sx q[3];
rz(-2.5194019) q[3];
sx q[3];
rz(0.49294642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2049415) q[2];
sx q[2];
rz(-3.1320269) q[2];
sx q[2];
rz(0.19634518) q[2];
rz(-2.8593072) q[3];
sx q[3];
rz(-2.0245656) q[3];
sx q[3];
rz(2.096368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.1333756) q[0];
sx q[0];
rz(-2.5947925) q[0];
sx q[0];
rz(-0.38544449) q[0];
rz(-1.4133981) q[1];
sx q[1];
rz(-1.2047267) q[1];
sx q[1];
rz(2.8916496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88943931) q[0];
sx q[0];
rz(-1.982054) q[0];
sx q[0];
rz(2.2338339) q[0];
rz(0.087935178) q[2];
sx q[2];
rz(-1.5107949) q[2];
sx q[2];
rz(2.0363925) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.085332) q[1];
sx q[1];
rz(-2.0592494) q[1];
sx q[1];
rz(-1.7621098) q[1];
rz(-pi) q[2];
x q[2];
rz(2.945845) q[3];
sx q[3];
rz(-2.3382814) q[3];
sx q[3];
rz(-2.3565528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6803153) q[2];
sx q[2];
rz(-1.8467434) q[2];
sx q[2];
rz(-2.9052022) q[2];
rz(1.2605028) q[3];
sx q[3];
rz(-2.8915296) q[3];
sx q[3];
rz(0.67729706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8054955) q[0];
sx q[0];
rz(-1.5376872) q[0];
sx q[0];
rz(-2.6564823) q[0];
rz(-1.1803892) q[1];
sx q[1];
rz(-1.5943269) q[1];
sx q[1];
rz(0.83650437) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1331951) q[0];
sx q[0];
rz(-2.6566165) q[0];
sx q[0];
rz(0.72686813) q[0];
rz(-0.91412394) q[2];
sx q[2];
rz(-0.49490041) q[2];
sx q[2];
rz(-0.47593853) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4662469) q[1];
sx q[1];
rz(-2.2398178) q[1];
sx q[1];
rz(0.061208486) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4007217) q[3];
sx q[3];
rz(-1.5127137) q[3];
sx q[3];
rz(-2.7207014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2430719) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(2.3720429) q[2];
rz(0.92929333) q[3];
sx q[3];
rz(-2.010689) q[3];
sx q[3];
rz(2.9088959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.1997851) q[0];
sx q[0];
rz(-2.6564044) q[0];
sx q[0];
rz(0.76882452) q[0];
rz(-1.3517514) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(-0.88968712) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.453196) q[0];
sx q[0];
rz(-1.9510117) q[0];
sx q[0];
rz(-2.9158457) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5403929) q[2];
sx q[2];
rz(-0.044375751) q[2];
sx q[2];
rz(2.6868683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8823679) q[1];
sx q[1];
rz(-0.15341694) q[1];
sx q[1];
rz(-1.2620366) q[1];
x q[2];
rz(1.6857288) q[3];
sx q[3];
rz(-2.0191666) q[3];
sx q[3];
rz(1.2796206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6512904) q[2];
sx q[2];
rz(-1.4375261) q[2];
sx q[2];
rz(-1.5865405) q[2];
rz(1.2706612) q[3];
sx q[3];
rz(-0.41897604) q[3];
sx q[3];
rz(3.0095625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9689869) q[0];
sx q[0];
rz(-1.1382599) q[0];
sx q[0];
rz(-1.1440811) q[0];
rz(0.83089685) q[1];
sx q[1];
rz(-1.7832719) q[1];
sx q[1];
rz(-2.3544618) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4517196) q[0];
sx q[0];
rz(-1.4362808) q[0];
sx q[0];
rz(-0.13264199) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2058555) q[2];
sx q[2];
rz(-2.3496186) q[2];
sx q[2];
rz(-2.9284649) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.28879189) q[1];
sx q[1];
rz(-1.1026942) q[1];
sx q[1];
rz(2.1062984) q[1];
rz(-2.5821286) q[3];
sx q[3];
rz(-1.9064404) q[3];
sx q[3];
rz(-1.0192724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0234915) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(-0.033163158) q[2];
rz(1.5432594) q[3];
sx q[3];
rz(-2.4256746) q[3];
sx q[3];
rz(1.5020812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24932662) q[0];
sx q[0];
rz(-1.5557657) q[0];
sx q[0];
rz(1.3931042) q[0];
rz(2.6847367) q[1];
sx q[1];
rz(-1.1445878) q[1];
sx q[1];
rz(1.1005864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7129214) q[0];
sx q[0];
rz(-1.572519) q[0];
sx q[0];
rz(-0.07600204) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69533657) q[2];
sx q[2];
rz(-1.0461263) q[2];
sx q[2];
rz(0.52548835) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17328429) q[1];
sx q[1];
rz(-1.5368764) q[1];
sx q[1];
rz(-0.90356253) q[1];
rz(-pi) q[2];
rz(-2.8131717) q[3];
sx q[3];
rz(-2.7304318) q[3];
sx q[3];
rz(2.823165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.042772375) q[2];
sx q[2];
rz(-0.49751147) q[2];
sx q[2];
rz(1.5607321) q[2];
rz(-0.76357311) q[3];
sx q[3];
rz(-1.5127425) q[3];
sx q[3];
rz(-1.0360576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40912691) q[0];
sx q[0];
rz(-2.3865073) q[0];
sx q[0];
rz(-0.29148802) q[0];
rz(-1.905722) q[1];
sx q[1];
rz(-1.9449077) q[1];
sx q[1];
rz(-3.018697) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50704573) q[0];
sx q[0];
rz(-1.6222008) q[0];
sx q[0];
rz(0.34684166) q[0];
x q[1];
rz(-2.2332623) q[2];
sx q[2];
rz(-0.35229063) q[2];
sx q[2];
rz(0.34991821) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5220241) q[1];
sx q[1];
rz(-1.777473) q[1];
sx q[1];
rz(-3.0289438) q[1];
x q[2];
rz(-3.1355643) q[3];
sx q[3];
rz(-1.42618) q[3];
sx q[3];
rz(0.58352375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5128532) q[2];
sx q[2];
rz(-1.7491919) q[2];
sx q[2];
rz(2.0780308) q[2];
rz(0.17744803) q[3];
sx q[3];
rz(-0.95824233) q[3];
sx q[3];
rz(1.0183081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7040831) q[0];
sx q[0];
rz(-2.6884485) q[0];
sx q[0];
rz(-1.4310687) q[0];
rz(0.60072947) q[1];
sx q[1];
rz(-2.6917916) q[1];
sx q[1];
rz(2.5240555) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3356614) q[0];
sx q[0];
rz(-1.7863723) q[0];
sx q[0];
rz(3.0144601) q[0];
x q[1];
rz(-2.3138758) q[2];
sx q[2];
rz(-1.1679782) q[2];
sx q[2];
rz(3.0981491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.70676196) q[1];
sx q[1];
rz(-1.037957) q[1];
sx q[1];
rz(-3.117242) q[1];
x q[2];
rz(-1.2631868) q[3];
sx q[3];
rz(-1.0196536) q[3];
sx q[3];
rz(1.8637509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2403229) q[2];
sx q[2];
rz(-1.0003041) q[2];
sx q[2];
rz(-2.1642302) q[2];
rz(1.0824925) q[3];
sx q[3];
rz(-0.87477028) q[3];
sx q[3];
rz(3.0860743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.33901535) q[0];
sx q[0];
rz(-2.3727198) q[0];
sx q[0];
rz(2.4045237) q[0];
rz(3.0958685) q[1];
sx q[1];
rz(-2.3241691) q[1];
sx q[1];
rz(-1.587422) q[1];
rz(-2.3297119) q[2];
sx q[2];
rz(-2.4627081) q[2];
sx q[2];
rz(-2.2344786) q[2];
rz(1.4010728) q[3];
sx q[3];
rz(-2.2206497) q[3];
sx q[3];
rz(2.7322265) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
