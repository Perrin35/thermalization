OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6558134) q[0];
sx q[0];
rz(-1.2278779) q[0];
sx q[0];
rz(1.2113843) q[0];
rz(-0.76531305) q[2];
sx q[2];
rz(-1.0368477) q[2];
sx q[2];
rz(0.050616654) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.624144) q[1];
sx q[1];
rz(-1.9074719) q[1];
sx q[1];
rz(1.3144073) q[1];
x q[2];
rz(-1.0899815) q[3];
sx q[3];
rz(-0.40502031) q[3];
sx q[3];
rz(2.4570176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87542614) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(-1.1323294) q[2];
rz(-1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(-2.1291389) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-0.18584132) q[0];
rz(-2.5813685) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(2.9247608) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2730392) q[0];
sx q[0];
rz(-0.93398636) q[0];
sx q[0];
rz(-0.36006948) q[0];
rz(-pi) q[1];
rz(-0.88044135) q[2];
sx q[2];
rz(-2.464622) q[2];
sx q[2];
rz(-1.134269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.195897) q[1];
sx q[1];
rz(-0.90598124) q[1];
sx q[1];
rz(-2.871454) q[1];
x q[2];
rz(0.64036815) q[3];
sx q[3];
rz(-2.159517) q[3];
sx q[3];
rz(2.0224188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.310114) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(1.2878093) q[2];
rz(-0.76256049) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-0.30502239) q[3];
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
rz(-1.4644311) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(0.60423869) q[0];
rz(-1.8151981) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(2.2089829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9661449) q[0];
sx q[0];
rz(-1.6245337) q[0];
sx q[0];
rz(1.8885814) q[0];
x q[1];
rz(0.18205299) q[2];
sx q[2];
rz(-1.7943873) q[2];
sx q[2];
rz(-1.455866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4438666) q[1];
sx q[1];
rz(-1.6892471) q[1];
sx q[1];
rz(0.5603793) q[1];
rz(-pi) q[2];
rz(2.1266537) q[3];
sx q[3];
rz(-1.1851289) q[3];
sx q[3];
rz(2.8876497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9937667) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(1.0926584) q[2];
rz(-0.5422194) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.7820691) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(0.50022593) q[0];
rz(-2.3362828) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(-1.6436228) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8805807) q[0];
sx q[0];
rz(-1.2769165) q[0];
sx q[0];
rz(2.0902993) q[0];
rz(-pi) q[1];
rz(3.0849886) q[2];
sx q[2];
rz(-1.3805693) q[2];
sx q[2];
rz(-2.9375926) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8343463) q[1];
sx q[1];
rz(-1.3386968) q[1];
sx q[1];
rz(0.26718617) q[1];
rz(-pi) q[2];
rz(0.76969947) q[3];
sx q[3];
rz(-0.61004988) q[3];
sx q[3];
rz(1.1806012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.74636373) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(0.70181075) q[2];
rz(2.3102405) q[3];
sx q[3];
rz(-0.96389198) q[3];
sx q[3];
rz(2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2410626) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(-0.81533122) q[0];
rz(1.5218081) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(1.048208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5883023) q[0];
sx q[0];
rz(-2.7748845) q[0];
sx q[0];
rz(-2.328863) q[0];
rz(3.1254966) q[2];
sx q[2];
rz(-0.65158366) q[2];
sx q[2];
rz(-0.96166699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8480307) q[1];
sx q[1];
rz(-1.8082431) q[1];
sx q[1];
rz(-2.8192987) q[1];
x q[2];
rz(-2.1225554) q[3];
sx q[3];
rz(-2.4487547) q[3];
sx q[3];
rz(0.90941959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52577019) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(2.0416416) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-2.0992978) q[3];
sx q[3];
rz(-0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0681756) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(-2.2391879) q[0];
rz(1.0166608) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(3.0117603) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3467305) q[0];
sx q[0];
rz(-0.71160337) q[0];
sx q[0];
rz(-2.5559588) q[0];
x q[1];
rz(-2.1075222) q[2];
sx q[2];
rz(-0.64646361) q[2];
sx q[2];
rz(1.549364) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1631158) q[1];
sx q[1];
rz(-2.4448131) q[1];
sx q[1];
rz(-0.58660581) q[1];
rz(-2.5927605) q[3];
sx q[3];
rz(-1.7389113) q[3];
sx q[3];
rz(2.7041534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8292024) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(2.9373346) q[2];
rz(1.2060818) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234574) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(-1.4468505) q[0];
rz(1.2591259) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(0.68626219) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670358) q[0];
sx q[0];
rz(-1.0970763) q[0];
sx q[0];
rz(2.0635598) q[0];
rz(-pi) q[1];
rz(-2.1967728) q[2];
sx q[2];
rz(-0.84569028) q[2];
sx q[2];
rz(-0.041989728) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8441019) q[1];
sx q[1];
rz(-2.9417848) q[1];
sx q[1];
rz(1.1872477) q[1];
rz(-pi) q[2];
rz(-1.2495743) q[3];
sx q[3];
rz(-1.9739082) q[3];
sx q[3];
rz(-1.5461127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4454322) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(2.5799675) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(1.5047489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5381662) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(1.4461393) q[0];
rz(2.360545) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(-1.5015645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7555996) q[0];
sx q[0];
rz(-2.0634683) q[0];
sx q[0];
rz(-3.1194869) q[0];
rz(-1.5090452) q[2];
sx q[2];
rz(-1.5973063) q[2];
sx q[2];
rz(2.3960631) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1135243) q[1];
sx q[1];
rz(-1.4816195) q[1];
sx q[1];
rz(0.32797565) q[1];
x q[2];
rz(-0.6286962) q[3];
sx q[3];
rz(-2.0572212) q[3];
sx q[3];
rz(-1.4592255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35187307) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(1.3191351) q[2];
rz(-1.9296648) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33655745) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(1.9375027) q[0];
rz(2.7583292) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(-0.35167545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4597804) q[0];
sx q[0];
rz(-0.60428719) q[0];
sx q[0];
rz(-0.87862815) q[0];
rz(-0.69182379) q[2];
sx q[2];
rz(-1.8080538) q[2];
sx q[2];
rz(-2.0123864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55191509) q[1];
sx q[1];
rz(-2.3304686) q[1];
sx q[1];
rz(0.07304904) q[1];
rz(-pi) q[2];
rz(-2.6978108) q[3];
sx q[3];
rz(-3.0203331) q[3];
sx q[3];
rz(-1.4951984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(-1.8593672) q[2];
rz(1.4964237) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(2.0857247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6431817) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(0.19432755) q[0];
rz(-1.0378029) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(-1.0338354) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4170096) q[0];
sx q[0];
rz(-1.7010744) q[0];
sx q[0];
rz(0.91086046) q[0];
rz(-pi) q[1];
rz(-2.1543703) q[2];
sx q[2];
rz(-0.7910896) q[2];
sx q[2];
rz(2.1949878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6453213) q[1];
sx q[1];
rz(-1.4693345) q[1];
sx q[1];
rz(-2.1437777) q[1];
x q[2];
rz(1.7887573) q[3];
sx q[3];
rz(-2.2721604) q[3];
sx q[3];
rz(2.4234114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0620492) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(-0.27030269) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.6132145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939659) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(1.7059965) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(-3.1096935) q[2];
sx q[2];
rz(-0.96822856) q[2];
sx q[2];
rz(-0.4005489) q[2];
rz(1.773949) q[3];
sx q[3];
rz(-0.33645867) q[3];
sx q[3];
rz(-2.4154739) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
