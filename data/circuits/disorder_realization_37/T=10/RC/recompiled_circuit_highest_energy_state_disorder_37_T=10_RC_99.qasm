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
rz(-1.3353424) q[0];
sx q[0];
rz(3.5080533) q[0];
sx q[0];
rz(8.9656497) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(-1.4067283) q[1];
sx q[1];
rz(-0.23174098) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.988294) q[0];
sx q[0];
rz(-0.35044119) q[0];
sx q[0];
rz(-0.20163433) q[0];
rz(-pi) q[1];
rz(-0.094005748) q[2];
sx q[2];
rz(-2.5931907) q[2];
sx q[2];
rz(1.7552055) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8084131) q[1];
sx q[1];
rz(-1.1929379) q[1];
sx q[1];
rz(0.66587944) q[1];
x q[2];
rz(-1.6123338) q[3];
sx q[3];
rz(-1.9088863) q[3];
sx q[3];
rz(1.5134144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0509402) q[2];
sx q[2];
rz(-2.6068164) q[2];
sx q[2];
rz(-1.862662) q[2];
rz(-2.2517962) q[3];
sx q[3];
rz(-1.5225478) q[3];
sx q[3];
rz(0.33862996) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58562529) q[0];
sx q[0];
rz(-2.2323759) q[0];
sx q[0];
rz(0.92581785) q[0];
rz(1.0649118) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(-3.056114) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31562284) q[0];
sx q[0];
rz(-1.1809826) q[0];
sx q[0];
rz(-1.5672383) q[0];
rz(-pi) q[1];
rz(-2.1824942) q[2];
sx q[2];
rz(-1.2818205) q[2];
sx q[2];
rz(1.1346357) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.036219941) q[1];
sx q[1];
rz(-2.5565127) q[1];
sx q[1];
rz(1.1452504) q[1];
x q[2];
rz(2.5898629) q[3];
sx q[3];
rz(-2.3750616) q[3];
sx q[3];
rz(2.6006002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.815879) q[2];
sx q[2];
rz(-1.6625983) q[2];
sx q[2];
rz(2.7549287) q[2];
rz(1.9246842) q[3];
sx q[3];
rz(-2.7972126) q[3];
sx q[3];
rz(0.2200505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81095186) q[0];
sx q[0];
rz(-2.7702259) q[0];
sx q[0];
rz(-2.6445342) q[0];
rz(-1.0423543) q[1];
sx q[1];
rz(-0.94756871) q[1];
sx q[1];
rz(0.29464468) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2089068) q[0];
sx q[0];
rz(-2.8021376) q[0];
sx q[0];
rz(0.20821388) q[0];
rz(-pi) q[1];
rz(-1.7435837) q[2];
sx q[2];
rz(-2.2026988) q[2];
sx q[2];
rz(-1.8586707) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.55185071) q[1];
sx q[1];
rz(-2.183118) q[1];
sx q[1];
rz(2.603884) q[1];
x q[2];
rz(-1.7792652) q[3];
sx q[3];
rz(-2.5304171) q[3];
sx q[3];
rz(-2.0198422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3543388) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(1.5617237) q[2];
rz(1.076237) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(0.92897433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8266066) q[0];
sx q[0];
rz(-1.6153233) q[0];
sx q[0];
rz(-0.22931799) q[0];
rz(1.5062821) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(-0.062072676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4631043) q[0];
sx q[0];
rz(-2.0651428) q[0];
sx q[0];
rz(-0.93017471) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1770085) q[2];
sx q[2];
rz(-2.6754489) q[2];
sx q[2];
rz(-0.073748253) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8530555) q[1];
sx q[1];
rz(-1.2934577) q[1];
sx q[1];
rz(-0.71373516) q[1];
rz(-pi) q[2];
rz(-0.43020474) q[3];
sx q[3];
rz(-2.9211126) q[3];
sx q[3];
rz(1.1663495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0486003) q[2];
sx q[2];
rz(-2.2255662) q[2];
sx q[2];
rz(-1.3108866) q[2];
rz(-3.1033031) q[3];
sx q[3];
rz(-1.3524651) q[3];
sx q[3];
rz(0.37461764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49486092) q[0];
sx q[0];
rz(-2.4202388) q[0];
sx q[0];
rz(-0.89299655) q[0];
rz(-1.0294186) q[1];
sx q[1];
rz(-2.4118377) q[1];
sx q[1];
rz(-0.095887862) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5928837) q[0];
sx q[0];
rz(-2.0796695) q[0];
sx q[0];
rz(2.6611317) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0711543) q[2];
sx q[2];
rz(-1.6595006) q[2];
sx q[2];
rz(-0.34102893) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3600625) q[1];
sx q[1];
rz(-1.9340098) q[1];
sx q[1];
rz(2.5468154) q[1];
rz(-1.1981691) q[3];
sx q[3];
rz(-1.0851139) q[3];
sx q[3];
rz(-1.612965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.15478495) q[2];
sx q[2];
rz(-1.3900577) q[2];
sx q[2];
rz(-2.7161157) q[2];
rz(-0.42243877) q[3];
sx q[3];
rz(-2.0368302) q[3];
sx q[3];
rz(-0.5031684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4130037) q[0];
sx q[0];
rz(-0.63218963) q[0];
sx q[0];
rz(-0.11716209) q[0];
rz(1.6475742) q[1];
sx q[1];
rz(-2.2393176) q[1];
sx q[1];
rz(-1.0711627) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8534096) q[0];
sx q[0];
rz(-0.57153915) q[0];
sx q[0];
rz(2.6009212) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2920678) q[2];
sx q[2];
rz(-1.0459002) q[2];
sx q[2];
rz(-0.43290813) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2177762) q[1];
sx q[1];
rz(-2.4144396) q[1];
sx q[1];
rz(-1.9132861) q[1];
rz(-pi) q[2];
rz(2.8240194) q[3];
sx q[3];
rz(-1.7264719) q[3];
sx q[3];
rz(-2.6359216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5693207) q[2];
sx q[2];
rz(-0.57491493) q[2];
sx q[2];
rz(0.039483698) q[2];
rz(0.076400541) q[3];
sx q[3];
rz(-1.1698086) q[3];
sx q[3];
rz(-0.45342818) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4385248) q[0];
sx q[0];
rz(-2.330133) q[0];
sx q[0];
rz(-0.95034289) q[0];
rz(1.1514661) q[1];
sx q[1];
rz(-1.5299503) q[1];
sx q[1];
rz(1.8340402) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261953) q[0];
sx q[0];
rz(-1.2457677) q[0];
sx q[0];
rz(2.9103878) q[0];
rz(-pi) q[1];
rz(-1.6436395) q[2];
sx q[2];
rz(-1.1955639) q[2];
sx q[2];
rz(-0.90922395) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2491409) q[1];
sx q[1];
rz(-1.080852) q[1];
sx q[1];
rz(0.23995121) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96486196) q[3];
sx q[3];
rz(-1.4601225) q[3];
sx q[3];
rz(-2.0697274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82039708) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(0.55366984) q[2];
rz(-0.6959483) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(1.6863916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221136) q[0];
sx q[0];
rz(-1.6895634) q[0];
sx q[0];
rz(-2.8346862) q[0];
rz(-1.8652929) q[1];
sx q[1];
rz(-0.46638322) q[1];
sx q[1];
rz(1.9006405) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0725482) q[0];
sx q[0];
rz(-0.47397754) q[0];
sx q[0];
rz(-1.9182253) q[0];
x q[1];
rz(0.5979171) q[2];
sx q[2];
rz(-1.2094524) q[2];
sx q[2];
rz(-1.8449699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.13682374) q[1];
sx q[1];
rz(-1.3826177) q[1];
sx q[1];
rz(0.31463639) q[1];
rz(-pi) q[2];
rz(1.9718902) q[3];
sx q[3];
rz(-0.99906427) q[3];
sx q[3];
rz(-0.16756646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79157311) q[2];
sx q[2];
rz(-2.3363523) q[2];
sx q[2];
rz(-1.8512858) q[2];
rz(-0.6811412) q[3];
sx q[3];
rz(-2.2414424) q[3];
sx q[3];
rz(-1.6397938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27555585) q[0];
sx q[0];
rz(-1.0365726) q[0];
sx q[0];
rz(-1.1736897) q[0];
rz(-0.58700079) q[1];
sx q[1];
rz(-0.78158164) q[1];
sx q[1];
rz(3.0618844) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82725924) q[0];
sx q[0];
rz(-1.200935) q[0];
sx q[0];
rz(-0.24263675) q[0];
rz(-pi) q[1];
rz(-2.9201512) q[2];
sx q[2];
rz(-1.6802854) q[2];
sx q[2];
rz(1.2564645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.798493) q[1];
sx q[1];
rz(-2.1329555) q[1];
sx q[1];
rz(-0.72748277) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3354638) q[3];
sx q[3];
rz(-1.7793831) q[3];
sx q[3];
rz(-0.94193469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6012663) q[2];
sx q[2];
rz(-1.2286295) q[2];
sx q[2];
rz(-1.3339174) q[2];
rz(1.5909083) q[3];
sx q[3];
rz(-1.0100789) q[3];
sx q[3];
rz(-0.18868119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(0.48225668) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(-1.2905066) q[0];
rz(-0.036529649) q[1];
sx q[1];
rz(-1.9772269) q[1];
sx q[1];
rz(-2.0972924) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14257763) q[0];
sx q[0];
rz(-1.8685088) q[0];
sx q[0];
rz(-2.0172869) q[0];
x q[1];
rz(2.2453111) q[2];
sx q[2];
rz(-2.7283629) q[2];
sx q[2];
rz(-0.21597029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56230169) q[1];
sx q[1];
rz(-2.3962483) q[1];
sx q[1];
rz(-1.3424804) q[1];
rz(-pi) q[2];
rz(-2.812522) q[3];
sx q[3];
rz(-1.6794723) q[3];
sx q[3];
rz(1.7038801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2269939) q[2];
sx q[2];
rz(-2.1103766) q[2];
sx q[2];
rz(1.3247789) q[2];
rz(0.75657183) q[3];
sx q[3];
rz(-2.2533267) q[3];
sx q[3];
rz(2.2911086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.666438) q[0];
sx q[0];
rz(-1.403724) q[0];
sx q[0];
rz(-2.4028461) q[0];
rz(1.4229763) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(-2.6773166) q[2];
sx q[2];
rz(-2.3515679) q[2];
sx q[2];
rz(-0.49754561) q[2];
rz(1.1005836) q[3];
sx q[3];
rz(-0.96405021) q[3];
sx q[3];
rz(-1.8586803) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
