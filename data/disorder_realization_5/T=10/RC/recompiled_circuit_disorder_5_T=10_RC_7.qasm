OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(2.6846057) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(-2.4609202) q[1];
sx q[1];
rz(-1.2759804) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9306643) q[0];
sx q[0];
rz(-1.352248) q[0];
sx q[0];
rz(1.4320089) q[0];
x q[1];
rz(-3.0511191) q[2];
sx q[2];
rz(-1.4960559) q[2];
sx q[2];
rz(0.59404101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4133271) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(3.0800746) q[1];
rz(-0.39235093) q[3];
sx q[3];
rz(-2.1808743) q[3];
sx q[3];
rz(-0.073079212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(0.16513744) q[2];
rz(-1.5103546) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(-2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5052658) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(0.045036137) q[0];
rz(-0.33915195) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(-0.80274686) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0720285) q[0];
sx q[0];
rz(-1.1840128) q[0];
sx q[0];
rz(0.2441498) q[0];
rz(-pi) q[1];
rz(1.3834125) q[2];
sx q[2];
rz(-0.51088453) q[2];
sx q[2];
rz(0.62142205) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3853559) q[1];
sx q[1];
rz(-1.6223575) q[1];
sx q[1];
rz(0.95183422) q[1];
rz(-2.2400168) q[3];
sx q[3];
rz(-1.672861) q[3];
sx q[3];
rz(-1.4249742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37614432) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.7155898) q[2];
rz(-2.5358893) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894492) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(2.8564575) q[0];
rz(0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(-1.3396938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08683603) q[0];
sx q[0];
rz(-1.9924379) q[0];
sx q[0];
rz(-0.7936759) q[0];
rz(-1.2013024) q[2];
sx q[2];
rz(-2.4174147) q[2];
sx q[2];
rz(2.5232814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7741751) q[1];
sx q[1];
rz(-1.4111102) q[1];
sx q[1];
rz(-1.3819441) q[1];
x q[2];
rz(-0.47311802) q[3];
sx q[3];
rz(-2.3948673) q[3];
sx q[3];
rz(2.8603922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4703935) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(-1.7335256) q[2];
rz(0.18313289) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-0.071921913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4801487) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(-2.9273422) q[0];
rz(0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(2.4750211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65581215) q[0];
sx q[0];
rz(-2.862145) q[0];
sx q[0];
rz(2.4026472) q[0];
rz(-pi) q[1];
rz(-2.4273657) q[2];
sx q[2];
rz(-1.3492609) q[2];
sx q[2];
rz(-2.5941338) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2506634) q[1];
sx q[1];
rz(-0.86859413) q[1];
sx q[1];
rz(0.82171792) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.781226) q[3];
sx q[3];
rz(-0.20877148) q[3];
sx q[3];
rz(-2.3695932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(2.4812223) q[2];
rz(-1.9541698) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(-0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289571) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-0.83706013) q[0];
rz(0.95056668) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(-1.1594835) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9966853) q[0];
sx q[0];
rz(-2.3527745) q[0];
sx q[0];
rz(-0.44992723) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3443089) q[2];
sx q[2];
rz(-0.61033568) q[2];
sx q[2];
rz(0.047705334) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3573208) q[1];
sx q[1];
rz(-1.1457774) q[1];
sx q[1];
rz(-0.97370633) q[1];
x q[2];
rz(-1.0347576) q[3];
sx q[3];
rz(-1.0039819) q[3];
sx q[3];
rz(-2.00622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7065457) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(-1.2716028) q[2];
rz(1.050625) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5669252) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(2.0264453) q[0];
rz(-0.043958157) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(0.10087068) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2642919) q[0];
sx q[0];
rz(-1.722324) q[0];
sx q[0];
rz(0.63440462) q[0];
rz(-pi) q[1];
rz(-1.2783373) q[2];
sx q[2];
rz(-2.0403701) q[2];
sx q[2];
rz(1.3172319) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85943356) q[1];
sx q[1];
rz(-1.7267425) q[1];
sx q[1];
rz(-1.366051) q[1];
rz(2.7363339) q[3];
sx q[3];
rz(-1.3558946) q[3];
sx q[3];
rz(-2.1586771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7373401) q[2];
sx q[2];
rz(-1.5161783) q[2];
sx q[2];
rz(-0.92528382) q[2];
rz(2.960079) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(-1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3826542) q[0];
sx q[0];
rz(-2.721334) q[0];
sx q[0];
rz(2.1784901) q[0];
rz(1.5513647) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(-2.4538453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9822385) q[0];
sx q[0];
rz(-1.237545) q[0];
sx q[0];
rz(0.23157816) q[0];
rz(-pi) q[1];
rz(-0.22600941) q[2];
sx q[2];
rz(-1.2603716) q[2];
sx q[2];
rz(-1.7058536) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2076599) q[1];
sx q[1];
rz(-1.4414409) q[1];
sx q[1];
rz(1.5473066) q[1];
rz(0.01597605) q[3];
sx q[3];
rz(-1.0025347) q[3];
sx q[3];
rz(-2.361134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(-2.2040099) q[2];
rz(-0.98085105) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5279609) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(-3.0623073) q[0];
rz(-1.1212564) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(-1.1987196) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51533651) q[0];
sx q[0];
rz(-1.5021828) q[0];
sx q[0];
rz(-2.7068044) q[0];
x q[1];
rz(-1.933868) q[2];
sx q[2];
rz(-1.6512401) q[2];
sx q[2];
rz(2.9094839) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.291154) q[1];
sx q[1];
rz(-1.6206695) q[1];
sx q[1];
rz(0.50317851) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3098605) q[3];
sx q[3];
rz(-1.121215) q[3];
sx q[3];
rz(-0.98228067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2856059) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(2.3525227) q[2];
rz(2.7834535) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(-1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.566074) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(3.1148124) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(2.2081597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8342455) q[0];
sx q[0];
rz(-1.5114307) q[0];
sx q[0];
rz(-3.0349915) q[0];
rz(2.2869296) q[2];
sx q[2];
rz(-0.48644201) q[2];
sx q[2];
rz(-1.9180627) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.583608) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(-3.0834553) q[1];
rz(-2.6870319) q[3];
sx q[3];
rz(-0.70485669) q[3];
sx q[3];
rz(2.9398033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.062722) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(-0.36994568) q[2];
rz(-0.89921078) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(2.1883011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(2.488234) q[0];
rz(1.6409142) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(-2.9577589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95589721) q[0];
sx q[0];
rz(-0.16793185) q[0];
sx q[0];
rz(-0.76783617) q[0];
rz(-2.3046938) q[2];
sx q[2];
rz(-2.1914334) q[2];
sx q[2];
rz(1.9377973) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4348236) q[1];
sx q[1];
rz(-2.4425047) q[1];
sx q[1];
rz(-1.8301151) q[1];
rz(0.58335431) q[3];
sx q[3];
rz(-1.5837985) q[3];
sx q[3];
rz(2.0875723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.9188149) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(-0.20239057) q[2];
rz(1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68207537) q[0];
sx q[0];
rz(-2.2820602) q[0];
sx q[0];
rz(2.5862502) q[0];
rz(-2.7783685) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(0.72485483) q[2];
sx q[2];
rz(-0.58314322) q[2];
sx q[2];
rz(-2.4287139) q[2];
rz(-0.17584569) q[3];
sx q[3];
rz(-0.98777117) q[3];
sx q[3];
rz(-2.5143757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];