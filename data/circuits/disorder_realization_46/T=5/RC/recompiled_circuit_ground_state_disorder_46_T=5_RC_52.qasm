OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98193869) q[0];
sx q[0];
rz(-1.5768134) q[0];
sx q[0];
rz(-1.2793581) q[0];
rz(0.98969412) q[1];
sx q[1];
rz(-1.1969748) q[1];
sx q[1];
rz(2.9762414) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26583126) q[0];
sx q[0];
rz(-0.94801694) q[0];
sx q[0];
rz(-0.39259194) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3897091) q[2];
sx q[2];
rz(-2.5640629) q[2];
sx q[2];
rz(-0.60345381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.27439427) q[1];
sx q[1];
rz(-1.9748678) q[1];
sx q[1];
rz(-1.0342802) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2591441) q[3];
sx q[3];
rz(-1.8168279) q[3];
sx q[3];
rz(-0.450799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2692261) q[2];
sx q[2];
rz(-1.9908315) q[2];
sx q[2];
rz(1.9744138) q[2];
rz(0.40927467) q[3];
sx q[3];
rz(-0.27547488) q[3];
sx q[3];
rz(-2.1854713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.069365) q[0];
sx q[0];
rz(-2.9874711) q[0];
sx q[0];
rz(-1.3519721) q[0];
rz(-1.6701472) q[1];
sx q[1];
rz(-0.92266005) q[1];
sx q[1];
rz(-2.2460489) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2118788) q[0];
sx q[0];
rz(-1.3037056) q[0];
sx q[0];
rz(-1.5729089) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47426736) q[2];
sx q[2];
rz(-1.5014428) q[2];
sx q[2];
rz(2.4276395) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2221007) q[1];
sx q[1];
rz(-2.2513012) q[1];
sx q[1];
rz(-0.51503985) q[1];
rz(-pi) q[2];
rz(2.0655153) q[3];
sx q[3];
rz(-0.46929324) q[3];
sx q[3];
rz(-2.9772907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46513778) q[2];
sx q[2];
rz(-0.82200161) q[2];
sx q[2];
rz(1.4466393) q[2];
rz(2.1220574) q[3];
sx q[3];
rz(-1.3186224) q[3];
sx q[3];
rz(2.8442966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4296221) q[0];
sx q[0];
rz(-1.977704) q[0];
sx q[0];
rz(-3.1373366) q[0];
rz(2.440522) q[1];
sx q[1];
rz(-2.6316167) q[1];
sx q[1];
rz(3.1140936) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0365899) q[0];
sx q[0];
rz(-1.9677015) q[0];
sx q[0];
rz(2.9734625) q[0];
x q[1];
rz(-1.4239935) q[2];
sx q[2];
rz(-0.73063421) q[2];
sx q[2];
rz(-0.50792009) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8054652) q[1];
sx q[1];
rz(-1.256811) q[1];
sx q[1];
rz(2.9538395) q[1];
x q[2];
rz(-0.013935907) q[3];
sx q[3];
rz(-1.2346141) q[3];
sx q[3];
rz(2.5232836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49715257) q[2];
sx q[2];
rz(-1.3500328) q[2];
sx q[2];
rz(2.7437317) q[2];
rz(-1.0128939) q[3];
sx q[3];
rz(-1.0254859) q[3];
sx q[3];
rz(2.7329172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560646) q[0];
sx q[0];
rz(-1.2701472) q[0];
sx q[0];
rz(2.2960466) q[0];
rz(-3.1327418) q[1];
sx q[1];
rz(-0.84819853) q[1];
sx q[1];
rz(-1.1525851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51868472) q[0];
sx q[0];
rz(-0.9810514) q[0];
sx q[0];
rz(0.49341644) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0716702) q[2];
sx q[2];
rz(-0.50412699) q[2];
sx q[2];
rz(1.1386516) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8804255) q[1];
sx q[1];
rz(-2.5486055) q[1];
sx q[1];
rz(-0.050131948) q[1];
rz(2.827876) q[3];
sx q[3];
rz(-1.6339746) q[3];
sx q[3];
rz(0.74387315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3919966) q[2];
sx q[2];
rz(-1.0276724) q[2];
sx q[2];
rz(1.9409625) q[2];
rz(-2.6642753) q[3];
sx q[3];
rz(-0.47003191) q[3];
sx q[3];
rz(-2.4376455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93477997) q[0];
sx q[0];
rz(-2.2581357) q[0];
sx q[0];
rz(-2.2373037) q[0];
rz(0.9710871) q[1];
sx q[1];
rz(-1.1957518) q[1];
sx q[1];
rz(-2.8579874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47434092) q[0];
sx q[0];
rz(-1.5548238) q[0];
sx q[0];
rz(-3.0706997) q[0];
x q[1];
rz(-2.0118158) q[2];
sx q[2];
rz(-1.9649817) q[2];
sx q[2];
rz(1.4798543) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0086669027) q[1];
sx q[1];
rz(-2.2382793) q[1];
sx q[1];
rz(-0.26393601) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0002076) q[3];
sx q[3];
rz(-2.3688487) q[3];
sx q[3];
rz(-0.34602133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35563719) q[2];
sx q[2];
rz(-1.9060308) q[2];
sx q[2];
rz(-2.9034485) q[2];
rz(0.98020482) q[3];
sx q[3];
rz(-1.1583867) q[3];
sx q[3];
rz(0.75077209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39764431) q[0];
sx q[0];
rz(-1.5357966) q[0];
sx q[0];
rz(-0.407298) q[0];
rz(0.40335718) q[1];
sx q[1];
rz(-0.74059161) q[1];
sx q[1];
rz(-0.86722803) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24466022) q[0];
sx q[0];
rz(-1.5112025) q[0];
sx q[0];
rz(-2.1367349) q[0];
rz(-pi) q[1];
rz(-0.32829702) q[2];
sx q[2];
rz(-0.5813501) q[2];
sx q[2];
rz(2.341908) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2867409) q[1];
sx q[1];
rz(-2.2959202) q[1];
sx q[1];
rz(1.185536) q[1];
x q[2];
rz(1.8676742) q[3];
sx q[3];
rz(-1.6921077) q[3];
sx q[3];
rz(-2.6685126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0200218) q[2];
sx q[2];
rz(-0.82841221) q[2];
sx q[2];
rz(1.7106445) q[2];
rz(-2.6089148) q[3];
sx q[3];
rz(-2.1055652) q[3];
sx q[3];
rz(-1.4177586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1581887) q[0];
sx q[0];
rz(-0.15022763) q[0];
sx q[0];
rz(-1.0839373) q[0];
rz(-0.051941959) q[1];
sx q[1];
rz(-0.40819326) q[1];
sx q[1];
rz(0.025024978) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0662224) q[0];
sx q[0];
rz(-2.0522723) q[0];
sx q[0];
rz(-1.6027149) q[0];
rz(0.81871732) q[2];
sx q[2];
rz(-1.229964) q[2];
sx q[2];
rz(2.9342295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.53326239) q[1];
sx q[1];
rz(-0.46583781) q[1];
sx q[1];
rz(1.9598258) q[1];
x q[2];
rz(-0.30625108) q[3];
sx q[3];
rz(-2.6970377) q[3];
sx q[3];
rz(-1.2769092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0742801) q[2];
sx q[2];
rz(-1.68579) q[2];
sx q[2];
rz(2.32453) q[2];
rz(0.95064154) q[3];
sx q[3];
rz(-1.0167511) q[3];
sx q[3];
rz(-1.1869717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.10512146) q[0];
sx q[0];
rz(-2.8558185) q[0];
sx q[0];
rz(2.5407319) q[0];
rz(-0.034424456) q[1];
sx q[1];
rz(-1.8079146) q[1];
sx q[1];
rz(1.5404125) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2297771) q[0];
sx q[0];
rz(-0.86965771) q[0];
sx q[0];
rz(0.83155737) q[0];
rz(0.95048423) q[2];
sx q[2];
rz(-0.99603876) q[2];
sx q[2];
rz(1.7983914) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2331761) q[1];
sx q[1];
rz(-1.5926835) q[1];
sx q[1];
rz(-2.7018956) q[1];
rz(0.65798379) q[3];
sx q[3];
rz(-2.4870928) q[3];
sx q[3];
rz(-0.94098781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2711266) q[2];
sx q[2];
rz(-1.8693962) q[2];
sx q[2];
rz(1.6477443) q[2];
rz(2.3624524) q[3];
sx q[3];
rz(-1.6611764) q[3];
sx q[3];
rz(-0.65472764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43799022) q[0];
sx q[0];
rz(-2.0605189) q[0];
sx q[0];
rz(-0.78824747) q[0];
rz(-1.2429271) q[1];
sx q[1];
rz(-1.5595167) q[1];
sx q[1];
rz(-2.8020614) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5632322) q[0];
sx q[0];
rz(-1.6867016) q[0];
sx q[0];
rz(-1.3770335) q[0];
rz(-pi) q[1];
rz(2.7044136) q[2];
sx q[2];
rz(-0.60761062) q[2];
sx q[2];
rz(-2.2504582) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.06188678) q[1];
sx q[1];
rz(-1.8822274) q[1];
sx q[1];
rz(2.7175886) q[1];
rz(-2.3655534) q[3];
sx q[3];
rz(-1.0747915) q[3];
sx q[3];
rz(-2.1797594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1303611) q[2];
sx q[2];
rz(-1.9731382) q[2];
sx q[2];
rz(1.2569024) q[2];
rz(-1.0907178) q[3];
sx q[3];
rz(-1.9382449) q[3];
sx q[3];
rz(-2.2244661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.286769) q[0];
sx q[0];
rz(-1.6463065) q[0];
sx q[0];
rz(2.0538034) q[0];
rz(0.65025672) q[1];
sx q[1];
rz(-1.6842664) q[1];
sx q[1];
rz(-3.1204209) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0567704) q[0];
sx q[0];
rz(-2.6426043) q[0];
sx q[0];
rz(-0.4468687) q[0];
rz(-pi) q[1];
rz(-1.3743873) q[2];
sx q[2];
rz(-2.0337262) q[2];
sx q[2];
rz(0.48505515) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.260878) q[1];
sx q[1];
rz(-1.78111) q[1];
sx q[1];
rz(2.3349663) q[1];
rz(-3.0983119) q[3];
sx q[3];
rz(-1.0387522) q[3];
sx q[3];
rz(-0.30902853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0040969) q[2];
sx q[2];
rz(-1.7986412) q[2];
sx q[2];
rz(0.98319483) q[2];
rz(2.148597) q[3];
sx q[3];
rz(-1.1119305) q[3];
sx q[3];
rz(0.23428169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8773593) q[0];
sx q[0];
rz(-0.79240427) q[0];
sx q[0];
rz(-2.0808676) q[0];
rz(-1.1119153) q[1];
sx q[1];
rz(-2.8363375) q[1];
sx q[1];
rz(-3.0189966) q[1];
rz(-1.1373043) q[2];
sx q[2];
rz(-0.52327427) q[2];
sx q[2];
rz(0.56868194) q[2];
rz(-3.1226782) q[3];
sx q[3];
rz(-0.52302845) q[3];
sx q[3];
rz(1.1136644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
