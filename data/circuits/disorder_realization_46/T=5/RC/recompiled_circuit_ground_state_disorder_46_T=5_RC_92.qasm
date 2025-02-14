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
rz(5.0862105) q[1];
sx q[1];
rz(12.401019) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5419614) q[0];
sx q[0];
rz(-1.2548057) q[0];
sx q[0];
rz(2.2315079) q[0];
rz(-pi) q[1];
rz(-1.1520391) q[2];
sx q[2];
rz(-1.9809696) q[2];
sx q[2];
rz(2.9048186) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.074881) q[1];
sx q[1];
rz(-1.0815269) q[1];
sx q[1];
rz(0.46164767) q[1];
rz(-1.2591441) q[3];
sx q[3];
rz(-1.3247648) q[3];
sx q[3];
rz(2.6907937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2692261) q[2];
sx q[2];
rz(-1.1507611) q[2];
sx q[2];
rz(1.1671789) q[2];
rz(0.40927467) q[3];
sx q[3];
rz(-2.8661178) q[3];
sx q[3];
rz(2.1854713) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.069365) q[0];
sx q[0];
rz(-0.15412155) q[0];
sx q[0];
rz(1.3519721) q[0];
rz(1.4714454) q[1];
sx q[1];
rz(-0.92266005) q[1];
sx q[1];
rz(-2.2460489) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35947511) q[0];
sx q[0];
rz(-1.572834) q[0];
sx q[0];
rz(-2.8745013) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4928713) q[2];
sx q[2];
rz(-2.04383) q[2];
sx q[2];
rz(-0.89240197) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9546763) q[1];
sx q[1];
rz(-2.3136931) q[1];
sx q[1];
rz(-2.1174341) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9053365) q[3];
sx q[3];
rz(-1.1614262) q[3];
sx q[3];
rz(-2.7619113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6764549) q[2];
sx q[2];
rz(-2.319591) q[2];
sx q[2];
rz(-1.6949534) q[2];
rz(-2.1220574) q[3];
sx q[3];
rz(-1.3186224) q[3];
sx q[3];
rz(-2.8442966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4296221) q[0];
sx q[0];
rz(-1.1638887) q[0];
sx q[0];
rz(-0.0042560552) q[0];
rz(-0.70107067) q[1];
sx q[1];
rz(-2.6316167) q[1];
sx q[1];
rz(3.1140936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6228559) q[0];
sx q[0];
rz(-0.42930005) q[0];
sx q[0];
rz(1.9506428) q[0];
x q[1];
rz(3.0112634) q[2];
sx q[2];
rz(-0.8497592) q[2];
sx q[2];
rz(-0.70391612) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.78464064) q[1];
sx q[1];
rz(-2.7773547) q[1];
sx q[1];
rz(1.0491153) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6106538) q[3];
sx q[3];
rz(-2.8051326) q[3];
sx q[3];
rz(0.57608673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.49715257) q[2];
sx q[2];
rz(-1.3500328) q[2];
sx q[2];
rz(-2.7437317) q[2];
rz(1.0128939) q[3];
sx q[3];
rz(-1.0254859) q[3];
sx q[3];
rz(0.40867543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48552805) q[0];
sx q[0];
rz(-1.2701472) q[0];
sx q[0];
rz(0.84554607) q[0];
rz(-3.1327418) q[1];
sx q[1];
rz(-0.84819853) q[1];
sx q[1];
rz(-1.1525851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3801013) q[0];
sx q[0];
rz(-1.9753755) q[0];
sx q[0];
rz(0.92104162) q[0];
rz(-pi) q[1];
rz(-0.25817008) q[2];
sx q[2];
rz(-2.008778) q[2];
sx q[2];
rz(-0.5817619) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3215969) q[1];
sx q[1];
rz(-2.1629371) q[1];
sx q[1];
rz(-1.5370398) q[1];
rz(-0.20219965) q[3];
sx q[3];
rz(-0.31980896) q[3];
sx q[3];
rz(0.6347189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3919966) q[2];
sx q[2];
rz(-2.1139202) q[2];
sx q[2];
rz(-1.9409625) q[2];
rz(2.6642753) q[3];
sx q[3];
rz(-2.6715607) q[3];
sx q[3];
rz(-2.4376455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93477997) q[0];
sx q[0];
rz(-2.2581357) q[0];
sx q[0];
rz(-2.2373037) q[0];
rz(-0.9710871) q[1];
sx q[1];
rz(-1.1957518) q[1];
sx q[1];
rz(-0.28360525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0975896) q[0];
sx q[0];
rz(-1.6416802) q[0];
sx q[0];
rz(1.5547836) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1297768) q[2];
sx q[2];
rz(-1.9649817) q[2];
sx q[2];
rz(-1.6617384) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.42029494) q[1];
sx q[1];
rz(-0.71027139) q[1];
sx q[1];
rz(1.8904449) q[1];
x q[2];
rz(1.0002076) q[3];
sx q[3];
rz(-0.772744) q[3];
sx q[3];
rz(-2.7955713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7859555) q[2];
sx q[2];
rz(-1.2355618) q[2];
sx q[2];
rz(0.23814417) q[2];
rz(2.1613878) q[3];
sx q[3];
rz(-1.1583867) q[3];
sx q[3];
rz(-0.75077209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7439483) q[0];
sx q[0];
rz(-1.605796) q[0];
sx q[0];
rz(-0.407298) q[0];
rz(-2.7382355) q[1];
sx q[1];
rz(-0.74059161) q[1];
sx q[1];
rz(2.2743646) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8969324) q[0];
sx q[0];
rz(-1.5112025) q[0];
sx q[0];
rz(-2.1367349) q[0];
rz(-pi) q[1];
rz(-2.5851493) q[2];
sx q[2];
rz(-1.3927937) q[2];
sx q[2];
rz(2.0931505) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6883863) q[1];
sx q[1];
rz(-1.8558985) q[1];
sx q[1];
rz(-0.76304719) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9656455) q[3];
sx q[3];
rz(-2.8215652) q[3];
sx q[3];
rz(1.4744028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12157089) q[2];
sx q[2];
rz(-2.3131804) q[2];
sx q[2];
rz(1.4309481) q[2];
rz(2.6089148) q[3];
sx q[3];
rz(-1.0360274) q[3];
sx q[3];
rz(1.7238341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1581887) q[0];
sx q[0];
rz(-0.15022763) q[0];
sx q[0];
rz(-2.0576553) q[0];
rz(0.051941959) q[1];
sx q[1];
rz(-2.7333994) q[1];
sx q[1];
rz(-3.1165677) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0753703) q[0];
sx q[0];
rz(-2.0522723) q[0];
sx q[0];
rz(-1.5388778) q[0];
rz(2.3228753) q[2];
sx q[2];
rz(-1.229964) q[2];
sx q[2];
rz(0.2073632) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6083303) q[1];
sx q[1];
rz(-0.46583781) q[1];
sx q[1];
rz(-1.9598258) q[1];
rz(0.30625108) q[3];
sx q[3];
rz(-0.44455498) q[3];
sx q[3];
rz(-1.2769092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0742801) q[2];
sx q[2];
rz(-1.68579) q[2];
sx q[2];
rz(2.32453) q[2];
rz(0.95064154) q[3];
sx q[3];
rz(-2.1248415) q[3];
sx q[3];
rz(-1.954621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10512146) q[0];
sx q[0];
rz(-2.8558185) q[0];
sx q[0];
rz(2.5407319) q[0];
rz(-3.1071682) q[1];
sx q[1];
rz(-1.8079146) q[1];
sx q[1];
rz(-1.5404125) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2297771) q[0];
sx q[0];
rz(-2.2719349) q[0];
sx q[0];
rz(-2.3100353) q[0];
rz(-pi) q[1];
rz(-0.67229597) q[2];
sx q[2];
rz(-1.0611254) q[2];
sx q[2];
rz(-0.59802645) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7936777) q[1];
sx q[1];
rz(-1.131212) q[1];
sx q[1];
rz(1.5949834) q[1];
rz(2.595946) q[3];
sx q[3];
rz(-1.9522523) q[3];
sx q[3];
rz(1.9617401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8704661) q[2];
sx q[2];
rz(-1.2721964) q[2];
sx q[2];
rz(1.4938483) q[2];
rz(-2.3624524) q[3];
sx q[3];
rz(-1.4804163) q[3];
sx q[3];
rz(-0.65472764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7036024) q[0];
sx q[0];
rz(-1.0810738) q[0];
sx q[0];
rz(-0.78824747) q[0];
rz(1.8986656) q[1];
sx q[1];
rz(-1.582076) q[1];
sx q[1];
rz(-16/(15*pi)) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96974715) q[0];
sx q[0];
rz(-1.7632428) q[0];
sx q[0];
rz(-3.0234973) q[0];
rz(1.8571157) q[2];
sx q[2];
rz(-2.114388) q[2];
sx q[2];
rz(1.7329777) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0362944) q[1];
sx q[1];
rz(-2.621161) q[1];
sx q[1];
rz(2.4776513) q[1];
x q[2];
rz(2.2195283) q[3];
sx q[3];
rz(-2.2345103) q[3];
sx q[3];
rz(-1.0459378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1303611) q[2];
sx q[2];
rz(-1.9731382) q[2];
sx q[2];
rz(1.2569024) q[2];
rz(2.0508749) q[3];
sx q[3];
rz(-1.2033477) q[3];
sx q[3];
rz(2.2244661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.286769) q[0];
sx q[0];
rz(-1.6463065) q[0];
sx q[0];
rz(-2.0538034) q[0];
rz(-0.65025672) q[1];
sx q[1];
rz(-1.4573263) q[1];
sx q[1];
rz(0.021171721) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084822239) q[0];
sx q[0];
rz(-0.49898832) q[0];
sx q[0];
rz(2.694724) q[0];
rz(-pi) q[1];
rz(-0.4707321) q[2];
sx q[2];
rz(-1.7463049) q[2];
sx q[2];
rz(2.144475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0291527) q[1];
sx q[1];
rz(-0.8275367) q[1];
sx q[1];
rz(-2.854101) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.043280799) q[3];
sx q[3];
rz(-1.0387522) q[3];
sx q[3];
rz(-2.8325641) q[3];
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
rz(-2.148597) q[3];
sx q[3];
rz(-2.0296622) q[3];
sx q[3];
rz(-2.907311) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2642333) q[0];
sx q[0];
rz(-2.3491884) q[0];
sx q[0];
rz(1.0607251) q[0];
rz(2.0296774) q[1];
sx q[1];
rz(-2.8363375) q[1];
sx q[1];
rz(-3.0189966) q[1];
rz(-0.23774679) q[2];
sx q[2];
rz(-1.100111) q[2];
sx q[2];
rz(-3.0636655) q[2];
rz(-0.018914472) q[3];
sx q[3];
rz(-2.6185642) q[3];
sx q[3];
rz(-2.0279283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
