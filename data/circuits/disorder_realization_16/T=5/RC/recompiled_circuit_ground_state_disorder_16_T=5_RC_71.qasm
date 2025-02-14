OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0082173) q[0];
sx q[0];
rz(5.0157265) q[0];
sx q[0];
rz(9.4118543) q[0];
rz(-2.456993) q[1];
sx q[1];
rz(-0.79889387) q[1];
sx q[1];
rz(-1.0577143) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2051281) q[0];
sx q[0];
rz(-1.797344) q[0];
sx q[0];
rz(2.4656057) q[0];
rz(-pi) q[1];
rz(-1.0619668) q[2];
sx q[2];
rz(-0.86039174) q[2];
sx q[2];
rz(1.8343385) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.94016852) q[1];
sx q[1];
rz(-0.88646171) q[1];
sx q[1];
rz(-2.0444872) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4679883) q[3];
sx q[3];
rz(-2.7479726) q[3];
sx q[3];
rz(3.032544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5554123) q[2];
sx q[2];
rz(-1.5427898) q[2];
sx q[2];
rz(3.0482698) q[2];
rz(-0.020545067) q[3];
sx q[3];
rz(-13/(16*pi)) q[3];
sx q[3];
rz(1.7790022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071863197) q[0];
sx q[0];
rz(-1.3988031) q[0];
sx q[0];
rz(0.82759696) q[0];
rz(-0.96356511) q[1];
sx q[1];
rz(-1.5255442) q[1];
sx q[1];
rz(2.4049984) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85918713) q[0];
sx q[0];
rz(-1.6466337) q[0];
sx q[0];
rz(-1.8453159) q[0];
x q[1];
rz(-2.170855) q[2];
sx q[2];
rz(-1.0781556) q[2];
sx q[2];
rz(-0.095183177) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.40562661) q[1];
sx q[1];
rz(-2.7752004) q[1];
sx q[1];
rz(-0.97586164) q[1];
rz(-pi) q[2];
rz(-2.9773744) q[3];
sx q[3];
rz(-0.13623691) q[3];
sx q[3];
rz(-0.0020023684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57614148) q[2];
sx q[2];
rz(-0.57500035) q[2];
sx q[2];
rz(-2.3973993) q[2];
rz(2.5326676) q[3];
sx q[3];
rz(-0.78151339) q[3];
sx q[3];
rz(-1.6412546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7805444) q[0];
sx q[0];
rz(-2.2760976) q[0];
sx q[0];
rz(-1.3737099) q[0];
rz(-0.79958493) q[1];
sx q[1];
rz(-2.1345963) q[1];
sx q[1];
rz(2.0094357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4904259) q[0];
sx q[0];
rz(-1.4803783) q[0];
sx q[0];
rz(-1.3549442) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5485053) q[2];
sx q[2];
rz(-1.4906851) q[2];
sx q[2];
rz(1.3673283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.76946041) q[1];
sx q[1];
rz(-2.4677708) q[1];
sx q[1];
rz(-3.0173183) q[1];
x q[2];
rz(0.062643073) q[3];
sx q[3];
rz(-2.5317305) q[3];
sx q[3];
rz(-0.26882879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.75480294) q[2];
sx q[2];
rz(-0.58480442) q[2];
sx q[2];
rz(1.7542138) q[2];
rz(-2.7925708) q[3];
sx q[3];
rz(-1.6882221) q[3];
sx q[3];
rz(-2.6121228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.741852) q[0];
sx q[0];
rz(-2.1735503) q[0];
sx q[0];
rz(-2.7834748) q[0];
rz(2.070836) q[1];
sx q[1];
rz(-0.64859575) q[1];
sx q[1];
rz(1.7020285) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1308129) q[0];
sx q[0];
rz(-0.81310105) q[0];
sx q[0];
rz(1.2050425) q[0];
rz(-pi) q[1];
rz(2.2622243) q[2];
sx q[2];
rz(-1.9676529) q[2];
sx q[2];
rz(-0.59890998) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.52709748) q[1];
sx q[1];
rz(-2.1367837) q[1];
sx q[1];
rz(0.66733349) q[1];
rz(-0.26243383) q[3];
sx q[3];
rz(-2.35733) q[3];
sx q[3];
rz(-1.5465496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7122571) q[2];
sx q[2];
rz(-1.2306932) q[2];
sx q[2];
rz(0.62475359) q[2];
rz(-1.5077) q[3];
sx q[3];
rz(-0.75993901) q[3];
sx q[3];
rz(1.1225351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7655012) q[0];
sx q[0];
rz(-2.2608345) q[0];
sx q[0];
rz(-1.1464024) q[0];
rz(1.435185) q[1];
sx q[1];
rz(-0.51992661) q[1];
sx q[1];
rz(2.3211839) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4389184) q[0];
sx q[0];
rz(-1.6587388) q[0];
sx q[0];
rz(1.2984718) q[0];
x q[1];
rz(1.7960288) q[2];
sx q[2];
rz(-0.67906717) q[2];
sx q[2];
rz(-2.1175543) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31215224) q[1];
sx q[1];
rz(-1.7746801) q[1];
sx q[1];
rz(-1.3711434) q[1];
rz(-pi) q[2];
rz(0.049656258) q[3];
sx q[3];
rz(-2.5512902) q[3];
sx q[3];
rz(1.4242537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4121805) q[2];
sx q[2];
rz(-2.086144) q[2];
sx q[2];
rz(0.5385651) q[2];
rz(-1.4696848) q[3];
sx q[3];
rz(-1.6939751) q[3];
sx q[3];
rz(3.0830429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84022123) q[0];
sx q[0];
rz(-3.0026307) q[0];
sx q[0];
rz(-0.090601623) q[0];
rz(2.7719356) q[1];
sx q[1];
rz(-1.548111) q[1];
sx q[1];
rz(2.7634117) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3017186) q[0];
sx q[0];
rz(-1.0953971) q[0];
sx q[0];
rz(1.8283707) q[0];
x q[1];
rz(-3.0514977) q[2];
sx q[2];
rz(-1.6125154) q[2];
sx q[2];
rz(0.99064529) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82434618) q[1];
sx q[1];
rz(-0.50889824) q[1];
sx q[1];
rz(-2.8864278) q[1];
rz(-3.0621594) q[3];
sx q[3];
rz(-1.380668) q[3];
sx q[3];
rz(0.058323764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52046627) q[2];
sx q[2];
rz(-1.4755321) q[2];
sx q[2];
rz(-3.0401518) q[2];
rz(-1.8488047) q[3];
sx q[3];
rz(-0.83270508) q[3];
sx q[3];
rz(1.4147991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8349649) q[0];
sx q[0];
rz(-1.2766301) q[0];
sx q[0];
rz(2.2391338) q[0];
rz(-0.2392256) q[1];
sx q[1];
rz(-1.3419515) q[1];
sx q[1];
rz(-0.36144027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28736605) q[0];
sx q[0];
rz(-2.0660095) q[0];
sx q[0];
rz(0.12676858) q[0];
x q[1];
rz(-2.7290384) q[2];
sx q[2];
rz(-3.0027886) q[2];
sx q[2];
rz(-0.96497646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16365227) q[1];
sx q[1];
rz(-1.1457902) q[1];
sx q[1];
rz(0.24042701) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42487259) q[3];
sx q[3];
rz(-0.36965224) q[3];
sx q[3];
rz(1.2583789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0685737) q[2];
sx q[2];
rz(-1.8014427) q[2];
sx q[2];
rz(-0.87693357) q[2];
rz(0.98313156) q[3];
sx q[3];
rz(-1.5777595) q[3];
sx q[3];
rz(0.94299281) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8125732) q[0];
sx q[0];
rz(-0.1952157) q[0];
sx q[0];
rz(0.83874291) q[0];
rz(-0.70873952) q[1];
sx q[1];
rz(-0.63116169) q[1];
sx q[1];
rz(-2.2355524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4562389) q[0];
sx q[0];
rz(-0.8461915) q[0];
sx q[0];
rz(-1.5426209) q[0];
x q[1];
rz(-2.9259813) q[2];
sx q[2];
rz(-1.9136568) q[2];
sx q[2];
rz(2.947399) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3248277) q[1];
sx q[1];
rz(-0.93597368) q[1];
sx q[1];
rz(0.88714182) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25211199) q[3];
sx q[3];
rz(-2.1641755) q[3];
sx q[3];
rz(-0.93170792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16443843) q[2];
sx q[2];
rz(-2.437037) q[2];
sx q[2];
rz(-2.0257115) q[2];
rz(2.5130533) q[3];
sx q[3];
rz(-1.081859) q[3];
sx q[3];
rz(0.61454296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81166613) q[0];
sx q[0];
rz(-2.3217432) q[0];
sx q[0];
rz(-0.16127583) q[0];
rz(-0.8423841) q[1];
sx q[1];
rz(-1.4295652) q[1];
sx q[1];
rz(-0.17002034) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-7.13037e-05) q[0];
sx q[0];
rz(-1.562644) q[0];
sx q[0];
rz(-3.1233112) q[0];
rz(-pi) q[1];
rz(0.34509322) q[2];
sx q[2];
rz(-2.6775914) q[2];
sx q[2];
rz(-2.4008) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0898897) q[1];
sx q[1];
rz(-1.1384299) q[1];
sx q[1];
rz(-0.21413762) q[1];
rz(-pi) q[2];
rz(0.028548553) q[3];
sx q[3];
rz(-1.953509) q[3];
sx q[3];
rz(0.0023502758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9876447) q[2];
sx q[2];
rz(-1.5445856) q[2];
sx q[2];
rz(-1.4863996) q[2];
rz(0.060128309) q[3];
sx q[3];
rz(-0.82258737) q[3];
sx q[3];
rz(1.7000343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5104367) q[0];
sx q[0];
rz(-3.1102409) q[0];
sx q[0];
rz(1.4309058) q[0];
rz(0.36262861) q[1];
sx q[1];
rz(-1.5321833) q[1];
sx q[1];
rz(-1.3154202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71821868) q[0];
sx q[0];
rz(-0.43180433) q[0];
sx q[0];
rz(0.81599094) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8706546) q[2];
sx q[2];
rz(-0.39405566) q[2];
sx q[2];
rz(-1.9410417) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3676995) q[1];
sx q[1];
rz(-0.65441346) q[1];
sx q[1];
rz(0.44632895) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0995448) q[3];
sx q[3];
rz(-2.3662851) q[3];
sx q[3];
rz(-0.14107832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.057283904) q[2];
sx q[2];
rz(-1.5466651) q[2];
sx q[2];
rz(-2.9810737) q[2];
rz(0.48216835) q[3];
sx q[3];
rz(-2.4653258) q[3];
sx q[3];
rz(-2.4619861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.01874825) q[0];
sx q[0];
rz(-1.7457122) q[0];
sx q[0];
rz(2.2298298) q[0];
rz(-1.1625166) q[1];
sx q[1];
rz(-1.5717506) q[1];
sx q[1];
rz(2.7350978) q[1];
rz(1.4760426) q[2];
sx q[2];
rz(-1.6577402) q[2];
sx q[2];
rz(0.18862373) q[2];
rz(-2.1382469) q[3];
sx q[3];
rz(-2.2570845) q[3];
sx q[3];
rz(1.7751638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
