OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5207198) q[0];
sx q[0];
rz(-1.7680661) q[0];
sx q[0];
rz(-1.5078478) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(2.642282) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7802785) q[0];
sx q[0];
rz(-1.8457992) q[0];
sx q[0];
rz(-1.5048024) q[0];
rz(-2.2917065) q[2];
sx q[2];
rz(-1.7011142) q[2];
sx q[2];
rz(-2.7352114) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62568134) q[1];
sx q[1];
rz(-1.4206919) q[1];
sx q[1];
rz(2.4019269) q[1];
rz(-pi) q[2];
rz(3.1148881) q[3];
sx q[3];
rz(-1.4674125) q[3];
sx q[3];
rz(-2.0863233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(1.0144368) q[2];
rz(2.9075918) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4734128) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(1.6218119) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(1.0027592) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5402055) q[0];
sx q[0];
rz(-1.1367072) q[0];
sx q[0];
rz(-2.4853404) q[0];
rz(-pi) q[1];
rz(-1.3833212) q[2];
sx q[2];
rz(-1.7406929) q[2];
sx q[2];
rz(1.0152917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1042852) q[1];
sx q[1];
rz(-0.98854317) q[1];
sx q[1];
rz(0.18584713) q[1];
rz(-0.49591222) q[3];
sx q[3];
rz(-2.4818588) q[3];
sx q[3];
rz(-2.0663313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.26943794) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(0.41444591) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(0.088236563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8931483) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(-0.77899581) q[0];
rz(0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-2.2580106) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8944091) q[0];
sx q[0];
rz(-2.7054225) q[0];
sx q[0];
rz(0.39788525) q[0];
rz(1.9687037) q[2];
sx q[2];
rz(-1.8430317) q[2];
sx q[2];
rz(-2.9583601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2391501) q[1];
sx q[1];
rz(-0.28581866) q[1];
sx q[1];
rz(1.6509389) q[1];
x q[2];
rz(2.9900842) q[3];
sx q[3];
rz(-0.66473367) q[3];
sx q[3];
rz(-2.1745149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.8847382) q[2];
sx q[2];
rz(-0.42923129) q[2];
rz(2.1515576) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(-2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7097968) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(2.0344095) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(2.591419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.600425) q[0];
sx q[0];
rz(-0.74099243) q[0];
sx q[0];
rz(3.0279972) q[0];
rz(-pi) q[1];
rz(0.29454622) q[2];
sx q[2];
rz(-1.4809594) q[2];
sx q[2];
rz(1.3603269) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78754567) q[1];
sx q[1];
rz(-0.58632942) q[1];
sx q[1];
rz(-1.7802618) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5189734) q[3];
sx q[3];
rz(-1.7997777) q[3];
sx q[3];
rz(2.1932972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6198373) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(-0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6977285) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(-0.87669796) q[0];
rz(2.450401) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(0.91526389) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1137431) q[0];
sx q[0];
rz(-2.5115974) q[0];
sx q[0];
rz(1.4907452) q[0];
rz(-pi) q[1];
rz(1.5259597) q[2];
sx q[2];
rz(-2.0955288) q[2];
sx q[2];
rz(2.5068138) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76449672) q[1];
sx q[1];
rz(-2.9829574) q[1];
sx q[1];
rz(-1.7736969) q[1];
rz(-0.74113412) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(-2.4853064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(2.6780224) q[2];
rz(-0.56435895) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(-2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-2.0625431) q[0];
rz(0.59533978) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(2.7450096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30124861) q[0];
sx q[0];
rz(-0.52919555) q[0];
sx q[0];
rz(1.0521786) q[0];
x q[1];
rz(2.6655212) q[2];
sx q[2];
rz(-1.2928315) q[2];
sx q[2];
rz(1.4897886) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.086416883) q[1];
sx q[1];
rz(-0.91311087) q[1];
sx q[1];
rz(-0.1187101) q[1];
rz(2.7820884) q[3];
sx q[3];
rz(-0.61509575) q[3];
sx q[3];
rz(-0.86715992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98809272) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(1.5863824) q[2];
rz(1.4554626) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(-2.356785) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(2.0152337) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2130177) q[0];
sx q[0];
rz(-1.2733766) q[0];
sx q[0];
rz(1.4777044) q[0];
x q[1];
rz(-1.8225841) q[2];
sx q[2];
rz(-2.6195824) q[2];
sx q[2];
rz(0.28282794) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.78325242) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(-1.3621484) q[1];
x q[2];
rz(1.9184291) q[3];
sx q[3];
rz(-1.6602483) q[3];
sx q[3];
rz(1.7355433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9324947) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(2.9879925) q[2];
rz(2.8364654) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(1.7677527) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(0.11288189) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(-3.1088366) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.207064) q[0];
sx q[0];
rz(-1.8257739) q[0];
sx q[0];
rz(-1.2364788) q[0];
x q[1];
rz(0.92808) q[2];
sx q[2];
rz(-0.66814458) q[2];
sx q[2];
rz(1.1754787) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0748464) q[1];
sx q[1];
rz(-0.42145887) q[1];
sx q[1];
rz(-1.4261817) q[1];
x q[2];
rz(-1.350358) q[3];
sx q[3];
rz(-2.7535097) q[3];
sx q[3];
rz(0.9006075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.96413606) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(0.62210554) q[2];
rz(1.1671676) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(2.7511403) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(1.5266248) q[0];
rz(-2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(-2.656235) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.582726) q[0];
sx q[0];
rz(-2.9968046) q[0];
sx q[0];
rz(2.7607714) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4114686) q[2];
sx q[2];
rz(-2.1630641) q[2];
sx q[2];
rz(-1.2194022) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7945054) q[1];
sx q[1];
rz(-2.3350888) q[1];
sx q[1];
rz(-1.1361213) q[1];
rz(2.5076809) q[3];
sx q[3];
rz(-0.68448193) q[3];
sx q[3];
rz(-0.67542911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0316524) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(0.1594485) q[2];
rz(-1.4032646) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195246) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(1.8822949) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(2.5433345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3529417) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.3593332) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7190476) q[2];
sx q[2];
rz(-1.3607927) q[2];
sx q[2];
rz(1.460618) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1640537) q[1];
sx q[1];
rz(-2.7643449) q[1];
sx q[1];
rz(-1.8305199) q[1];
rz(-pi) q[2];
rz(-0.15828295) q[3];
sx q[3];
rz(-1.4083107) q[3];
sx q[3];
rz(0.76009258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.082211994) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(2.4895978) q[2];
rz(0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(-1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.839529) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(-1.2364173) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(2.9700206) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(2.0486352) q[3];
sx q[3];
rz(-0.7328877) q[3];
sx q[3];
rz(2.4536798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
