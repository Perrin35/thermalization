OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(1.5078478) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(2.642282) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2274269) q[0];
sx q[0];
rz(-1.6343071) q[0];
sx q[0];
rz(-2.8660197) q[0];
x q[1];
rz(2.9688641) q[2];
sx q[2];
rz(-2.2842801) q[2];
sx q[2];
rz(1.2781065) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5159113) q[1];
sx q[1];
rz(-1.7209007) q[1];
sx q[1];
rz(0.73966571) q[1];
rz(-pi) q[2];
rz(1.8226837) q[3];
sx q[3];
rz(-0.10676521) q[3];
sx q[3];
rz(1.3085384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9177861) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(2.1271558) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4734128) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(1.5197808) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(1.0027592) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013872) q[0];
sx q[0];
rz(-2.0048855) q[0];
sx q[0];
rz(-0.65625221) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3833212) q[2];
sx q[2];
rz(-1.4008998) q[2];
sx q[2];
rz(-1.0152917) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1042852) q[1];
sx q[1];
rz(-2.1530495) q[1];
sx q[1];
rz(2.9557455) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2172132) q[3];
sx q[3];
rz(-1.0014605) q[3];
sx q[3];
rz(-2.6667037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26943794) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(-2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2484444) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(2.3625968) q[0];
rz(-0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-0.88358203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6813864) q[0];
sx q[0];
rz(-1.1707414) q[0];
sx q[0];
rz(1.3921188) q[0];
x q[1];
rz(-1.1728889) q[2];
sx q[2];
rz(-1.298561) q[2];
sx q[2];
rz(2.9583601) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2391501) q[1];
sx q[1];
rz(-2.855774) q[1];
sx q[1];
rz(-1.4906537) q[1];
x q[2];
rz(0.15150841) q[3];
sx q[3];
rz(-0.66473367) q[3];
sx q[3];
rz(-0.96707771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.8847382) q[2];
sx q[2];
rz(2.7123614) q[2];
rz(0.99003506) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7097968) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(-2.0344095) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(2.591419) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69458354) q[0];
sx q[0];
rz(-0.83568474) q[0];
sx q[0];
rz(1.6741333) q[0];
rz(-pi) q[1];
rz(1.4769396) q[2];
sx q[2];
rz(-1.2774733) q[2];
sx q[2];
rz(2.9039127) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.95851129) q[1];
sx q[1];
rz(-1.6861048) q[1];
sx q[1];
rz(2.1469841) q[1];
rz(-pi) q[2];
rz(0.21869603) q[3];
sx q[3];
rz(-2.9069206) q[3];
sx q[3];
rz(-2.4179539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.52175534) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(-2.8660529) q[2];
rz(-3.0299305) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(-0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4438641) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-0.87669796) q[0];
rz(2.450401) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(2.2263288) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6633776) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(-2.1992654) q[0];
x q[1];
rz(-1.5259597) q[2];
sx q[2];
rz(-2.0955288) q[2];
sx q[2];
rz(-2.5068138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60587864) q[1];
sx q[1];
rz(-1.5389581) q[1];
sx q[1];
rz(1.4153626) q[1];
x q[2];
rz(-2.3239273) q[3];
sx q[3];
rz(-2.2125803) q[3];
sx q[3];
rz(-1.4828009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(-0.4635703) q[2];
rz(-0.56435895) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(1.0790496) q[0];
rz(0.59533978) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(2.7450096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2561589) q[0];
sx q[0];
rz(-2.0246756) q[0];
sx q[0];
rz(-0.28215779) q[0];
x q[1];
rz(-0.47607143) q[2];
sx q[2];
rz(-1.8487612) q[2];
sx q[2];
rz(-1.4897886) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10627667) q[1];
sx q[1];
rz(-2.4748487) q[1];
sx q[1];
rz(-1.4186526) q[1];
rz(-2.7820884) q[3];
sx q[3];
rz(-2.5264969) q[3];
sx q[3];
rz(-0.86715992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1534999) q[2];
sx q[2];
rz(-2.2160539) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(-0.78480762) q[0];
rz(1.2706884) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(1.126359) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2369909) q[0];
sx q[0];
rz(-0.31123529) q[0];
sx q[0];
rz(-2.8471332) q[0];
x q[1];
rz(2.9992505) q[2];
sx q[2];
rz(-2.0747613) q[2];
sx q[2];
rz(-3.1359283) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3583402) q[1];
sx q[1];
rz(-1.1249152) q[1];
sx q[1];
rz(-1.7794442) q[1];
rz(1.9184291) q[3];
sx q[3];
rz(-1.4813444) q[3];
sx q[3];
rz(1.4060494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.209098) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(2.9879925) q[2];
rz(-2.8364654) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37041935) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(-0.11288189) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(0.032756068) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93452867) q[0];
sx q[0];
rz(-1.8257739) q[0];
sx q[0];
rz(-1.9051139) q[0];
x q[1];
rz(2.6997386) q[2];
sx q[2];
rz(-1.0519069) q[2];
sx q[2];
rz(-1.9372802) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7777562) q[1];
sx q[1];
rz(-1.5118074) q[1];
sx q[1];
rz(-1.9883518) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7912346) q[3];
sx q[3];
rz(-0.38808295) q[3];
sx q[3];
rz(-0.9006075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1774566) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(-0.73293066) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(2.656235) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.582726) q[0];
sx q[0];
rz(-0.14478806) q[0];
sx q[0];
rz(-2.7607714) q[0];
x q[1];
rz(1.7301241) q[2];
sx q[2];
rz(-2.1630641) q[2];
sx q[2];
rz(1.2194022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7945054) q[1];
sx q[1];
rz(-0.80650389) q[1];
sx q[1];
rz(-2.0054714) q[1];
x q[2];
rz(0.63391179) q[3];
sx q[3];
rz(-0.68448193) q[3];
sx q[3];
rz(-2.4661635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1099403) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(-0.1594485) q[2];
rz(-1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.52206802) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(1.4341226) q[0];
rz(-1.2592978) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(2.5433345) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3132968) q[0];
sx q[0];
rz(-1.7772355) q[0];
sx q[0];
rz(2.9199827) q[0];
rz(-pi) q[1];
rz(-2.5355849) q[2];
sx q[2];
rz(-2.8851644) q[2];
sx q[2];
rz(-2.3026349) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.64893374) q[1];
sx q[1];
rz(-1.6655386) q[1];
sx q[1];
rz(-1.2050864) q[1];
rz(1.7353021) q[3];
sx q[3];
rz(-1.4146155) q[3];
sx q[3];
rz(2.3567049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.082211994) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(0.65199488) q[2];
rz(-2.5478798) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(-1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3020637) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(1.9051753) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(-2.5219953) q[2];
sx q[2];
rz(-1.6710812) q[2];
sx q[2];
rz(-2.8166213) q[2];
rz(1.0929575) q[3];
sx q[3];
rz(-2.408705) q[3];
sx q[3];
rz(-0.68791289) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
