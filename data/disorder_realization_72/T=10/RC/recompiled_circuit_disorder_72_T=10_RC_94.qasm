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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9141657) q[0];
sx q[0];
rz(-1.6343071) q[0];
sx q[0];
rz(-0.27557296) q[0];
x q[1];
rz(0.17272858) q[2];
sx q[2];
rz(-2.2842801) q[2];
sx q[2];
rz(1.8634862) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62568134) q[1];
sx q[1];
rz(-1.7209007) q[1];
sx q[1];
rz(-0.73966571) q[1];
rz(-pi) q[2];
x q[2];
rz(0.026704549) q[3];
sx q[3];
rz(-1.4674125) q[3];
sx q[3];
rz(-1.0552693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9177861) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(-1.0144368) q[2];
rz(0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(-0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.4734128) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(-1.0027592) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6013872) q[0];
sx q[0];
rz(-1.1367072) q[0];
sx q[0];
rz(-0.65625221) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9687256) q[2];
sx q[2];
rz(-1.7555408) q[2];
sx q[2];
rz(0.52344054) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7111295) q[1];
sx q[1];
rz(-1.7257479) q[1];
sx q[1];
rz(0.98053812) q[1];
rz(0.59870436) q[3];
sx q[3];
rz(-1.8667392) q[3];
sx q[3];
rz(0.89950365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8721547) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(-0.41444591) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(-0.088236563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.2484444) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(2.3625968) q[0];
rz(-0.39930725) q[1];
sx q[1];
rz(-1.2483968) q[1];
sx q[1];
rz(0.88358203) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012264) q[0];
sx q[0];
rz(-1.7352312) q[0];
sx q[0];
rz(0.40584392) q[0];
rz(0.29404624) q[2];
sx q[2];
rz(-1.9532734) q[2];
sx q[2];
rz(1.2750212) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2391501) q[1];
sx q[1];
rz(-0.28581866) q[1];
sx q[1];
rz(-1.6509389) q[1];
x q[2];
rz(-1.6885353) q[3];
sx q[3];
rz(-2.2265834) q[3];
sx q[3];
rz(1.1586787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3016004) q[2];
sx q[2];
rz(-1.8847382) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(-0.99003506) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(-2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4317959) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(2.0344095) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(-2.591419) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69458354) q[0];
sx q[0];
rz(-0.83568474) q[0];
sx q[0];
rz(1.4674594) q[0];
x q[1];
rz(1.6646531) q[2];
sx q[2];
rz(-1.8641194) q[2];
sx q[2];
rz(2.9039127) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.354047) q[1];
sx q[1];
rz(-2.5552632) q[1];
sx q[1];
rz(1.7802618) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21869603) q[3];
sx q[3];
rz(-0.23467206) q[3];
sx q[3];
rz(-2.4179539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6198373) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(-0.27553976) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(-0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4438641) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(-0.69119167) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(-0.91526389) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0147858) q[0];
sx q[0];
rz(-0.94313398) q[0];
sx q[0];
rz(-0.058237596) q[0];
x q[1];
rz(-1.5259597) q[2];
sx q[2];
rz(-2.0955288) q[2];
sx q[2];
rz(-2.5068138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3770959) q[1];
sx q[1];
rz(-0.15863523) q[1];
sx q[1];
rz(-1.7736969) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3239273) q[3];
sx q[3];
rz(-2.2125803) q[3];
sx q[3];
rz(1.6587917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-2.129014) q[2];
sx q[2];
rz(-2.6780224) q[2];
rz(-0.56435895) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(-2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0267462) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(-2.0625431) q[0];
rz(-0.59533978) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(-0.39658305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30124861) q[0];
sx q[0];
rz(-2.6123971) q[0];
sx q[0];
rz(2.089414) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55690342) q[2];
sx q[2];
rz(-0.54585005) q[2];
sx q[2];
rz(2.5715373) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.086416883) q[1];
sx q[1];
rz(-0.91311087) q[1];
sx q[1];
rz(3.0228826) q[1];
rz(-pi) q[2];
rz(2.7820884) q[3];
sx q[3];
rz(-0.61509575) q[3];
sx q[3];
rz(2.2744327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1534999) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(-1.5863824) q[2];
rz(1.4554626) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(-0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7431188) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(2.356785) q[0];
rz(-1.2706884) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(-2.0152337) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61486812) q[0];
sx q[0];
rz(-1.4818026) q[0];
sx q[0];
rz(0.29863775) q[0];
x q[1];
rz(1.0625213) q[2];
sx q[2];
rz(-1.4462496) q[2];
sx q[2];
rz(1.5073656) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3583402) q[1];
sx q[1];
rz(-1.1249152) q[1];
sx q[1];
rz(-1.7794442) q[1];
rz(-pi) q[2];
rz(1.8282312) q[3];
sx q[3];
rz(-2.7830887) q[3];
sx q[3];
rz(-0.076971019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.209098) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(-1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(-0.11288189) q[0];
rz(-1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(3.1088366) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8772802) q[0];
sx q[0];
rz(-0.41752975) q[0];
sx q[0];
rz(0.89950048) q[0];
x q[1];
rz(-2.1342437) q[2];
sx q[2];
rz(-1.1903561) q[2];
sx q[2];
rz(0.13605875) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2331088) q[1];
sx q[1];
rz(-1.1540124) q[1];
sx q[1];
rz(-3.0770739) q[1];
rz(1.7912346) q[3];
sx q[3];
rz(-2.7535097) q[3];
sx q[3];
rz(0.9006075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1774566) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(-2.5194871) q[2];
rz(1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0416097) q[0];
sx q[0];
rz(-2.6384625) q[0];
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
rz(-0.14478806) q[0];
sx q[0];
rz(-2.7607714) q[0];
rz(-pi) q[1];
x q[1];
rz(2.910026) q[2];
sx q[2];
rz(-0.61083691) q[2];
sx q[2];
rz(-1.4996741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.24385142) q[1];
sx q[1];
rz(-0.85695367) q[1];
sx q[1];
rz(-2.7276911) q[1];
rz(-pi) q[2];
rz(2.5599307) q[3];
sx q[3];
rz(-1.186944) q[3];
sx q[3];
rz(-1.7285085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0316524) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(-2.9821441) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6195246) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(-1.7074701) q[0];
rz(1.8822949) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(-2.5433345) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3132968) q[0];
sx q[0];
rz(-1.7772355) q[0];
sx q[0];
rz(0.22160991) q[0];
x q[1];
rz(-2.9293289) q[2];
sx q[2];
rz(-1.4258254) q[2];
sx q[2];
rz(-3.0002909) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64893374) q[1];
sx q[1];
rz(-1.4760541) q[1];
sx q[1];
rz(1.9365063) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15828295) q[3];
sx q[3];
rz(-1.7332819) q[3];
sx q[3];
rz(2.3815001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.082211994) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(2.5478798) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.9700206) q[2];
sx q[2];
rz(-2.5149859) q[2];
sx q[2];
rz(1.7563216) q[2];
rz(-0.39246172) q[3];
sx q[3];
rz(-0.93467181) q[3];
sx q[3];
rz(-1.296464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
