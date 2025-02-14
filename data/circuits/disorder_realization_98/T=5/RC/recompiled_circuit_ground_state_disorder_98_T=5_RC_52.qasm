OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7336361) q[0];
sx q[0];
rz(-0.1853369) q[0];
sx q[0];
rz(-1.6428525) q[0];
rz(-0.19417956) q[1];
sx q[1];
rz(-0.3436389) q[1];
sx q[1];
rz(-0.37771168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35750264) q[0];
sx q[0];
rz(-1.8121077) q[0];
sx q[0];
rz(1.7989798) q[0];
rz(-pi) q[1];
rz(0.80773662) q[2];
sx q[2];
rz(-2.5662072) q[2];
sx q[2];
rz(-2.6859716) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1520253) q[1];
sx q[1];
rz(-1.5366239) q[1];
sx q[1];
rz(-2.4844643) q[1];
rz(-pi) q[2];
rz(1.8355719) q[3];
sx q[3];
rz(-1.5989948) q[3];
sx q[3];
rz(-1.4169852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7517884) q[2];
sx q[2];
rz(-1.3499539) q[2];
sx q[2];
rz(0.26503116) q[2];
rz(-0.42936471) q[3];
sx q[3];
rz(-1.8931484) q[3];
sx q[3];
rz(-3.140669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65381831) q[0];
sx q[0];
rz(-0.14350292) q[0];
sx q[0];
rz(1.300746) q[0];
rz(-0.067151345) q[1];
sx q[1];
rz(-0.90922272) q[1];
sx q[1];
rz(-0.86004177) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3305682) q[0];
sx q[0];
rz(-2.1386792) q[0];
sx q[0];
rz(-0.83356838) q[0];
x q[1];
rz(2.4323787) q[2];
sx q[2];
rz(-2.1265891) q[2];
sx q[2];
rz(-1.9947987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7262267) q[1];
sx q[1];
rz(-0.72390717) q[1];
sx q[1];
rz(-0.93127802) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1116756) q[3];
sx q[3];
rz(-2.7949998) q[3];
sx q[3];
rz(2.5427713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.63090515) q[2];
sx q[2];
rz(-2.5248933) q[2];
sx q[2];
rz(-0.56646937) q[2];
rz(2.412292) q[3];
sx q[3];
rz(-0.84435487) q[3];
sx q[3];
rz(-2.9384889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1931964) q[0];
sx q[0];
rz(-1.0862792) q[0];
sx q[0];
rz(-0.36619827) q[0];
rz(-1.5858448) q[1];
sx q[1];
rz(-1.0428753) q[1];
sx q[1];
rz(-0.050447024) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90441762) q[0];
sx q[0];
rz(-1.5290998) q[0];
sx q[0];
rz(1.4916184) q[0];
rz(1.0773727) q[2];
sx q[2];
rz(-0.9443379) q[2];
sx q[2];
rz(2.3979417) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5844684) q[1];
sx q[1];
rz(-1.7377751) q[1];
sx q[1];
rz(-2.9191769) q[1];
x q[2];
rz(1.8837711) q[3];
sx q[3];
rz(-0.44383263) q[3];
sx q[3];
rz(-2.7601506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0744276) q[2];
sx q[2];
rz(-1.0600435) q[2];
sx q[2];
rz(1.3107497) q[2];
rz(-1.7835167) q[3];
sx q[3];
rz(-0.57099968) q[3];
sx q[3];
rz(1.3091492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.8686789) q[0];
sx q[0];
rz(-2.9153115) q[0];
sx q[0];
rz(-1.7568461) q[0];
rz(1.2387431) q[1];
sx q[1];
rz(-1.9107995) q[1];
sx q[1];
rz(1.967427) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096701972) q[0];
sx q[0];
rz(-0.47098038) q[0];
sx q[0];
rz(1.9975918) q[0];
rz(-1.8205244) q[2];
sx q[2];
rz(-1.3018697) q[2];
sx q[2];
rz(-0.099322546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0020743) q[1];
sx q[1];
rz(-1.0660831) q[1];
sx q[1];
rz(-2.7667749) q[1];
rz(2.2635095) q[3];
sx q[3];
rz(-1.2125848) q[3];
sx q[3];
rz(-0.75853759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7063286) q[2];
sx q[2];
rz(-1.3202983) q[2];
sx q[2];
rz(-0.038979385) q[2];
rz(-0.6428166) q[3];
sx q[3];
rz(-2.1233852) q[3];
sx q[3];
rz(-2.5207998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83127999) q[0];
sx q[0];
rz(-1.0142925) q[0];
sx q[0];
rz(0.020462791) q[0];
rz(-2.9488355) q[1];
sx q[1];
rz(-0.61734504) q[1];
sx q[1];
rz(-1.1654759) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7503238) q[0];
sx q[0];
rz(-1.170715) q[0];
sx q[0];
rz(0.66949908) q[0];
rz(-pi) q[1];
rz(-2.4679409) q[2];
sx q[2];
rz(-1.9314226) q[2];
sx q[2];
rz(-1.1471105) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.36042696) q[1];
sx q[1];
rz(-1.5910223) q[1];
sx q[1];
rz(-2.0745049) q[1];
rz(-pi) q[2];
rz(-1.3806512) q[3];
sx q[3];
rz(-1.4961637) q[3];
sx q[3];
rz(0.41342218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.3141025) q[2];
sx q[2];
rz(-0.9919439) q[2];
sx q[2];
rz(0.52725434) q[2];
rz(2.5984247) q[3];
sx q[3];
rz(-2.4542464) q[3];
sx q[3];
rz(-2.7988722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0813893) q[0];
sx q[0];
rz(-0.15601604) q[0];
sx q[0];
rz(0.40670893) q[0];
rz(1.9806865) q[1];
sx q[1];
rz(-0.91861594) q[1];
sx q[1];
rz(2.1312174) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.062101) q[0];
sx q[0];
rz(-1.6915503) q[0];
sx q[0];
rz(1.0213721) q[0];
rz(-pi) q[1];
rz(-0.40596227) q[2];
sx q[2];
rz(-0.34892198) q[2];
sx q[2];
rz(0.59970784) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6076489) q[1];
sx q[1];
rz(-1.2022809) q[1];
sx q[1];
rz(-1.3435958) q[1];
rz(-pi) q[2];
rz(0.25185092) q[3];
sx q[3];
rz(-2.3201019) q[3];
sx q[3];
rz(-1.7918394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0668209) q[2];
sx q[2];
rz(-1.8596884) q[2];
sx q[2];
rz(-2.1066966) q[2];
rz(-1.7847938) q[3];
sx q[3];
rz(-0.59485888) q[3];
sx q[3];
rz(-2.518173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9676301) q[0];
sx q[0];
rz(-1.6169463) q[0];
sx q[0];
rz(0.3279283) q[0];
rz(3.0294042) q[1];
sx q[1];
rz(-1.2022377) q[1];
sx q[1];
rz(-0.97253886) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93720651) q[0];
sx q[0];
rz(-0.68336672) q[0];
sx q[0];
rz(2.1478189) q[0];
rz(-pi) q[1];
rz(1.1292636) q[2];
sx q[2];
rz(-1.2988833) q[2];
sx q[2];
rz(-0.11833469) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9303846) q[1];
sx q[1];
rz(-0.58961419) q[1];
sx q[1];
rz(-3.1077887) q[1];
rz(-pi) q[2];
rz(2.5150033) q[3];
sx q[3];
rz(-1.2703151) q[3];
sx q[3];
rz(1.8412875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.605725) q[2];
sx q[2];
rz(-2.2650227) q[2];
sx q[2];
rz(-2.4884339) q[2];
rz(-0.43290916) q[3];
sx q[3];
rz(-0.50986367) q[3];
sx q[3];
rz(-2.9292817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2046278) q[0];
sx q[0];
rz(-2.300394) q[0];
sx q[0];
rz(0.2247819) q[0];
rz(0.79554355) q[1];
sx q[1];
rz(-1.3139775) q[1];
sx q[1];
rz(2.0733817) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8248937) q[0];
sx q[0];
rz(-2.2980656) q[0];
sx q[0];
rz(0.69490142) q[0];
x q[1];
rz(-2.178763) q[2];
sx q[2];
rz(-0.89917937) q[2];
sx q[2];
rz(1.2149908) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28212122) q[1];
sx q[1];
rz(-2.2582114) q[1];
sx q[1];
rz(1.6038546) q[1];
rz(-1.3146888) q[3];
sx q[3];
rz(-0.22897069) q[3];
sx q[3];
rz(-0.010502149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3855359) q[2];
sx q[2];
rz(-1.7588561) q[2];
sx q[2];
rz(-0.63560152) q[2];
rz(0.33291891) q[3];
sx q[3];
rz(-1.0115441) q[3];
sx q[3];
rz(0.46271589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3391089) q[0];
sx q[0];
rz(-3.1153296) q[0];
sx q[0];
rz(0.12301692) q[0];
rz(1.3975551) q[1];
sx q[1];
rz(-1.3879644) q[1];
sx q[1];
rz(0.77973286) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5132039) q[0];
sx q[0];
rz(-1.3003948) q[0];
sx q[0];
rz(2.2487075) q[0];
rz(-pi) q[1];
rz(2.9626368) q[2];
sx q[2];
rz(-0.6935941) q[2];
sx q[2];
rz(-0.041698448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.349682) q[1];
sx q[1];
rz(-1.1983219) q[1];
sx q[1];
rz(0.21362408) q[1];
rz(-pi) q[2];
rz(-1.443583) q[3];
sx q[3];
rz(-1.4641342) q[3];
sx q[3];
rz(-0.17388177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66703779) q[2];
sx q[2];
rz(-1.8755308) q[2];
sx q[2];
rz(0.090465821) q[2];
rz(2.7247834) q[3];
sx q[3];
rz(-0.79206812) q[3];
sx q[3];
rz(-0.99378234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-1.6592634) q[0];
sx q[0];
rz(-1.8195131) q[0];
sx q[0];
rz(-3.0521159) q[0];
rz(-0.80884519) q[1];
sx q[1];
rz(-0.50061148) q[1];
sx q[1];
rz(1.0577177) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5070478) q[0];
sx q[0];
rz(-1.9814361) q[0];
sx q[0];
rz(1.9837186) q[0];
rz(-pi) q[1];
rz(1.5799149) q[2];
sx q[2];
rz(-2.5895725) q[2];
sx q[2];
rz(0.43285757) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.48185086) q[1];
sx q[1];
rz(-0.54164499) q[1];
sx q[1];
rz(-2.2996344) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.026151733) q[3];
sx q[3];
rz(-0.75687486) q[3];
sx q[3];
rz(2.7084086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7384501) q[2];
sx q[2];
rz(-0.53692997) q[2];
sx q[2];
rz(-0.37330791) q[2];
rz(-0.25660723) q[3];
sx q[3];
rz(-1.4213057) q[3];
sx q[3];
rz(0.074450113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8858717) q[0];
sx q[0];
rz(-1.5626361) q[0];
sx q[0];
rz(1.7241021) q[0];
rz(-1.7120842) q[1];
sx q[1];
rz(-0.34965873) q[1];
sx q[1];
rz(0.8082334) q[1];
rz(-1.5803554) q[2];
sx q[2];
rz(-1.111996) q[2];
sx q[2];
rz(-1.5954191) q[2];
rz(3.1028845) q[3];
sx q[3];
rz(-0.9699655) q[3];
sx q[3];
rz(-1.7326596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
