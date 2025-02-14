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
rz(0.32956707) q[0];
sx q[0];
rz(0.94533935) q[0];
sx q[0];
rz(9.4234484) q[0];
rz(3.3450491) q[1];
sx q[1];
rz(3.3766881) q[1];
sx q[1];
rz(8.8535218) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0482727) q[0];
sx q[0];
rz(-2.0287345) q[0];
sx q[0];
rz(-2.4620374) q[0];
rz(2.4680092) q[2];
sx q[2];
rz(-2.9004708) q[2];
sx q[2];
rz(-0.6553638) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3498342) q[1];
sx q[1];
rz(-1.5877643) q[1];
sx q[1];
rz(-0.49032536) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68526973) q[3];
sx q[3];
rz(-0.56495404) q[3];
sx q[3];
rz(2.9034815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0433537) q[2];
sx q[2];
rz(-1.9184155) q[2];
sx q[2];
rz(2.5377048) q[2];
rz(1.0107001) q[3];
sx q[3];
rz(-2.5778008) q[3];
sx q[3];
rz(0.92196661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0589703) q[0];
sx q[0];
rz(-2.048546) q[0];
sx q[0];
rz(-3.1329204) q[0];
rz(1.1588833) q[1];
sx q[1];
rz(-2.4544139) q[1];
sx q[1];
rz(-3.029356) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2896947) q[0];
sx q[0];
rz(-1.0852975) q[0];
sx q[0];
rz(-1.8957036) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5331763) q[2];
sx q[2];
rz(-2.2171187) q[2];
sx q[2];
rz(-0.96137709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1232154) q[1];
sx q[1];
rz(-1.0800414) q[1];
sx q[1];
rz(1.823455) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26437624) q[3];
sx q[3];
rz(-1.8473139) q[3];
sx q[3];
rz(0.087527601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9021641) q[2];
sx q[2];
rz(-1.8205234) q[2];
sx q[2];
rz(0.93977896) q[2];
rz(-2.2043665) q[3];
sx q[3];
rz(-1.7669433) q[3];
sx q[3];
rz(0.37041131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9409907) q[0];
sx q[0];
rz(-0.90018278) q[0];
sx q[0];
rz(-0.58315939) q[0];
rz(-2.0896301) q[1];
sx q[1];
rz(-2.5556892) q[1];
sx q[1];
rz(0.016187035) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72618851) q[0];
sx q[0];
rz(-2.4016018) q[0];
sx q[0];
rz(-1.9053024) q[0];
x q[1];
rz(1.8772241) q[2];
sx q[2];
rz(-1.2867905) q[2];
sx q[2];
rz(-0.71137448) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5213373) q[1];
sx q[1];
rz(-2.9456821) q[1];
sx q[1];
rz(-0.98448648) q[1];
rz(-pi) q[2];
rz(-2.2233769) q[3];
sx q[3];
rz(-0.86091061) q[3];
sx q[3];
rz(0.48274279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6215324) q[2];
sx q[2];
rz(-1.3935057) q[2];
sx q[2];
rz(3.0690466) q[2];
rz(1.0091311) q[3];
sx q[3];
rz(-2.3435209) q[3];
sx q[3];
rz(0.23536853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16044727) q[0];
sx q[0];
rz(-2.7837842) q[0];
sx q[0];
rz(-0.88448802) q[0];
rz(-1.5011939) q[1];
sx q[1];
rz(-1.6694992) q[1];
sx q[1];
rz(-0.59008682) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9035066) q[0];
sx q[0];
rz(-2.0683388) q[0];
sx q[0];
rz(0.4273129) q[0];
x q[1];
rz(0.22802148) q[2];
sx q[2];
rz(-1.3022026) q[2];
sx q[2];
rz(3.0326469) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2324625) q[1];
sx q[1];
rz(-1.2386444) q[1];
sx q[1];
rz(1.0494768) q[1];
rz(-0.43729525) q[3];
sx q[3];
rz(-1.4623433) q[3];
sx q[3];
rz(-0.57009283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0044452) q[2];
sx q[2];
rz(-1.659617) q[2];
sx q[2];
rz(1.8972634) q[2];
rz(1.7826049) q[3];
sx q[3];
rz(-2.6748896) q[3];
sx q[3];
rz(2.9984503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.8878079) q[0];
sx q[0];
rz(-2.4847327) q[0];
sx q[0];
rz(-0.3669056) q[0];
rz(-1.3430355) q[1];
sx q[1];
rz(-2.8470706) q[1];
sx q[1];
rz(-1.3142745) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2283233) q[0];
sx q[0];
rz(-1.1342888) q[0];
sx q[0];
rz(2.4455347) q[0];
x q[1];
rz(-2.6571006) q[2];
sx q[2];
rz(-1.3132261) q[2];
sx q[2];
rz(0.07312873) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6180743) q[1];
sx q[1];
rz(-0.95361751) q[1];
sx q[1];
rz(-2.7568222) q[1];
x q[2];
rz(-0.35700123) q[3];
sx q[3];
rz(-1.9244266) q[3];
sx q[3];
rz(-0.90333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5985976) q[2];
sx q[2];
rz(-2.6948805) q[2];
sx q[2];
rz(0.57999769) q[2];
rz(-0.18481208) q[3];
sx q[3];
rz(-1.5744659) q[3];
sx q[3];
rz(0.68661657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0762894) q[0];
sx q[0];
rz(-0.46173254) q[0];
sx q[0];
rz(-0.28934685) q[0];
rz(0.099960001) q[1];
sx q[1];
rz(-2.3048765) q[1];
sx q[1];
rz(2.3438556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8160523) q[0];
sx q[0];
rz(-2.1959248) q[0];
sx q[0];
rz(0.16640551) q[0];
x q[1];
rz(-2.6358068) q[2];
sx q[2];
rz(-0.2295851) q[2];
sx q[2];
rz(-1.2248271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9790837) q[1];
sx q[1];
rz(-1.519527) q[1];
sx q[1];
rz(2.8460345) q[1];
x q[2];
rz(-3.0403071) q[3];
sx q[3];
rz(-2.2982979) q[3];
sx q[3];
rz(-1.5366713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.37819698) q[2];
sx q[2];
rz(-1.244647) q[2];
sx q[2];
rz(-2.3217616) q[2];
rz(-1.648692) q[3];
sx q[3];
rz(-2.7870352) q[3];
sx q[3];
rz(0.41560391) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019526871) q[0];
sx q[0];
rz(-2.8391835) q[0];
sx q[0];
rz(0.31461883) q[0];
rz(0.64612359) q[1];
sx q[1];
rz(-1.6982634) q[1];
sx q[1];
rz(-3.019928) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3356161) q[0];
sx q[0];
rz(-1.2947963) q[0];
sx q[0];
rz(1.1742623) q[0];
rz(-pi) q[1];
rz(0.88603772) q[2];
sx q[2];
rz(-0.90038634) q[2];
sx q[2];
rz(-1.3051906) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.52895228) q[1];
sx q[1];
rz(-1.4390872) q[1];
sx q[1];
rz(-2.6963833) q[1];
rz(-pi) q[2];
rz(-0.96243919) q[3];
sx q[3];
rz(-1.7311058) q[3];
sx q[3];
rz(0.95534217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6803711) q[2];
sx q[2];
rz(-1.9663726) q[2];
sx q[2];
rz(-0.97869527) q[2];
rz(2.3693502) q[3];
sx q[3];
rz(-2.6148836) q[3];
sx q[3];
rz(-1.0129119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0270281) q[0];
sx q[0];
rz(-1.9975198) q[0];
sx q[0];
rz(-1.8373328) q[0];
rz(-3.0005241) q[1];
sx q[1];
rz(-2.8619659) q[1];
sx q[1];
rz(-1.28349) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.210189) q[0];
sx q[0];
rz(-0.21912665) q[0];
sx q[0];
rz(-1.7396084) q[0];
x q[1];
rz(-1.8780297) q[2];
sx q[2];
rz(-2.8275785) q[2];
sx q[2];
rz(0.28672126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.56399689) q[1];
sx q[1];
rz(-1.8297927) q[1];
sx q[1];
rz(-1.4513272) q[1];
x q[2];
rz(3.0937681) q[3];
sx q[3];
rz(-2.3615814) q[3];
sx q[3];
rz(0.081950233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.86113247) q[2];
sx q[2];
rz(-2.2304163) q[2];
sx q[2];
rz(-2.148441) q[2];
rz(1.447621) q[3];
sx q[3];
rz(-1.9341035) q[3];
sx q[3];
rz(1.8705468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.3250378) q[0];
sx q[0];
rz(-0.10295454) q[0];
sx q[0];
rz(-2.2651941) q[0];
rz(1.5986298) q[1];
sx q[1];
rz(-1.2878659) q[1];
sx q[1];
rz(-1.1184982) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0792313) q[0];
sx q[0];
rz(-2.1665299) q[0];
sx q[0];
rz(-2.0078288) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1886979) q[2];
sx q[2];
rz(-2.2006773) q[2];
sx q[2];
rz(3.0840616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1269605) q[1];
sx q[1];
rz(-1.6829957) q[1];
sx q[1];
rz(-0.87027624) q[1];
rz(1.275609) q[3];
sx q[3];
rz(-0.94691197) q[3];
sx q[3];
rz(2.5518887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8886275) q[2];
sx q[2];
rz(-0.94433633) q[2];
sx q[2];
rz(-1.2523119) q[2];
rz(2.0400932) q[3];
sx q[3];
rz(-1.3809729) q[3];
sx q[3];
rz(-1.8496877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0970704) q[0];
sx q[0];
rz(-3.0945393) q[0];
sx q[0];
rz(2.4142921) q[0];
rz(-2.3749088) q[1];
sx q[1];
rz(-2.2946281) q[1];
sx q[1];
rz(-1.9511706) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56404468) q[0];
sx q[0];
rz(-0.8506895) q[0];
sx q[0];
rz(1.7386566) q[0];
x q[1];
rz(-1.6929564) q[2];
sx q[2];
rz(-0.66389951) q[2];
sx q[2];
rz(-0.66421504) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1228635) q[1];
sx q[1];
rz(-1.8179389) q[1];
sx q[1];
rz(-2.8309532) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6383357) q[3];
sx q[3];
rz(-1.0531593) q[3];
sx q[3];
rz(-2.518017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7850354) q[2];
sx q[2];
rz(-2.0544453) q[2];
sx q[2];
rz(-1.1615151) q[2];
rz(-1.131743) q[3];
sx q[3];
rz(-0.83860207) q[3];
sx q[3];
rz(1.0034126) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7575191) q[0];
sx q[0];
rz(-0.90255559) q[0];
sx q[0];
rz(2.3046816) q[0];
rz(1.8867672) q[1];
sx q[1];
rz(-1.8950987) q[1];
sx q[1];
rz(-1.7991039) q[1];
rz(-0.26623434) q[2];
sx q[2];
rz(-1.0221368) q[2];
sx q[2];
rz(-0.72697097) q[2];
rz(-1.499621) q[3];
sx q[3];
rz(-0.19656062) q[3];
sx q[3];
rz(1.82609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
