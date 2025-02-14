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
rz(-2.1962533) q[0];
sx q[0];
rz(0.0013295833) q[0];
rz(0.20345649) q[1];
sx q[1];
rz(-0.23509547) q[1];
sx q[1];
rz(-2.5703365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8205367) q[0];
sx q[0];
rz(-0.97194304) q[0];
sx q[0];
rz(-1.0060034) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18988256) q[2];
sx q[2];
rz(-1.4212835) q[2];
sx q[2];
rz(2.8854795) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9524051) q[1];
sx q[1];
rz(-2.6509977) q[1];
sx q[1];
rz(0.036019939) q[1];
x q[2];
rz(-0.68526973) q[3];
sx q[3];
rz(-0.56495404) q[3];
sx q[3];
rz(0.23811114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0433537) q[2];
sx q[2];
rz(-1.9184155) q[2];
sx q[2];
rz(2.5377048) q[2];
rz(-1.0107001) q[3];
sx q[3];
rz(-0.5637919) q[3];
sx q[3];
rz(0.92196661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.0826223) q[0];
sx q[0];
rz(-1.0930467) q[0];
sx q[0];
rz(0.0086722886) q[0];
rz(1.1588833) q[1];
sx q[1];
rz(-0.68717879) q[1];
sx q[1];
rz(-0.11223665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85189795) q[0];
sx q[0];
rz(-1.0852975) q[0];
sx q[0];
rz(1.8957036) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64666266) q[2];
sx q[2];
rz(-1.5407667) q[2];
sx q[2];
rz(-2.554837) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48277897) q[1];
sx q[1];
rz(-0.54721745) q[1];
sx q[1];
rz(0.43756715) q[1];
rz(2.8772164) q[3];
sx q[3];
rz(-1.2942787) q[3];
sx q[3];
rz(-3.0540651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23942854) q[2];
sx q[2];
rz(-1.3210693) q[2];
sx q[2];
rz(-2.2018137) q[2];
rz(0.93722614) q[3];
sx q[3];
rz(-1.7669433) q[3];
sx q[3];
rz(0.37041131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20060191) q[0];
sx q[0];
rz(-0.90018278) q[0];
sx q[0];
rz(-0.58315939) q[0];
rz(2.0896301) q[1];
sx q[1];
rz(-0.58590341) q[1];
sx q[1];
rz(0.016187035) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.095854) q[0];
sx q[0];
rz(-1.7940137) q[0];
sx q[0];
rz(-2.2824817) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80183713) q[2];
sx q[2];
rz(-2.7268598) q[2];
sx q[2];
rz(-1.5842337) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1165818) q[1];
sx q[1];
rz(-1.7336646) q[1];
sx q[1];
rz(-3.0322269) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2233769) q[3];
sx q[3];
rz(-0.86091061) q[3];
sx q[3];
rz(2.6588499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6215324) q[2];
sx q[2];
rz(-1.748087) q[2];
sx q[2];
rz(-0.072546093) q[2];
rz(-1.0091311) q[3];
sx q[3];
rz(-0.79807177) q[3];
sx q[3];
rz(0.23536853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9811454) q[0];
sx q[0];
rz(-0.35780847) q[0];
sx q[0];
rz(2.2571046) q[0];
rz(-1.5011939) q[1];
sx q[1];
rz(-1.6694992) q[1];
sx q[1];
rz(-0.59008682) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6654331) q[0];
sx q[0];
rz(-0.64393015) q[0];
sx q[0];
rz(-0.91899271) q[0];
rz(1.2954166) q[2];
sx q[2];
rz(-1.3510873) q[2];
sx q[2];
rz(1.4003458) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6650538) q[1];
sx q[1];
rz(-1.0805942) q[1];
sx q[1];
rz(-0.3785822) q[1];
x q[2];
rz(0.25165148) q[3];
sx q[3];
rz(-2.6918937) q[3];
sx q[3];
rz(-2.3684199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0044452) q[2];
sx q[2];
rz(-1.4819757) q[2];
sx q[2];
rz(-1.2443292) q[2];
rz(-1.3589877) q[3];
sx q[3];
rz(-2.6748896) q[3];
sx q[3];
rz(2.9984503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2537848) q[0];
sx q[0];
rz(-2.4847327) q[0];
sx q[0];
rz(2.7746871) q[0];
rz(-1.3430355) q[1];
sx q[1];
rz(-0.29452205) q[1];
sx q[1];
rz(-1.8273182) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2283233) q[0];
sx q[0];
rz(-1.1342888) q[0];
sx q[0];
rz(-2.4455347) q[0];
rz(2.6571006) q[2];
sx q[2];
rz(-1.3132261) q[2];
sx q[2];
rz(-0.07312873) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3245222) q[1];
sx q[1];
rz(-1.2597023) q[1];
sx q[1];
rz(2.2242311) q[1];
rz(0.81277992) q[3];
sx q[3];
rz(-0.4970937) q[3];
sx q[3];
rz(1.415783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5985976) q[2];
sx q[2];
rz(-0.44671217) q[2];
sx q[2];
rz(-0.57999769) q[2];
rz(-0.18481208) q[3];
sx q[3];
rz(-1.5744659) q[3];
sx q[3];
rz(0.68661657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0762894) q[0];
sx q[0];
rz(-2.6798601) q[0];
sx q[0];
rz(-0.28934685) q[0];
rz(-0.099960001) q[1];
sx q[1];
rz(-2.3048765) q[1];
sx q[1];
rz(0.79773703) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7983618) q[0];
sx q[0];
rz(-1.7055178) q[0];
sx q[0];
rz(-0.93905296) q[0];
rz(-pi) q[1];
rz(2.6358068) q[2];
sx q[2];
rz(-0.2295851) q[2];
sx q[2];
rz(-1.9167655) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5665566) q[1];
sx q[1];
rz(-2.8417491) q[1];
sx q[1];
rz(2.9672102) q[1];
x q[2];
rz(-3.0403071) q[3];
sx q[3];
rz(-2.2982979) q[3];
sx q[3];
rz(1.6049214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7633957) q[2];
sx q[2];
rz(-1.244647) q[2];
sx q[2];
rz(2.3217616) q[2];
rz(1.4929006) q[3];
sx q[3];
rz(-2.7870352) q[3];
sx q[3];
rz(0.41560391) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1220658) q[0];
sx q[0];
rz(-2.8391835) q[0];
sx q[0];
rz(0.31461883) q[0];
rz(-0.64612359) q[1];
sx q[1];
rz(-1.4433292) q[1];
sx q[1];
rz(0.12166469) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87843175) q[0];
sx q[0];
rz(-1.9515349) q[0];
sx q[0];
rz(2.8436766) q[0];
rz(0.88603772) q[2];
sx q[2];
rz(-2.2412063) q[2];
sx q[2];
rz(-1.836402) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6126404) q[1];
sx q[1];
rz(-1.7025055) q[1];
sx q[1];
rz(0.44520933) q[1];
x q[2];
rz(-2.9470351) q[3];
sx q[3];
rz(-2.1702507) q[3];
sx q[3];
rz(0.50473467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6803711) q[2];
sx q[2];
rz(-1.9663726) q[2];
sx q[2];
rz(0.97869527) q[2];
rz(0.77224246) q[3];
sx q[3];
rz(-0.52670908) q[3];
sx q[3];
rz(-1.0129119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0270281) q[0];
sx q[0];
rz(-1.1440729) q[0];
sx q[0];
rz(1.3042599) q[0];
rz(-3.0005241) q[1];
sx q[1];
rz(-2.8619659) q[1];
sx q[1];
rz(-1.28349) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52545628) q[0];
sx q[0];
rz(-1.5342664) q[0];
sx q[0];
rz(-1.7869048) q[0];
rz(-1.2705994) q[2];
sx q[2];
rz(-1.4772479) q[2];
sx q[2];
rz(-1.5771505) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0375335) q[1];
sx q[1];
rz(-1.4553299) q[1];
sx q[1];
rz(2.8808198) q[1];
rz(-2.3621534) q[3];
sx q[3];
rz(-1.5371684) q[3];
sx q[3];
rz(-1.5228576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2804602) q[2];
sx q[2];
rz(-0.91117636) q[2];
sx q[2];
rz(0.99315161) q[2];
rz(-1.6939717) q[3];
sx q[3];
rz(-1.2074892) q[3];
sx q[3];
rz(1.2710458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8165548) q[0];
sx q[0];
rz(-3.0386381) q[0];
sx q[0];
rz(0.8763985) q[0];
rz(1.5986298) q[1];
sx q[1];
rz(-1.2878659) q[1];
sx q[1];
rz(-1.1184982) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0792313) q[0];
sx q[0];
rz(-2.1665299) q[0];
sx q[0];
rz(-2.0078288) q[0];
rz(-pi) q[1];
rz(-2.1886979) q[2];
sx q[2];
rz(-0.94091533) q[2];
sx q[2];
rz(3.0840616) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7174445) q[1];
sx q[1];
rz(-2.433648) q[1];
sx q[1];
rz(-1.7438386) q[1];
rz(-pi) q[2];
rz(1.275609) q[3];
sx q[3];
rz(-2.1946807) q[3];
sx q[3];
rz(-2.5518887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8886275) q[2];
sx q[2];
rz(-2.1972563) q[2];
sx q[2];
rz(-1.2523119) q[2];
rz(-2.0400932) q[3];
sx q[3];
rz(-1.7606198) q[3];
sx q[3];
rz(1.291905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0445223) q[0];
sx q[0];
rz(-0.047053311) q[0];
sx q[0];
rz(-2.4142921) q[0];
rz(2.3749088) q[1];
sx q[1];
rz(-0.84696451) q[1];
sx q[1];
rz(-1.9511706) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3260314) q[0];
sx q[0];
rz(-2.4056068) q[0];
sx q[0];
rz(-0.18819564) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2310689) q[2];
sx q[2];
rz(-1.6459542) q[2];
sx q[2];
rz(2.1386168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3695581) q[1];
sx q[1];
rz(-1.8716964) q[1];
sx q[1];
rz(1.8298261) q[1];
rz(2.6383357) q[3];
sx q[3];
rz(-1.0531593) q[3];
sx q[3];
rz(2.518017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.35655725) q[2];
sx q[2];
rz(-2.0544453) q[2];
sx q[2];
rz(1.1615151) q[2];
rz(1.131743) q[3];
sx q[3];
rz(-0.83860207) q[3];
sx q[3];
rz(-1.0034126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3840735) q[0];
sx q[0];
rz(-0.90255559) q[0];
sx q[0];
rz(2.3046816) q[0];
rz(-1.8867672) q[1];
sx q[1];
rz(-1.246494) q[1];
sx q[1];
rz(1.3424887) q[1];
rz(2.135545) q[2];
sx q[2];
rz(-1.3443832) q[2];
sx q[2];
rz(0.98510712) q[2];
rz(-1.3747207) q[3];
sx q[3];
rz(-1.5569073) q[3];
sx q[3];
rz(0.18548439) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
