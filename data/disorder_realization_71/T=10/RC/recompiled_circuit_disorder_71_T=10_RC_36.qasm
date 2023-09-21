OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(0.82984501) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(-0.87632626) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012955879) q[0];
sx q[0];
rz(-0.46015938) q[0];
sx q[0];
rz(2.602306) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1562188) q[2];
sx q[2];
rz(-1.7314163) q[2];
sx q[2];
rz(-1.0759682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6688924) q[1];
sx q[1];
rz(-1.2408222) q[1];
sx q[1];
rz(-0.015923576) q[1];
rz(1.417744) q[3];
sx q[3];
rz(-1.9733841) q[3];
sx q[3];
rz(-0.11868417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17065389) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(2.412964) q[2];
rz(-0.5209926) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(-0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-2.0200502) q[0];
rz(-2.8858378) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(2.2671525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5912936) q[0];
sx q[0];
rz(-1.8739788) q[0];
sx q[0];
rz(1.1217872) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1463257) q[2];
sx q[2];
rz(-1.7003635) q[2];
sx q[2];
rz(2.1339983) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.892131) q[1];
sx q[1];
rz(-2.3453237) q[1];
sx q[1];
rz(-2.4888121) q[1];
x q[2];
rz(-1.3014684) q[3];
sx q[3];
rz(-1.6225617) q[3];
sx q[3];
rz(-0.89382899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4008537) q[2];
sx q[2];
rz(-2.6541371) q[2];
sx q[2];
rz(2.7056616) q[2];
rz(0.68108264) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(-2.7387103) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9044559) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(-1.0748192) q[0];
rz(2.3020321) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(0.39594617) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2303282) q[0];
sx q[0];
rz(-2.0718144) q[0];
sx q[0];
rz(2.688174) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6171574) q[2];
sx q[2];
rz(-0.29435396) q[2];
sx q[2];
rz(-2.4183395) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0242651) q[1];
sx q[1];
rz(-1.5803442) q[1];
sx q[1];
rz(-0.68318232) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.08926908) q[3];
sx q[3];
rz(-2.8850728) q[3];
sx q[3];
rz(1.4005816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5376771) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(-2.9023857) q[2];
rz(3.0662597) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(1.8384365) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72162119) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(-2.3216632) q[0];
rz(2.6539102) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(-2.908169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6596286) q[0];
sx q[0];
rz(-3.0273962) q[0];
sx q[0];
rz(-1.0114848) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31077023) q[2];
sx q[2];
rz(-0.66590532) q[2];
sx q[2];
rz(-3.0096465) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1366795) q[1];
sx q[1];
rz(-2.3466952) q[1];
sx q[1];
rz(-2.6585048) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7531284) q[3];
sx q[3];
rz(-0.82287517) q[3];
sx q[3];
rz(2.8188761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0908115) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(0.74742571) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(-2.8994765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6915879) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(0.064237021) q[0];
rz(-2.1977987) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(-0.76104004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9984263) q[0];
sx q[0];
rz(-0.26253065) q[0];
sx q[0];
rz(1.4545928) q[0];
rz(1.7663076) q[2];
sx q[2];
rz(-1.7695145) q[2];
sx q[2];
rz(-2.9544601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0423454) q[1];
sx q[1];
rz(-1.595101) q[1];
sx q[1];
rz(-2.821032) q[1];
rz(1.8499591) q[3];
sx q[3];
rz(-1.8431292) q[3];
sx q[3];
rz(2.2945822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.80660194) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(-0.09207329) q[2];
rz(0.66172415) q[3];
sx q[3];
rz(-0.79143733) q[3];
sx q[3];
rz(1.3823284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6732366) q[0];
sx q[0];
rz(-1.8259003) q[0];
sx q[0];
rz(-0.026542149) q[0];
rz(2.2684855) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(0.32593265) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78544261) q[0];
sx q[0];
rz(-0.66473648) q[0];
sx q[0];
rz(2.507693) q[0];
rz(-0.44547703) q[2];
sx q[2];
rz(-1.2947086) q[2];
sx q[2];
rz(2.7518227) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0210452) q[1];
sx q[1];
rz(-1.435934) q[1];
sx q[1];
rz(0.30121505) q[1];
x q[2];
rz(1.5876706) q[3];
sx q[3];
rz(-1.6718739) q[3];
sx q[3];
rz(0.24731393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5439593) q[2];
sx q[2];
rz(-1.3175069) q[2];
sx q[2];
rz(-0.65417543) q[2];
rz(1.7116961) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-0.54106075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27286801) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-2.4196999) q[0];
rz(1.7294653) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(0.11925764) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8598547) q[0];
sx q[0];
rz(-2.2681232) q[0];
sx q[0];
rz(1.2233234) q[0];
x q[1];
rz(3.0316333) q[2];
sx q[2];
rz(-1.7320776) q[2];
sx q[2];
rz(1.1683299) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99265656) q[1];
sx q[1];
rz(-1.7033556) q[1];
sx q[1];
rz(2.0723144) q[1];
x q[2];
rz(-0.70763208) q[3];
sx q[3];
rz(-1.8579357) q[3];
sx q[3];
rz(1.6361145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(1.3593486) q[2];
rz(-0.75891495) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83157241) q[0];
sx q[0];
rz(-0.65905237) q[0];
sx q[0];
rz(-2.0781562) q[0];
rz(-2.8670782) q[1];
sx q[1];
rz(-1.9332705) q[1];
sx q[1];
rz(-0.88561052) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0793593) q[0];
sx q[0];
rz(-1.0842807) q[0];
sx q[0];
rz(-1.3079206) q[0];
rz(-pi) q[1];
rz(3.1259414) q[2];
sx q[2];
rz(-0.99132292) q[2];
sx q[2];
rz(1.5660812) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.13247989) q[1];
sx q[1];
rz(-1.1785893) q[1];
sx q[1];
rz(2.1177887) q[1];
x q[2];
rz(-2.3732244) q[3];
sx q[3];
rz(-1.0230912) q[3];
sx q[3];
rz(-0.87793575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.72835913) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(1.1317066) q[2];
rz(1.0845832) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8383012) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(2.0595179) q[0];
rz(1.8661631) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(-1.1358322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71400434) q[0];
sx q[0];
rz(-1.4052608) q[0];
sx q[0];
rz(-1.0323314) q[0];
rz(-pi) q[1];
rz(-1.5485974) q[2];
sx q[2];
rz(-2.0054842) q[2];
sx q[2];
rz(-1.7322025) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.218704) q[1];
sx q[1];
rz(-1.522038) q[1];
sx q[1];
rz(2.7152275) q[1];
rz(-pi) q[2];
rz(2.7157852) q[3];
sx q[3];
rz(-2.1771181) q[3];
sx q[3];
rz(3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11848005) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(0.87289587) q[2];
rz(0.84351271) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(-0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2492367) q[0];
sx q[0];
rz(-1.7974239) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(-1.0247914) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(-1.1970253) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7322757) q[0];
sx q[0];
rz(-1.1336375) q[0];
sx q[0];
rz(0.93426312) q[0];
x q[1];
rz(0.12038259) q[2];
sx q[2];
rz(-1.3345846) q[2];
sx q[2];
rz(1.2506968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4492606) q[1];
sx q[1];
rz(-2.3093866) q[1];
sx q[1];
rz(2.2663121) q[1];
rz(-pi) q[2];
rz(-1.193612) q[3];
sx q[3];
rz(-1.0918655) q[3];
sx q[3];
rz(2.9951028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3060351) q[2];
sx q[2];
rz(-0.71283895) q[2];
sx q[2];
rz(-2.3416134) q[2];
rz(1.1768613) q[3];
sx q[3];
rz(-1.4326982) q[3];
sx q[3];
rz(-2.1879788) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70893127) q[0];
sx q[0];
rz(-2.9928757) q[0];
sx q[0];
rz(-2.3401674) q[0];
rz(0.52195436) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(-2.3847053) q[2];
sx q[2];
rz(-2.6464528) q[2];
sx q[2];
rz(1.1337627) q[2];
rz(-2.2189191) q[3];
sx q[3];
rz(-2.1019983) q[3];
sx q[3];
rz(-2.0207873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
