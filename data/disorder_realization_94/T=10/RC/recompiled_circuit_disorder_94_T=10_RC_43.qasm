OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.86413971) q[0];
sx q[0];
rz(-1.5530518) q[0];
sx q[0];
rz(1.6341524) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(0.61520666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078755137) q[0];
sx q[0];
rz(-2.3427999) q[0];
sx q[0];
rz(-0.90057217) q[0];
rz(-pi) q[1];
rz(-0.46703672) q[2];
sx q[2];
rz(-2.7454498) q[2];
sx q[2];
rz(3.0774088) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7429744) q[1];
sx q[1];
rz(-2.3906039) q[1];
sx q[1];
rz(-2.2258334) q[1];
x q[2];
rz(-2.6184222) q[3];
sx q[3];
rz(-2.7912931) q[3];
sx q[3];
rz(-1.3594128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.4404099) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(0.33828503) q[2];
rz(1.4398549) q[3];
sx q[3];
rz(-0.9153291) q[3];
sx q[3];
rz(2.2556944) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97025362) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(-3.1112444) q[0];
rz(0.066210315) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(1.617584) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.108792) q[0];
sx q[0];
rz(-1.5691225) q[0];
sx q[0];
rz(-1.770442) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5288058) q[2];
sx q[2];
rz(-2.6868372) q[2];
sx q[2];
rz(3.0595879) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.46088947) q[1];
sx q[1];
rz(-1.5017121) q[1];
sx q[1];
rz(2.1417888) q[1];
rz(-pi) q[2];
x q[2];
rz(0.094035427) q[3];
sx q[3];
rz(-1.7574851) q[3];
sx q[3];
rz(1.3446913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7559738) q[2];
sx q[2];
rz(-2.1531343) q[2];
sx q[2];
rz(1.9937817) q[2];
rz(1.8148445) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(-2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1266992) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(2.8258064) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(2.8895203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7092428) q[0];
sx q[0];
rz(-2.3253257) q[0];
sx q[0];
rz(-0.66246756) q[0];
rz(-2.2601068) q[2];
sx q[2];
rz(-0.93886095) q[2];
sx q[2];
rz(1.880868) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7604916) q[1];
sx q[1];
rz(-0.64844202) q[1];
sx q[1];
rz(-1.2566503) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4113703) q[3];
sx q[3];
rz(-2.632004) q[3];
sx q[3];
rz(-2.2024221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0198274) q[2];
sx q[2];
rz(-2.6941507) q[2];
sx q[2];
rz(3.1075409) q[2];
rz(0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-1.0954558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5144192) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(-1.5267641) q[0];
rz(1.0871672) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(2.4345051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1336466) q[0];
sx q[0];
rz(-2.5292853) q[0];
sx q[0];
rz(-0.72511073) q[0];
x q[1];
rz(-0.25755067) q[2];
sx q[2];
rz(-0.45986816) q[2];
sx q[2];
rz(0.79007733) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.9634339) q[1];
sx q[1];
rz(-1.4512832) q[1];
sx q[1];
rz(-2.2258289) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2624192) q[3];
sx q[3];
rz(-1.3652507) q[3];
sx q[3];
rz(-1.8431078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2924071) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(-2.6814931) q[2];
rz(-1.7442616) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.3209155) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(-1.7472349) q[0];
rz(1.0955411) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-0.25462338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6760315) q[0];
sx q[0];
rz(-1.6499632) q[0];
sx q[0];
rz(0.013750793) q[0];
x q[1];
rz(1.5484372) q[2];
sx q[2];
rz(-0.5776814) q[2];
sx q[2];
rz(-2.3235842) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8498807) q[1];
sx q[1];
rz(-2.2176718) q[1];
sx q[1];
rz(1.2667659) q[1];
rz(-pi) q[2];
rz(-1.0765431) q[3];
sx q[3];
rz(-2.3688865) q[3];
sx q[3];
rz(-2.5172174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1191117) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(-1.3767892) q[2];
rz(1.6453751) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(-2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6901533) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-2.8421463) q[0];
rz(-2.1014138) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(-0.20656955) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20139192) q[0];
sx q[0];
rz(-0.79027806) q[0];
sx q[0];
rz(-1.9076365) q[0];
rz(-2.0986404) q[2];
sx q[2];
rz(-1.9163418) q[2];
sx q[2];
rz(2.8258459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8787074) q[1];
sx q[1];
rz(-0.72808121) q[1];
sx q[1];
rz(-1.3036149) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9566831) q[3];
sx q[3];
rz(-0.69291249) q[3];
sx q[3];
rz(3.0544359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2999337) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(-2.858813) q[2];
rz(2.3287866) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2151826) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(2.7600631) q[0];
rz(-2.5577257) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(1.3279703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8570003) q[0];
sx q[0];
rz(-2.6866331) q[0];
sx q[0];
rz(-1.8332464) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96372693) q[2];
sx q[2];
rz(-0.71787314) q[2];
sx q[2];
rz(-1.8075862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92762891) q[1];
sx q[1];
rz(-2.2447526) q[1];
sx q[1];
rz(-0.2143292) q[1];
rz(-2.6286969) q[3];
sx q[3];
rz(-1.1937871) q[3];
sx q[3];
rz(0.60037724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61838377) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(1.3605114) q[2];
rz(-1.4303738) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(0.18187901) q[0];
rz(0.47422844) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(-0.95091933) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5247105) q[0];
sx q[0];
rz(-1.8572154) q[0];
sx q[0];
rz(-2.8911203) q[0];
rz(-pi) q[1];
rz(3.0792564) q[2];
sx q[2];
rz(-1.8364292) q[2];
sx q[2];
rz(-0.072364256) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7656895) q[1];
sx q[1];
rz(-1.5058869) q[1];
sx q[1];
rz(-0.20903559) q[1];
rz(0.37787921) q[3];
sx q[3];
rz(-1.1933019) q[3];
sx q[3];
rz(2.8187403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2237079) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(-3.0656832) q[2];
rz(-2.5935796) q[3];
sx q[3];
rz(-1.8173822) q[3];
sx q[3];
rz(0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.6178745) q[0];
sx q[0];
rz(-1.0567559) q[0];
sx q[0];
rz(1.3762208) q[0];
rz(-2.7245522) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(0.65972796) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2840246) q[0];
sx q[0];
rz(-2.095247) q[0];
sx q[0];
rz(0.88733034) q[0];
rz(-0.78166878) q[2];
sx q[2];
rz(-1.726892) q[2];
sx q[2];
rz(-2.8761656) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8068741) q[1];
sx q[1];
rz(-1.8611307) q[1];
sx q[1];
rz(-2.2335386) q[1];
rz(-0.057007313) q[3];
sx q[3];
rz(-2.1087286) q[3];
sx q[3];
rz(1.9407335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.518121) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(-1.0260322) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-0.025645105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24348564) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(-1.8359258) q[0];
rz(-1.1765515) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(-1.0356888) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.068533) q[0];
sx q[0];
rz(-1.8951299) q[0];
sx q[0];
rz(2.7102094) q[0];
rz(-2.5752441) q[2];
sx q[2];
rz(-0.95745917) q[2];
sx q[2];
rz(-1.4373506) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2346238) q[1];
sx q[1];
rz(-1.8112438) q[1];
sx q[1];
rz(-2.8333227) q[1];
rz(-pi) q[2];
rz(0.74929897) q[3];
sx q[3];
rz(-1.370508) q[3];
sx q[3];
rz(-2.0405958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62853652) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(-1.9899842) q[2];
rz(2.628905) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(-0.36469665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37968996) q[0];
sx q[0];
rz(-1.3544449) q[0];
sx q[0];
rz(0.83308573) q[0];
rz(-1.6336541) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(1.6133616) q[2];
sx q[2];
rz(-1.1602989) q[2];
sx q[2];
rz(1.5699671) q[2];
rz(-1.5796173) q[3];
sx q[3];
rz(-2.1891441) q[3];
sx q[3];
rz(1.7563663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];