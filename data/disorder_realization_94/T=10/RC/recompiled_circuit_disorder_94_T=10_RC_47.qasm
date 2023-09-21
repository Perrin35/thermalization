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
rz(-1.5074402) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(2.5453321) q[1];
sx q[1];
rz(8.8095713) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9278487) q[0];
sx q[0];
rz(-0.97457492) q[0];
sx q[0];
rz(2.5736789) q[0];
x q[1];
rz(-1.3846606) q[2];
sx q[2];
rz(-1.2190483) q[2];
sx q[2];
rz(-2.577201) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.93278904) q[1];
sx q[1];
rz(-2.1425769) q[1];
sx q[1];
rz(-0.51704452) q[1];
x q[2];
rz(2.6184222) q[3];
sx q[3];
rz(-0.35029951) q[3];
sx q[3];
rz(1.7821799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7011828) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(0.33828503) q[2];
rz(-1.7017378) q[3];
sx q[3];
rz(-0.9153291) q[3];
sx q[3];
rz(-0.88589823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171339) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(3.1112444) q[0];
rz(-0.066210315) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(-1.617584) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47027662) q[0];
sx q[0];
rz(-2.94194) q[0];
sx q[0];
rz(1.5792363) q[0];
rz(1.5288058) q[2];
sx q[2];
rz(-0.4547555) q[2];
sx q[2];
rz(-3.0595879) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.46088947) q[1];
sx q[1];
rz(-1.5017121) q[1];
sx q[1];
rz(0.99980385) q[1];
rz(0.094035427) q[3];
sx q[3];
rz(-1.3841076) q[3];
sx q[3];
rz(1.7969014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38561884) q[2];
sx q[2];
rz(-2.1531343) q[2];
sx q[2];
rz(-1.9937817) q[2];
rz(-1.8148445) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(-0.23708788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1266992) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(0.31578627) q[0];
rz(-0.93859998) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(-0.25207239) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35226563) q[0];
sx q[0];
rz(-2.0354712) q[0];
sx q[0];
rz(-2.4436823) q[0];
rz(-2.3825233) q[2];
sx q[2];
rz(-1.0319064) q[2];
sx q[2];
rz(-2.3784504) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1356126) q[1];
sx q[1];
rz(-0.95893919) q[1];
sx q[1];
rz(-0.23000418) q[1];
x q[2];
rz(0.39450816) q[3];
sx q[3];
rz(-1.9022226) q[3];
sx q[3];
rz(0.031771914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1217653) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(0.034051731) q[2];
rz(0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5144192) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(1.5267641) q[0];
rz(-1.0871672) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(-0.70708752) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83311659) q[0];
sx q[0];
rz(-2.0154698) q[0];
sx q[0];
rz(1.1348669) q[0];
rz(-pi) q[1];
rz(0.25755067) q[2];
sx q[2];
rz(-0.45986816) q[2];
sx q[2];
rz(-0.79007733) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9634339) q[1];
sx q[1];
rz(-1.4512832) q[1];
sx q[1];
rz(0.91576373) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96890038) q[3];
sx q[3];
rz(-0.36877353) q[3];
sx q[3];
rz(0.2975279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84918555) q[2];
sx q[2];
rz(-1.9533998) q[2];
sx q[2];
rz(-2.6814931) q[2];
rz(1.397331) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8206772) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(1.3943577) q[0];
rz(2.0460515) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(0.25462338) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8481962) q[0];
sx q[0];
rz(-0.080349803) q[0];
sx q[0];
rz(1.3991762) q[0];
rz(-1.5931555) q[2];
sx q[2];
rz(-2.5639113) q[2];
sx q[2];
rz(-0.81800848) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8498807) q[1];
sx q[1];
rz(-0.92392081) q[1];
sx q[1];
rz(1.2667659) q[1];
rz(-pi) q[2];
rz(0.43318627) q[3];
sx q[3];
rz(-2.232589) q[3];
sx q[3];
rz(-1.2695241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.022481) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(-1.3767892) q[2];
rz(-1.6453751) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45143932) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-2.8421463) q[0];
rz(-2.1014138) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(-2.9350231) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26039133) q[0];
sx q[0];
rz(-0.83575373) q[0];
sx q[0];
rz(-0.3221237) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7467696) q[2];
sx q[2];
rz(-1.0770814) q[2];
sx q[2];
rz(-2.0815108) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5094604) q[1];
sx q[1];
rz(-1.7473979) q[1];
sx q[1];
rz(-2.280974) q[1];
rz(1.9566831) q[3];
sx q[3];
rz(-0.69291249) q[3];
sx q[3];
rz(-0.087156765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(2.858813) q[2];
rz(-0.81280604) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(-2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2151826) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(-2.7600631) q[0];
rz(-0.58386699) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(1.8136224) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0921811) q[0];
sx q[0];
rz(-1.4565399) q[0];
sx q[0];
rz(2.0121502) q[0];
rz(-2.1930201) q[2];
sx q[2];
rz(-1.186139) q[2];
sx q[2];
rz(-2.4228061) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3634062) q[1];
sx q[1];
rz(-1.4038329) q[1];
sx q[1];
rz(-0.88552514) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6286969) q[3];
sx q[3];
rz(-1.9478056) q[3];
sx q[3];
rz(2.5412154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5232089) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(-1.7810812) q[2];
rz(-1.7112188) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-2.2935304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1469864) q[0];
sx q[0];
rz(-1.9817579) q[0];
sx q[0];
rz(-2.9597136) q[0];
rz(-2.6673642) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(-2.1906733) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3527746) q[0];
sx q[0];
rz(-0.37820658) q[0];
sx q[0];
rz(-2.2703855) q[0];
rz(0.062336246) q[2];
sx q[2];
rz(-1.3051635) q[2];
sx q[2];
rz(3.0692284) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9329405) q[1];
sx q[1];
rz(-1.7793852) q[1];
sx q[1];
rz(-1.5044466) q[1];
x q[2];
rz(0.37787921) q[3];
sx q[3];
rz(-1.9482908) q[3];
sx q[3];
rz(-2.8187403) q[3];
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
rz(2.5935796) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(-2.5089335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.52371812) q[0];
sx q[0];
rz(-1.0567559) q[0];
sx q[0];
rz(-1.3762208) q[0];
rz(2.7245522) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(2.4818647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8575681) q[0];
sx q[0];
rz(-2.095247) q[0];
sx q[0];
rz(0.88733034) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3599239) q[2];
sx q[2];
rz(-1.4147007) q[2];
sx q[2];
rz(0.26542703) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.58762) q[1];
sx q[1];
rz(-0.71466043) q[1];
sx q[1];
rz(-1.1187394) q[1];
rz(-pi) q[2];
rz(3.0845853) q[3];
sx q[3];
rz(-1.032864) q[3];
sx q[3];
rz(1.2008592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.518121) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(1.0260322) q[2];
rz(-0.14885151) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-3.1159475) q[3];
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
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24348564) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(1.8359258) q[0];
rz(-1.1765515) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(-1.0356888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4982088) q[0];
sx q[0];
rz(-1.9783101) q[0];
sx q[0];
rz(-1.9252752) q[0];
rz(-pi) q[1];
rz(-2.2220988) q[2];
sx q[2];
rz(-2.3323625) q[2];
sx q[2];
rz(-0.86916718) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2346238) q[1];
sx q[1];
rz(-1.3303489) q[1];
sx q[1];
rz(0.30827) q[1];
x q[2];
rz(-2.3922937) q[3];
sx q[3];
rz(-1.370508) q[3];
sx q[3];
rz(1.1009969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5130561) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(-1.9899842) q[2];
rz(2.628905) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(-1.528231) q[2];
sx q[2];
rz(-1.1602989) q[2];
sx q[2];
rz(1.5699671) q[2];
rz(3.1291943) q[3];
sx q[3];
rz(-0.61840246) q[3];
sx q[3];
rz(1.7411504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];