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
rz(0.54654044) q[0];
sx q[0];
rz(-2.9822783) q[0];
sx q[0];
rz(1.8848609) q[0];
rz(-2.9873084) q[1];
sx q[1];
rz(-1.0839387) q[1];
sx q[1];
rz(-1.3561603) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0480373) q[0];
sx q[0];
rz(-2.0394271) q[0];
sx q[0];
rz(-2.8891536) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5370813) q[2];
sx q[2];
rz(-1.4178992) q[2];
sx q[2];
rz(2.4375985) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6572842) q[1];
sx q[1];
rz(-2.2806232) q[1];
sx q[1];
rz(-1.4906989) q[1];
rz(-1.1514444) q[3];
sx q[3];
rz(-2.0128801) q[3];
sx q[3];
rz(2.6579554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95947444) q[2];
sx q[2];
rz(-1.139816) q[2];
sx q[2];
rz(0.095890447) q[2];
rz(-0.31283665) q[3];
sx q[3];
rz(-1.0597119) q[3];
sx q[3];
rz(-2.634341) q[3];
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
rz(-1.7392015) q[0];
sx q[0];
rz(-1.5132138) q[0];
sx q[0];
rz(1.8550523) q[0];
rz(1.8402428) q[1];
sx q[1];
rz(-1.8845314) q[1];
sx q[1];
rz(-1.4305065) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98393082) q[0];
sx q[0];
rz(-1.0583911) q[0];
sx q[0];
rz(-0.061467193) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3545121) q[2];
sx q[2];
rz(-0.59709096) q[2];
sx q[2];
rz(-2.7379089) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1455866) q[1];
sx q[1];
rz(-1.4544444) q[1];
sx q[1];
rz(0.12980588) q[1];
rz(2.4336171) q[3];
sx q[3];
rz(-1.7574508) q[3];
sx q[3];
rz(0.1992697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8282738) q[2];
sx q[2];
rz(-0.88130772) q[2];
sx q[2];
rz(0.71448294) q[2];
rz(-1.03164) q[3];
sx q[3];
rz(-2.1561626) q[3];
sx q[3];
rz(-0.45903444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76630509) q[0];
sx q[0];
rz(-2.3630688) q[0];
sx q[0];
rz(2.9189723) q[0];
rz(-0.34046945) q[1];
sx q[1];
rz(-1.6727996) q[1];
sx q[1];
rz(0.31390831) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1134046) q[0];
sx q[0];
rz(-0.88498275) q[0];
sx q[0];
rz(0.96341204) q[0];
rz(-2.6572793) q[2];
sx q[2];
rz(-2.2398758) q[2];
sx q[2];
rz(-2.1824257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7063278) q[1];
sx q[1];
rz(-1.8774596) q[1];
sx q[1];
rz(0.72172647) q[1];
x q[2];
rz(-2.7438358) q[3];
sx q[3];
rz(-0.53103775) q[3];
sx q[3];
rz(-0.49285313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24806222) q[2];
sx q[2];
rz(-1.6174822) q[2];
sx q[2];
rz(-1.0244055) q[2];
rz(-1.8223358) q[3];
sx q[3];
rz(-2.8596467) q[3];
sx q[3];
rz(-0.2969186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9803186) q[0];
sx q[0];
rz(-1.4264822) q[0];
sx q[0];
rz(2.1153765) q[0];
rz(-0.16464344) q[1];
sx q[1];
rz(-2.7963729) q[1];
sx q[1];
rz(-0.77273291) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5948295) q[0];
sx q[0];
rz(-1.208361) q[0];
sx q[0];
rz(2.4700463) q[0];
rz(-pi) q[1];
x q[1];
rz(0.010920694) q[2];
sx q[2];
rz(-0.35775634) q[2];
sx q[2];
rz(-0.27586473) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0028006) q[1];
sx q[1];
rz(-1.7874663) q[1];
sx q[1];
rz(-2.2598221) q[1];
x q[2];
rz(-1.4882632) q[3];
sx q[3];
rz(-3.0898819) q[3];
sx q[3];
rz(-2.5178096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.905895) q[2];
sx q[2];
rz(-0.57196456) q[2];
sx q[2];
rz(0.31615654) q[2];
rz(0.21137992) q[3];
sx q[3];
rz(-1.6930765) q[3];
sx q[3];
rz(2.0785418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0036156) q[0];
sx q[0];
rz(-1.5253541) q[0];
sx q[0];
rz(-0.77950087) q[0];
rz(1.4315073) q[1];
sx q[1];
rz(-1.5305488) q[1];
sx q[1];
rz(0.12758189) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12870991) q[0];
sx q[0];
rz(-1.5361551) q[0];
sx q[0];
rz(-1.4564464) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63795264) q[2];
sx q[2];
rz(-2.0918796) q[2];
sx q[2];
rz(-0.47547728) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7320805) q[1];
sx q[1];
rz(-1.3139551) q[1];
sx q[1];
rz(1.9496296) q[1];
rz(0.74381103) q[3];
sx q[3];
rz(-2.3530745) q[3];
sx q[3];
rz(0.34495993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7369507) q[2];
sx q[2];
rz(-1.5664132) q[2];
sx q[2];
rz(-2.289782) q[2];
rz(-0.24438721) q[3];
sx q[3];
rz(-0.71072018) q[3];
sx q[3];
rz(-1.3062668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39563018) q[0];
sx q[0];
rz(-1.4816477) q[0];
sx q[0];
rz(-0.10598824) q[0];
rz(-1.7164187) q[1];
sx q[1];
rz(-1.1499848) q[1];
sx q[1];
rz(-2.0521767) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7558173) q[0];
sx q[0];
rz(-1.5013278) q[0];
sx q[0];
rz(-1.8088732) q[0];
rz(-pi) q[1];
rz(2.0136497) q[2];
sx q[2];
rz(-0.93292716) q[2];
sx q[2];
rz(1.4952212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.32192595) q[1];
sx q[1];
rz(-1.8882671) q[1];
sx q[1];
rz(3.0078933) q[1];
rz(-pi) q[2];
rz(1.7884729) q[3];
sx q[3];
rz(-2.5158415) q[3];
sx q[3];
rz(1.1661287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.085422903) q[2];
sx q[2];
rz(-1.1916173) q[2];
sx q[2];
rz(0.58758152) q[2];
rz(-1.2299296) q[3];
sx q[3];
rz(-2.7432224) q[3];
sx q[3];
rz(-0.57268322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9375482) q[0];
sx q[0];
rz(-1.2296822) q[0];
sx q[0];
rz(0.44573927) q[0];
rz(2.2988689) q[1];
sx q[1];
rz(-0.80519599) q[1];
sx q[1];
rz(-2.3420948) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2989395) q[0];
sx q[0];
rz(-1.1976722) q[0];
sx q[0];
rz(-1.1375269) q[0];
rz(-2.4386625) q[2];
sx q[2];
rz(-2.0893731) q[2];
sx q[2];
rz(1.2805243) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2098391) q[1];
sx q[1];
rz(-2.5343905) q[1];
sx q[1];
rz(0.74895133) q[1];
rz(-1.6715253) q[3];
sx q[3];
rz(-2.5593649) q[3];
sx q[3];
rz(-1.738172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2004956) q[2];
sx q[2];
rz(-2.7636187) q[2];
sx q[2];
rz(-0.7824347) q[2];
rz(-2.2109219) q[3];
sx q[3];
rz(-1.2717671) q[3];
sx q[3];
rz(-0.29606393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3928669) q[0];
sx q[0];
rz(-0.94806945) q[0];
sx q[0];
rz(1.123708) q[0];
rz(-0.8404845) q[1];
sx q[1];
rz(-1.1437462) q[1];
sx q[1];
rz(-1.7488272) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8776317) q[0];
sx q[0];
rz(-0.11271206) q[0];
sx q[0];
rz(2.0196223) q[0];
rz(-1.7818591) q[2];
sx q[2];
rz(-1.2540134) q[2];
sx q[2];
rz(2.359129) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1906793) q[1];
sx q[1];
rz(-1.2251405) q[1];
sx q[1];
rz(2.9596364) q[1];
rz(2.5244765) q[3];
sx q[3];
rz(-1.502862) q[3];
sx q[3];
rz(-1.2319433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.30354083) q[2];
sx q[2];
rz(-0.76798648) q[2];
sx q[2];
rz(1.5194019) q[2];
rz(-2.0766025) q[3];
sx q[3];
rz(-2.0678554) q[3];
sx q[3];
rz(-0.85552335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4143455) q[0];
sx q[0];
rz(-0.94731826) q[0];
sx q[0];
rz(-1.974768) q[0];
rz(-0.30544063) q[1];
sx q[1];
rz(-2.0716397) q[1];
sx q[1];
rz(2.6397612) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4763757) q[0];
sx q[0];
rz(-1.7656275) q[0];
sx q[0];
rz(-1.4754773) q[0];
rz(-2.3911693) q[2];
sx q[2];
rz(-1.0142676) q[2];
sx q[2];
rz(-1.5518853) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8055685) q[1];
sx q[1];
rz(-0.37266392) q[1];
sx q[1];
rz(1.7737081) q[1];
x q[2];
rz(2.1360399) q[3];
sx q[3];
rz(-1.0985507) q[3];
sx q[3];
rz(2.9173571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66494232) q[2];
sx q[2];
rz(-0.75935894) q[2];
sx q[2];
rz(-0.13258983) q[2];
rz(-2.8442123) q[3];
sx q[3];
rz(-1.5404276) q[3];
sx q[3];
rz(-0.69445777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36841682) q[0];
sx q[0];
rz(-1.4765803) q[0];
sx q[0];
rz(1.6126527) q[0];
rz(-1.348314) q[1];
sx q[1];
rz(-1.3765843) q[1];
sx q[1];
rz(2.7682159) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2397176) q[0];
sx q[0];
rz(-1.4148226) q[0];
sx q[0];
rz(1.9330935) q[0];
rz(-1.8739204) q[2];
sx q[2];
rz(-2.0344007) q[2];
sx q[2];
rz(1.3410717) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5525145) q[1];
sx q[1];
rz(-0.92570549) q[1];
sx q[1];
rz(-2.6721775) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3530497) q[3];
sx q[3];
rz(-1.3065803) q[3];
sx q[3];
rz(-1.4261725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1030964) q[2];
sx q[2];
rz(-3.0615443) q[2];
sx q[2];
rz(-0.45269629) q[2];
rz(-2.2321841) q[3];
sx q[3];
rz(-2.1287287) q[3];
sx q[3];
rz(-2.4694841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.4368923) q[0];
sx q[0];
rz(-1.9743275) q[0];
sx q[0];
rz(-2.0335017) q[0];
rz(-2.3469901) q[1];
sx q[1];
rz(-1.4789076) q[1];
sx q[1];
rz(2.2750003) q[1];
rz(3.1187155) q[2];
sx q[2];
rz(-0.87999095) q[2];
sx q[2];
rz(-2.8578514) q[2];
rz(-2.361479) q[3];
sx q[3];
rz(-1.2419392) q[3];
sx q[3];
rz(1.8942647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
