OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0680024) q[0];
sx q[0];
rz(-1.468714) q[0];
sx q[0];
rz(0.52748632) q[0];
rz(0.5440076) q[1];
sx q[1];
rz(3.5825621) q[1];
sx q[1];
rz(9.4234186) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5088592) q[0];
sx q[0];
rz(-1.5799045) q[0];
sx q[0];
rz(-1.5004116) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5293526) q[2];
sx q[2];
rz(-1.3748825) q[2];
sx q[2];
rz(1.4961835) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2058918) q[1];
sx q[1];
rz(-2.08826) q[1];
sx q[1];
rz(-2.0175354) q[1];
x q[2];
rz(-1.6) q[3];
sx q[3];
rz(-1.7701704) q[3];
sx q[3];
rz(-2.3674611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.15343836) q[2];
sx q[2];
rz(-0.64186382) q[2];
sx q[2];
rz(1.362907) q[2];
rz(-3.0617132) q[3];
sx q[3];
rz(-2.2814543) q[3];
sx q[3];
rz(3.1352654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7774571) q[0];
sx q[0];
rz(-1.2292925) q[0];
sx q[0];
rz(0.51032132) q[0];
rz(1.8136884) q[1];
sx q[1];
rz(-2.4102305) q[1];
sx q[1];
rz(-2.7046943) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53117673) q[0];
sx q[0];
rz(-2.057909) q[0];
sx q[0];
rz(0.25419828) q[0];
rz(-pi) q[1];
rz(-2.6742009) q[2];
sx q[2];
rz(-0.49023489) q[2];
sx q[2];
rz(0.64783123) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4480312) q[1];
sx q[1];
rz(-2.3400819) q[1];
sx q[1];
rz(1.526725) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0776229) q[3];
sx q[3];
rz(-1.973458) q[3];
sx q[3];
rz(-0.10884604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3029311) q[2];
sx q[2];
rz(-0.51753664) q[2];
sx q[2];
rz(-0.17520629) q[2];
rz(-0.13904275) q[3];
sx q[3];
rz(-1.6486721) q[3];
sx q[3];
rz(2.2429332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43056968) q[0];
sx q[0];
rz(-0.88453203) q[0];
sx q[0];
rz(2.7810466) q[0];
rz(1.644246) q[1];
sx q[1];
rz(-1.2253573) q[1];
sx q[1];
rz(2.5490733) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07537341) q[0];
sx q[0];
rz(-2.1751715) q[0];
sx q[0];
rz(-2.3492667) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63590598) q[2];
sx q[2];
rz(-1.4472258) q[2];
sx q[2];
rz(-1.4452782) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4289817) q[1];
sx q[1];
rz(-1.4498222) q[1];
sx q[1];
rz(1.9209204) q[1];
rz(-pi) q[2];
rz(2.6733562) q[3];
sx q[3];
rz(-1.1854725) q[3];
sx q[3];
rz(1.6637529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3561919) q[2];
sx q[2];
rz(-2.5742026) q[2];
sx q[2];
rz(0.94245911) q[2];
rz(2.5475907) q[3];
sx q[3];
rz(-2.3022251) q[3];
sx q[3];
rz(0.12282898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.9451611) q[0];
sx q[0];
rz(-0.62569797) q[0];
sx q[0];
rz(-0.64836597) q[0];
rz(2.2658589) q[1];
sx q[1];
rz(-1.6301194) q[1];
sx q[1];
rz(-1.7902364) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195456) q[0];
sx q[0];
rz(-0.81455671) q[0];
sx q[0];
rz(-3.087939) q[0];
rz(-pi) q[1];
rz(0.8114795) q[2];
sx q[2];
rz(-1.4373903) q[2];
sx q[2];
rz(1.8638944) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9571805) q[1];
sx q[1];
rz(-2.4676128) q[1];
sx q[1];
rz(1.5247702) q[1];
x q[2];
rz(-1.625678) q[3];
sx q[3];
rz(-1.1883231) q[3];
sx q[3];
rz(-0.11688133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15328345) q[2];
sx q[2];
rz(-0.51258665) q[2];
sx q[2];
rz(-1.886606) q[2];
rz(-2.5399688) q[3];
sx q[3];
rz(-2.2646077) q[3];
sx q[3];
rz(-1.4582483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0772142) q[0];
sx q[0];
rz(-2.0266396) q[0];
sx q[0];
rz(-0.13378046) q[0];
rz(1.6556219) q[1];
sx q[1];
rz(-2.8120815) q[1];
sx q[1];
rz(-2.9697184) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3456589) q[0];
sx q[0];
rz(-2.5194114) q[0];
sx q[0];
rz(0.90779974) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36891035) q[2];
sx q[2];
rz(-2.2930573) q[2];
sx q[2];
rz(1.5992355) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50087029) q[1];
sx q[1];
rz(-1.7224887) q[1];
sx q[1];
rz(1.2176563) q[1];
rz(-pi) q[2];
rz(-1.9952568) q[3];
sx q[3];
rz(-2.0133284) q[3];
sx q[3];
rz(-2.8660021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6493426) q[2];
sx q[2];
rz(-2.616373) q[2];
sx q[2];
rz(1.4225175) q[2];
rz(-0.20497841) q[3];
sx q[3];
rz(-0.96115464) q[3];
sx q[3];
rz(-2.3375296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.499046) q[0];
sx q[0];
rz(-2.5048984) q[0];
sx q[0];
rz(-0.194304) q[0];
rz(2.4957472) q[1];
sx q[1];
rz(-1.9648809) q[1];
sx q[1];
rz(-0.35020721) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98726942) q[0];
sx q[0];
rz(-1.7065757) q[0];
sx q[0];
rz(2.6063347) q[0];
rz(-0.62927134) q[2];
sx q[2];
rz(-2.2221178) q[2];
sx q[2];
rz(0.68683147) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3694075) q[1];
sx q[1];
rz(-2.6168452) q[1];
sx q[1];
rz(2.9725542) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.086669162) q[3];
sx q[3];
rz(-2.4695167) q[3];
sx q[3];
rz(-1.8928526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.084006) q[2];
sx q[2];
rz(-0.39613327) q[2];
sx q[2];
rz(0.32220379) q[2];
rz(2.3816439) q[3];
sx q[3];
rz(-2.6252169) q[3];
sx q[3];
rz(3.0907536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8754804) q[0];
sx q[0];
rz(-2.4732944) q[0];
sx q[0];
rz(-0.71568263) q[0];
rz(-3.074805) q[1];
sx q[1];
rz(-2.5804434) q[1];
sx q[1];
rz(-0.36398789) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2192236) q[0];
sx q[0];
rz(-0.69747335) q[0];
sx q[0];
rz(-2.1261901) q[0];
rz(-pi) q[1];
rz(-2.5189713) q[2];
sx q[2];
rz(-1.004815) q[2];
sx q[2];
rz(-0.37812585) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5753117) q[1];
sx q[1];
rz(-0.76427459) q[1];
sx q[1];
rz(-3.1180095) q[1];
rz(-pi) q[2];
rz(-0.19192275) q[3];
sx q[3];
rz(-2.4867184) q[3];
sx q[3];
rz(0.45134967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62465847) q[2];
sx q[2];
rz(-0.90183574) q[2];
sx q[2];
rz(1.2348403) q[2];
rz(0.45461795) q[3];
sx q[3];
rz(-1.8574628) q[3];
sx q[3];
rz(3.0668018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5816077) q[0];
sx q[0];
rz(-0.11542628) q[0];
sx q[0];
rz(0.64062947) q[0];
rz(0.36264125) q[1];
sx q[1];
rz(-2.8022712) q[1];
sx q[1];
rz(-1.6222662) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59121038) q[0];
sx q[0];
rz(-1.2347295) q[0];
sx q[0];
rz(2.1318011) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2799543) q[2];
sx q[2];
rz(-2.2802326) q[2];
sx q[2];
rz(2.4091313) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2920529) q[1];
sx q[1];
rz(-0.89703945) q[1];
sx q[1];
rz(-0.11815355) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0985116) q[3];
sx q[3];
rz(-1.660822) q[3];
sx q[3];
rz(1.2028026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92056876) q[2];
sx q[2];
rz(-1.2732882) q[2];
sx q[2];
rz(-1.1576687) q[2];
rz(-0.32933346) q[3];
sx q[3];
rz(-0.18940997) q[3];
sx q[3];
rz(-2.2277189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90410239) q[0];
sx q[0];
rz(-2.3512023) q[0];
sx q[0];
rz(2.7303586) q[0];
rz(-2.5143738) q[1];
sx q[1];
rz(-2.5721481) q[1];
sx q[1];
rz(-1.2818744) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3498308) q[0];
sx q[0];
rz(-1.786953) q[0];
sx q[0];
rz(2.7695275) q[0];
x q[1];
rz(0.55686624) q[2];
sx q[2];
rz(-1.0039731) q[2];
sx q[2];
rz(-1.8447529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5465849) q[1];
sx q[1];
rz(-0.91239415) q[1];
sx q[1];
rz(-2.9987572) q[1];
rz(-pi) q[2];
rz(-2.0342213) q[3];
sx q[3];
rz(-2.7181566) q[3];
sx q[3];
rz(0.23611072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3605448) q[2];
sx q[2];
rz(-1.5534399) q[2];
sx q[2];
rz(0.9786728) q[2];
rz(0.22748889) q[3];
sx q[3];
rz(-0.74930185) q[3];
sx q[3];
rz(2.6849875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6351673) q[0];
sx q[0];
rz(-0.058110617) q[0];
sx q[0];
rz(2.3458922) q[0];
rz(-2.6129163) q[1];
sx q[1];
rz(-2.1541336) q[1];
sx q[1];
rz(-2.9472369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07819019) q[0];
sx q[0];
rz(-1.9424129) q[0];
sx q[0];
rz(0.097784575) q[0];
x q[1];
rz(-1.8647381) q[2];
sx q[2];
rz(-1.1938098) q[2];
sx q[2];
rz(1.3177208) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4116507) q[1];
sx q[1];
rz(-0.46272181) q[1];
sx q[1];
rz(-0.28280854) q[1];
rz(0.50823516) q[3];
sx q[3];
rz(-2.0401032) q[3];
sx q[3];
rz(-1.5526916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.84870321) q[2];
sx q[2];
rz(-1.6699426) q[2];
sx q[2];
rz(0.29937747) q[2];
rz(2.7401466) q[3];
sx q[3];
rz(-0.58599389) q[3];
sx q[3];
rz(-2.6011023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.9837579) q[0];
sx q[0];
rz(-2.3102289) q[0];
sx q[0];
rz(1.8752444) q[0];
rz(-0.10706317) q[1];
sx q[1];
rz(-2.3657847) q[1];
sx q[1];
rz(1.8484144) q[1];
rz(-0.54218311) q[2];
sx q[2];
rz(-1.6169006) q[2];
sx q[2];
rz(-2.9666881) q[2];
rz(-0.43184256) q[3];
sx q[3];
rz(-2.8819537) q[3];
sx q[3];
rz(1.0859539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
