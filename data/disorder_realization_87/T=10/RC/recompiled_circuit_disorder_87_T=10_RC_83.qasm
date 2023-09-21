OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23624578) q[0];
sx q[0];
rz(-2.4155004) q[0];
sx q[0];
rz(0.2015764) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(2.6013241) q[1];
sx q[1];
rz(7.2202914) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21298458) q[0];
sx q[0];
rz(-1.0092508) q[0];
sx q[0];
rz(1.0667849) q[0];
rz(-pi) q[1];
rz(-2.5295528) q[2];
sx q[2];
rz(-1.9844374) q[2];
sx q[2];
rz(-0.94758247) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5697437) q[1];
sx q[1];
rz(-1.7535926) q[1];
sx q[1];
rz(-0.079449541) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9107781) q[3];
sx q[3];
rz(-2.8570606) q[3];
sx q[3];
rz(1.9576548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15930882) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(3.0554331) q[2];
rz(-2.384095) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(-2.0479726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57698292) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(-1.8151059) q[0];
rz(1.8857229) q[1];
sx q[1];
rz(-1.565226) q[1];
sx q[1];
rz(2.870141) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4685681) q[0];
sx q[0];
rz(-2.4789171) q[0];
sx q[0];
rz(2.2492692) q[0];
rz(-pi) q[1];
rz(-0.51606744) q[2];
sx q[2];
rz(-1.8415673) q[2];
sx q[2];
rz(-0.64112907) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2084864) q[1];
sx q[1];
rz(-0.83175627) q[1];
sx q[1];
rz(-2.0887124) q[1];
rz(-pi) q[2];
rz(-2.7852887) q[3];
sx q[3];
rz(-0.84173991) q[3];
sx q[3];
rz(-0.19405288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0469971) q[2];
sx q[2];
rz(-1.3826028) q[2];
sx q[2];
rz(-0.57717741) q[2];
rz(-0.92352891) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(-1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975824) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(-0.14257167) q[0];
rz(1.7890731) q[1];
sx q[1];
rz(-1.0357772) q[1];
sx q[1];
rz(-2.9325063) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3554879) q[0];
sx q[0];
rz(-0.84083637) q[0];
sx q[0];
rz(0.79746042) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4019037) q[2];
sx q[2];
rz(-1.2581173) q[2];
sx q[2];
rz(0.82175053) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8612954) q[1];
sx q[1];
rz(-1.1251083) q[1];
sx q[1];
rz(0.46511005) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61061065) q[3];
sx q[3];
rz(-2.7770677) q[3];
sx q[3];
rz(0.025346905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3645939) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(0.40412942) q[2];
rz(1.2858307) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053112) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(2.1787815) q[0];
rz(0.46936938) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(-0.00096360047) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9199333) q[0];
sx q[0];
rz(-2.1533305) q[0];
sx q[0];
rz(-2.4525053) q[0];
rz(-pi) q[1];
rz(0.12400603) q[2];
sx q[2];
rz(-1.5713072) q[2];
sx q[2];
rz(2.309547) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8648659) q[1];
sx q[1];
rz(-2.556567) q[1];
sx q[1];
rz(1.3328972) q[1];
rz(-pi) q[2];
rz(2.5721781) q[3];
sx q[3];
rz(-1.6462407) q[3];
sx q[3];
rz(-0.30573341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0659539) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(0.11165079) q[2];
rz(-0.81104898) q[3];
sx q[3];
rz(-0.45447293) q[3];
sx q[3];
rz(-0.013899175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1902996) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(-2.7850889) q[0];
rz(-2.6351392) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(-0.26062632) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60359943) q[0];
sx q[0];
rz(-1.7126181) q[0];
sx q[0];
rz(-1.6130559) q[0];
x q[1];
rz(-1.8243276) q[2];
sx q[2];
rz(-1.3273444) q[2];
sx q[2];
rz(1.7340811) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.645748) q[1];
sx q[1];
rz(-2.145088) q[1];
sx q[1];
rz(-1.5917642) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7138249) q[3];
sx q[3];
rz(-1.7149394) q[3];
sx q[3];
rz(-0.10728697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2300718) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(2.9122706) q[2];
rz(0.54245943) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(-1.0361766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77263537) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(-1.649958) q[0];
rz(2.1024599) q[1];
sx q[1];
rz(-1.8455448) q[1];
sx q[1];
rz(-1.414149) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0522239) q[0];
sx q[0];
rz(-1.5234689) q[0];
sx q[0];
rz(2.8663859) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51775198) q[2];
sx q[2];
rz(-1.3377829) q[2];
sx q[2];
rz(2.0875967) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27281877) q[1];
sx q[1];
rz(-2.0588074) q[1];
sx q[1];
rz(-0.52656071) q[1];
x q[2];
rz(-1.3544975) q[3];
sx q[3];
rz(-1.4970475) q[3];
sx q[3];
rz(-2.4019965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1386537) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(0.49267832) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(-1.7475351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7572927) q[0];
sx q[0];
rz(-1.5889656) q[0];
sx q[0];
rz(-0.12167715) q[0];
rz(-1.9901468) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(-2.7391403) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50169045) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(0.34696607) q[0];
x q[1];
rz(0.17860883) q[2];
sx q[2];
rz(-3.0948234) q[2];
sx q[2];
rz(-0.59431078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7546141) q[1];
sx q[1];
rz(-1.3941892) q[1];
sx q[1];
rz(2.8547924) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6472858) q[3];
sx q[3];
rz(-0.34919958) q[3];
sx q[3];
rz(-0.66222144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1272614) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(0.86223117) q[2];
rz(-0.47752738) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(2.2235218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35571337) q[0];
sx q[0];
rz(-2.9861351) q[0];
sx q[0];
rz(2.4086337) q[0];
rz(3.0015302) q[1];
sx q[1];
rz(-0.99761325) q[1];
sx q[1];
rz(-1.0345116) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2037172) q[0];
sx q[0];
rz(-0.28521842) q[0];
sx q[0];
rz(2.2144149) q[0];
rz(-pi) q[1];
rz(-0.60136749) q[2];
sx q[2];
rz(-3.0172918) q[2];
sx q[2];
rz(-2.3483495) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5518783) q[1];
sx q[1];
rz(-1.0052181) q[1];
sx q[1];
rz(0.043885529) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41820742) q[3];
sx q[3];
rz(-1.4718664) q[3];
sx q[3];
rz(-1.6855155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90116477) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(2.365716) q[2];
rz(-2.4173229) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(1.4340713) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6999321) q[0];
sx q[0];
rz(-0.76403809) q[0];
sx q[0];
rz(3.124776) q[0];
rz(0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-0.7787849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9433702) q[0];
sx q[0];
rz(-1.3226042) q[0];
sx q[0];
rz(2.7944837) q[0];
x q[1];
rz(-2.8674556) q[2];
sx q[2];
rz(-2.0275653) q[2];
sx q[2];
rz(2.0704839) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2494431) q[1];
sx q[1];
rz(-1.3378694) q[1];
sx q[1];
rz(-1.1942785) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0770244) q[3];
sx q[3];
rz(-1.350292) q[3];
sx q[3];
rz(-1.2793503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(-2.4251535) q[2];
rz(1.6843494) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0338106) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-0.15144908) q[0];
rz(-2.9653213) q[1];
sx q[1];
rz(-1.1947894) q[1];
sx q[1];
rz(0.72296468) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13847362) q[0];
sx q[0];
rz(-1.6771294) q[0];
sx q[0];
rz(1.4195725) q[0];
rz(-pi) q[1];
rz(0.16913551) q[2];
sx q[2];
rz(-1.198487) q[2];
sx q[2];
rz(-1.7793836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42230095) q[1];
sx q[1];
rz(-2.0740168) q[1];
sx q[1];
rz(-0.06312381) q[1];
x q[2];
rz(-0.39977269) q[3];
sx q[3];
rz(-2.3120566) q[3];
sx q[3];
rz(1.6403891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4621949) q[2];
sx q[2];
rz(-1.2827736) q[2];
sx q[2];
rz(-0.52552137) q[2];
rz(0.28371352) q[3];
sx q[3];
rz(-2.0339537) q[3];
sx q[3];
rz(2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491966) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(2.0422968) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(0.47808403) q[2];
sx q[2];
rz(-2.1838084) q[2];
sx q[2];
rz(2.9391391) q[2];
rz(-1.7952193) q[3];
sx q[3];
rz(-1.8891469) q[3];
sx q[3];
rz(-0.3286152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];