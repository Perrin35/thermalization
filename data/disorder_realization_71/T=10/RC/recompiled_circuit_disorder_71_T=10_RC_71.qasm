OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61385566) q[0];
sx q[0];
rz(-1.6439438) q[0];
sx q[0];
rz(-0.82984501) q[0];
rz(3.9217477) q[1];
sx q[1];
rz(5.2182066) q[1];
sx q[1];
rz(10.301104) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012955879) q[0];
sx q[0];
rz(-0.46015938) q[0];
sx q[0];
rz(-0.53928661) q[0];
rz(-pi) q[1];
rz(2.949602) q[2];
sx q[2];
rz(-2.1477094) q[2];
sx q[2];
rz(-0.38919762) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6688924) q[1];
sx q[1];
rz(-1.2408222) q[1];
sx q[1];
rz(-3.1256691) q[1];
rz(-pi) q[2];
rz(1.7238486) q[3];
sx q[3];
rz(-1.9733841) q[3];
sx q[3];
rz(0.11868417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17065389) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(2.412964) q[2];
rz(-2.6206) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(2.9339824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.306863) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(2.0200502) q[0];
rz(0.25575486) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(-2.2671525) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16337285) q[0];
sx q[0];
rz(-1.997943) q[0];
sx q[0];
rz(0.33421974) q[0];
rz(0.99526694) q[2];
sx q[2];
rz(-1.4412291) q[2];
sx q[2];
rz(2.1339983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9718711) q[1];
sx q[1];
rz(-1.1217146) q[1];
sx q[1];
rz(2.4596618) q[1];
x q[2];
rz(-1.8401243) q[3];
sx q[3];
rz(-1.519031) q[3];
sx q[3];
rz(-0.89382899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4008537) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(-2.7056616) q[2];
rz(-2.46051) q[3];
sx q[3];
rz(-0.77107945) q[3];
sx q[3];
rz(2.7387103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9044559) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(-2.0667734) q[0];
rz(0.83956051) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(2.7456465) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9112644) q[0];
sx q[0];
rz(-2.0718144) q[0];
sx q[0];
rz(2.688174) q[0];
rz(1.6171574) q[2];
sx q[2];
rz(-2.8472387) q[2];
sx q[2];
rz(0.72325318) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44569762) q[1];
sx q[1];
rz(-2.2539415) q[1];
sx q[1];
rz(1.5831069) q[1];
x q[2];
rz(1.5941761) q[3];
sx q[3];
rz(-1.8262719) q[3];
sx q[3];
rz(1.3083096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5376771) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(2.9023857) q[2];
rz(-0.075332969) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199715) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(2.3216632) q[0];
rz(-0.48768249) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(0.23342361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4963835) q[0];
sx q[0];
rz(-1.6312946) q[0];
sx q[0];
rz(1.4738826) q[0];
x q[1];
rz(2.8308224) q[2];
sx q[2];
rz(-0.66590532) q[2];
sx q[2];
rz(-3.0096465) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.64775002) q[1];
sx q[1];
rz(-0.88652459) q[1];
sx q[1];
rz(-1.1286331) q[1];
x q[2];
rz(-0.38846429) q[3];
sx q[3];
rz(-2.3187175) q[3];
sx q[3];
rz(2.8188761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0908115) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(-2.3941669) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(-0.24211611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500047) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(3.0773556) q[0];
rz(2.1977987) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(-2.3805526) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8262186) q[0];
sx q[0];
rz(-1.5407019) q[0];
sx q[0];
rz(-1.8316359) q[0];
rz(1.7663076) q[2];
sx q[2];
rz(-1.3720781) q[2];
sx q[2];
rz(-0.18713258) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4796175) q[1];
sx q[1];
rz(-1.8912589) q[1];
sx q[1];
rz(-1.596405) q[1];
rz(-1.8499591) q[3];
sx q[3];
rz(-1.8431292) q[3];
sx q[3];
rz(0.84701049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.80660194) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(-0.09207329) q[2];
rz(2.4798685) q[3];
sx q[3];
rz(-0.79143733) q[3];
sx q[3];
rz(1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-2.81566) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.35615) q[0];
sx q[0];
rz(-2.4768562) q[0];
sx q[0];
rz(0.63389969) q[0];
rz(-2.5599401) q[2];
sx q[2];
rz(-0.51917167) q[2];
sx q[2];
rz(2.4793064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6335878) q[1];
sx q[1];
rz(-1.272401) q[1];
sx q[1];
rz(1.7119346) q[1];
x q[2];
rz(-1.5876706) q[3];
sx q[3];
rz(-1.6718739) q[3];
sx q[3];
rz(-0.24731393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5439593) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(-2.4874172) q[2];
rz(-1.7116961) q[3];
sx q[3];
rz(-1.1681898) q[3];
sx q[3];
rz(2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8687246) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-2.4196999) q[0];
rz(1.4121274) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(3.022335) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3463466) q[0];
sx q[0];
rz(-2.3756785) q[0];
sx q[0];
rz(2.7555097) q[0];
rz(3.0316333) q[2];
sx q[2];
rz(-1.409515) q[2];
sx q[2];
rz(1.9732628) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1489361) q[1];
sx q[1];
rz(-1.7033556) q[1];
sx q[1];
rz(1.0692783) q[1];
x q[2];
rz(1.9414385) q[3];
sx q[3];
rz(-2.2439085) q[3];
sx q[3];
rz(-0.17236575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44935903) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(-1.7822441) q[2];
rz(-0.75891495) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(-0.53708491) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3100202) q[0];
sx q[0];
rz(-0.65905237) q[0];
sx q[0];
rz(-2.0781562) q[0];
rz(0.27451441) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(0.88561052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6016156) q[0];
sx q[0];
rz(-0.5479387) q[0];
sx q[0];
rz(0.45666306) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1503259) q[2];
sx q[2];
rz(-1.5838924) q[2];
sx q[2];
rz(-3.1377369) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.209621) q[1];
sx q[1];
rz(-2.0721657) q[1];
sx q[1];
rz(-2.6905836) q[1];
x q[2];
rz(0.76836821) q[3];
sx q[3];
rz(-2.1185015) q[3];
sx q[3];
rz(0.87793575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4132335) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(1.1317066) q[2];
rz(-1.0845832) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(-1.2148946) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8383012) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(-1.0820748) q[0];
rz(1.2754296) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(1.1358322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95490676) q[0];
sx q[0];
rz(-1.0404772) q[0];
sx q[0];
rz(-2.9493939) q[0];
x q[1];
rz(0.047770569) q[2];
sx q[2];
rz(-2.7063745) q[2];
sx q[2];
rz(-1.4620632) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.51582) q[1];
sx q[1];
rz(-1.9966218) q[1];
sx q[1];
rz(-1.5172525) q[1];
rz(-2.2215861) q[3];
sx q[3];
rz(-1.2244867) q[3];
sx q[3];
rz(-1.8370093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0231126) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(-2.2686968) q[2];
rz(-2.2980799) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(2.9516454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2492367) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(1.2783485) q[0];
rz(-1.0247914) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(-1.1970253) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68162936) q[0];
sx q[0];
rz(-2.3869793) q[0];
sx q[0];
rz(2.2370536) q[0];
rz(-3.0212101) q[2];
sx q[2];
rz(-1.8070081) q[2];
sx q[2];
rz(-1.2506968) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.69233209) q[1];
sx q[1];
rz(-0.83220607) q[1];
sx q[1];
rz(0.87528054) q[1];
rz(-2.6322369) q[3];
sx q[3];
rz(-1.2378113) q[3];
sx q[3];
rz(1.8978564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3060351) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(-0.79997921) q[2];
rz(1.9647313) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(-2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.4326614) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(-2.6196383) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(1.21571) q[2];
sx q[2];
rz(-1.9234895) q[2];
sx q[2];
rz(-1.1870155) q[2];
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