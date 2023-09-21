OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(2.2709742) q[0];
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38116954) q[0];
sx q[0];
rz(-1.086735) q[0];
sx q[0];
rz(2.3056187) q[0];
rz(0.50617354) q[2];
sx q[2];
rz(-2.0548327) q[2];
sx q[2];
rz(1.6137705) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95853165) q[1];
sx q[1];
rz(-1.6089988) q[1];
sx q[1];
rz(1.3903635) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9505694) q[3];
sx q[3];
rz(-1.028562) q[3];
sx q[3];
rz(-2.1477826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8831138) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(2.1885833) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54884058) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(0.63252226) q[0];
rz(2.6951492) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(2.4893563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3632293) q[0];
sx q[0];
rz(-3.1104381) q[0];
sx q[0];
rz(-0.40701436) q[0];
rz(1.863443) q[2];
sx q[2];
rz(-1.7990944) q[2];
sx q[2];
rz(-1.4676859) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3995041) q[1];
sx q[1];
rz(-1.3474476) q[1];
sx q[1];
rz(2.8469267) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90918031) q[3];
sx q[3];
rz(-1.0717857) q[3];
sx q[3];
rz(0.0024851174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5941045) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(1.1616421) q[2];
rz(1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(-1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(2.2913349) q[0];
rz(2.6440874) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(1.3495548) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.772086) q[0];
sx q[0];
rz(-1.3319351) q[0];
sx q[0];
rz(0.68349616) q[0];
rz(-1.6367958) q[2];
sx q[2];
rz(-1.5335576) q[2];
sx q[2];
rz(-0.50022349) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9393443) q[1];
sx q[1];
rz(-1.0679686) q[1];
sx q[1];
rz(1.3015675) q[1];
rz(-0.0098185929) q[3];
sx q[3];
rz(-1.9968642) q[3];
sx q[3];
rz(0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0321908) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(0.90399495) q[2];
rz(2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(1.3999456) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49144739) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(1.7310463) q[0];
rz(-0.63181216) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(0.036380336) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0050126652) q[0];
sx q[0];
rz(-0.67626017) q[0];
sx q[0];
rz(-1.2269292) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58156275) q[2];
sx q[2];
rz(-1.2956603) q[2];
sx q[2];
rz(3.0951701) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1069146) q[1];
sx q[1];
rz(-0.87768302) q[1];
sx q[1];
rz(0.17130674) q[1];
x q[2];
rz(-0.39156885) q[3];
sx q[3];
rz(-0.25300004) q[3];
sx q[3];
rz(-2.5391425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2146384) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(1.9533763) q[2];
rz(0.67048091) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(-2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(0.16648509) q[0];
rz(2.3855551) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(2.9072445) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74041022) q[0];
sx q[0];
rz(-2.4966842) q[0];
sx q[0];
rz(2.5572204) q[0];
rz(-pi) q[1];
rz(-0.095586153) q[2];
sx q[2];
rz(-1.0786973) q[2];
sx q[2];
rz(-2.0795859) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.25917945) q[1];
sx q[1];
rz(-2.3098574) q[1];
sx q[1];
rz(-1.3684567) q[1];
rz(-pi) q[2];
rz(-1.2469532) q[3];
sx q[3];
rz(-1.4959335) q[3];
sx q[3];
rz(-0.74136855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70871893) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(-0.43705127) q[3];
sx q[3];
rz(-2.1190937) q[3];
sx q[3];
rz(-2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41546145) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-0.81714001) q[0];
rz(2.5754886) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(-1.9979427) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4139347) q[0];
sx q[0];
rz(-1.2000788) q[0];
sx q[0];
rz(1.6137705) q[0];
x q[1];
rz(0.67310682) q[2];
sx q[2];
rz(-1.7761201) q[2];
sx q[2];
rz(-0.38773195) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5864582) q[1];
sx q[1];
rz(-2.0124334) q[1];
sx q[1];
rz(-0.099979062) q[1];
rz(-pi) q[2];
rz(-2.6236344) q[3];
sx q[3];
rz(-1.3266139) q[3];
sx q[3];
rz(1.1740008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-0.4894408) q[2];
rz(-2.9135381) q[3];
sx q[3];
rz(-1.8830048) q[3];
sx q[3];
rz(-2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5722826) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.8547159) q[0];
rz(2.4781748) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(-1.2333262) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9702643) q[0];
sx q[0];
rz(-1.5597222) q[0];
sx q[0];
rz(-1.1867255) q[0];
rz(-pi) q[1];
rz(2.5160518) q[2];
sx q[2];
rz(-0.95226804) q[2];
sx q[2];
rz(-2.2369838) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8804633) q[1];
sx q[1];
rz(-1.5848586) q[1];
sx q[1];
rz(-3.0304099) q[1];
rz(2.535378) q[3];
sx q[3];
rz(-0.57176916) q[3];
sx q[3];
rz(-1.9619463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.104091) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(2.1255778) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0163517) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(-1.4455147) q[0];
rz(-2.9267172) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(-1.258237) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94851516) q[0];
sx q[0];
rz(-2.6876039) q[0];
sx q[0];
rz(1.8817188) q[0];
rz(-pi) q[1];
rz(-2.4969205) q[2];
sx q[2];
rz(-1.03627) q[2];
sx q[2];
rz(2.1246186) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82842365) q[1];
sx q[1];
rz(-2.1143882) q[1];
sx q[1];
rz(-1.4491175) q[1];
rz(-pi) q[2];
rz(2.6120841) q[3];
sx q[3];
rz(-2.2141075) q[3];
sx q[3];
rz(0.95638004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4593279) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-0.56274596) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(-2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257618) q[0];
sx q[0];
rz(-2.0985726) q[0];
sx q[0];
rz(-0.2510221) q[0];
rz(0.42731467) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(-0.02773157) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6513034) q[0];
sx q[0];
rz(-2.2130744) q[0];
sx q[0];
rz(0.62255967) q[0];
x q[1];
rz(2.2698195) q[2];
sx q[2];
rz(-1.7773526) q[2];
sx q[2];
rz(0.86824647) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9412649) q[1];
sx q[1];
rz(-0.97201921) q[1];
sx q[1];
rz(1.6914678) q[1];
rz(-1.2003044) q[3];
sx q[3];
rz(-3.070188) q[3];
sx q[3];
rz(2.2927473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1256844) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(-1.998418) q[2];
rz(-0.14287359) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(-2.4699396) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(-2.8840816) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2430902) q[0];
sx q[0];
rz(-1.1872963) q[0];
sx q[0];
rz(2.6570508) q[0];
rz(1.0742513) q[2];
sx q[2];
rz(-2.0601344) q[2];
sx q[2];
rz(0.41000965) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32138667) q[1];
sx q[1];
rz(-2.4452129) q[1];
sx q[1];
rz(-0.022298261) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0912283) q[3];
sx q[3];
rz(-1.3551095) q[3];
sx q[3];
rz(1.2558162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8132849) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(2.5349687) q[2];
rz(2.666752) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(0.84038466) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7941147) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(2.9150302) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(1.1889585) q[2];
sx q[2];
rz(-1.1321862) q[2];
sx q[2];
rz(1.95375) q[2];
rz(2.0541035) q[3];
sx q[3];
rz(-1.3369505) q[3];
sx q[3];
rz(-1.7575775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
