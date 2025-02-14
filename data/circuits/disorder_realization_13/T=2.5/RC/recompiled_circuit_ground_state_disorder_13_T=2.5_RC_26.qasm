OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90091997) q[0];
sx q[0];
rz(2.9958041) q[0];
sx q[0];
rz(11.267405) q[0];
rz(2.9453912) q[1];
sx q[1];
rz(-1.8008404) q[1];
sx q[1];
rz(0.10032108) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31831384) q[0];
sx q[0];
rz(-0.66368503) q[0];
sx q[0];
rz(0.084245988) q[0];
x q[1];
rz(1.2353363) q[2];
sx q[2];
rz(-2.0192207) q[2];
sx q[2];
rz(-1.264073) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92857526) q[1];
sx q[1];
rz(-1.8452106) q[1];
sx q[1];
rz(-1.7875415) q[1];
x q[2];
rz(-0.45083188) q[3];
sx q[3];
rz(-0.63967645) q[3];
sx q[3];
rz(-1.0420639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7413062) q[2];
sx q[2];
rz(-2.9897959) q[2];
sx q[2];
rz(-2.3347704) q[2];
rz(2.412879) q[3];
sx q[3];
rz(-2.3829134) q[3];
sx q[3];
rz(-1.2046643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5173986) q[0];
sx q[0];
rz(-2.875858) q[0];
sx q[0];
rz(0.30701315) q[0];
rz(-1.864805) q[1];
sx q[1];
rz(-2.0025608) q[1];
sx q[1];
rz(2.0397287) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66070304) q[0];
sx q[0];
rz(-1.0202978) q[0];
sx q[0];
rz(0.20396978) q[0];
x q[1];
rz(-0.17942682) q[2];
sx q[2];
rz(-1.5080001) q[2];
sx q[2];
rz(-1.6321317) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.700775) q[1];
sx q[1];
rz(-1.3937181) q[1];
sx q[1];
rz(3.0040279) q[1];
rz(0.78897055) q[3];
sx q[3];
rz(-1.0010825) q[3];
sx q[3];
rz(-0.33984847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.028194204) q[2];
sx q[2];
rz(-2.5598309) q[2];
sx q[2];
rz(3.0451575) q[2];
rz(0.15245572) q[3];
sx q[3];
rz(-1.5072482) q[3];
sx q[3];
rz(0.69795394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089712791) q[0];
sx q[0];
rz(-1.1002325) q[0];
sx q[0];
rz(2.7413947) q[0];
rz(-1.6290889) q[1];
sx q[1];
rz(-0.17833231) q[1];
sx q[1];
rz(-0.2581183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5921051) q[0];
sx q[0];
rz(-1.7096115) q[0];
sx q[0];
rz(-1.8068683) q[0];
rz(-pi) q[1];
rz(-2.4674008) q[2];
sx q[2];
rz(-3.0140244) q[2];
sx q[2];
rz(-1.3889165) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5096881) q[1];
sx q[1];
rz(-0.013876112) q[1];
sx q[1];
rz(-0.3664151) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24127023) q[3];
sx q[3];
rz(-1.4033969) q[3];
sx q[3];
rz(2.3679855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1199946) q[2];
sx q[2];
rz(-1.9436516) q[2];
sx q[2];
rz(3.1094587) q[2];
rz(2.8739127) q[3];
sx q[3];
rz(-1.6240424) q[3];
sx q[3];
rz(1.5599498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.5716008) q[0];
sx q[0];
rz(-1.5084234) q[0];
sx q[0];
rz(-3.0573523) q[0];
rz(-0.035942297) q[1];
sx q[1];
rz(-0.033999559) q[1];
sx q[1];
rz(2.8004004) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0628898) q[0];
sx q[0];
rz(-0.69403946) q[0];
sx q[0];
rz(-2.7258123) q[0];
x q[1];
rz(-2.555055) q[2];
sx q[2];
rz(-1.7866724) q[2];
sx q[2];
rz(3.0160273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2349639) q[1];
sx q[1];
rz(-2.2412657) q[1];
sx q[1];
rz(-2.5690298) q[1];
rz(-pi) q[2];
rz(0.61670035) q[3];
sx q[3];
rz(-2.4444207) q[3];
sx q[3];
rz(0.32000772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42698947) q[2];
sx q[2];
rz(-1.0726856) q[2];
sx q[2];
rz(1.6061456) q[2];
rz(-2.8793907) q[3];
sx q[3];
rz(-1.6180792) q[3];
sx q[3];
rz(-2.1867627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75761211) q[0];
sx q[0];
rz(-2.7181427) q[0];
sx q[0];
rz(1.0773995) q[0];
rz(2.7491838) q[1];
sx q[1];
rz(-3.0632186) q[1];
sx q[1];
rz(-1.0307301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01973132) q[0];
sx q[0];
rz(-1.1543589) q[0];
sx q[0];
rz(2.5520578) q[0];
rz(-pi) q[1];
rz(-1.0723128) q[2];
sx q[2];
rz(-1.4118115) q[2];
sx q[2];
rz(-2.0778401) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2490152) q[1];
sx q[1];
rz(-1.7356145) q[1];
sx q[1];
rz(1.3457764) q[1];
x q[2];
rz(2.7983973) q[3];
sx q[3];
rz(-2.0569909) q[3];
sx q[3];
rz(-2.9064158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3173759) q[2];
sx q[2];
rz(-0.65418303) q[2];
sx q[2];
rz(0.83924323) q[2];
rz(-0.9134891) q[3];
sx q[3];
rz(-1.313611) q[3];
sx q[3];
rz(0.14348468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6827253) q[0];
sx q[0];
rz(-0.23696466) q[0];
sx q[0];
rz(1.6656026) q[0];
rz(2.7492375) q[1];
sx q[1];
rz(-1.0959492) q[1];
sx q[1];
rz(2.5501693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77350835) q[0];
sx q[0];
rz(-2.4249833) q[0];
sx q[0];
rz(-2.1518097) q[0];
x q[1];
rz(1.3900969) q[2];
sx q[2];
rz(-1.2998253) q[2];
sx q[2];
rz(-0.17943741) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6900031) q[1];
sx q[1];
rz(-1.065548) q[1];
sx q[1];
rz(3.0662159) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2675681) q[3];
sx q[3];
rz(-1.8667779) q[3];
sx q[3];
rz(3.0324453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.300294) q[2];
sx q[2];
rz(-2.4755307) q[2];
sx q[2];
rz(0.57233125) q[2];
rz(-2.9193997) q[3];
sx q[3];
rz(-2.7100345) q[3];
sx q[3];
rz(0.64479327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38729024) q[0];
sx q[0];
rz(-3.001725) q[0];
sx q[0];
rz(0.40827665) q[0];
rz(2.4093742) q[1];
sx q[1];
rz(-3.0156942) q[1];
sx q[1];
rz(-0.29762038) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0669516) q[0];
sx q[0];
rz(-0.56109259) q[0];
sx q[0];
rz(-0.69343167) q[0];
rz(-pi) q[1];
rz(-3.1213644) q[2];
sx q[2];
rz(-2.0654404) q[2];
sx q[2];
rz(-1.9488364) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1302766) q[1];
sx q[1];
rz(-1.374048) q[1];
sx q[1];
rz(-1.9185683) q[1];
x q[2];
rz(1.7209531) q[3];
sx q[3];
rz(-0.87267002) q[3];
sx q[3];
rz(-2.0719178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.97062651) q[2];
sx q[2];
rz(-1.8793224) q[2];
sx q[2];
rz(-0.69620281) q[2];
rz(-0.89037406) q[3];
sx q[3];
rz(-1.1768769) q[3];
sx q[3];
rz(-1.6740929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9625229) q[0];
sx q[0];
rz(-3.1130377) q[0];
sx q[0];
rz(-0.21275511) q[0];
rz(-0.46956024) q[1];
sx q[1];
rz(-0.9649562) q[1];
sx q[1];
rz(-0.75417095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51242669) q[0];
sx q[0];
rz(-0.34722024) q[0];
sx q[0];
rz(-2.1512593) q[0];
rz(0.55108549) q[2];
sx q[2];
rz(-0.37691016) q[2];
sx q[2];
rz(-0.92632471) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2060125) q[1];
sx q[1];
rz(-1.4569062) q[1];
sx q[1];
rz(-2.4920032) q[1];
rz(2.2175104) q[3];
sx q[3];
rz(-1.0093186) q[3];
sx q[3];
rz(-1.1738861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0065877) q[2];
sx q[2];
rz(-2.2392515) q[2];
sx q[2];
rz(0.78224409) q[2];
rz(1.4462645) q[3];
sx q[3];
rz(-2.5952314) q[3];
sx q[3];
rz(2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049147216) q[0];
sx q[0];
rz(-2.6706084) q[0];
sx q[0];
rz(-0.96442047) q[0];
rz(-1.2760705) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(-1.5020348) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3827189) q[0];
sx q[0];
rz(-1.8028462) q[0];
sx q[0];
rz(-2.9358572) q[0];
rz(2.658074) q[2];
sx q[2];
rz(-2.3502825) q[2];
sx q[2];
rz(-0.4936337) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6148541) q[1];
sx q[1];
rz(-0.51797082) q[1];
sx q[1];
rz(1.5244085) q[1];
rz(-1.4235017) q[3];
sx q[3];
rz(-1.5039267) q[3];
sx q[3];
rz(-1.221958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8883349) q[2];
sx q[2];
rz(-1.2829245) q[2];
sx q[2];
rz(-2.1962732) q[2];
rz(-0.82593289) q[3];
sx q[3];
rz(-1.5086987) q[3];
sx q[3];
rz(-0.53502423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0655521) q[0];
sx q[0];
rz(-1.1689508) q[0];
sx q[0];
rz(2.3384576) q[0];
rz(1.5638634) q[1];
sx q[1];
rz(-1.4811265) q[1];
sx q[1];
rz(0.28958431) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824376) q[0];
sx q[0];
rz(-1.4652325) q[0];
sx q[0];
rz(3.1292874) q[0];
x q[1];
rz(-1.2425735) q[2];
sx q[2];
rz(-1.8716836) q[2];
sx q[2];
rz(-2.2102578) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9152569) q[1];
sx q[1];
rz(-0.42789547) q[1];
sx q[1];
rz(1.0029327) q[1];
x q[2];
rz(-2.2749994) q[3];
sx q[3];
rz(-2.7150053) q[3];
sx q[3];
rz(1.0556737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3760066) q[2];
sx q[2];
rz(-0.12066081) q[2];
sx q[2];
rz(0.96013367) q[2];
rz(-0.58297408) q[3];
sx q[3];
rz(-0.65810242) q[3];
sx q[3];
rz(2.2768903) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4509907) q[0];
sx q[0];
rz(-1.7091746) q[0];
sx q[0];
rz(-1.5102392) q[0];
rz(-3.1008537) q[1];
sx q[1];
rz(-0.67650411) q[1];
sx q[1];
rz(0.13112851) q[1];
rz(2.1292674) q[2];
sx q[2];
rz(-2.6016015) q[2];
sx q[2];
rz(0.38228959) q[2];
rz(-0.05091713) q[3];
sx q[3];
rz(-1.49972) q[3];
sx q[3];
rz(-1.6928991) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
