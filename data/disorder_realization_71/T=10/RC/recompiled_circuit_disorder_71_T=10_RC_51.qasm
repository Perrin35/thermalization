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
rz(-2.3117476) q[0];
rz(-2.3614376) q[1];
sx q[1];
rz(-1.0649788) q[1];
sx q[1];
rz(0.87632626) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5397545) q[0];
sx q[0];
rz(-1.1798501) q[0];
sx q[0];
rz(-1.8200309) q[0];
x q[1];
rz(0.19199065) q[2];
sx q[2];
rz(-2.1477094) q[2];
sx q[2];
rz(0.38919762) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4235916) q[1];
sx q[1];
rz(-2.8112486) q[1];
sx q[1];
rz(1.6172536) q[1];
rz(-pi) q[2];
rz(-0.40684367) q[3];
sx q[3];
rz(-1.7115271) q[3];
sx q[3];
rz(-1.391747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9709388) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(0.7286287) q[2];
rz(-2.6206) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347297) q[0];
sx q[0];
rz(-1.9711718) q[0];
sx q[0];
rz(1.1215425) q[0];
rz(2.8858378) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(-0.87444011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16337285) q[0];
sx q[0];
rz(-1.1436497) q[0];
sx q[0];
rz(2.8073729) q[0];
x q[1];
rz(-0.15408709) q[2];
sx q[2];
rz(-1.0006957) q[2];
sx q[2];
rz(-0.64683435) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.892131) q[1];
sx q[1];
rz(-2.3453237) q[1];
sx q[1];
rz(-0.65278058) q[1];
rz(-0.053697649) q[3];
sx q[3];
rz(-1.8397545) q[3];
sx q[3];
rz(0.69124903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.740739) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(0.43593105) q[2];
rz(-2.46051) q[3];
sx q[3];
rz(-0.77107945) q[3];
sx q[3];
rz(2.7387103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9044559) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(2.0667734) q[0];
rz(0.83956051) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(-2.7456465) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0229189) q[0];
sx q[0];
rz(-0.6624822) q[0];
sx q[0];
rz(0.89612095) q[0];
x q[1];
rz(-1.6171574) q[2];
sx q[2];
rz(-0.29435396) q[2];
sx q[2];
rz(-2.4183395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.695895) q[1];
sx q[1];
rz(-2.2539415) q[1];
sx q[1];
rz(1.5831069) q[1];
rz(0.25554244) q[3];
sx q[3];
rz(-1.5481755) q[3];
sx q[3];
rz(2.8850151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6039156) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(-0.23920693) q[2];
rz(3.0662597) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(-1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4199715) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(0.81992942) q[0];
rz(2.6539102) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(-2.908169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4963835) q[0];
sx q[0];
rz(-1.5102981) q[0];
sx q[0];
rz(-1.66771) q[0];
rz(-1.3350305) q[2];
sx q[2];
rz(-2.1996017) q[2];
sx q[2];
rz(-0.51970383) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1366795) q[1];
sx q[1];
rz(-2.3466952) q[1];
sx q[1];
rz(-2.6585048) q[1];
rz(0.38846429) q[3];
sx q[3];
rz(-2.3187175) q[3];
sx q[3];
rz(0.32271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0908115) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(2.3941669) q[2];
rz(-2.9181972) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(0.24211611) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500047) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(-0.064237021) q[0];
rz(-0.94379395) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(2.3805526) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8262186) q[0];
sx q[0];
rz(-1.5407019) q[0];
sx q[0];
rz(-1.3099567) q[0];
rz(0.76743482) q[2];
sx q[2];
rz(-0.27786294) q[2];
sx q[2];
rz(-2.1674736) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4796175) q[1];
sx q[1];
rz(-1.2503337) q[1];
sx q[1];
rz(1.596405) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8499591) q[3];
sx q[3];
rz(-1.8431292) q[3];
sx q[3];
rz(0.84701049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.80660194) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(-3.0495194) q[2];
rz(2.4798685) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(-1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6732366) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(0.026542149) q[0];
rz(2.2684855) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(2.81566) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.034176055) q[0];
sx q[0];
rz(-2.0909485) q[0];
sx q[0];
rz(-2.0053894) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8750149) q[2];
sx q[2];
rz(-1.1433257) q[2];
sx q[2];
rz(1.8311335) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95820108) q[1];
sx q[1];
rz(-2.8124053) q[1];
sx q[1];
rz(0.42894657) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10109191) q[3];
sx q[3];
rz(-1.5540082) q[3];
sx q[3];
rz(1.8198131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59763336) q[2];
sx q[2];
rz(-1.3175069) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(-1.4298965) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-0.54106075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8687246) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(2.4196999) q[0];
rz(-1.7294653) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(-3.022335) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8598547) q[0];
sx q[0];
rz(-2.2681232) q[0];
sx q[0];
rz(-1.2233234) q[0];
rz(0.97738816) q[2];
sx q[2];
rz(-0.19492976) q[2];
sx q[2];
rz(-2.5755663) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6357928) q[1];
sx q[1];
rz(-2.0675106) q[1];
sx q[1];
rz(0.15091166) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70763208) q[3];
sx q[3];
rz(-1.8579357) q[3];
sx q[3];
rz(1.6361145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44935903) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(1.3593486) q[2];
rz(-0.75891495) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(-0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83157241) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(-1.0634364) q[0];
rz(2.8670782) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(-0.88561052) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0622334) q[0];
sx q[0];
rz(-2.057312) q[0];
sx q[0];
rz(-1.3079206) q[0];
x q[1];
rz(-2.1503259) q[2];
sx q[2];
rz(-1.5577003) q[2];
sx q[2];
rz(-0.0038557204) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.13247989) q[1];
sx q[1];
rz(-1.9630034) q[1];
sx q[1];
rz(1.023804) q[1];
x q[2];
rz(2.4212491) q[3];
sx q[3];
rz(-2.2317436) q[3];
sx q[3];
rz(-1.9545481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4132335) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(-2.0098861) q[2];
rz(1.0845832) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(1.926698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.30329147) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(2.0595179) q[0];
rz(1.2754296) q[1];
sx q[1];
rz(-1.0042896) q[1];
sx q[1];
rz(2.0057604) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5539615) q[0];
sx q[0];
rz(-2.580664) q[0];
sx q[0];
rz(-1.8857303) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43478195) q[2];
sx q[2];
rz(-1.5506622) q[2];
sx q[2];
rz(-0.15205631) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.218704) q[1];
sx q[1];
rz(-1.522038) q[1];
sx q[1];
rz(0.4263652) q[1];
x q[2];
rz(-1.0335835) q[3];
sx q[3];
rz(-2.4164003) q[3];
sx q[3];
rz(2.4560526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0231126) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(-2.2686968) q[2];
rz(-2.2980799) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(0.18994722) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2492367) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(1.0247914) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(1.9445673) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.409317) q[0];
sx q[0];
rz(-1.1336375) q[0];
sx q[0];
rz(0.93426312) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1080152) q[2];
sx q[2];
rz(-2.8769851) q[2];
sx q[2];
rz(-2.3679784) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7750818) q[1];
sx q[1];
rz(-1.0771891) q[1];
sx q[1];
rz(-2.2713186) q[1];
x q[2];
rz(-2.6322369) q[3];
sx q[3];
rz(-1.2378113) q[3];
sx q[3];
rz(1.8978564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8355576) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(-0.79997921) q[2];
rz(-1.1768613) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(-2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4326614) q[0];
sx q[0];
rz(-2.9928757) q[0];
sx q[0];
rz(-2.3401674) q[0];
rz(2.6196383) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(2.3847053) q[2];
sx q[2];
rz(-0.49513985) q[2];
sx q[2];
rz(-2.0078299) q[2];
rz(2.5064777) q[3];
sx q[3];
rz(-1.0233581) q[3];
sx q[3];
rz(-0.81627853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];