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
rz(-1.3136366) q[0];
sx q[0];
rz(-3.0744636) q[0];
sx q[0];
rz(3.1188174) q[0];
rz(1.0442806) q[1];
sx q[1];
rz(-2.2106946) q[1];
sx q[1];
rz(0.013962362) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305796) q[0];
sx q[0];
rz(-1.530181) q[0];
sx q[0];
rz(-0.75675772) q[0];
rz(-pi) q[1];
rz(1.1544873) q[2];
sx q[2];
rz(-1.1616366) q[2];
sx q[2];
rz(0.43658999) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1779132) q[1];
sx q[1];
rz(-1.3909893) q[1];
sx q[1];
rz(-2.9397029) q[1];
rz(-pi) q[2];
rz(-2.4926643) q[3];
sx q[3];
rz(-2.5960659) q[3];
sx q[3];
rz(2.3758171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.969101) q[2];
sx q[2];
rz(-1.8164941) q[2];
sx q[2];
rz(1.3410404) q[2];
rz(1.4260882) q[3];
sx q[3];
rz(-1.1172349) q[3];
sx q[3];
rz(-2.8360227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7476927) q[0];
sx q[0];
rz(-1.1423528) q[0];
sx q[0];
rz(-0.72545141) q[0];
rz(-2.2254288) q[1];
sx q[1];
rz(-2.3326645) q[1];
sx q[1];
rz(2.5221672) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9206132) q[0];
sx q[0];
rz(-2.0720464) q[0];
sx q[0];
rz(0.46904012) q[0];
rz(-pi) q[1];
rz(-2.393707) q[2];
sx q[2];
rz(-1.1074142) q[2];
sx q[2];
rz(-0.057967535) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6065154) q[1];
sx q[1];
rz(-1.1454574) q[1];
sx q[1];
rz(2.7414447) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1257915) q[3];
sx q[3];
rz(-1.3135305) q[3];
sx q[3];
rz(-0.98350888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3118423) q[2];
sx q[2];
rz(-2.1205015) q[2];
sx q[2];
rz(2.1799083) q[2];
rz(-2.035615) q[3];
sx q[3];
rz(-2.1092236) q[3];
sx q[3];
rz(-0.34524125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.2894534) q[0];
sx q[0];
rz(-2.0134605) q[0];
sx q[0];
rz(-1.3683251) q[0];
rz(3.1276357) q[1];
sx q[1];
rz(-1.114926) q[1];
sx q[1];
rz(-1.7272635) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089979261) q[0];
sx q[0];
rz(-0.3246626) q[0];
sx q[0];
rz(1.8568296) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0711381) q[2];
sx q[2];
rz(-1.2025598) q[2];
sx q[2];
rz(0.28063831) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3625379) q[1];
sx q[1];
rz(-1.0666872) q[1];
sx q[1];
rz(-2.6466531) q[1];
rz(-pi) q[2];
rz(-1.3120804) q[3];
sx q[3];
rz(-1.1911827) q[3];
sx q[3];
rz(0.9791475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3783375) q[2];
sx q[2];
rz(-1.9941284) q[2];
sx q[2];
rz(-2.9249127) q[2];
rz(0.8832776) q[3];
sx q[3];
rz(-2.7933385) q[3];
sx q[3];
rz(-1.1368375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.20275177) q[0];
sx q[0];
rz(-2.1348248) q[0];
sx q[0];
rz(-2.3034565) q[0];
rz(2.5996767) q[1];
sx q[1];
rz(-2.7653265) q[1];
sx q[1];
rz(-0.045104973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7823001) q[0];
sx q[0];
rz(-1.4385537) q[0];
sx q[0];
rz(-2.537893) q[0];
x q[1];
rz(0.25928478) q[2];
sx q[2];
rz(-1.102299) q[2];
sx q[2];
rz(-1.1069654) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.20238961) q[1];
sx q[1];
rz(-1.09859) q[1];
sx q[1];
rz(-3.0665728) q[1];
x q[2];
rz(-0.18121775) q[3];
sx q[3];
rz(-2.6406796) q[3];
sx q[3];
rz(1.8883226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0185467) q[2];
sx q[2];
rz(-0.16877731) q[2];
sx q[2];
rz(-0.65400845) q[2];
rz(0.6319913) q[3];
sx q[3];
rz(-1.4258823) q[3];
sx q[3];
rz(-0.33647195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80394799) q[0];
sx q[0];
rz(-2.4027282) q[0];
sx q[0];
rz(1.7919354) q[0];
rz(2.3140287) q[1];
sx q[1];
rz(-2.5065828) q[1];
sx q[1];
rz(0.48447022) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2506566) q[0];
sx q[0];
rz(-1.93205) q[0];
sx q[0];
rz(2.7879232) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1250317) q[2];
sx q[2];
rz(-1.4942188) q[2];
sx q[2];
rz(-0.34293338) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9974629) q[1];
sx q[1];
rz(-1.1511377) q[1];
sx q[1];
rz(-1.4317572) q[1];
x q[2];
rz(-2.6788533) q[3];
sx q[3];
rz(-0.92904186) q[3];
sx q[3];
rz(-0.79824191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95939246) q[2];
sx q[2];
rz(-1.0687989) q[2];
sx q[2];
rz(-2.5051795) q[2];
rz(-0.098879769) q[3];
sx q[3];
rz(-1.5560047) q[3];
sx q[3];
rz(-1.4720565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8168378) q[0];
sx q[0];
rz(-0.88961283) q[0];
sx q[0];
rz(-1.9689993) q[0];
rz(-1.6088387) q[1];
sx q[1];
rz(-1.0120665) q[1];
sx q[1];
rz(3.1408659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6470203) q[0];
sx q[0];
rz(-1.9969517) q[0];
sx q[0];
rz(-2.1151353) q[0];
rz(-pi) q[1];
rz(-0.22063271) q[2];
sx q[2];
rz(-1.4743525) q[2];
sx q[2];
rz(0.19410832) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.128617) q[1];
sx q[1];
rz(-0.3518577) q[1];
sx q[1];
rz(-0.33280876) q[1];
x q[2];
rz(-0.29720593) q[3];
sx q[3];
rz(-2.5416221) q[3];
sx q[3];
rz(-1.6855437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.26649228) q[2];
sx q[2];
rz(-2.2537587) q[2];
sx q[2];
rz(-0.20981851) q[2];
rz(2.3328414) q[3];
sx q[3];
rz(-2.5918312) q[3];
sx q[3];
rz(-2.2426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784191) q[0];
sx q[0];
rz(-1.2043948) q[0];
sx q[0];
rz(-3.0294982) q[0];
rz(-0.049292715) q[1];
sx q[1];
rz(-2.6105328) q[1];
sx q[1];
rz(-1.8045527) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5264915) q[0];
sx q[0];
rz(-1.3630629) q[0];
sx q[0];
rz(-2.2206942) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3434197) q[2];
sx q[2];
rz(-2.6223287) q[2];
sx q[2];
rz(2.8092172) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35838764) q[1];
sx q[1];
rz(-1.9175944) q[1];
sx q[1];
rz(1.5601182) q[1];
rz(1.2483761) q[3];
sx q[3];
rz(-0.60203505) q[3];
sx q[3];
rz(0.69794929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8423395) q[2];
sx q[2];
rz(-0.82483333) q[2];
sx q[2];
rz(2.4315368) q[2];
rz(0.0013141343) q[3];
sx q[3];
rz(-0.90130663) q[3];
sx q[3];
rz(-2.5406751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8173219) q[0];
sx q[0];
rz(-1.7621499) q[0];
sx q[0];
rz(-2.5589909) q[0];
rz(-2.461589) q[1];
sx q[1];
rz(-2.1425207) q[1];
sx q[1];
rz(-1.8780139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1909785) q[0];
sx q[0];
rz(-2.593747) q[0];
sx q[0];
rz(-0.56444278) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.787705) q[2];
sx q[2];
rz(-1.3706651) q[2];
sx q[2];
rz(-1.4544124) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3547825) q[1];
sx q[1];
rz(-0.92325961) q[1];
sx q[1];
rz(1.7969653) q[1];
x q[2];
rz(-0.42383343) q[3];
sx q[3];
rz(-1.5409711) q[3];
sx q[3];
rz(0.78889293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6363643) q[2];
sx q[2];
rz(-1.4437081) q[2];
sx q[2];
rz(3.0155724) q[2];
rz(-1.5382918) q[3];
sx q[3];
rz(-2.3647629) q[3];
sx q[3];
rz(-2.7974424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6525604) q[0];
sx q[0];
rz(-1.5131938) q[0];
sx q[0];
rz(-1.162758) q[0];
rz(1.8310422) q[1];
sx q[1];
rz(-2.4030011) q[1];
sx q[1];
rz(1.7576677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3726295) q[0];
sx q[0];
rz(-1.2314737) q[0];
sx q[0];
rz(0.047554544) q[0];
rz(-pi) q[1];
rz(-3.018749) q[2];
sx q[2];
rz(-1.4990988) q[2];
sx q[2];
rz(-0.724585) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0392929) q[1];
sx q[1];
rz(-2.1792534) q[1];
sx q[1];
rz(0.14157544) q[1];
rz(-pi) q[2];
x q[2];
rz(1.27326) q[3];
sx q[3];
rz(-0.14500824) q[3];
sx q[3];
rz(0.48806897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80013529) q[2];
sx q[2];
rz(-1.8693482) q[2];
sx q[2];
rz(-2.6696491) q[2];
rz(-2.3927472) q[3];
sx q[3];
rz(-2.2180836) q[3];
sx q[3];
rz(0.81764618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5188468) q[0];
sx q[0];
rz(-0.8569583) q[0];
sx q[0];
rz(2.6527606) q[0];
rz(-2.7986616) q[1];
sx q[1];
rz(-0.88645005) q[1];
sx q[1];
rz(2.0874646) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80741548) q[0];
sx q[0];
rz(-1.6082967) q[0];
sx q[0];
rz(-1.6223909) q[0];
x q[1];
rz(-1.2230049) q[2];
sx q[2];
rz(-0.985983) q[2];
sx q[2];
rz(1.0725759) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17133896) q[1];
sx q[1];
rz(-1.043817) q[1];
sx q[1];
rz(-2.8026366) q[1];
rz(-3.0929186) q[3];
sx q[3];
rz(-1.0992388) q[3];
sx q[3];
rz(-2.8296628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1663345) q[2];
sx q[2];
rz(-2.9322093) q[2];
sx q[2];
rz(0.091910467) q[2];
rz(2.6722243) q[3];
sx q[3];
rz(-2.2769603) q[3];
sx q[3];
rz(-0.11274591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5975006) q[0];
sx q[0];
rz(-1.8202029) q[0];
sx q[0];
rz(1.2910917) q[0];
rz(0.76580936) q[1];
sx q[1];
rz(-1.3594834) q[1];
sx q[1];
rz(1.8654738) q[1];
rz(-1.932939) q[2];
sx q[2];
rz(-1.3315143) q[2];
sx q[2];
rz(-1.3523921) q[2];
rz(1.3113931) q[3];
sx q[3];
rz(-1.048199) q[3];
sx q[3];
rz(-0.38182624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
