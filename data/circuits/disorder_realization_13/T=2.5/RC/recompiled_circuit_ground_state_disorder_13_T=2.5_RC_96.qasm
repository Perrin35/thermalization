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
rz(-0.19620148) q[1];
sx q[1];
rz(-1.3407522) q[1];
sx q[1];
rz(3.0412716) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8226877) q[0];
sx q[0];
rz(-1.5189369) q[0];
sx q[0];
rz(2.479631) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60007976) q[2];
sx q[2];
rz(-0.55309764) q[2];
sx q[2];
rz(1.2002522) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70182204) q[1];
sx q[1];
rz(-1.3622829) q[1];
sx q[1];
rz(-2.8609402) q[1];
x q[2];
rz(-2.5514929) q[3];
sx q[3];
rz(-1.8339155) q[3];
sx q[3];
rz(0.89917574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4002865) q[2];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62419409) q[0];
sx q[0];
rz(-0.26573467) q[0];
sx q[0];
rz(2.8345795) q[0];
rz(1.864805) q[1];
sx q[1];
rz(-1.1390319) q[1];
sx q[1];
rz(-1.101864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3392838) q[0];
sx q[0];
rz(-1.3972939) q[0];
sx q[0];
rz(2.1306778) q[0];
rz(1.5069783) q[2];
sx q[2];
rz(-1.391727) q[2];
sx q[2];
rz(0.049953559) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.15436048) q[1];
sx q[1];
rz(-1.4353961) q[1];
sx q[1];
rz(1.7495278) q[1];
x q[2];
rz(0.78897055) q[3];
sx q[3];
rz(-1.0010825) q[3];
sx q[3];
rz(2.8017442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1133984) q[2];
sx q[2];
rz(-0.58176175) q[2];
sx q[2];
rz(0.096435189) q[2];
rz(-2.9891369) q[3];
sx q[3];
rz(-1.5072482) q[3];
sx q[3];
rz(-2.4436387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089712791) q[0];
sx q[0];
rz(-2.0413601) q[0];
sx q[0];
rz(-2.7413947) q[0];
rz(-1.5125037) q[1];
sx q[1];
rz(-0.17833231) q[1];
sx q[1];
rz(-2.8834744) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1296279) q[0];
sx q[0];
rz(-1.3370378) q[0];
sx q[0];
rz(2.9988704) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4908954) q[2];
sx q[2];
rz(-1.4712442) q[2];
sx q[2];
rz(2.4308506) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3052747) q[1];
sx q[1];
rz(-1.5757676) q[1];
sx q[1];
rz(-3.1286376) q[1];
rz(-0.6155562) q[3];
sx q[3];
rz(-0.29272348) q[3];
sx q[3];
rz(0.20197257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.021598024) q[2];
sx q[2];
rz(-1.9436516) q[2];
sx q[2];
rz(-3.1094587) q[2];
rz(2.8739127) q[3];
sx q[3];
rz(-1.6240424) q[3];
sx q[3];
rz(1.5599498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(1.5716008) q[0];
sx q[0];
rz(-1.5084234) q[0];
sx q[0];
rz(-3.0573523) q[0];
rz(-3.1056504) q[1];
sx q[1];
rz(-3.1075931) q[1];
sx q[1];
rz(-0.34119225) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81927903) q[0];
sx q[0];
rz(-1.8321165) q[0];
sx q[0];
rz(-0.65066353) q[0];
rz(-pi) q[1];
rz(1.313339) q[2];
sx q[2];
rz(-0.99960589) q[2];
sx q[2];
rz(-1.5866304) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1028776) q[1];
sx q[1];
rz(-0.85188085) q[1];
sx q[1];
rz(2.170156) q[1];
x q[2];
rz(-2.5248923) q[3];
sx q[3];
rz(-2.4444207) q[3];
sx q[3];
rz(-2.8215849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7146032) q[2];
sx q[2];
rz(-1.0726856) q[2];
sx q[2];
rz(1.6061456) q[2];
rz(0.26220194) q[3];
sx q[3];
rz(-1.5235135) q[3];
sx q[3];
rz(2.1867627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3839805) q[0];
sx q[0];
rz(-0.42344991) q[0];
sx q[0];
rz(-1.0773995) q[0];
rz(0.39240882) q[1];
sx q[1];
rz(-3.0632186) q[1];
sx q[1];
rz(-2.1108625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0466246) q[0];
sx q[0];
rz(-2.4343581) q[0];
sx q[0];
rz(0.67202248) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0723128) q[2];
sx q[2];
rz(-1.7297812) q[2];
sx q[2];
rz(1.0637525) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3001316) q[1];
sx q[1];
rz(-0.2780973) q[1];
sx q[1];
rz(-0.93021955) q[1];
rz(2.7983973) q[3];
sx q[3];
rz(-2.0569909) q[3];
sx q[3];
rz(0.23517683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82421676) q[2];
sx q[2];
rz(-0.65418303) q[2];
sx q[2];
rz(2.3023494) q[2];
rz(2.2281036) q[3];
sx q[3];
rz(-1.313611) q[3];
sx q[3];
rz(0.14348468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4588673) q[0];
sx q[0];
rz(-0.23696466) q[0];
sx q[0];
rz(1.4759901) q[0];
rz(0.39235517) q[1];
sx q[1];
rz(-1.0959492) q[1];
sx q[1];
rz(0.59142339) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8845733) q[0];
sx q[0];
rz(-1.2019751) q[0];
sx q[0];
rz(-0.9414282) q[0];
rz(-pi) q[1];
rz(-2.5673836) q[2];
sx q[2];
rz(-2.817135) q[2];
sx q[2];
rz(-0.41947075) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5352262) q[1];
sx q[1];
rz(-2.6312345) q[1];
sx q[1];
rz(-1.4354857) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37844946) q[3];
sx q[3];
rz(-0.90988084) q[3];
sx q[3];
rz(-1.4405516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8412987) q[2];
sx q[2];
rz(-2.4755307) q[2];
sx q[2];
rz(2.5692614) q[2];
rz(-2.9193997) q[3];
sx q[3];
rz(-0.43155813) q[3];
sx q[3];
rz(-0.64479327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7543024) q[0];
sx q[0];
rz(-3.001725) q[0];
sx q[0];
rz(-2.733316) q[0];
rz(2.4093742) q[1];
sx q[1];
rz(-3.0156942) q[1];
sx q[1];
rz(2.8439723) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0746411) q[0];
sx q[0];
rz(-0.56109259) q[0];
sx q[0];
rz(2.448161) q[0];
x q[1];
rz(1.5333129) q[2];
sx q[2];
rz(-2.6465694) q[2];
sx q[2];
rz(1.2353473) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.054033259) q[1];
sx q[1];
rz(-2.7440024) q[1];
sx q[1];
rz(-2.1000186) q[1];
x q[2];
rz(2.4378889) q[3];
sx q[3];
rz(-1.6856442) q[3];
sx q[3];
rz(0.59806693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97062651) q[2];
sx q[2];
rz(-1.8793224) q[2];
sx q[2];
rz(-2.4453898) q[2];
rz(0.89037406) q[3];
sx q[3];
rz(-1.9647157) q[3];
sx q[3];
rz(1.4674998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9625229) q[0];
sx q[0];
rz(-0.028554976) q[0];
sx q[0];
rz(-2.9288375) q[0];
rz(0.46956024) q[1];
sx q[1];
rz(-2.1766365) q[1];
sx q[1];
rz(-0.75417095) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6109723) q[0];
sx q[0];
rz(-1.3830796) q[0];
sx q[0];
rz(1.8646445) q[0];
rz(-pi) q[1];
rz(-2.8163359) q[2];
sx q[2];
rz(-1.3768679) q[2];
sx q[2];
rz(-3.0162899) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9355802) q[1];
sx q[1];
rz(-1.4569062) q[1];
sx q[1];
rz(0.64958944) q[1];
rz(-2.3776618) q[3];
sx q[3];
rz(-0.82909938) q[3];
sx q[3];
rz(-1.0111077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0065877) q[2];
sx q[2];
rz(-2.2392515) q[2];
sx q[2];
rz(-0.78224409) q[2];
rz(1.6953281) q[3];
sx q[3];
rz(-2.5952314) q[3];
sx q[3];
rz(-2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924454) q[0];
sx q[0];
rz(-0.47098422) q[0];
sx q[0];
rz(-2.1771722) q[0];
rz(1.2760705) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(-1.6395578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0014711) q[0];
sx q[0];
rz(-1.3706511) q[0];
sx q[0];
rz(1.8076623) q[0];
x q[1];
rz(-2.0105069) q[2];
sx q[2];
rz(-2.2520116) q[2];
sx q[2];
rz(1.1352254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5267386) q[1];
sx q[1];
rz(-2.6236218) q[1];
sx q[1];
rz(1.5244085) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1427059) q[3];
sx q[3];
rz(-2.9799298) q[3];
sx q[3];
rz(2.3695994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2532578) q[2];
sx q[2];
rz(-1.2829245) q[2];
sx q[2];
rz(-0.9453195) q[2];
rz(0.82593289) q[3];
sx q[3];
rz(-1.5086987) q[3];
sx q[3];
rz(-2.6065684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076040529) q[0];
sx q[0];
rz(-1.1689508) q[0];
sx q[0];
rz(0.8031351) q[0];
rz(1.5638634) q[1];
sx q[1];
rz(-1.6604661) q[1];
sx q[1];
rz(-0.28958431) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824376) q[0];
sx q[0];
rz(-1.4652325) q[0];
sx q[0];
rz(-0.012305208) q[0];
rz(-pi) q[1];
rz(-1.2425735) q[2];
sx q[2];
rz(-1.8716836) q[2];
sx q[2];
rz(0.93133486) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9601396) q[1];
sx q[1];
rz(-1.3457237) q[1];
sx q[1];
rz(-1.9378661) q[1];
rz(-1.9042468) q[3];
sx q[3];
rz(-1.2996009) q[3];
sx q[3];
rz(-1.9681794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76558602) q[2];
sx q[2];
rz(-3.0209318) q[2];
sx q[2];
rz(-0.96013367) q[2];
rz(-2.5586186) q[3];
sx q[3];
rz(-2.4834902) q[3];
sx q[3];
rz(2.2768903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4509907) q[0];
sx q[0];
rz(-1.4324181) q[0];
sx q[0];
rz(1.6313534) q[0];
rz(-3.1008537) q[1];
sx q[1];
rz(-0.67650411) q[1];
sx q[1];
rz(0.13112851) q[1];
rz(-0.30754752) q[2];
sx q[2];
rz(-2.0219621) q[2];
sx q[2];
rz(-2.1297217) q[2];
rz(-1.4996281) q[3];
sx q[3];
rz(-1.5200079) q[3];
sx q[3];
rz(-0.11848371) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
