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
rz(1.963653) q[0];
sx q[0];
rz(-0.093129245) q[0];
sx q[0];
rz(7.0153305) q[0];
rz(1.572345) q[1];
sx q[1];
rz(0.88580004) q[1];
sx q[1];
rz(9.8635397) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6308545) q[0];
sx q[0];
rz(-2.8048212) q[0];
sx q[0];
rz(-1.4022409) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2874424) q[2];
sx q[2];
rz(-1.3204976) q[2];
sx q[2];
rz(1.3586501) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4784412) q[1];
sx q[1];
rz(-1.8290797) q[1];
sx q[1];
rz(-1.6186991) q[1];
x q[2];
rz(3.0436727) q[3];
sx q[3];
rz(-0.46028462) q[3];
sx q[3];
rz(1.3336045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0438805) q[2];
sx q[2];
rz(-0.040412929) q[2];
sx q[2];
rz(-2.2229693) q[2];
rz(2.0888603) q[3];
sx q[3];
rz(-0.91679263) q[3];
sx q[3];
rz(-0.91744939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.265428) q[0];
sx q[0];
rz(-2.3227203) q[0];
sx q[0];
rz(0.87232605) q[0];
rz(1.0126975) q[1];
sx q[1];
rz(-1.6302949) q[1];
sx q[1];
rz(-0.33908078) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6144086) q[0];
sx q[0];
rz(-2.0054584) q[0];
sx q[0];
rz(1.4285157) q[0];
x q[1];
rz(-2.2834899) q[2];
sx q[2];
rz(-0.65026186) q[2];
sx q[2];
rz(1.5257143) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6244858) q[1];
sx q[1];
rz(-0.042306012) q[1];
sx q[1];
rz(-2.9632764) q[1];
rz(1.2463739) q[3];
sx q[3];
rz(-1.4528465) q[3];
sx q[3];
rz(-1.1583386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33112153) q[2];
sx q[2];
rz(-1.1473742) q[2];
sx q[2];
rz(-2.5362711) q[2];
rz(-2.807054) q[3];
sx q[3];
rz(-2.7693222) q[3];
sx q[3];
rz(-0.076586671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.2452068) q[0];
sx q[0];
rz(-2.4672282) q[0];
sx q[0];
rz(-1.4728004) q[0];
rz(-0.32707602) q[1];
sx q[1];
rz(-1.3027124) q[1];
sx q[1];
rz(-0.27007857) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.208858) q[0];
sx q[0];
rz(-1.7678102) q[0];
sx q[0];
rz(-0.52893649) q[0];
rz(-1.4065341) q[2];
sx q[2];
rz(-2.7077419) q[2];
sx q[2];
rz(2.1151154) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9880319) q[1];
sx q[1];
rz(-0.84756339) q[1];
sx q[1];
rz(-0.60232298) q[1];
rz(-pi) q[2];
rz(-0.080167218) q[3];
sx q[3];
rz(-1.5746377) q[3];
sx q[3];
rz(0.14417917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0373056) q[2];
sx q[2];
rz(-1.1423926) q[2];
sx q[2];
rz(-1.8541065) q[2];
rz(1.9262975) q[3];
sx q[3];
rz(-1.3853962) q[3];
sx q[3];
rz(-2.9638885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52008587) q[0];
sx q[0];
rz(-2.6984213) q[0];
sx q[0];
rz(0.31975123) q[0];
rz(-2.7301835) q[1];
sx q[1];
rz(-1.5450059) q[1];
sx q[1];
rz(0.080332669) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4302931) q[0];
sx q[0];
rz(-1.7691433) q[0];
sx q[0];
rz(1.8848778) q[0];
rz(-pi) q[1];
rz(2.3351203) q[2];
sx q[2];
rz(-1.6739168) q[2];
sx q[2];
rz(-1.1191302) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.44583933) q[1];
sx q[1];
rz(-1.2994804) q[1];
sx q[1];
rz(-2.325227) q[1];
rz(0.61088722) q[3];
sx q[3];
rz(-1.7912994) q[3];
sx q[3];
rz(2.084651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2692928) q[2];
sx q[2];
rz(-2.9041957) q[2];
sx q[2];
rz(-3.0966975) q[2];
rz(-2.6063555) q[3];
sx q[3];
rz(-1.6343296) q[3];
sx q[3];
rz(-0.11307344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.991268) q[0];
sx q[0];
rz(-0.70882216) q[0];
sx q[0];
rz(2.2659361) q[0];
rz(0.21370299) q[1];
sx q[1];
rz(-0.49476606) q[1];
sx q[1];
rz(-1.5692086) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0914577) q[0];
sx q[0];
rz(-2.7612855) q[0];
sx q[0];
rz(0.97635348) q[0];
rz(2.7046811) q[2];
sx q[2];
rz(-1.3624117) q[2];
sx q[2];
rz(-2.2885099) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.24092785) q[1];
sx q[1];
rz(-0.87605175) q[1];
sx q[1];
rz(1.7068295) q[1];
rz(-pi) q[2];
rz(1.6650136) q[3];
sx q[3];
rz(-1.251128) q[3];
sx q[3];
rz(1.6188542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.956942) q[2];
sx q[2];
rz(-0.1447548) q[2];
sx q[2];
rz(2.0733898) q[2];
rz(-0.60643658) q[3];
sx q[3];
rz(-1.2111827) q[3];
sx q[3];
rz(-3.1312805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.40474263) q[0];
sx q[0];
rz(-0.24979845) q[0];
sx q[0];
rz(1.7506208) q[0];
rz(2.6178316) q[1];
sx q[1];
rz(-2.5660089) q[1];
sx q[1];
rz(1.1837122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8075631) q[0];
sx q[0];
rz(-1.5153756) q[0];
sx q[0];
rz(1.8593273) q[0];
x q[1];
rz(-3.0880586) q[2];
sx q[2];
rz(-1.6463037) q[2];
sx q[2];
rz(-1.8078505) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27734807) q[1];
sx q[1];
rz(-1.7226761) q[1];
sx q[1];
rz(-0.95452301) q[1];
rz(-1.7051172) q[3];
sx q[3];
rz(-1.4616579) q[3];
sx q[3];
rz(1.2827831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4411053) q[2];
sx q[2];
rz(-2.5429071) q[2];
sx q[2];
rz(2.109745) q[2];
rz(1.6796238) q[3];
sx q[3];
rz(-0.69457355) q[3];
sx q[3];
rz(-1.3612548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66330355) q[0];
sx q[0];
rz(-0.42600584) q[0];
sx q[0];
rz(1.4248832) q[0];
rz(0.58319631) q[1];
sx q[1];
rz(-1.2251264) q[1];
sx q[1];
rz(0.057417631) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9641477) q[0];
sx q[0];
rz(-0.76745196) q[0];
sx q[0];
rz(-1.5758118) q[0];
rz(-2.4671747) q[2];
sx q[2];
rz(-0.73235529) q[2];
sx q[2];
rz(-0.44841097) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.26099953) q[1];
sx q[1];
rz(-2.2003502) q[1];
sx q[1];
rz(1.7884607) q[1];
x q[2];
rz(-2.613861) q[3];
sx q[3];
rz(-0.70042983) q[3];
sx q[3];
rz(2.79984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3164369) q[2];
sx q[2];
rz(-1.8437443) q[2];
sx q[2];
rz(-1.8334897) q[2];
rz(1.3206652) q[3];
sx q[3];
rz(-0.84669176) q[3];
sx q[3];
rz(2.2804885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2945781) q[0];
sx q[0];
rz(-2.4569643) q[0];
sx q[0];
rz(-1.2680049) q[0];
rz(1.1216724) q[1];
sx q[1];
rz(-1.8662165) q[1];
sx q[1];
rz(0.65006382) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92597678) q[0];
sx q[0];
rz(-0.6840261) q[0];
sx q[0];
rz(2.2355272) q[0];
rz(-pi) q[1];
rz(-0.68393884) q[2];
sx q[2];
rz(-1.9994447) q[2];
sx q[2];
rz(-0.49698439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6396811) q[1];
sx q[1];
rz(-0.6957275) q[1];
sx q[1];
rz(0.63668107) q[1];
x q[2];
rz(2.6542957) q[3];
sx q[3];
rz(-1.0505884) q[3];
sx q[3];
rz(2.3765124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0012297) q[2];
sx q[2];
rz(-1.4823806) q[2];
sx q[2];
rz(-2.7492827) q[2];
rz(-2.9366734) q[3];
sx q[3];
rz(-2.6409918) q[3];
sx q[3];
rz(-2.4216381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.1392764) q[0];
sx q[0];
rz(-1.5895695) q[0];
sx q[0];
rz(-1.4578777) q[0];
rz(0.32683364) q[1];
sx q[1];
rz(-1.9953597) q[1];
sx q[1];
rz(-2.0976417) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69125873) q[0];
sx q[0];
rz(-2.674461) q[0];
sx q[0];
rz(-1.1410261) q[0];
x q[1];
rz(0.27803266) q[2];
sx q[2];
rz(-2.6671609) q[2];
sx q[2];
rz(2.5086046) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9324016) q[1];
sx q[1];
rz(-1.4889917) q[1];
sx q[1];
rz(1.9664759) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92071988) q[3];
sx q[3];
rz(-1.3082005) q[3];
sx q[3];
rz(-2.289734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30847183) q[2];
sx q[2];
rz(-1.0288419) q[2];
sx q[2];
rz(2.2157045) q[2];
rz(-1.067767) q[3];
sx q[3];
rz(-2.0177149) q[3];
sx q[3];
rz(-1.7759751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2866216) q[0];
sx q[0];
rz(-1.6678896) q[0];
sx q[0];
rz(0.58706748) q[0];
rz(-0.58285561) q[1];
sx q[1];
rz(-0.28935495) q[1];
sx q[1];
rz(0.48097441) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0868549) q[0];
sx q[0];
rz(-1.361892) q[0];
sx q[0];
rz(0.41685327) q[0];
rz(-pi) q[1];
rz(-2.9513861) q[2];
sx q[2];
rz(-1.2163194) q[2];
sx q[2];
rz(0.5917393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1744057) q[1];
sx q[1];
rz(-2.4250826) q[1];
sx q[1];
rz(-2.4661676) q[1];
rz(-pi) q[2];
rz(1.770788) q[3];
sx q[3];
rz(-1.963435) q[3];
sx q[3];
rz(-0.61905608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19758548) q[2];
sx q[2];
rz(-2.6052167) q[2];
sx q[2];
rz(1.3295004) q[2];
rz(0.78209376) q[3];
sx q[3];
rz(-0.32556459) q[3];
sx q[3];
rz(-1.1480924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6999577) q[0];
sx q[0];
rz(-1.4186207) q[0];
sx q[0];
rz(1.6765208) q[0];
rz(-1.5326473) q[1];
sx q[1];
rz(-1.1679222) q[1];
sx q[1];
rz(2.4293778) q[1];
rz(2.4423626) q[2];
sx q[2];
rz(-0.86730994) q[2];
sx q[2];
rz(2.2326474) q[2];
rz(-0.28859517) q[3];
sx q[3];
rz(-0.51641083) q[3];
sx q[3];
rz(-0.33752059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
