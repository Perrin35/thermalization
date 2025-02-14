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
rz(0.022194447) q[0];
sx q[0];
rz(3.7137346) q[0];
sx q[0];
rz(9.619286) q[0];
rz(-1.0021078) q[1];
sx q[1];
rz(-0.78295308) q[1];
sx q[1];
rz(3.1257544) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59796625) q[0];
sx q[0];
rz(-0.84613325) q[0];
sx q[0];
rz(-1.0084125) q[0];
rz(-2.3443647) q[2];
sx q[2];
rz(-1.8626894) q[2];
sx q[2];
rz(2.6780617) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8324278) q[1];
sx q[1];
rz(-2.2057071) q[1];
sx q[1];
rz(3.0135703) q[1];
rz(-1.7884618) q[3];
sx q[3];
rz(-2.7870745) q[3];
sx q[3];
rz(0.61787546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1281434) q[2];
sx q[2];
rz(-1.4712508) q[2];
sx q[2];
rz(0.61689287) q[2];
rz(-1.8938176) q[3];
sx q[3];
rz(-0.82472491) q[3];
sx q[3];
rz(-2.8491546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1013252) q[0];
sx q[0];
rz(-1.9693002) q[0];
sx q[0];
rz(-0.030666703) q[0];
rz(1.6568294) q[1];
sx q[1];
rz(-0.29850423) q[1];
sx q[1];
rz(-2.8224077) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6761665) q[0];
sx q[0];
rz(-2.5990708) q[0];
sx q[0];
rz(1.8886306) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99700751) q[2];
sx q[2];
rz(-2.6827723) q[2];
sx q[2];
rz(2.8850537) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.044925) q[1];
sx q[1];
rz(-0.54122347) q[1];
sx q[1];
rz(-2.2071532) q[1];
x q[2];
rz(-1.033554) q[3];
sx q[3];
rz(-1.1059009) q[3];
sx q[3];
rz(-3.1102095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.14898758) q[2];
sx q[2];
rz(-2.7379898) q[2];
sx q[2];
rz(-0.65275025) q[2];
rz(-2.2720689) q[3];
sx q[3];
rz(-2.0448989) q[3];
sx q[3];
rz(-2.7147527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94153786) q[0];
sx q[0];
rz(-2.8917199) q[0];
sx q[0];
rz(0.13724929) q[0];
rz(-2.6380154) q[1];
sx q[1];
rz(-2.5027687) q[1];
sx q[1];
rz(-2.8403179) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5001016) q[0];
sx q[0];
rz(-1.2610241) q[0];
sx q[0];
rz(0.70112438) q[0];
x q[1];
rz(-1.7687479) q[2];
sx q[2];
rz(-0.29445172) q[2];
sx q[2];
rz(2.1535923) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5168686) q[1];
sx q[1];
rz(-1.4249603) q[1];
sx q[1];
rz(-2.9353492) q[1];
rz(-pi) q[2];
rz(2.4739315) q[3];
sx q[3];
rz(-1.6958837) q[3];
sx q[3];
rz(1.1043105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0990024) q[2];
sx q[2];
rz(-2.4799535) q[2];
sx q[2];
rz(-0.46690565) q[2];
rz(-3.1344938) q[3];
sx q[3];
rz(-3.0341798) q[3];
sx q[3];
rz(-2.1527009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0840833) q[0];
sx q[0];
rz(-3.1355661) q[0];
sx q[0];
rz(2.4868593) q[0];
rz(-1.5701125) q[1];
sx q[1];
rz(-1.7127769) q[1];
sx q[1];
rz(0.44081259) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27106287) q[0];
sx q[0];
rz(-1.54633) q[0];
sx q[0];
rz(1.5464252) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.019182656) q[2];
sx q[2];
rz(-0.47878107) q[2];
sx q[2];
rz(-0.60730241) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7880169) q[1];
sx q[1];
rz(-1.792479) q[1];
sx q[1];
rz(0.30239971) q[1];
rz(-2.58617) q[3];
sx q[3];
rz(-0.54612386) q[3];
sx q[3];
rz(-2.6046942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.64678994) q[2];
sx q[2];
rz(-1.9014634) q[2];
sx q[2];
rz(0.49171641) q[2];
rz(1.803319) q[3];
sx q[3];
rz(-2.0526142) q[3];
sx q[3];
rz(1.4664388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76984513) q[0];
sx q[0];
rz(-1.0862229) q[0];
sx q[0];
rz(-1.0787429) q[0];
rz(0.046401333) q[1];
sx q[1];
rz(-1.6175783) q[1];
sx q[1];
rz(2.0414415) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4217401) q[0];
sx q[0];
rz(-1.5800086) q[0];
sx q[0];
rz(-3.0250027) q[0];
rz(0.82450878) q[2];
sx q[2];
rz(-0.84133278) q[2];
sx q[2];
rz(1.0223323) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0667923) q[1];
sx q[1];
rz(-0.88035816) q[1];
sx q[1];
rz(-1.2509173) q[1];
rz(-pi) q[2];
rz(1.6836105) q[3];
sx q[3];
rz(-1.8387018) q[3];
sx q[3];
rz(-0.15427854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2940353) q[2];
sx q[2];
rz(-0.39176771) q[2];
sx q[2];
rz(2.7993287) q[2];
rz(-1.3058454) q[3];
sx q[3];
rz(-1.2021474) q[3];
sx q[3];
rz(-1.1024629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24868988) q[0];
sx q[0];
rz(-1.1868287) q[0];
sx q[0];
rz(-2.7500395) q[0];
rz(0.64741099) q[1];
sx q[1];
rz(-2.7029111) q[1];
sx q[1];
rz(-2.2364) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0150512) q[0];
sx q[0];
rz(-2.371696) q[0];
sx q[0];
rz(-2.4036744) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91475418) q[2];
sx q[2];
rz(-2.0080655) q[2];
sx q[2];
rz(0.8503051) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.47899095) q[1];
sx q[1];
rz(-2.4756065) q[1];
sx q[1];
rz(2.7061092) q[1];
rz(-0.16512434) q[3];
sx q[3];
rz(-0.29224631) q[3];
sx q[3];
rz(-1.0632893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4539926) q[2];
sx q[2];
rz(-0.31012055) q[2];
sx q[2];
rz(2.0886759) q[2];
rz(-0.29371253) q[3];
sx q[3];
rz(-0.83455825) q[3];
sx q[3];
rz(-0.53148758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3783136) q[0];
sx q[0];
rz(-1.2792307) q[0];
sx q[0];
rz(0.026750201) q[0];
rz(-1.1862952) q[1];
sx q[1];
rz(-0.18272884) q[1];
sx q[1];
rz(0.15348405) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54220056) q[0];
sx q[0];
rz(-1.3017103) q[0];
sx q[0];
rz(2.0890636) q[0];
rz(0.085603733) q[2];
sx q[2];
rz(-2.589648) q[2];
sx q[2];
rz(-2.7526223) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.80375043) q[1];
sx q[1];
rz(-2.1598247) q[1];
sx q[1];
rz(0.71658012) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1128902) q[3];
sx q[3];
rz(-1.6316292) q[3];
sx q[3];
rz(-0.83567747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.89580506) q[2];
sx q[2];
rz(-0.20496002) q[2];
sx q[2];
rz(-1.2650355) q[2];
rz(-0.17299077) q[3];
sx q[3];
rz(-1.3909307) q[3];
sx q[3];
rz(2.0632108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10385253) q[0];
sx q[0];
rz(-0.2094035) q[0];
sx q[0];
rz(0.10845342) q[0];
rz(-1.7697889) q[1];
sx q[1];
rz(-1.3550974) q[1];
sx q[1];
rz(3.1350737) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2222424) q[0];
sx q[0];
rz(-1.1092289) q[0];
sx q[0];
rz(0.95141769) q[0];
rz(-pi) q[1];
rz(0.44640572) q[2];
sx q[2];
rz(-2.8978928) q[2];
sx q[2];
rz(-2.6350465) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7830917) q[1];
sx q[1];
rz(-2.3066562) q[1];
sx q[1];
rz(2.3651644) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7192332) q[3];
sx q[3];
rz(-0.82882753) q[3];
sx q[3];
rz(-3.1008848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.45695496) q[2];
sx q[2];
rz(-0.64720398) q[2];
sx q[2];
rz(-1.4400488) q[2];
rz(0.35259926) q[3];
sx q[3];
rz(-0.86360258) q[3];
sx q[3];
rz(0.45261228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0532992) q[0];
sx q[0];
rz(-1.091401) q[0];
sx q[0];
rz(1.1540867) q[0];
rz(-1.6905009) q[1];
sx q[1];
rz(-2.4595551) q[1];
sx q[1];
rz(2.3114204) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3367171) q[0];
sx q[0];
rz(-1.5815539) q[0];
sx q[0];
rz(-1.5740001) q[0];
rz(-pi) q[1];
rz(-1.7042034) q[2];
sx q[2];
rz(-2.0683943) q[2];
sx q[2];
rz(-1.069151) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7651319) q[1];
sx q[1];
rz(-1.9674976) q[1];
sx q[1];
rz(2.6668616) q[1];
x q[2];
rz(-0.19751246) q[3];
sx q[3];
rz(-2.2220051) q[3];
sx q[3];
rz(-0.75266389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.95149583) q[2];
sx q[2];
rz(-0.87679619) q[2];
sx q[2];
rz(2.4693176) q[2];
rz(-0.43664524) q[3];
sx q[3];
rz(-1.3029706) q[3];
sx q[3];
rz(-2.8521027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.745568) q[0];
sx q[0];
rz(-2.4085299) q[0];
sx q[0];
rz(-0.42098862) q[0];
rz(1.1517395) q[1];
sx q[1];
rz(-0.20944171) q[1];
sx q[1];
rz(2.4236325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67506266) q[0];
sx q[0];
rz(-2.8631239) q[0];
sx q[0];
rz(-1.6519) q[0];
rz(-2.6804982) q[2];
sx q[2];
rz(-2.3873219) q[2];
sx q[2];
rz(2.8476771) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3612671) q[1];
sx q[1];
rz(-0.33814374) q[1];
sx q[1];
rz(-1.9664155) q[1];
x q[2];
rz(2.8676492) q[3];
sx q[3];
rz(-2.0060325) q[3];
sx q[3];
rz(2.1147605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.026160508) q[2];
sx q[2];
rz(-2.5371976) q[2];
sx q[2];
rz(-1.3755414) q[2];
rz(-0.17624217) q[3];
sx q[3];
rz(-0.79564375) q[3];
sx q[3];
rz(0.073504329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7701223) q[0];
sx q[0];
rz(-1.6902516) q[0];
sx q[0];
rz(-1.1023735) q[0];
rz(-2.2813028) q[1];
sx q[1];
rz(-1.5033036) q[1];
sx q[1];
rz(-1.1183429) q[1];
rz(-3.0768298) q[2];
sx q[2];
rz(-0.76640609) q[2];
sx q[2];
rz(0.7647209) q[2];
rz(2.7793814) q[3];
sx q[3];
rz(-1.9345761) q[3];
sx q[3];
rz(-3.0349845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
