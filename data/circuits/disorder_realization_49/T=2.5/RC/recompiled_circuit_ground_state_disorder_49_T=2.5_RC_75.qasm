OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1385652) q[0];
sx q[0];
rz(-0.87831098) q[0];
sx q[0];
rz(-0.83100975) q[0];
rz(-0.57549685) q[1];
sx q[1];
rz(3.8776445) q[1];
sx q[1];
rz(10.164645) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54300302) q[0];
sx q[0];
rz(-2.9393688) q[0];
sx q[0];
rz(-0.36665066) q[0];
x q[1];
rz(-2.1388571) q[2];
sx q[2];
rz(-1.6913927) q[2];
sx q[2];
rz(1.2029778) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4783096) q[1];
sx q[1];
rz(-1.7170719) q[1];
sx q[1];
rz(1.0895132) q[1];
rz(-0.44828592) q[3];
sx q[3];
rz(-2.7196251) q[3];
sx q[3];
rz(-1.937605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2089219) q[2];
sx q[2];
rz(-0.17503665) q[2];
sx q[2];
rz(0.33898655) q[2];
rz(-2.637376) q[3];
sx q[3];
rz(-2.1239069) q[3];
sx q[3];
rz(-1.1241815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9826688) q[0];
sx q[0];
rz(-0.6944387) q[0];
sx q[0];
rz(-2.6320631) q[0];
rz(1.5024028) q[1];
sx q[1];
rz(-0.27703151) q[1];
sx q[1];
rz(0.94430077) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751289) q[0];
sx q[0];
rz(-0.24953609) q[0];
sx q[0];
rz(-0.93231045) q[0];
rz(1.239907) q[2];
sx q[2];
rz(-2.9941786) q[2];
sx q[2];
rz(-0.40357631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1183853) q[1];
sx q[1];
rz(-1.002335) q[1];
sx q[1];
rz(0.92250198) q[1];
rz(0.34667947) q[3];
sx q[3];
rz(-2.2190385) q[3];
sx q[3];
rz(-2.7559506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3462191) q[2];
sx q[2];
rz(-1.5495164) q[2];
sx q[2];
rz(-0.53595558) q[2];
rz(-1.0653982) q[3];
sx q[3];
rz(-0.25748101) q[3];
sx q[3];
rz(-2.4901701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1123493) q[0];
sx q[0];
rz(-1.9922682) q[0];
sx q[0];
rz(-2.405622) q[0];
rz(2.60587) q[1];
sx q[1];
rz(-1.842272) q[1];
sx q[1];
rz(-0.89964286) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9743703) q[0];
sx q[0];
rz(-1.456454) q[0];
sx q[0];
rz(-3.0024282) q[0];
rz(-2.1754335) q[2];
sx q[2];
rz(-1.7636429) q[2];
sx q[2];
rz(1.218534) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.31832987) q[1];
sx q[1];
rz(-1.0771275) q[1];
sx q[1];
rz(0.089463316) q[1];
rz(-2.8923762) q[3];
sx q[3];
rz(-1.936541) q[3];
sx q[3];
rz(-2.5178227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6003517) q[2];
sx q[2];
rz(-1.9466126) q[2];
sx q[2];
rz(0.80292732) q[2];
rz(-2.3927355) q[3];
sx q[3];
rz(-1.3978981) q[3];
sx q[3];
rz(-2.9822541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6244741) q[0];
sx q[0];
rz(-1.4117389) q[0];
sx q[0];
rz(-3.0112322) q[0];
rz(-1.8761926) q[1];
sx q[1];
rz(-1.1771026) q[1];
sx q[1];
rz(3.0317543) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8396436) q[0];
sx q[0];
rz(-0.57900864) q[0];
sx q[0];
rz(-2.2173475) q[0];
rz(-0.420094) q[2];
sx q[2];
rz(-0.60852988) q[2];
sx q[2];
rz(2.9836754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7848282) q[1];
sx q[1];
rz(-2.3127529) q[1];
sx q[1];
rz(-1.0395365) q[1];
rz(-pi) q[2];
rz(0.62105824) q[3];
sx q[3];
rz(-0.73535669) q[3];
sx q[3];
rz(1.3842954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.639223) q[2];
sx q[2];
rz(-1.1797735) q[2];
sx q[2];
rz(2.7184674) q[2];
rz(2.1028178) q[3];
sx q[3];
rz(-2.7602502) q[3];
sx q[3];
rz(1.0264621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0005242) q[0];
sx q[0];
rz(-1.2491913) q[0];
sx q[0];
rz(2.8619859) q[0];
rz(0.33310834) q[1];
sx q[1];
rz(-2.6393642) q[1];
sx q[1];
rz(2.8531029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97388291) q[0];
sx q[0];
rz(-1.5076734) q[0];
sx q[0];
rz(-1.0991251) q[0];
rz(1.8158378) q[2];
sx q[2];
rz(-1.6769655) q[2];
sx q[2];
rz(-2.9965056) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7964871) q[1];
sx q[1];
rz(-3.0287841) q[1];
sx q[1];
rz(-2.2276001) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1909423) q[3];
sx q[3];
rz(-2.4831746) q[3];
sx q[3];
rz(-1.0059716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.592411) q[2];
sx q[2];
rz(-1.2191399) q[2];
sx q[2];
rz(-2.9782817) q[2];
rz(-1.9773989) q[3];
sx q[3];
rz(-2.4653698) q[3];
sx q[3];
rz(-1.3588847) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0802245) q[0];
sx q[0];
rz(-3.1243262) q[0];
sx q[0];
rz(1.226271) q[0];
rz(1.8105043) q[1];
sx q[1];
rz(-2.2115579) q[1];
sx q[1];
rz(2.5221882) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50766599) q[0];
sx q[0];
rz(-1.0867501) q[0];
sx q[0];
rz(2.593301) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7707497) q[2];
sx q[2];
rz(-1.8189578) q[2];
sx q[2];
rz(-1.8979817) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.631397) q[1];
sx q[1];
rz(-1.7185154) q[1];
sx q[1];
rz(0.057675511) q[1];
rz(-pi) q[2];
rz(1.1167239) q[3];
sx q[3];
rz(-2.7259856) q[3];
sx q[3];
rz(-2.6799283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9161393) q[2];
sx q[2];
rz(-1.4832387) q[2];
sx q[2];
rz(-2.7062866) q[2];
rz(-0.12886038) q[3];
sx q[3];
rz(-1.83056) q[3];
sx q[3];
rz(0.35663566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.72047609) q[0];
sx q[0];
rz(-0.55835503) q[0];
sx q[0];
rz(2.5893353) q[0];
rz(2.5241959) q[1];
sx q[1];
rz(-1.9300902) q[1];
sx q[1];
rz(1.6389821) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5867859) q[0];
sx q[0];
rz(-1.9585591) q[0];
sx q[0];
rz(-0.047090637) q[0];
rz(-pi) q[1];
rz(-2.9947979) q[2];
sx q[2];
rz(-2.0710315) q[2];
sx q[2];
rz(-3.0584832) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0956816) q[1];
sx q[1];
rz(-1.0675745) q[1];
sx q[1];
rz(-2.9614425) q[1];
rz(2.8822717) q[3];
sx q[3];
rz(-1.590456) q[3];
sx q[3];
rz(-0.033084083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5804533) q[2];
sx q[2];
rz(-2.9353607) q[2];
sx q[2];
rz(-2.0533766) q[2];
rz(0.020180833) q[3];
sx q[3];
rz(-1.2729278) q[3];
sx q[3];
rz(-1.5816241) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35767558) q[0];
sx q[0];
rz(-2.7516784) q[0];
sx q[0];
rz(-1.0274603) q[0];
rz(-1.1032392) q[1];
sx q[1];
rz(-1.5682861) q[1];
sx q[1];
rz(-1.2148414) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88893696) q[0];
sx q[0];
rz(-1.2644469) q[0];
sx q[0];
rz(-3.0333748) q[0];
x q[1];
rz(1.096345) q[2];
sx q[2];
rz(-1.9812968) q[2];
sx q[2];
rz(0.83555789) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5826792) q[1];
sx q[1];
rz(-1.2238992) q[1];
sx q[1];
rz(-0.87474645) q[1];
rz(2.5276466) q[3];
sx q[3];
rz(-1.9910195) q[3];
sx q[3];
rz(-0.16638923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8814298) q[2];
sx q[2];
rz(-1.9378928) q[2];
sx q[2];
rz(2.0737958) q[2];
rz(-2.2504375) q[3];
sx q[3];
rz(-2.6813337) q[3];
sx q[3];
rz(1.1394181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41336173) q[0];
sx q[0];
rz(-0.75465337) q[0];
sx q[0];
rz(0.2555787) q[0];
rz(0.74343395) q[1];
sx q[1];
rz(-1.5579222) q[1];
sx q[1];
rz(1.5240634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96252464) q[0];
sx q[0];
rz(-2.0605378) q[0];
sx q[0];
rz(-2.6262002) q[0];
rz(0.41478283) q[2];
sx q[2];
rz(-2.4786502) q[2];
sx q[2];
rz(-2.489733) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0164106) q[1];
sx q[1];
rz(-2.4577854) q[1];
sx q[1];
rz(1.3586302) q[1];
x q[2];
rz(-2.2117046) q[3];
sx q[3];
rz(-1.7853338) q[3];
sx q[3];
rz(0.39232871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0301547) q[2];
sx q[2];
rz(-1.43575) q[2];
sx q[2];
rz(-1.6072404) q[2];
rz(-1.0715019) q[3];
sx q[3];
rz(-1.6000308) q[3];
sx q[3];
rz(-1.1736386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8906422) q[0];
sx q[0];
rz(-2.8548456) q[0];
sx q[0];
rz(-0.51837921) q[0];
rz(-2.3333343) q[1];
sx q[1];
rz(-1.4603442) q[1];
sx q[1];
rz(1.6917276) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1113365) q[0];
sx q[0];
rz(-1.5606631) q[0];
sx q[0];
rz(1.2819321) q[0];
x q[1];
rz(1.1495744) q[2];
sx q[2];
rz(-1.0738157) q[2];
sx q[2];
rz(0.62825655) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5603578) q[1];
sx q[1];
rz(-0.61040813) q[1];
sx q[1];
rz(-2.6950652) q[1];
rz(0.73734452) q[3];
sx q[3];
rz(-0.72190815) q[3];
sx q[3];
rz(-1.9487716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1125674) q[2];
sx q[2];
rz(-1.9064648) q[2];
sx q[2];
rz(0.9355363) q[2];
rz(-1.6701291) q[3];
sx q[3];
rz(-1.6664489) q[3];
sx q[3];
rz(-1.0940301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4569296) q[0];
sx q[0];
rz(-2.5419432) q[0];
sx q[0];
rz(2.3210617) q[0];
rz(-2.6928071) q[1];
sx q[1];
rz(-1.9955336) q[1];
sx q[1];
rz(1.21036) q[1];
rz(-1.6937428) q[2];
sx q[2];
rz(-0.37005432) q[2];
sx q[2];
rz(2.4126929) q[2];
rz(-2.5815677) q[3];
sx q[3];
rz(-1.8631794) q[3];
sx q[3];
rz(1.2003492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
