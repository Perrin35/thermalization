OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5175944) q[0];
sx q[0];
rz(-2.6156293) q[0];
sx q[0];
rz(2.9232803) q[0];
rz(-1.6649618) q[1];
sx q[1];
rz(2.6638439) q[1];
sx q[1];
rz(12.01241) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71611878) q[0];
sx q[0];
rz(-1.1083352) q[0];
sx q[0];
rz(1.0008971) q[0];
rz(-pi) q[1];
rz(0.10138114) q[2];
sx q[2];
rz(-0.53115618) q[2];
sx q[2];
rz(0.44975933) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.40383717) q[1];
sx q[1];
rz(-0.74450508) q[1];
sx q[1];
rz(-1.8831403) q[1];
rz(0.64726909) q[3];
sx q[3];
rz(-0.64220253) q[3];
sx q[3];
rz(0.49261478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.37609425) q[2];
sx q[2];
rz(-2.8480397) q[2];
sx q[2];
rz(-1.2676839) q[2];
rz(-2.3119161) q[3];
sx q[3];
rz(-1.6206348) q[3];
sx q[3];
rz(-0.42384306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053452881) q[0];
sx q[0];
rz(-1.8064878) q[0];
sx q[0];
rz(-1.394519) q[0];
rz(-2.246619) q[1];
sx q[1];
rz(-2.112969) q[1];
sx q[1];
rz(1.5825533) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6000309) q[0];
sx q[0];
rz(-2.2293334) q[0];
sx q[0];
rz(-0.70777871) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3560654) q[2];
sx q[2];
rz(-1.967397) q[2];
sx q[2];
rz(0.44677904) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6239418) q[1];
sx q[1];
rz(-1.883916) q[1];
sx q[1];
rz(2.4409358) q[1];
rz(-pi) q[2];
rz(0.30541916) q[3];
sx q[3];
rz(-1.1608184) q[3];
sx q[3];
rz(0.78193203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6643657) q[2];
sx q[2];
rz(-1.5218647) q[2];
sx q[2];
rz(0.030755432) q[2];
rz(0.45298806) q[3];
sx q[3];
rz(-0.22946295) q[3];
sx q[3];
rz(-0.17791137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7396616) q[0];
sx q[0];
rz(-3.0406096) q[0];
sx q[0];
rz(-2.3133551) q[0];
rz(3.0939057) q[1];
sx q[1];
rz(-2.2779155) q[1];
sx q[1];
rz(1.9140859) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21624529) q[0];
sx q[0];
rz(-2.0106533) q[0];
sx q[0];
rz(0.60103215) q[0];
x q[1];
rz(2.9256608) q[2];
sx q[2];
rz(-1.1732444) q[2];
sx q[2];
rz(0.88653195) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32900287) q[1];
sx q[1];
rz(-2.0101476) q[1];
sx q[1];
rz(-2.8529608) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2994386) q[3];
sx q[3];
rz(-1.531736) q[3];
sx q[3];
rz(-1.0674764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0487655) q[2];
sx q[2];
rz(-2.2266677) q[2];
sx q[2];
rz(1.1009334) q[2];
rz(-3.1277505) q[3];
sx q[3];
rz(-1.3661386) q[3];
sx q[3];
rz(-2.3069416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5562627) q[0];
sx q[0];
rz(-1.4034554) q[0];
sx q[0];
rz(-2.8072667) q[0];
rz(-2.4090134) q[1];
sx q[1];
rz(-2.1926447) q[1];
sx q[1];
rz(0.11071959) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77483141) q[0];
sx q[0];
rz(-2.3322421) q[0];
sx q[0];
rz(-1.5917718) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16558318) q[2];
sx q[2];
rz(-1.4787758) q[2];
sx q[2];
rz(-1.657287) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9690937) q[1];
sx q[1];
rz(-2.2374472) q[1];
sx q[1];
rz(2.4039925) q[1];
rz(-0.7752876) q[3];
sx q[3];
rz(-1.9248157) q[3];
sx q[3];
rz(-2.275685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4577786) q[2];
sx q[2];
rz(-0.28202287) q[2];
sx q[2];
rz(1.4078183) q[2];
rz(-1.1553361) q[3];
sx q[3];
rz(-1.1563533) q[3];
sx q[3];
rz(0.19481625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4487576) q[0];
sx q[0];
rz(-0.91788569) q[0];
sx q[0];
rz(2.1441929) q[0];
rz(1.6150486) q[1];
sx q[1];
rz(-2.5020182) q[1];
sx q[1];
rz(-1.4195199) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90516312) q[0];
sx q[0];
rz(-1.4078724) q[0];
sx q[0];
rz(1.7927367) q[0];
rz(-pi) q[1];
rz(-1.4155343) q[2];
sx q[2];
rz(-2.3387032) q[2];
sx q[2];
rz(2.2096095) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.43806009) q[1];
sx q[1];
rz(-2.4171481) q[1];
sx q[1];
rz(-1.0781592) q[1];
rz(0.12357281) q[3];
sx q[3];
rz(-1.6032748) q[3];
sx q[3];
rz(-2.6952254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.851696) q[2];
sx q[2];
rz(-0.09446129) q[2];
sx q[2];
rz(-0.32290253) q[2];
rz(-2.0276535) q[3];
sx q[3];
rz(-2.0506004) q[3];
sx q[3];
rz(-2.7285301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7776529) q[0];
sx q[0];
rz(-2.8033065) q[0];
sx q[0];
rz(-0.80663484) q[0];
rz(2.5576162) q[1];
sx q[1];
rz(-1.1176502) q[1];
sx q[1];
rz(-1.1899828) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2428044) q[0];
sx q[0];
rz(-1.4477647) q[0];
sx q[0];
rz(-1.9210299) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7940668) q[2];
sx q[2];
rz(-2.579877) q[2];
sx q[2];
rz(-2.7507741) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.62184238) q[1];
sx q[1];
rz(-1.6222266) q[1];
sx q[1];
rz(-1.7396881) q[1];
rz(2.9526918) q[3];
sx q[3];
rz(-2.402225) q[3];
sx q[3];
rz(-2.1738571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7334062) q[2];
sx q[2];
rz(-1.6383645) q[2];
sx q[2];
rz(-0.1850941) q[2];
rz(1.5271651) q[3];
sx q[3];
rz(-1.7388758) q[3];
sx q[3];
rz(-0.68814284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.9516893) q[0];
sx q[0];
rz(-2.1898495) q[0];
sx q[0];
rz(2.9290747) q[0];
rz(2.3566133) q[1];
sx q[1];
rz(-1.6467983) q[1];
sx q[1];
rz(-0.77883887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23379743) q[0];
sx q[0];
rz(-2.1229232) q[0];
sx q[0];
rz(0.88498451) q[0];
rz(-pi) q[1];
rz(-2.5539407) q[2];
sx q[2];
rz(-1.2562804) q[2];
sx q[2];
rz(1.8314198) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1800715) q[1];
sx q[1];
rz(-2.8043282) q[1];
sx q[1];
rz(0.86581655) q[1];
x q[2];
rz(-2.6565348) q[3];
sx q[3];
rz(-2.2211255) q[3];
sx q[3];
rz(1.2984087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6276041) q[2];
sx q[2];
rz(-0.49391654) q[2];
sx q[2];
rz(-2.7331875) q[2];
rz(0.70185316) q[3];
sx q[3];
rz(-0.93188325) q[3];
sx q[3];
rz(-1.2592038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34847611) q[0];
sx q[0];
rz(-1.2154673) q[0];
sx q[0];
rz(-0.066019639) q[0];
rz(1.6325715) q[1];
sx q[1];
rz(-1.0450109) q[1];
sx q[1];
rz(2.1836233) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97236605) q[0];
sx q[0];
rz(-1.1145076) q[0];
sx q[0];
rz(1.2537987) q[0];
x q[1];
rz(0.81664576) q[2];
sx q[2];
rz(-1.1743675) q[2];
sx q[2];
rz(-2.3289837) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4730395) q[1];
sx q[1];
rz(-0.91382342) q[1];
sx q[1];
rz(-1.9690352) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.041298) q[3];
sx q[3];
rz(-1.417932) q[3];
sx q[3];
rz(2.4748442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3079188) q[2];
sx q[2];
rz(-1.5235528) q[2];
sx q[2];
rz(-1.7306805) q[2];
rz(2.446567) q[3];
sx q[3];
rz(-1.6618988) q[3];
sx q[3];
rz(3.1237349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0628292) q[0];
sx q[0];
rz(-1.48209) q[0];
sx q[0];
rz(-0.98989809) q[0];
rz(0.46317378) q[1];
sx q[1];
rz(-1.4521234) q[1];
sx q[1];
rz(1.2407726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.038886) q[0];
sx q[0];
rz(-0.30065824) q[0];
sx q[0];
rz(-1.7666398) q[0];
rz(-1.1282721) q[2];
sx q[2];
rz(-1.754369) q[2];
sx q[2];
rz(-1.6115481) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5338143) q[1];
sx q[1];
rz(-1.3097714) q[1];
sx q[1];
rz(2.8487474) q[1];
rz(-pi) q[2];
rz(-2.1911591) q[3];
sx q[3];
rz(-1.1945121) q[3];
sx q[3];
rz(-2.2971414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3339633) q[2];
sx q[2];
rz(-1.3796076) q[2];
sx q[2];
rz(2.4228952) q[2];
rz(1.1527609) q[3];
sx q[3];
rz(-1.4216239) q[3];
sx q[3];
rz(-1.0144455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54866791) q[0];
sx q[0];
rz(-2.7935226) q[0];
sx q[0];
rz(-2.3396709) q[0];
rz(-2.0536664) q[1];
sx q[1];
rz(-1.8049003) q[1];
sx q[1];
rz(0.71802872) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3182326) q[0];
sx q[0];
rz(-1.8917068) q[0];
sx q[0];
rz(-0.72707392) q[0];
x q[1];
rz(-0.52387107) q[2];
sx q[2];
rz(-2.0064559) q[2];
sx q[2];
rz(2.4031529) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1433318) q[1];
sx q[1];
rz(-1.623053) q[1];
sx q[1];
rz(1.2388171) q[1];
rz(-pi) q[2];
rz(-2.2949785) q[3];
sx q[3];
rz(-0.9808971) q[3];
sx q[3];
rz(-2.5835832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3512909) q[2];
sx q[2];
rz(-2.9238034) q[2];
sx q[2];
rz(2.6289319) q[2];
rz(-0.74448186) q[3];
sx q[3];
rz(-2.1148041) q[3];
sx q[3];
rz(2.3775533) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22513334) q[0];
sx q[0];
rz(-1.2703348) q[0];
sx q[0];
rz(-1.2402007) q[0];
rz(0.71612877) q[1];
sx q[1];
rz(-2.5934673) q[1];
sx q[1];
rz(0.60636884) q[1];
rz(-3.0307583) q[2];
sx q[2];
rz(-1.7946984) q[2];
sx q[2];
rz(2.6738965) q[2];
rz(-1.2342831) q[3];
sx q[3];
rz(-1.6675622) q[3];
sx q[3];
rz(0.29940816) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
