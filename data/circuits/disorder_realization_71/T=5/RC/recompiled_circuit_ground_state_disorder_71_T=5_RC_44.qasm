OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7856287) q[0];
sx q[0];
rz(-0.2956737) q[0];
sx q[0];
rz(9.8386968) q[0];
rz(-2.8726481) q[1];
sx q[1];
rz(-2.7719331) q[1];
sx q[1];
rz(1.2764021) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7463425) q[0];
sx q[0];
rz(-1.479404) q[0];
sx q[0];
rz(-1.6618098) q[0];
rz(-pi) q[1];
rz(-0.062687473) q[2];
sx q[2];
rz(-0.90342316) q[2];
sx q[2];
rz(1.1476715) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.40117404) q[1];
sx q[1];
rz(-0.83073069) q[1];
sx q[1];
rz(0.38895815) q[1];
x q[2];
rz(-1.7539361) q[3];
sx q[3];
rz(-2.1570342) q[3];
sx q[3];
rz(2.8318015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3966177) q[2];
sx q[2];
rz(-1.9048385) q[2];
sx q[2];
rz(1.2006987) q[2];
rz(-2.4602304) q[3];
sx q[3];
rz(-1.5228289) q[3];
sx q[3];
rz(-2.7549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.50614) q[0];
sx q[0];
rz(-2.8400087) q[0];
sx q[0];
rz(-2.8801081) q[0];
rz(-1.5226978) q[1];
sx q[1];
rz(-0.50917429) q[1];
sx q[1];
rz(-1.8656628) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8573924) q[0];
sx q[0];
rz(-2.5821857) q[0];
sx q[0];
rz(0.18226923) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66565565) q[2];
sx q[2];
rz(-2.3710069) q[2];
sx q[2];
rz(1.4879172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6258982) q[1];
sx q[1];
rz(-2.1014903) q[1];
sx q[1];
rz(1.1926535) q[1];
rz(-pi) q[2];
rz(-1.3456144) q[3];
sx q[3];
rz(-2.9888267) q[3];
sx q[3];
rz(1.6880715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2168938) q[2];
sx q[2];
rz(-1.0328707) q[2];
sx q[2];
rz(2.4135598) q[2];
rz(-1.3701471) q[3];
sx q[3];
rz(-2.5244505) q[3];
sx q[3];
rz(-3.106015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1074693) q[0];
sx q[0];
rz(-1.1119482) q[0];
sx q[0];
rz(0.6849826) q[0];
rz(1.1032387) q[1];
sx q[1];
rz(-0.61956844) q[1];
sx q[1];
rz(-1.0122976) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89268003) q[0];
sx q[0];
rz(-1.2272226) q[0];
sx q[0];
rz(1.1569886) q[0];
rz(-pi) q[1];
rz(0.84350297) q[2];
sx q[2];
rz(-2.8828104) q[2];
sx q[2];
rz(2.7121787) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.60437459) q[1];
sx q[1];
rz(-2.8979725) q[1];
sx q[1];
rz(0.55018665) q[1];
rz(-0.52148444) q[3];
sx q[3];
rz(-2.3064724) q[3];
sx q[3];
rz(-2.2948752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95969069) q[2];
sx q[2];
rz(-1.5896553) q[2];
sx q[2];
rz(0.8054136) q[2];
rz(-0.81625932) q[3];
sx q[3];
rz(-2.5160242) q[3];
sx q[3];
rz(3.0474385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-3.1212946) q[0];
sx q[0];
rz(-2.2327406) q[0];
sx q[0];
rz(2.5138309) q[0];
rz(0.10920564) q[1];
sx q[1];
rz(-0.86995482) q[1];
sx q[1];
rz(2.3039718) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19310424) q[0];
sx q[0];
rz(-2.645506) q[0];
sx q[0];
rz(-0.14004616) q[0];
rz(2.1659746) q[2];
sx q[2];
rz(-1.383923) q[2];
sx q[2];
rz(2.486791) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8404897) q[1];
sx q[1];
rz(-0.87574848) q[1];
sx q[1];
rz(1.7453421) q[1];
rz(-pi) q[2];
rz(-0.37552278) q[3];
sx q[3];
rz(-1.8584455) q[3];
sx q[3];
rz(2.875553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8236905) q[2];
sx q[2];
rz(-2.6573942) q[2];
sx q[2];
rz(-2.79706) q[2];
rz(1.4065929) q[3];
sx q[3];
rz(-0.35511261) q[3];
sx q[3];
rz(-0.047903456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.542881) q[0];
sx q[0];
rz(-0.6364091) q[0];
sx q[0];
rz(2.6357546) q[0];
rz(2.089031) q[1];
sx q[1];
rz(-0.86741766) q[1];
sx q[1];
rz(2.1956086) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81954581) q[0];
sx q[0];
rz(-1.2975386) q[0];
sx q[0];
rz(-1.8332493) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7713114) q[2];
sx q[2];
rz(-0.38114377) q[2];
sx q[2];
rz(-0.76424341) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8648423) q[1];
sx q[1];
rz(-0.71886834) q[1];
sx q[1];
rz(-2.5548773) q[1];
x q[2];
rz(-1.1422154) q[3];
sx q[3];
rz(-0.16332291) q[3];
sx q[3];
rz(1.4455011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.08789173) q[2];
sx q[2];
rz(-1.0975857) q[2];
sx q[2];
rz(1.5778479) q[2];
rz(2.0569233) q[3];
sx q[3];
rz(-0.8756777) q[3];
sx q[3];
rz(0.55115551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56187335) q[0];
sx q[0];
rz(-0.60573524) q[0];
sx q[0];
rz(-2.9015923) q[0];
rz(2.2198246) q[1];
sx q[1];
rz(-1.0488989) q[1];
sx q[1];
rz(1.1704495) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3055127) q[0];
sx q[0];
rz(-1.5275947) q[0];
sx q[0];
rz(-0.57010285) q[0];
rz(-pi) q[1];
rz(0.35640772) q[2];
sx q[2];
rz(-1.6898167) q[2];
sx q[2];
rz(-1.0072034) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7266143) q[1];
sx q[1];
rz(-1.1220334) q[1];
sx q[1];
rz(-2.1244551) q[1];
rz(-pi) q[2];
rz(-1.9715973) q[3];
sx q[3];
rz(-1.6605536) q[3];
sx q[3];
rz(-2.2211162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74448284) q[2];
sx q[2];
rz(-2.3220671) q[2];
sx q[2];
rz(-2.5605555) q[2];
rz(-0.93823141) q[3];
sx q[3];
rz(-2.5751028) q[3];
sx q[3];
rz(0.70403045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13643232) q[0];
sx q[0];
rz(-2.4649824) q[0];
sx q[0];
rz(2.6410979) q[0];
rz(2.9035134) q[1];
sx q[1];
rz(-1.6903189) q[1];
sx q[1];
rz(2.4949825) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3377285) q[0];
sx q[0];
rz(-3.1067305) q[0];
sx q[0];
rz(0.050970391) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4562723) q[2];
sx q[2];
rz(-1.6045286) q[2];
sx q[2];
rz(-2.9303868) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.889786) q[1];
sx q[1];
rz(-1.4507035) q[1];
sx q[1];
rz(1.0838064) q[1];
x q[2];
rz(-3.079861) q[3];
sx q[3];
rz(-1.2128856) q[3];
sx q[3];
rz(0.45321143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2639192) q[2];
sx q[2];
rz(-2.4223902) q[2];
sx q[2];
rz(0.70362299) q[2];
rz(-1.8536812) q[3];
sx q[3];
rz(-2.0813007) q[3];
sx q[3];
rz(-0.50584403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25381655) q[0];
sx q[0];
rz(-2.0075338) q[0];
sx q[0];
rz(-1.8164841) q[0];
rz(-2.3690986) q[1];
sx q[1];
rz(-1.1224116) q[1];
sx q[1];
rz(-2.966029) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85925519) q[0];
sx q[0];
rz(-1.3842177) q[0];
sx q[0];
rz(-0.76754359) q[0];
x q[1];
rz(0.80412038) q[2];
sx q[2];
rz(-1.0808699) q[2];
sx q[2];
rz(-2.8049289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7372706) q[1];
sx q[1];
rz(-2.0312738) q[1];
sx q[1];
rz(-0.15562017) q[1];
rz(-pi) q[2];
rz(1.6962849) q[3];
sx q[3];
rz(-0.85832046) q[3];
sx q[3];
rz(0.10601302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10302155) q[2];
sx q[2];
rz(-2.0671637) q[2];
sx q[2];
rz(-2.0477022) q[2];
rz(-1.1755747) q[3];
sx q[3];
rz(-2.3837377) q[3];
sx q[3];
rz(2.8567543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4232101) q[0];
sx q[0];
rz(-0.045561401) q[0];
sx q[0];
rz(1.8157995) q[0];
rz(-2.6153053) q[1];
sx q[1];
rz(-2.278639) q[1];
sx q[1];
rz(2.3525499) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1068971) q[0];
sx q[0];
rz(-1.7321087) q[0];
sx q[0];
rz(-3.0437896) q[0];
rz(-pi) q[1];
rz(1.3157719) q[2];
sx q[2];
rz(-1.2256283) q[2];
sx q[2];
rz(-1.6127301) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.10875094) q[1];
sx q[1];
rz(-2.1819253) q[1];
sx q[1];
rz(2.717589) q[1];
rz(-pi) q[2];
rz(2.3655543) q[3];
sx q[3];
rz(-2.0774842) q[3];
sx q[3];
rz(-0.53467487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.21251707) q[2];
sx q[2];
rz(-2.2079461) q[2];
sx q[2];
rz(2.2361501) q[2];
rz(-0.91600156) q[3];
sx q[3];
rz(-1.6986877) q[3];
sx q[3];
rz(-2.0974832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7020029) q[0];
sx q[0];
rz(-0.8332533) q[0];
sx q[0];
rz(-2.2138017) q[0];
rz(-2.0480428) q[1];
sx q[1];
rz(-2.5826726) q[1];
sx q[1];
rz(-0.13339001) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0463842) q[0];
sx q[0];
rz(-0.36050561) q[0];
sx q[0];
rz(0.62490873) q[0];
x q[1];
rz(0.20602823) q[2];
sx q[2];
rz(-1.0090172) q[2];
sx q[2];
rz(0.59625193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9815753) q[1];
sx q[1];
rz(-0.80343738) q[1];
sx q[1];
rz(0.74635069) q[1];
rz(-pi) q[2];
rz(2.9319256) q[3];
sx q[3];
rz(-1.121796) q[3];
sx q[3];
rz(1.0202788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6371969) q[2];
sx q[2];
rz(-1.7375676) q[2];
sx q[2];
rz(-1.8782328) q[2];
rz(-2.0718306) q[3];
sx q[3];
rz(-1.0354778) q[3];
sx q[3];
rz(-3.0647762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99061154) q[0];
sx q[0];
rz(-2.4438416) q[0];
sx q[0];
rz(-1.7846815) q[0];
rz(1.6928584) q[1];
sx q[1];
rz(-0.81304638) q[1];
sx q[1];
rz(-0.1437694) q[1];
rz(-2.0781381) q[2];
sx q[2];
rz(-0.58488556) q[2];
sx q[2];
rz(1.9716138) q[2];
rz(1.9142022) q[3];
sx q[3];
rz(-1.3385942) q[3];
sx q[3];
rz(-0.75029324) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
