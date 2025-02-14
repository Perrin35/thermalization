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
rz(0.12240527) q[0];
sx q[0];
rz(-0.91949099) q[0];
sx q[0];
rz(1.9424196) q[0];
rz(0.1872669) q[1];
sx q[1];
rz(-2.5993102) q[1];
sx q[1];
rz(1.6220925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9334979) q[0];
sx q[0];
rz(-1.6483432) q[0];
sx q[0];
rz(2.5258875) q[0];
rz(-pi) q[1];
rz(1.1178944) q[2];
sx q[2];
rz(-1.2317961) q[2];
sx q[2];
rz(1.9725368) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.853316) q[1];
sx q[1];
rz(-1.9717934) q[1];
sx q[1];
rz(2.6021491) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19272645) q[3];
sx q[3];
rz(-2.8929919) q[3];
sx q[3];
rz(-0.981244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2970994) q[2];
sx q[2];
rz(-0.62701925) q[2];
sx q[2];
rz(1.373488) q[2];
rz(0.80396906) q[3];
sx q[3];
rz(-1.6425902) q[3];
sx q[3];
rz(-1.9544301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3961769) q[0];
sx q[0];
rz(-2.4276623) q[0];
sx q[0];
rz(-2.7667238) q[0];
rz(-0.51775852) q[1];
sx q[1];
rz(-1.2034143) q[1];
sx q[1];
rz(1.3688603) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2320011) q[0];
sx q[0];
rz(-1.5717432) q[0];
sx q[0];
rz(-1.5707004) q[0];
x q[1];
rz(-1.5438186) q[2];
sx q[2];
rz(-1.3744053) q[2];
sx q[2];
rz(2.2324004) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79530935) q[1];
sx q[1];
rz(-2.4716271) q[1];
sx q[1];
rz(-2.4479292) q[1];
rz(-1.0090461) q[3];
sx q[3];
rz(-0.25190946) q[3];
sx q[3];
rz(-2.589165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.123473) q[2];
sx q[2];
rz(-2.4407083) q[2];
sx q[2];
rz(2.4269721) q[2];
rz(2.9361652) q[3];
sx q[3];
rz(-1.3222062) q[3];
sx q[3];
rz(1.2986758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6767204) q[0];
sx q[0];
rz(-1.0370075) q[0];
sx q[0];
rz(2.6439457) q[0];
rz(2.5045555) q[1];
sx q[1];
rz(-0.62890816) q[1];
sx q[1];
rz(-0.10674891) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80947948) q[0];
sx q[0];
rz(-0.3705855) q[0];
sx q[0];
rz(-0.95032145) q[0];
rz(-1.5791513) q[2];
sx q[2];
rz(-0.66973493) q[2];
sx q[2];
rz(1.005203) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5046204) q[1];
sx q[1];
rz(-2.2832738) q[1];
sx q[1];
rz(3.0709355) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.867093) q[3];
sx q[3];
rz(-1.955964) q[3];
sx q[3];
rz(2.9134458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.60483021) q[2];
sx q[2];
rz(-2.3894775) q[2];
sx q[2];
rz(2.9467648) q[2];
rz(-0.005006494) q[3];
sx q[3];
rz(-2.5275793) q[3];
sx q[3];
rz(2.3495638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0358589) q[0];
sx q[0];
rz(-2.6073313) q[0];
sx q[0];
rz(2.1015097) q[0];
rz(0.19373521) q[1];
sx q[1];
rz(-1.6400784) q[1];
sx q[1];
rz(-0.74660444) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0464942) q[0];
sx q[0];
rz(-1.455716) q[0];
sx q[0];
rz(1.1287354) q[0];
x q[1];
rz(-2.2708621) q[2];
sx q[2];
rz(-1.4615897) q[2];
sx q[2];
rz(-2.2094215) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0185701) q[1];
sx q[1];
rz(-2.445652) q[1];
sx q[1];
rz(-0.080841149) q[1];
rz(-pi) q[2];
rz(-1.6799404) q[3];
sx q[3];
rz(-1.8092844) q[3];
sx q[3];
rz(-1.7643238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99183434) q[2];
sx q[2];
rz(-1.6712302) q[2];
sx q[2];
rz(-0.20450083) q[2];
rz(1.9483942) q[3];
sx q[3];
rz(-2.2843993) q[3];
sx q[3];
rz(0.96113718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1292773) q[0];
sx q[0];
rz(-2.9670872) q[0];
sx q[0];
rz(1.0804863) q[0];
rz(-0.37733817) q[1];
sx q[1];
rz(-1.6308547) q[1];
sx q[1];
rz(-0.89471716) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6166234) q[0];
sx q[0];
rz(-1.3436514) q[0];
sx q[0];
rz(2.7779915) q[0];
rz(-0.9099877) q[2];
sx q[2];
rz(-2.2468781) q[2];
sx q[2];
rz(-2.8582339) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1929568) q[1];
sx q[1];
rz(-2.0564449) q[1];
sx q[1];
rz(0.57106496) q[1];
x q[2];
rz(-0.36578806) q[3];
sx q[3];
rz(-0.96082965) q[3];
sx q[3];
rz(2.6637885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.47199029) q[2];
sx q[2];
rz(-0.44595465) q[2];
sx q[2];
rz(-0.10776821) q[2];
rz(0.17555155) q[3];
sx q[3];
rz(-1.8864417) q[3];
sx q[3];
rz(0.56041437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0948697) q[0];
sx q[0];
rz(-2.6304498) q[0];
sx q[0];
rz(-1.0119337) q[0];
rz(1.5049505) q[1];
sx q[1];
rz(-0.71957809) q[1];
sx q[1];
rz(1.0171657) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95421529) q[0];
sx q[0];
rz(-0.34437505) q[0];
sx q[0];
rz(-1.0402753) q[0];
rz(-pi) q[1];
rz(-2.1332729) q[2];
sx q[2];
rz(-0.89300821) q[2];
sx q[2];
rz(0.00046367292) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0067082891) q[1];
sx q[1];
rz(-1.5794288) q[1];
sx q[1];
rz(-2.951202) q[1];
rz(-pi) q[2];
rz(-1.8886376) q[3];
sx q[3];
rz(-0.46964619) q[3];
sx q[3];
rz(1.8152678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73198685) q[2];
sx q[2];
rz(-0.66183949) q[2];
sx q[2];
rz(-0.99676639) q[2];
rz(3.0766727) q[3];
sx q[3];
rz(-1.4109979) q[3];
sx q[3];
rz(-1.1499278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0885334) q[0];
sx q[0];
rz(-2.2011338) q[0];
sx q[0];
rz(-0.9084107) q[0];
rz(-2.9216595) q[1];
sx q[1];
rz(-1.5426153) q[1];
sx q[1];
rz(-1.251108) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044949637) q[0];
sx q[0];
rz(-1.310955) q[0];
sx q[0];
rz(-1.0837062) q[0];
rz(-2.8562137) q[2];
sx q[2];
rz(-1.4754461) q[2];
sx q[2];
rz(-0.91172632) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0617705) q[1];
sx q[1];
rz(-2.1822189) q[1];
sx q[1];
rz(-2.3794214) q[1];
rz(-0.12507579) q[3];
sx q[3];
rz(-0.8440869) q[3];
sx q[3];
rz(-0.39966003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1088341) q[2];
sx q[2];
rz(-2.3075576) q[2];
sx q[2];
rz(-1.6081107) q[2];
rz(1.1533302) q[3];
sx q[3];
rz(-1.0554375) q[3];
sx q[3];
rz(-1.6495033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9913427) q[0];
sx q[0];
rz(-0.38337502) q[0];
sx q[0];
rz(-0.94069329) q[0];
rz(-3.0062145) q[1];
sx q[1];
rz(-1.4043413) q[1];
sx q[1];
rz(-2.1167596) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0198284) q[0];
sx q[0];
rz(-1.8786977) q[0];
sx q[0];
rz(0.93754365) q[0];
rz(-pi) q[1];
rz(1.631279) q[2];
sx q[2];
rz(-0.96394682) q[2];
sx q[2];
rz(-1.902193) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26014806) q[1];
sx q[1];
rz(-1.3760412) q[1];
sx q[1];
rz(1.4347726) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7303194) q[3];
sx q[3];
rz(-1.6182634) q[3];
sx q[3];
rz(-0.91225831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3465053) q[2];
sx q[2];
rz(-2.2515209) q[2];
sx q[2];
rz(-0.65004641) q[2];
rz(2.5551689) q[3];
sx q[3];
rz(-1.6610049) q[3];
sx q[3];
rz(0.75064269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2670249) q[0];
sx q[0];
rz(-2.6331007) q[0];
sx q[0];
rz(1.1453999) q[0];
rz(-1.8264495) q[1];
sx q[1];
rz(-1.9086842) q[1];
sx q[1];
rz(-2.4536536) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75632771) q[0];
sx q[0];
rz(-1.9419799) q[0];
sx q[0];
rz(0.20275499) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97369377) q[2];
sx q[2];
rz(-1.8988238) q[2];
sx q[2];
rz(-2.6231678) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8677788) q[1];
sx q[1];
rz(-1.4090589) q[1];
sx q[1];
rz(-0.75839569) q[1];
rz(-2.2086772) q[3];
sx q[3];
rz(-2.5543) q[3];
sx q[3];
rz(-1.4838532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.73100662) q[2];
sx q[2];
rz(-2.035391) q[2];
sx q[2];
rz(0.86622396) q[2];
rz(-0.90803641) q[3];
sx q[3];
rz(-2.3767545) q[3];
sx q[3];
rz(1.0959371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1753801) q[0];
sx q[0];
rz(-1.1840273) q[0];
sx q[0];
rz(-0.41561919) q[0];
rz(-2.3497154) q[1];
sx q[1];
rz(-1.2245347) q[1];
sx q[1];
rz(0.22722879) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5281866) q[0];
sx q[0];
rz(-0.019917073) q[0];
sx q[0];
rz(-2.5100757) q[0];
rz(-pi) q[1];
rz(-0.077543295) q[2];
sx q[2];
rz(-0.80813404) q[2];
sx q[2];
rz(-1.5002973) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1243678) q[1];
sx q[1];
rz(-1.827888) q[1];
sx q[1];
rz(-0.59558792) q[1];
rz(2.6260904) q[3];
sx q[3];
rz(-1.1519264) q[3];
sx q[3];
rz(0.91332289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8772584) q[2];
sx q[2];
rz(-2.4206968) q[2];
sx q[2];
rz(-2.7016675) q[2];
rz(0.59709966) q[3];
sx q[3];
rz(-2.4720981) q[3];
sx q[3];
rz(1.7841024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6482342) q[0];
sx q[0];
rz(-1.5879205) q[0];
sx q[0];
rz(1.2863202) q[0];
rz(-1.413912) q[1];
sx q[1];
rz(-1.021011) q[1];
sx q[1];
rz(0.11722142) q[1];
rz(3.0138409) q[2];
sx q[2];
rz(-0.41018872) q[2];
sx q[2];
rz(-1.142516) q[2];
rz(-2.8942378) q[3];
sx q[3];
rz(-0.83360278) q[3];
sx q[3];
rz(-3.0616888) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
