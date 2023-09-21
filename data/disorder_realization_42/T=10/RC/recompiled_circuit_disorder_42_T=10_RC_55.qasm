OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5126098) q[0];
sx q[0];
rz(-1.0273758) q[0];
sx q[0];
rz(-2.7789814) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(-0.4904823) q[1];
sx q[1];
rz(-0.34163707) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1385961) q[0];
sx q[0];
rz(-1.9435104) q[0];
sx q[0];
rz(1.7025823) q[0];
rz(2.6356959) q[2];
sx q[2];
rz(-2.4873173) q[2];
sx q[2];
rz(-1.3327395) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4075549) q[1];
sx q[1];
rz(-1.2734379) q[1];
sx q[1];
rz(2.107891) q[1];
rz(-2.9176641) q[3];
sx q[3];
rz(-1.180598) q[3];
sx q[3];
rz(-1.3003365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.59387702) q[2];
sx q[2];
rz(-1.9888069) q[2];
sx q[2];
rz(-0.95648742) q[2];
rz(2.9603413) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(-1.6199002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-3.0759838) q[0];
sx q[0];
rz(1.99615) q[0];
rz(1.068813) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(-2.4157445) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84859914) q[0];
sx q[0];
rz(-1.2147875) q[0];
sx q[0];
rz(2.5550935) q[0];
rz(2.0306573) q[2];
sx q[2];
rz(-2.1150555) q[2];
sx q[2];
rz(2.0366675) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1314288) q[1];
sx q[1];
rz(-1.5958852) q[1];
sx q[1];
rz(-1.2618834) q[1];
rz(2.8511091) q[3];
sx q[3];
rz(-0.64290128) q[3];
sx q[3];
rz(3.052352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(0.99728161) q[2];
rz(1.7287792) q[3];
sx q[3];
rz(-1.0935067) q[3];
sx q[3];
rz(-0.93878448) q[3];
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
rz(-0.41855758) q[0];
sx q[0];
rz(-2.1756873) q[0];
sx q[0];
rz(-1.1285271) q[0];
rz(1.3035125) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(-0.9054786) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83429503) q[0];
sx q[0];
rz(-0.30954888) q[0];
sx q[0];
rz(2.4718667) q[0];
x q[1];
rz(2.3891874) q[2];
sx q[2];
rz(-1.2453326) q[2];
sx q[2];
rz(2.4706555) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9943774) q[1];
sx q[1];
rz(-0.71502393) q[1];
sx q[1];
rz(-1.6836402) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6204897) q[3];
sx q[3];
rz(-2.4274785) q[3];
sx q[3];
rz(-0.36220887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3811615) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(-2.1543489) q[2];
rz(0.045771249) q[3];
sx q[3];
rz(-1.8493303) q[3];
sx q[3];
rz(-2.9614017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0320597) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(-3.0155244) q[0];
rz(-0.18445045) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(-1.9741845) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4520893) q[0];
sx q[0];
rz(-0.79012442) q[0];
sx q[0];
rz(0.28496123) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6378176) q[2];
sx q[2];
rz(-1.4132032) q[2];
sx q[2];
rz(0.7047082) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.57707225) q[1];
sx q[1];
rz(-2.2153691) q[1];
sx q[1];
rz(-1.2381899) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4950124) q[3];
sx q[3];
rz(-1.0160035) q[3];
sx q[3];
rz(-0.085689714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2771153) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(-0.63151044) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(-0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1642078) q[0];
sx q[0];
rz(-3.003484) q[0];
sx q[0];
rz(-0.50267977) q[0];
rz(2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(0.46708333) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53411667) q[0];
sx q[0];
rz(-2.7152938) q[0];
sx q[0];
rz(2.0712453) q[0];
rz(-pi) q[1];
rz(-2.5555243) q[2];
sx q[2];
rz(-1.7826826) q[2];
sx q[2];
rz(0.52473247) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5284164) q[1];
sx q[1];
rz(-1.3532552) q[1];
sx q[1];
rz(1.1385285) q[1];
x q[2];
rz(-2.2561982) q[3];
sx q[3];
rz(-1.6939031) q[3];
sx q[3];
rz(0.76243329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46665329) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(0.49218991) q[2];
rz(-2.8594033) q[3];
sx q[3];
rz(-1.4331093) q[3];
sx q[3];
rz(-1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.870134) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(-3.0555994) q[0];
rz(1.7135235) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(-1.6451947) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39280415) q[0];
sx q[0];
rz(-1.7668056) q[0];
sx q[0];
rz(-1.0402354) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1804579) q[2];
sx q[2];
rz(-2.5679553) q[2];
sx q[2];
rz(0.2804642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9073346) q[1];
sx q[1];
rz(-1.9252216) q[1];
sx q[1];
rz(0.025269421) q[1];
rz(-pi) q[2];
rz(0.71293998) q[3];
sx q[3];
rz(-1.561165) q[3];
sx q[3];
rz(-0.18129098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4795586) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(-2.7453444) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(-1.5007639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2632161) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(-2.0492045) q[0];
rz(-1.3573525) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(-0.011627442) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8788293) q[0];
sx q[0];
rz(-1.3063523) q[0];
sx q[0];
rz(0.39475616) q[0];
rz(-pi) q[1];
rz(0.83519148) q[2];
sx q[2];
rz(-3.1036658) q[2];
sx q[2];
rz(-1.4284301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7203622) q[1];
sx q[1];
rz(-1.0248529) q[1];
sx q[1];
rz(-0.027688428) q[1];
rz(-pi) q[2];
rz(-1.8947951) q[3];
sx q[3];
rz(-2.218459) q[3];
sx q[3];
rz(0.049829359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.45550436) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(2.8536076) q[2];
rz(-1.1503495) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(2.5434125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312254) q[0];
sx q[0];
rz(-0.53124017) q[0];
sx q[0];
rz(-1.9294552) q[0];
rz(-0.37777004) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(-1.6479962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15570116) q[0];
sx q[0];
rz(-2.771286) q[0];
sx q[0];
rz(2.5168688) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4385857) q[2];
sx q[2];
rz(-1.6301179) q[2];
sx q[2];
rz(1.6249663) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89477506) q[1];
sx q[1];
rz(-0.90118876) q[1];
sx q[1];
rz(-1.0424022) q[1];
x q[2];
rz(-0.11843189) q[3];
sx q[3];
rz(-1.0716972) q[3];
sx q[3];
rz(-2.0407608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7887855) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(1.8898233) q[2];
rz(-2.2552323) q[3];
sx q[3];
rz(-0.76615196) q[3];
sx q[3];
rz(-0.10929508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66824526) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(-0.1782724) q[0];
rz(0.072862236) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(-1.1791621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060886325) q[0];
sx q[0];
rz(-2.076425) q[0];
sx q[0];
rz(-1.2348742) q[0];
rz(-pi) q[1];
rz(-0.945325) q[2];
sx q[2];
rz(-1.9918459) q[2];
sx q[2];
rz(2.9025214) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33680962) q[1];
sx q[1];
rz(-2.0840624) q[1];
sx q[1];
rz(-1.8670765) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7618802) q[3];
sx q[3];
rz(-2.0413997) q[3];
sx q[3];
rz(-3.0203366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.69958413) q[2];
sx q[2];
rz(-2.1506491) q[2];
sx q[2];
rz(1.0317831) q[2];
rz(1.7992841) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(-2.9163196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853302) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(-0.061696079) q[0];
rz(-1.306698) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(-2.0526989) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13578829) q[0];
sx q[0];
rz(-3.0514768) q[0];
sx q[0];
rz(0.80149217) q[0];
rz(-2.225869) q[2];
sx q[2];
rz(-1.7172126) q[2];
sx q[2];
rz(0.54619782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1354462) q[1];
sx q[1];
rz(-1.0137614) q[1];
sx q[1];
rz(2.1250493) q[1];
rz(-pi) q[2];
rz(-1.2286931) q[3];
sx q[3];
rz(-1.7497239) q[3];
sx q[3];
rz(0.10781328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5237727) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(3.0715959) q[2];
rz(-2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(2.783412) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(-1.6246673) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(2.2286227) q[2];
sx q[2];
rz(-1.2428478) q[2];
sx q[2];
rz(-1.7075677) q[2];
rz(0.041257507) q[3];
sx q[3];
rz(-1.2357124) q[3];
sx q[3];
rz(0.99832051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];