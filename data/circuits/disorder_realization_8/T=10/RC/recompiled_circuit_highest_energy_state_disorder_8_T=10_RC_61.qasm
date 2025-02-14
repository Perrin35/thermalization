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
rz(-0.8166135) q[0];
sx q[0];
rz(-0.54083523) q[0];
sx q[0];
rz(1.1582561) q[0];
rz(6.3732014) q[1];
sx q[1];
rz(6.7572588) q[1];
sx q[1];
rz(13.401539) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0220563) q[0];
sx q[0];
rz(-1.2588663) q[0];
sx q[0];
rz(2.3989912) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64106854) q[2];
sx q[2];
rz(-1.8279229) q[2];
sx q[2];
rz(-2.7066305) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1013284) q[1];
sx q[1];
rz(-0.44875408) q[1];
sx q[1];
rz(0.59228102) q[1];
x q[2];
rz(0.99906875) q[3];
sx q[3];
rz(-0.99410996) q[3];
sx q[3];
rz(-0.028923361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0917255) q[2];
sx q[2];
rz(-0.89605248) q[2];
sx q[2];
rz(-2.8095165) q[2];
rz(2.8927228) q[3];
sx q[3];
rz(-1.9460461) q[3];
sx q[3];
rz(0.88768774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
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
rz(-2.2700972) q[0];
sx q[0];
rz(-0.29093727) q[0];
sx q[0];
rz(2.6089597) q[0];
rz(0.10781413) q[1];
sx q[1];
rz(-1.1153406) q[1];
sx q[1];
rz(-0.32049387) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5607213) q[0];
sx q[0];
rz(-0.91976368) q[0];
sx q[0];
rz(3.0614775) q[0];
rz(1.3391431) q[2];
sx q[2];
rz(-1.4024874) q[2];
sx q[2];
rz(1.2141808) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1793666) q[1];
sx q[1];
rz(-1.4186064) q[1];
sx q[1];
rz(2.9801912) q[1];
rz(-pi) q[2];
rz(-1.3555834) q[3];
sx q[3];
rz(-1.2485412) q[3];
sx q[3];
rz(-2.2029869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33114854) q[2];
sx q[2];
rz(-0.30511567) q[2];
sx q[2];
rz(1.5020874) q[2];
rz(0.33809996) q[3];
sx q[3];
rz(-0.88518849) q[3];
sx q[3];
rz(2.6147208) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51760393) q[0];
sx q[0];
rz(-1.9094587) q[0];
sx q[0];
rz(2.7841618) q[0];
rz(2.7492211) q[1];
sx q[1];
rz(-0.79634276) q[1];
sx q[1];
rz(-1.2145112) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93335184) q[0];
sx q[0];
rz(-0.14259556) q[0];
sx q[0];
rz(1.6928133) q[0];
rz(-pi) q[1];
rz(-2.8896595) q[2];
sx q[2];
rz(-2.5802543) q[2];
sx q[2];
rz(0.47719819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4068661) q[1];
sx q[1];
rz(-2.7541783) q[1];
sx q[1];
rz(1.563036) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9830832) q[3];
sx q[3];
rz(-1.5395035) q[3];
sx q[3];
rz(-2.8060117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7963205) q[2];
sx q[2];
rz(-2.1638963) q[2];
sx q[2];
rz(0.39682445) q[2];
rz(-1.2567629) q[3];
sx q[3];
rz(-1.4182914) q[3];
sx q[3];
rz(2.5313012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20763718) q[0];
sx q[0];
rz(-1.3960681) q[0];
sx q[0];
rz(-1.9130094) q[0];
rz(-1.2127016) q[1];
sx q[1];
rz(-2.1804501) q[1];
sx q[1];
rz(-0.62087762) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8437794) q[0];
sx q[0];
rz(-0.50619805) q[0];
sx q[0];
rz(1.9202581) q[0];
rz(-0.16571705) q[2];
sx q[2];
rz(-1.724337) q[2];
sx q[2];
rz(2.4927947) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.32737) q[1];
sx q[1];
rz(-1.3124221) q[1];
sx q[1];
rz(-1.3034921) q[1];
rz(1.4520763) q[3];
sx q[3];
rz(-1.9510036) q[3];
sx q[3];
rz(-0.79352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88317251) q[2];
sx q[2];
rz(-1.1383388) q[2];
sx q[2];
rz(2.9849198) q[2];
rz(-0.91642085) q[3];
sx q[3];
rz(-2.1082924) q[3];
sx q[3];
rz(-0.55287439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857392) q[0];
sx q[0];
rz(-1.1612949) q[0];
sx q[0];
rz(-0.39749417) q[0];
rz(3.1187348) q[1];
sx q[1];
rz(-2.6409179) q[1];
sx q[1];
rz(-0.674725) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8431339) q[0];
sx q[0];
rz(-1.805691) q[0];
sx q[0];
rz(0.21505925) q[0];
rz(-1.6653453) q[2];
sx q[2];
rz(-1.4473905) q[2];
sx q[2];
rz(-2.9837554) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6276363) q[1];
sx q[1];
rz(-1.7412724) q[1];
sx q[1];
rz(1.4161311) q[1];
x q[2];
rz(0.85002331) q[3];
sx q[3];
rz(-2.7394419) q[3];
sx q[3];
rz(-1.7274203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.61458331) q[2];
sx q[2];
rz(-0.60636568) q[2];
sx q[2];
rz(-0.69925365) q[2];
rz(-1.3462542) q[3];
sx q[3];
rz(-2.5739539) q[3];
sx q[3];
rz(-2.4301372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.72494495) q[0];
sx q[0];
rz(-2.7236433) q[0];
sx q[0];
rz(-2.7432192) q[0];
rz(-0.97459546) q[1];
sx q[1];
rz(-1.3037126) q[1];
sx q[1];
rz(-1.7030565) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8038626) q[0];
sx q[0];
rz(-2.8210905) q[0];
sx q[0];
rz(-0.33945531) q[0];
x q[1];
rz(-0.010377093) q[2];
sx q[2];
rz(-1.3149909) q[2];
sx q[2];
rz(2.1880045) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7185091) q[1];
sx q[1];
rz(-0.72186493) q[1];
sx q[1];
rz(-1.2150498) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8277728) q[3];
sx q[3];
rz(-2.0517618) q[3];
sx q[3];
rz(-2.0723267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.67821104) q[2];
sx q[2];
rz(-1.3621829) q[2];
sx q[2];
rz(1.3055275) q[2];
rz(-1.8259004) q[3];
sx q[3];
rz(-1.7782327) q[3];
sx q[3];
rz(2.2731884) q[3];
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
rz(1.0618133) q[0];
sx q[0];
rz(-0.4158622) q[0];
sx q[0];
rz(1.6449991) q[0];
rz(-0.98622259) q[1];
sx q[1];
rz(-1.5602427) q[1];
sx q[1];
rz(0.49759069) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4301028) q[0];
sx q[0];
rz(-0.27207366) q[0];
sx q[0];
rz(-1.2077232) q[0];
rz(2.6025101) q[2];
sx q[2];
rz(-1.0175127) q[2];
sx q[2];
rz(0.48083559) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6141245) q[1];
sx q[1];
rz(-1.9640199) q[1];
sx q[1];
rz(2.9725463) q[1];
rz(1.6646181) q[3];
sx q[3];
rz(-1.2872496) q[3];
sx q[3];
rz(-0.16167262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8387973) q[2];
sx q[2];
rz(-1.5747728) q[2];
sx q[2];
rz(-0.29067579) q[2];
rz(-2.931328) q[3];
sx q[3];
rz(-2.1472774) q[3];
sx q[3];
rz(1.0342342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5656089) q[0];
sx q[0];
rz(-1.9714332) q[0];
sx q[0];
rz(-1.959311) q[0];
rz(-0.15444175) q[1];
sx q[1];
rz(-1.4192105) q[1];
sx q[1];
rz(-0.90528893) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2700896) q[0];
sx q[0];
rz(-0.2729899) q[0];
sx q[0];
rz(1.5259302) q[0];
x q[1];
rz(-1.8087093) q[2];
sx q[2];
rz(-2.015997) q[2];
sx q[2];
rz(-0.92724909) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33069995) q[1];
sx q[1];
rz(-2.3726683) q[1];
sx q[1];
rz(-2.9309209) q[1];
rz(-pi) q[2];
rz(-2.1240299) q[3];
sx q[3];
rz(-0.39504566) q[3];
sx q[3];
rz(-2.1086804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0514544) q[2];
sx q[2];
rz(-1.5608414) q[2];
sx q[2];
rz(-0.95019379) q[2];
rz(-0.072362445) q[3];
sx q[3];
rz(-1.3772929) q[3];
sx q[3];
rz(-2.4322521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224834) q[0];
sx q[0];
rz(-0.95471946) q[0];
sx q[0];
rz(2.0461653) q[0];
rz(2.0504045) q[1];
sx q[1];
rz(-0.90091101) q[1];
sx q[1];
rz(0.99517623) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93928775) q[0];
sx q[0];
rz(-1.6350593) q[0];
sx q[0];
rz(0.66575428) q[0];
x q[1];
rz(0.27492304) q[2];
sx q[2];
rz(-1.2464393) q[2];
sx q[2];
rz(2.7664326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.51982626) q[1];
sx q[1];
rz(-1.6556181) q[1];
sx q[1];
rz(1.5250526) q[1];
x q[2];
rz(-1.3714637) q[3];
sx q[3];
rz(-1.3889894) q[3];
sx q[3];
rz(-1.1246641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6301849) q[2];
sx q[2];
rz(-0.46773043) q[2];
sx q[2];
rz(0.67997813) q[2];
rz(1.4558815) q[3];
sx q[3];
rz(-0.89434353) q[3];
sx q[3];
rz(-1.9623914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7162914) q[0];
sx q[0];
rz(-2.20708) q[0];
sx q[0];
rz(-0.61538482) q[0];
rz(1.3353434) q[1];
sx q[1];
rz(-1.6067303) q[1];
sx q[1];
rz(-1.7937484) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6123209) q[0];
sx q[0];
rz(-2.7833496) q[0];
sx q[0];
rz(0.1310346) q[0];
rz(-pi) q[1];
rz(0.81735264) q[2];
sx q[2];
rz(-2.0869227) q[2];
sx q[2];
rz(-1.1990666) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.19619689) q[1];
sx q[1];
rz(-2.7295503) q[1];
sx q[1];
rz(-2.7315188) q[1];
rz(-1.6928635) q[3];
sx q[3];
rz(-2.3231835) q[3];
sx q[3];
rz(-2.9451452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24796692) q[2];
sx q[2];
rz(-1.3004356) q[2];
sx q[2];
rz(2.628053) q[2];
rz(-2.9945471) q[3];
sx q[3];
rz(-2.5976318) q[3];
sx q[3];
rz(1.8769544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2697987) q[0];
sx q[0];
rz(-2.2525621) q[0];
sx q[0];
rz(2.9004108) q[0];
rz(-1.3006032) q[1];
sx q[1];
rz(-1.069297) q[1];
sx q[1];
rz(-0.020513608) q[1];
rz(-0.14221556) q[2];
sx q[2];
rz(-1.1995865) q[2];
sx q[2];
rz(-2.8778278) q[2];
rz(0.62773283) q[3];
sx q[3];
rz(-1.2017602) q[3];
sx q[3];
rz(-1.1084262) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
