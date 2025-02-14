OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73622048) q[0];
sx q[0];
rz(3.1625746) q[0];
sx q[0];
rz(7.4641135) q[0];
rz(-0.48038545) q[1];
sx q[1];
rz(-1.9863702) q[1];
sx q[1];
rz(2.6748599) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48343998) q[0];
sx q[0];
rz(-0.83388336) q[0];
sx q[0];
rz(1.1841274) q[0];
rz(0.67961116) q[2];
sx q[2];
rz(-0.70290297) q[2];
sx q[2];
rz(1.694569) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8267104) q[1];
sx q[1];
rz(-0.2383543) q[1];
sx q[1];
rz(-1.5520067) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0084413) q[3];
sx q[3];
rz(-0.97259023) q[3];
sx q[3];
rz(1.7266718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0878318) q[2];
sx q[2];
rz(-1.6873282) q[2];
sx q[2];
rz(1.9161179) q[2];
rz(-2.7671704) q[3];
sx q[3];
rz(-2.008308) q[3];
sx q[3];
rz(0.55330127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9541009) q[0];
sx q[0];
rz(-2.2611389) q[0];
sx q[0];
rz(-0.53400293) q[0];
rz(-2.7502637) q[1];
sx q[1];
rz(-0.44328872) q[1];
sx q[1];
rz(-2.2543529) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0011883) q[0];
sx q[0];
rz(-0.47022223) q[0];
sx q[0];
rz(-0.73064877) q[0];
rz(-2.5975304) q[2];
sx q[2];
rz(-2.3070222) q[2];
sx q[2];
rz(-0.53235934) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56480592) q[1];
sx q[1];
rz(-0.24165711) q[1];
sx q[1];
rz(-2.5741626) q[1];
x q[2];
rz(2.128936) q[3];
sx q[3];
rz(-0.59132708) q[3];
sx q[3];
rz(1.8048353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5947764) q[2];
sx q[2];
rz(-1.8211326) q[2];
sx q[2];
rz(0.36396626) q[2];
rz(-1.4173896) q[3];
sx q[3];
rz(-2.851749) q[3];
sx q[3];
rz(-0.83436203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2367547) q[0];
sx q[0];
rz(-3.0433488) q[0];
sx q[0];
rz(-2.8114317) q[0];
rz(-2.5579021) q[1];
sx q[1];
rz(-2.3287562) q[1];
sx q[1];
rz(0.49461734) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1903822) q[0];
sx q[0];
rz(-1.8415762) q[0];
sx q[0];
rz(-0.27373224) q[0];
x q[1];
rz(-2.0614659) q[2];
sx q[2];
rz(-1.3441836) q[2];
sx q[2];
rz(-1.8267711) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3266625) q[1];
sx q[1];
rz(-2.0771793) q[1];
sx q[1];
rz(2.6145007) q[1];
rz(-0.98683896) q[3];
sx q[3];
rz(-1.6641326) q[3];
sx q[3];
rz(3.0350041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.79564774) q[2];
sx q[2];
rz(-2.4730885) q[2];
sx q[2];
rz(0.12269679) q[2];
rz(-3.0986541) q[3];
sx q[3];
rz(-2.0624845) q[3];
sx q[3];
rz(0.19449657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4075459) q[0];
sx q[0];
rz(-0.50739822) q[0];
sx q[0];
rz(-0.88974446) q[0];
rz(0.45097688) q[1];
sx q[1];
rz(-0.86800066) q[1];
sx q[1];
rz(-2.6873592) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20393237) q[0];
sx q[0];
rz(-1.9113921) q[0];
sx q[0];
rz(-1.6837337) q[0];
rz(-pi) q[1];
rz(-1.5726015) q[2];
sx q[2];
rz(-1.1434525) q[2];
sx q[2];
rz(-1.7693365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0955086) q[1];
sx q[1];
rz(-2.3254776) q[1];
sx q[1];
rz(-2.7296394) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6481712) q[3];
sx q[3];
rz(-1.2168125) q[3];
sx q[3];
rz(0.0019687584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98895994) q[2];
sx q[2];
rz(-2.5071414) q[2];
sx q[2];
rz(2.3839942) q[2];
rz(2.9518413) q[3];
sx q[3];
rz(-1.3187783) q[3];
sx q[3];
rz(-0.54722133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2106638) q[0];
sx q[0];
rz(-0.78721109) q[0];
sx q[0];
rz(1.2934562) q[0];
rz(2.5594607) q[1];
sx q[1];
rz(-1.7030092) q[1];
sx q[1];
rz(2.9621946) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85953259) q[0];
sx q[0];
rz(-2.8683887) q[0];
sx q[0];
rz(0.84635107) q[0];
rz(-2.8371528) q[2];
sx q[2];
rz(-2.4477262) q[2];
sx q[2];
rz(2.2562698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7475954) q[1];
sx q[1];
rz(-2.4863805) q[1];
sx q[1];
rz(1.2607695) q[1];
x q[2];
rz(0.50197451) q[3];
sx q[3];
rz(-2.2036512) q[3];
sx q[3];
rz(-0.33251921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3680215) q[2];
sx q[2];
rz(-0.37006912) q[2];
sx q[2];
rz(2.5717112) q[2];
rz(-0.44024769) q[3];
sx q[3];
rz(-1.2443685) q[3];
sx q[3];
rz(-2.7888035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2920947) q[0];
sx q[0];
rz(-2.8773913) q[0];
sx q[0];
rz(1.1415035) q[0];
rz(1.3453206) q[1];
sx q[1];
rz(-2.1514386) q[1];
sx q[1];
rz(-2.9703531) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9123403) q[0];
sx q[0];
rz(-1.0919337) q[0];
sx q[0];
rz(-1.6845869) q[0];
rz(-pi) q[1];
rz(2.0948579) q[2];
sx q[2];
rz(-1.1772794) q[2];
sx q[2];
rz(-1.6167906) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.332883) q[1];
sx q[1];
rz(-0.58296219) q[1];
sx q[1];
rz(-2.8364532) q[1];
rz(-1.0140795) q[3];
sx q[3];
rz(-1.5601592) q[3];
sx q[3];
rz(-0.47350804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3693927) q[2];
sx q[2];
rz(-0.93595305) q[2];
sx q[2];
rz(-0.13747036) q[2];
rz(-1.2482268) q[3];
sx q[3];
rz(-0.56709254) q[3];
sx q[3];
rz(1.9198157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0565979) q[0];
sx q[0];
rz(-2.5039112) q[0];
sx q[0];
rz(-1.6531264) q[0];
rz(2.3562863) q[1];
sx q[1];
rz(-1.4005125) q[1];
sx q[1];
rz(-1.8280425) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0037601622) q[0];
sx q[0];
rz(-1.6562471) q[0];
sx q[0];
rz(-1.6944044) q[0];
rz(0.19124548) q[2];
sx q[2];
rz(-1.4594263) q[2];
sx q[2];
rz(-2.0759866) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.056712778) q[1];
sx q[1];
rz(-3.0482349) q[1];
sx q[1];
rz(0.30594029) q[1];
rz(-pi) q[2];
x q[2];
rz(0.051387473) q[3];
sx q[3];
rz(-0.87552658) q[3];
sx q[3];
rz(-1.072708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.58275783) q[2];
sx q[2];
rz(-2.1119327) q[2];
sx q[2];
rz(-2.6666857) q[2];
rz(2.9389935) q[3];
sx q[3];
rz(-0.83297268) q[3];
sx q[3];
rz(1.1420265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9433893) q[0];
sx q[0];
rz(-0.13309637) q[0];
sx q[0];
rz(0.30580795) q[0];
rz(0.43931475) q[1];
sx q[1];
rz(-0.82249928) q[1];
sx q[1];
rz(-1.3154715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1179349) q[0];
sx q[0];
rz(-0.88738897) q[0];
sx q[0];
rz(-1.5853496) q[0];
rz(0.28355174) q[2];
sx q[2];
rz(-2.1075948) q[2];
sx q[2];
rz(-0.62284046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2478359) q[1];
sx q[1];
rz(-2.3512321) q[1];
sx q[1];
rz(-1.7414581) q[1];
x q[2];
rz(-0.9813811) q[3];
sx q[3];
rz(-1.9187315) q[3];
sx q[3];
rz(0.98435005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4149912) q[2];
sx q[2];
rz(-1.5982268) q[2];
sx q[2];
rz(-0.42465633) q[2];
rz(-0.031938227) q[3];
sx q[3];
rz(-1.6274803) q[3];
sx q[3];
rz(2.4922075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6040819) q[0];
sx q[0];
rz(-0.6518971) q[0];
sx q[0];
rz(-2.5867468) q[0];
rz(0.26508731) q[1];
sx q[1];
rz(-1.6731429) q[1];
sx q[1];
rz(1.2010942) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6938615) q[0];
sx q[0];
rz(-2.434512) q[0];
sx q[0];
rz(-2.5140544) q[0];
rz(0.55055203) q[2];
sx q[2];
rz(-0.32062396) q[2];
sx q[2];
rz(-2.6132646) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20165138) q[1];
sx q[1];
rz(-1.2229511) q[1];
sx q[1];
rz(-1.6099206) q[1];
x q[2];
rz(1.5952121) q[3];
sx q[3];
rz(-2.3204969) q[3];
sx q[3];
rz(-1.2670307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25192866) q[2];
sx q[2];
rz(-2.1519075) q[2];
sx q[2];
rz(0.84452334) q[2];
rz(2.0097513) q[3];
sx q[3];
rz(-2.2364538) q[3];
sx q[3];
rz(1.7822781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97245401) q[0];
sx q[0];
rz(-0.3903946) q[0];
sx q[0];
rz(1.1727232) q[0];
rz(0.68924618) q[1];
sx q[1];
rz(-2.0447562) q[1];
sx q[1];
rz(0.26395878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64484064) q[0];
sx q[0];
rz(-0.94366108) q[0];
sx q[0];
rz(-2.1261313) q[0];
rz(-0.9969292) q[2];
sx q[2];
rz(-0.80748122) q[2];
sx q[2];
rz(-1.6107744) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8431153) q[1];
sx q[1];
rz(-2.7727371) q[1];
sx q[1];
rz(-2.1082676) q[1];
x q[2];
rz(-0.12216025) q[3];
sx q[3];
rz(-2.2111243) q[3];
sx q[3];
rz(1.8807008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76137233) q[2];
sx q[2];
rz(-1.1470046) q[2];
sx q[2];
rz(-2.0210361) q[2];
rz(-1.606696) q[3];
sx q[3];
rz(-2.1963019) q[3];
sx q[3];
rz(1.631261) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0705538) q[0];
sx q[0];
rz(-2.4865535) q[0];
sx q[0];
rz(3.0336663) q[0];
rz(0.72036605) q[1];
sx q[1];
rz(-1.9279059) q[1];
sx q[1];
rz(2.7005213) q[1];
rz(2.2137523) q[2];
sx q[2];
rz(-0.8561047) q[2];
sx q[2];
rz(2.718953) q[2];
rz(1.4955487) q[3];
sx q[3];
rz(-0.9552707) q[3];
sx q[3];
rz(0.09494119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
