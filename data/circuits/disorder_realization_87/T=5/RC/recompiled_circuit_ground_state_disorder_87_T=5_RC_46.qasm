OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4053722) q[0];
sx q[0];
rz(-0.020981941) q[0];
sx q[0];
rz(1.9606645) q[0];
rz(-0.48038545) q[1];
sx q[1];
rz(-1.9863702) q[1];
sx q[1];
rz(-0.46673271) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8202756) q[0];
sx q[0];
rz(-1.853824) q[0];
sx q[0];
rz(0.77518605) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67961116) q[2];
sx q[2];
rz(-2.4386897) q[2];
sx q[2];
rz(1.4470237) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8267104) q[1];
sx q[1];
rz(-2.9032384) q[1];
sx q[1];
rz(-1.5895859) q[1];
rz(-pi) q[2];
rz(0.67812293) q[3];
sx q[3];
rz(-1.1145385) q[3];
sx q[3];
rz(0.18517906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0537609) q[2];
sx q[2];
rz(-1.6873282) q[2];
sx q[2];
rz(-1.9161179) q[2];
rz(2.7671704) q[3];
sx q[3];
rz(-2.008308) q[3];
sx q[3];
rz(-0.55330127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.18749172) q[0];
sx q[0];
rz(-0.88045374) q[0];
sx q[0];
rz(0.53400293) q[0];
rz(-0.39132896) q[1];
sx q[1];
rz(-2.6983039) q[1];
sx q[1];
rz(0.88723976) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0011883) q[0];
sx q[0];
rz(-2.6713704) q[0];
sx q[0];
rz(-0.73064877) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3848955) q[2];
sx q[2];
rz(-1.9644418) q[2];
sx q[2];
rz(-1.7172888) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6897517) q[1];
sx q[1];
rz(-1.4418169) q[1];
sx q[1];
rz(2.9366628) q[1];
x q[2];
rz(-0.34167413) q[3];
sx q[3];
rz(-1.0782584) q[3];
sx q[3];
rz(1.1600174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54681626) q[2];
sx q[2];
rz(-1.3204601) q[2];
sx q[2];
rz(2.7776264) q[2];
rz(1.724203) q[3];
sx q[3];
rz(-2.851749) q[3];
sx q[3];
rz(-0.83436203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.904838) q[0];
sx q[0];
rz(-3.0433488) q[0];
sx q[0];
rz(-0.33016095) q[0];
rz(2.5579021) q[1];
sx q[1];
rz(-2.3287562) q[1];
sx q[1];
rz(-0.49461734) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8361417) q[0];
sx q[0];
rz(-1.8343121) q[0];
sx q[0];
rz(1.8515153) q[0];
x q[1];
rz(-0.25569368) q[2];
sx q[2];
rz(-1.0937368) q[2];
sx q[2];
rz(0.13651234) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51926326) q[1];
sx q[1];
rz(-1.1153145) q[1];
sx q[1];
rz(1.0002656) q[1];
rz(-pi) q[2];
rz(1.7389789) q[3];
sx q[3];
rz(-2.5510812) q[3];
sx q[3];
rz(1.3240977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3459449) q[2];
sx q[2];
rz(-0.66850418) q[2];
sx q[2];
rz(0.12269679) q[2];
rz(-3.0986541) q[3];
sx q[3];
rz(-1.0791082) q[3];
sx q[3];
rz(2.9470961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73404679) q[0];
sx q[0];
rz(-2.6341944) q[0];
sx q[0];
rz(2.2518482) q[0];
rz(-0.45097688) q[1];
sx q[1];
rz(-0.86800066) q[1];
sx q[1];
rz(-0.45423347) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20393237) q[0];
sx q[0];
rz(-1.2302006) q[0];
sx q[0];
rz(1.457859) q[0];
rz(-pi) q[1];
rz(0.003963917) q[2];
sx q[2];
rz(-0.42734738) q[2];
sx q[2];
rz(1.3766118) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.046084) q[1];
sx q[1];
rz(-2.3254776) q[1];
sx q[1];
rz(-0.41195324) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4934214) q[3];
sx q[3];
rz(-1.2168125) q[3];
sx q[3];
rz(3.1396239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1526327) q[2];
sx q[2];
rz(-2.5071414) q[2];
sx q[2];
rz(0.75759849) q[2];
rz(-2.9518413) q[3];
sx q[3];
rz(-1.8228143) q[3];
sx q[3];
rz(2.5943713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2106638) q[0];
sx q[0];
rz(-2.3543816) q[0];
sx q[0];
rz(-1.8481365) q[0];
rz(2.5594607) q[1];
sx q[1];
rz(-1.4385834) q[1];
sx q[1];
rz(-2.9621946) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5388193) q[0];
sx q[0];
rz(-1.3673377) q[0];
sx q[0];
rz(0.18360965) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67082884) q[2];
sx q[2];
rz(-1.7636904) q[2];
sx q[2];
rz(-2.2191032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3939973) q[1];
sx q[1];
rz(-2.4863805) q[1];
sx q[1];
rz(1.2607695) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1513624) q[3];
sx q[3];
rz(-0.78573686) q[3];
sx q[3];
rz(-0.41555575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.77357117) q[2];
sx q[2];
rz(-0.37006912) q[2];
sx q[2];
rz(0.5698815) q[2];
rz(0.44024769) q[3];
sx q[3];
rz(-1.2443685) q[3];
sx q[3];
rz(2.7888035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84949795) q[0];
sx q[0];
rz(-2.8773913) q[0];
sx q[0];
rz(1.1415035) q[0];
rz(1.796272) q[1];
sx q[1];
rz(-2.1514386) q[1];
sx q[1];
rz(2.9703531) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2889338) q[0];
sx q[0];
rz(-1.6717413) q[0];
sx q[0];
rz(-2.6600718) q[0];
rz(-pi) q[1];
rz(2.6944567) q[2];
sx q[2];
rz(-2.0512181) q[2];
sx q[2];
rz(0.26408476) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6936612) q[1];
sx q[1];
rz(-1.0180001) q[1];
sx q[1];
rz(-1.7663683) q[1];
rz(-pi) q[2];
x q[2];
rz(0.01252894) q[3];
sx q[3];
rz(-2.1274779) q[3];
sx q[3];
rz(-2.0509246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3693927) q[2];
sx q[2];
rz(-2.2056396) q[2];
sx q[2];
rz(-3.0041223) q[2];
rz(1.8933659) q[3];
sx q[3];
rz(-2.5745001) q[3];
sx q[3];
rz(-1.9198157) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084994706) q[0];
sx q[0];
rz(-2.5039112) q[0];
sx q[0];
rz(1.4884663) q[0];
rz(0.78530637) q[1];
sx q[1];
rz(-1.7410802) q[1];
sx q[1];
rz(1.3135501) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5851597) q[0];
sx q[0];
rz(-1.4476416) q[0];
sx q[0];
rz(-0.086104546) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4573757) q[2];
sx q[2];
rz(-1.7608425) q[2];
sx q[2];
rz(0.52670497) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.32282) q[1];
sx q[1];
rz(-1.5427151) q[1];
sx q[1];
rz(3.0525467) q[1];
rz(-pi) q[2];
rz(-0.87487674) q[3];
sx q[3];
rz(-1.6102487) q[3];
sx q[3];
rz(-2.6764398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5588348) q[2];
sx q[2];
rz(-2.1119327) q[2];
sx q[2];
rz(-2.6666857) q[2];
rz(-0.20259914) q[3];
sx q[3];
rz(-2.30862) q[3];
sx q[3];
rz(-1.1420265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1982034) q[0];
sx q[0];
rz(-3.0084963) q[0];
sx q[0];
rz(0.30580795) q[0];
rz(-2.7022779) q[1];
sx q[1];
rz(-0.82249928) q[1];
sx q[1];
rz(1.8261212) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1179349) q[0];
sx q[0];
rz(-2.2542037) q[0];
sx q[0];
rz(-1.5562431) q[0];
rz(2.1256746) q[2];
sx q[2];
rz(-1.8135929) q[2];
sx q[2];
rz(-1.0958874) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.89375673) q[1];
sx q[1];
rz(-2.3512321) q[1];
sx q[1];
rz(-1.4001346) q[1];
rz(-pi) q[2];
rz(0.41141182) q[3];
sx q[3];
rz(-1.0209382) q[3];
sx q[3];
rz(-0.362277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4149912) q[2];
sx q[2];
rz(-1.5982268) q[2];
sx q[2];
rz(-0.42465633) q[2];
rz(3.1096544) q[3];
sx q[3];
rz(-1.5141124) q[3];
sx q[3];
rz(0.64938515) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53751078) q[0];
sx q[0];
rz(-2.4896955) q[0];
sx q[0];
rz(0.55484581) q[0];
rz(0.26508731) q[1];
sx q[1];
rz(-1.4684497) q[1];
sx q[1];
rz(1.9404985) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6938615) q[0];
sx q[0];
rz(-0.70708067) q[0];
sx q[0];
rz(0.62753824) q[0];
x q[1];
rz(-0.55055203) q[2];
sx q[2];
rz(-2.8209687) q[2];
sx q[2];
rz(0.52832802) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9399413) q[1];
sx q[1];
rz(-1.2229511) q[1];
sx q[1];
rz(1.531672) q[1];
rz(-pi) q[2];
x q[2];
rz(3.115377) q[3];
sx q[3];
rz(-2.391572) q[3];
sx q[3];
rz(-1.2312082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.889664) q[2];
sx q[2];
rz(-0.98968518) q[2];
sx q[2];
rz(2.2970693) q[2];
rz(2.0097513) q[3];
sx q[3];
rz(-0.90513888) q[3];
sx q[3];
rz(-1.7822781) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1691386) q[0];
sx q[0];
rz(-0.3903946) q[0];
sx q[0];
rz(1.1727232) q[0];
rz(0.68924618) q[1];
sx q[1];
rz(-1.0968364) q[1];
sx q[1];
rz(2.8776339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4580973) q[0];
sx q[0];
rz(-0.8120233) q[0];
sx q[0];
rz(2.5126636) q[0];
x q[1];
rz(-2.1446635) q[2];
sx q[2];
rz(-0.80748122) q[2];
sx q[2];
rz(1.6107744) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8431153) q[1];
sx q[1];
rz(-2.7727371) q[1];
sx q[1];
rz(-1.0333251) q[1];
rz(-1.4086776) q[3];
sx q[3];
rz(-2.4913284) q[3];
sx q[3];
rz(1.6780323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76137233) q[2];
sx q[2];
rz(-1.1470046) q[2];
sx q[2];
rz(-2.0210361) q[2];
rz(-1.606696) q[3];
sx q[3];
rz(-0.94529072) q[3];
sx q[3];
rz(-1.631261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071038889) q[0];
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
rz(-0.10590774) q[3];
sx q[3];
rz(-2.5220768) q[3];
sx q[3];
rz(0.22477748) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
