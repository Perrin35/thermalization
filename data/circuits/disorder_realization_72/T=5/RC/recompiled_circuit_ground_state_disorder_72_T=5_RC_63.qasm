OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.44719625) q[0];
sx q[0];
rz(-1.5153272) q[0];
sx q[0];
rz(0.84994999) q[0];
rz(-1.7605468) q[1];
sx q[1];
rz(-2.2989506) q[1];
sx q[1];
rz(3.091264) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.91431) q[0];
sx q[0];
rz(-1.7512055) q[0];
sx q[0];
rz(1.4854027) q[0];
rz(-pi) q[1];
rz(-0.16762208) q[2];
sx q[2];
rz(-2.1838038) q[2];
sx q[2];
rz(-0.64152628) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1427372) q[1];
sx q[1];
rz(-1.6405447) q[1];
sx q[1];
rz(-0.19474366) q[1];
x q[2];
rz(1.6292455) q[3];
sx q[3];
rz(-0.93114595) q[3];
sx q[3];
rz(2.9906038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3927268) q[2];
sx q[2];
rz(-1.2458845) q[2];
sx q[2];
rz(2.5457814) q[2];
rz(-0.81123224) q[3];
sx q[3];
rz(-1.2657284) q[3];
sx q[3];
rz(0.59734145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3101462) q[0];
sx q[0];
rz(-0.50066384) q[0];
sx q[0];
rz(-1.3432107) q[0];
rz(1.385618) q[1];
sx q[1];
rz(-1.132553) q[1];
sx q[1];
rz(-1.064942) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32414054) q[0];
sx q[0];
rz(-1.5498501) q[0];
sx q[0];
rz(3.0684242) q[0];
rz(-pi) q[1];
rz(2.166719) q[2];
sx q[2];
rz(-0.6069912) q[2];
sx q[2];
rz(-0.18894228) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2328223) q[1];
sx q[1];
rz(-1.5651287) q[1];
sx q[1];
rz(3.0012111) q[1];
rz(-1.2805544) q[3];
sx q[3];
rz(-0.10435552) q[3];
sx q[3];
rz(-0.39784986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.500835) q[2];
sx q[2];
rz(-1.0229599) q[2];
sx q[2];
rz(0.41066059) q[2];
rz(1.143035) q[3];
sx q[3];
rz(-2.394702) q[3];
sx q[3];
rz(0.66543287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(1.2086585) q[0];
sx q[0];
rz(-1.3422796) q[0];
sx q[0];
rz(-0.68516532) q[0];
rz(-0.67600983) q[1];
sx q[1];
rz(-1.6513377) q[1];
sx q[1];
rz(-2.4032059) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15085193) q[0];
sx q[0];
rz(-0.60078492) q[0];
sx q[0];
rz(2.5505801) q[0];
rz(0.73263158) q[2];
sx q[2];
rz(-1.5937798) q[2];
sx q[2];
rz(-2.2642291) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8540875) q[1];
sx q[1];
rz(-1.7190134) q[1];
sx q[1];
rz(0.052637859) q[1];
x q[2];
rz(2.4285391) q[3];
sx q[3];
rz(-1.0950297) q[3];
sx q[3];
rz(0.66407628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43397778) q[2];
sx q[2];
rz(-1.6699764) q[2];
sx q[2];
rz(-2.1493256) q[2];
rz(2.9662002) q[3];
sx q[3];
rz(-0.49687353) q[3];
sx q[3];
rz(2.8309256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53395143) q[0];
sx q[0];
rz(-2.0168004) q[0];
sx q[0];
rz(-2.9006309) q[0];
rz(-2.2843649) q[1];
sx q[1];
rz(-0.78397426) q[1];
sx q[1];
rz(1.2711752) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1403811) q[0];
sx q[0];
rz(-2.0095946) q[0];
sx q[0];
rz(-2.9649404) q[0];
rz(-2.6189162) q[2];
sx q[2];
rz(-1.9999522) q[2];
sx q[2];
rz(1.9945838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7422116) q[1];
sx q[1];
rz(-1.9866418) q[1];
sx q[1];
rz(2.0539961) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.70636) q[3];
sx q[3];
rz(-2.3381544) q[3];
sx q[3];
rz(-2.9627707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6406389) q[2];
sx q[2];
rz(-0.95879889) q[2];
sx q[2];
rz(2.1232088) q[2];
rz(-1.1650677) q[3];
sx q[3];
rz(-2.1924721) q[3];
sx q[3];
rz(0.44786662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1027706) q[0];
sx q[0];
rz(-0.7907246) q[0];
sx q[0];
rz(0.08654174) q[0];
rz(-1.4718919) q[1];
sx q[1];
rz(-2.4457928) q[1];
sx q[1];
rz(-2.5873628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80866586) q[0];
sx q[0];
rz(-0.73360591) q[0];
sx q[0];
rz(-0.33302078) q[0];
x q[1];
rz(1.4325028) q[2];
sx q[2];
rz(-0.90511647) q[2];
sx q[2];
rz(0.67065698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8468626) q[1];
sx q[1];
rz(-2.1536441) q[1];
sx q[1];
rz(2.1267049) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0114927) q[3];
sx q[3];
rz(-1.3050858) q[3];
sx q[3];
rz(1.2540224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8467796) q[2];
sx q[2];
rz(-2.7816935) q[2];
sx q[2];
rz(0.93152085) q[2];
rz(0.31342634) q[3];
sx q[3];
rz(-0.64160186) q[3];
sx q[3];
rz(0.48345598) q[3];
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
rz(2.038092) q[0];
sx q[0];
rz(-0.8041389) q[0];
sx q[0];
rz(-1.3793797) q[0];
rz(0.74229678) q[1];
sx q[1];
rz(-1.392044) q[1];
sx q[1];
rz(2.0753863) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.659935) q[0];
sx q[0];
rz(-1.602442) q[0];
sx q[0];
rz(-0.13755799) q[0];
rz(-pi) q[1];
rz(-0.13958884) q[2];
sx q[2];
rz(-2.0878279) q[2];
sx q[2];
rz(-1.6986183) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.226245) q[1];
sx q[1];
rz(-0.97999014) q[1];
sx q[1];
rz(0.3027447) q[1];
x q[2];
rz(2.8442758) q[3];
sx q[3];
rz(-2.8585189) q[3];
sx q[3];
rz(-2.360614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1833056) q[2];
sx q[2];
rz(-2.1669407) q[2];
sx q[2];
rz(0.24678123) q[2];
rz(-1.0684446) q[3];
sx q[3];
rz(-1.1698134) q[3];
sx q[3];
rz(-1.0698414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5088155) q[0];
sx q[0];
rz(-0.95938534) q[0];
sx q[0];
rz(-2.7799613) q[0];
rz(1.5536785) q[1];
sx q[1];
rz(-1.5767153) q[1];
sx q[1];
rz(-1.9379001) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7772352) q[0];
sx q[0];
rz(-0.89221749) q[0];
sx q[0];
rz(-0.70058544) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3839528) q[2];
sx q[2];
rz(-2.0645541) q[2];
sx q[2];
rz(0.46390033) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0140397) q[1];
sx q[1];
rz(-1.7726092) q[1];
sx q[1];
rz(2.1743808) q[1];
x q[2];
rz(-0.43872469) q[3];
sx q[3];
rz(-1.8774021) q[3];
sx q[3];
rz(2.1244989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2162073) q[2];
sx q[2];
rz(-1.46571) q[2];
sx q[2];
rz(-2.1290131) q[2];
rz(0.25977627) q[3];
sx q[3];
rz(-2.8949819) q[3];
sx q[3];
rz(2.7369505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15553661) q[0];
sx q[0];
rz(-1.8590834) q[0];
sx q[0];
rz(-2.3543661) q[0];
rz(0.8440482) q[1];
sx q[1];
rz(-2.7828352) q[1];
sx q[1];
rz(-1.2740096) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9106261) q[0];
sx q[0];
rz(-1.2342802) q[0];
sx q[0];
rz(-1.5084159) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34113555) q[2];
sx q[2];
rz(-0.76934177) q[2];
sx q[2];
rz(-3.0770707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9032517) q[1];
sx q[1];
rz(-2.1834186) q[1];
sx q[1];
rz(-0.56507052) q[1];
x q[2];
rz(1.6350217) q[3];
sx q[3];
rz(-1.9750486) q[3];
sx q[3];
rz(-0.4413213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0969703) q[2];
sx q[2];
rz(-1.5771834) q[2];
sx q[2];
rz(2.7596562) q[2];
rz(1.1218128) q[3];
sx q[3];
rz(-2.8048281) q[3];
sx q[3];
rz(-0.063966123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0762416) q[0];
sx q[0];
rz(-0.92472804) q[0];
sx q[0];
rz(-1.4724154) q[0];
rz(2.8352101) q[1];
sx q[1];
rz(-0.52426052) q[1];
sx q[1];
rz(2.870097) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6097019) q[0];
sx q[0];
rz(-1.176495) q[0];
sx q[0];
rz(2.0608725) q[0];
x q[1];
rz(-0.36868544) q[2];
sx q[2];
rz(-2.2767608) q[2];
sx q[2];
rz(2.5384704) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68743905) q[1];
sx q[1];
rz(-2.1104129) q[1];
sx q[1];
rz(1.0034877) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.027609) q[3];
sx q[3];
rz(-1.9577259) q[3];
sx q[3];
rz(-0.99150211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86758119) q[2];
sx q[2];
rz(-2.0143955) q[2];
sx q[2];
rz(-2.4809044) q[2];
rz(0.88179669) q[3];
sx q[3];
rz(-2.5776358) q[3];
sx q[3];
rz(-0.86625117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8771186) q[0];
sx q[0];
rz(-1.6082123) q[0];
sx q[0];
rz(1.8033173) q[0];
rz(-1.0461461) q[1];
sx q[1];
rz(-1.0271065) q[1];
sx q[1];
rz(2.963692) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.698211) q[0];
sx q[0];
rz(-2.215235) q[0];
sx q[0];
rz(2.3257735) q[0];
rz(-2.1523802) q[2];
sx q[2];
rz(-1.7187558) q[2];
sx q[2];
rz(-1.6616247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2819351) q[1];
sx q[1];
rz(-1.2392231) q[1];
sx q[1];
rz(2.8012707) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4334861) q[3];
sx q[3];
rz(-1.2009283) q[3];
sx q[3];
rz(0.22320492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.88343128) q[2];
sx q[2];
rz(-2.0475755) q[2];
sx q[2];
rz(0.9359614) q[2];
rz(0.59518138) q[3];
sx q[3];
rz(-1.7116356) q[3];
sx q[3];
rz(-0.87966758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7137939) q[0];
sx q[0];
rz(-1.6811163) q[0];
sx q[0];
rz(1.2291193) q[0];
rz(1.3270558) q[1];
sx q[1];
rz(-1.6358903) q[1];
sx q[1];
rz(-1.9930175) q[1];
rz(0.11360609) q[2];
sx q[2];
rz(-1.3945711) q[2];
sx q[2];
rz(-0.61749189) q[2];
rz(-0.26575539) q[3];
sx q[3];
rz(-1.8674217) q[3];
sx q[3];
rz(-2.7056497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
