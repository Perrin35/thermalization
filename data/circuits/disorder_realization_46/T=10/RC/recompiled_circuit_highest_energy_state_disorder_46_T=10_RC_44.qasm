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
rz(2.1815648) q[0];
sx q[0];
rz(-0.61909827) q[0];
sx q[0];
rz(1.4033138) q[0];
rz(2.4461441) q[1];
sx q[1];
rz(-0.56013501) q[1];
sx q[1];
rz(-2.4897895) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2346828) q[0];
sx q[0];
rz(-1.2359338) q[0];
sx q[0];
rz(-0.85483179) q[0];
rz(-pi) q[1];
rz(2.8601134) q[2];
sx q[2];
rz(-1.5764766) q[2];
sx q[2];
rz(-1.3303632) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6095602) q[1];
sx q[1];
rz(-1.3390307) q[1];
sx q[1];
rz(1.5941117) q[1];
rz(-pi) q[2];
rz(0.82020369) q[3];
sx q[3];
rz(-1.3705882) q[3];
sx q[3];
rz(0.28991308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2191676) q[2];
sx q[2];
rz(-1.3895915) q[2];
sx q[2];
rz(2.8233042) q[2];
rz(0.90218941) q[3];
sx q[3];
rz(-2.2212494) q[3];
sx q[3];
rz(-2.2778146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1282463) q[0];
sx q[0];
rz(-2.8885169) q[0];
sx q[0];
rz(2.6213562) q[0];
rz(3.0546313) q[1];
sx q[1];
rz(-2.088701) q[1];
sx q[1];
rz(0.797995) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3083852) q[0];
sx q[0];
rz(-1.7705581) q[0];
sx q[0];
rz(2.8248592) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3981061) q[2];
sx q[2];
rz(-1.6891589) q[2];
sx q[2];
rz(1.7497045) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2472331) q[1];
sx q[1];
rz(-1.9590568) q[1];
sx q[1];
rz(-1.7227933) q[1];
rz(-pi) q[2];
rz(-0.4466855) q[3];
sx q[3];
rz(-1.2036714) q[3];
sx q[3];
rz(-1.1243096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8331376) q[2];
sx q[2];
rz(-1.986958) q[2];
sx q[2];
rz(-2.7412565) q[2];
rz(-0.18276246) q[3];
sx q[3];
rz(-0.57923135) q[3];
sx q[3];
rz(1.8496752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12277814) q[0];
sx q[0];
rz(-0.084538758) q[0];
sx q[0];
rz(-1.2070745) q[0];
rz(1.6385551) q[1];
sx q[1];
rz(-1.259558) q[1];
sx q[1];
rz(2.9685793) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91813501) q[0];
sx q[0];
rz(-1.4854781) q[0];
sx q[0];
rz(1.1847909) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0315057) q[2];
sx q[2];
rz(-0.99186691) q[2];
sx q[2];
rz(-2.8814709) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8743778) q[1];
sx q[1];
rz(-2.3372531) q[1];
sx q[1];
rz(2.4206603) q[1];
rz(-pi) q[2];
rz(1.7812114) q[3];
sx q[3];
rz(-0.83411479) q[3];
sx q[3];
rz(0.93005108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.59306899) q[2];
sx q[2];
rz(-2.5544281) q[2];
sx q[2];
rz(2.1068088) q[2];
rz(-2.2506524) q[3];
sx q[3];
rz(-2.2726629) q[3];
sx q[3];
rz(1.2263891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064875038) q[0];
sx q[0];
rz(-0.95219505) q[0];
sx q[0];
rz(-0.12884831) q[0];
rz(-2.7128362) q[1];
sx q[1];
rz(-1.7196451) q[1];
sx q[1];
rz(1.4770003) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7744666) q[0];
sx q[0];
rz(-1.2259577) q[0];
sx q[0];
rz(1.4982629) q[0];
x q[1];
rz(-1.1294133) q[2];
sx q[2];
rz(-1.1611106) q[2];
sx q[2];
rz(-0.022155174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.99812388) q[1];
sx q[1];
rz(-1.3593771) q[1];
sx q[1];
rz(1.9969673) q[1];
rz(-3.1404895) q[3];
sx q[3];
rz(-0.97954053) q[3];
sx q[3];
rz(-1.7756893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9513272) q[2];
sx q[2];
rz(-1.724388) q[2];
sx q[2];
rz(-1.8180656) q[2];
rz(-1.6773112) q[3];
sx q[3];
rz(-2.5530294) q[3];
sx q[3];
rz(-0.67414362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73439634) q[0];
sx q[0];
rz(-2.7792327) q[0];
sx q[0];
rz(-1.4417484) q[0];
rz(2.7184519) q[1];
sx q[1];
rz(-0.95760456) q[1];
sx q[1];
rz(0.29022455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5231756) q[0];
sx q[0];
rz(-1.8632065) q[0];
sx q[0];
rz(-0.66177701) q[0];
rz(-pi) q[1];
rz(-2.5686479) q[2];
sx q[2];
rz(-2.101295) q[2];
sx q[2];
rz(-2.6823519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2331302) q[1];
sx q[1];
rz(-1.4428348) q[1];
sx q[1];
rz(-1.5217693) q[1];
x q[2];
rz(0.30778389) q[3];
sx q[3];
rz(-1.3213385) q[3];
sx q[3];
rz(-1.4741999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4238802) q[2];
sx q[2];
rz(-2.3545357) q[2];
sx q[2];
rz(0.94669739) q[2];
rz(-0.26149073) q[3];
sx q[3];
rz(-1.7777092) q[3];
sx q[3];
rz(2.7546895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70447266) q[0];
sx q[0];
rz(-1.3334714) q[0];
sx q[0];
rz(-0.76960808) q[0];
rz(-0.17598027) q[1];
sx q[1];
rz(-1.3409921) q[1];
sx q[1];
rz(-1.5716858) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0933857) q[0];
sx q[0];
rz(-2.1820659) q[0];
sx q[0];
rz(2.3552191) q[0];
rz(-2.2899175) q[2];
sx q[2];
rz(-1.6234509) q[2];
sx q[2];
rz(2.5824314) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4665199) q[1];
sx q[1];
rz(-2.7278455) q[1];
sx q[1];
rz(-2.6690261) q[1];
rz(-pi) q[2];
rz(0.50521781) q[3];
sx q[3];
rz(-1.1105624) q[3];
sx q[3];
rz(-1.203618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1311329) q[2];
sx q[2];
rz(-2.185148) q[2];
sx q[2];
rz(1.825911) q[2];
rz(-1.9192421) q[3];
sx q[3];
rz(-2.0171916) q[3];
sx q[3];
rz(0.30092064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2216126) q[0];
sx q[0];
rz(-1.6918809) q[0];
sx q[0];
rz(1.8020887) q[0];
rz(1.1969396) q[1];
sx q[1];
rz(-1.0547799) q[1];
sx q[1];
rz(-2.1766677) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20112637) q[0];
sx q[0];
rz(-2.622704) q[0];
sx q[0];
rz(0.19176439) q[0];
x q[1];
rz(2.8469725) q[2];
sx q[2];
rz(-1.2333394) q[2];
sx q[2];
rz(1.9532415) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8037012) q[1];
sx q[1];
rz(-1.6543904) q[1];
sx q[1];
rz(0.22050942) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32437848) q[3];
sx q[3];
rz(-2.6357963) q[3];
sx q[3];
rz(-1.5816816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7397466) q[2];
sx q[2];
rz(-1.2082938) q[2];
sx q[2];
rz(-1.9604663) q[2];
rz(3.127023) q[3];
sx q[3];
rz(-0.40112344) q[3];
sx q[3];
rz(1.1723664) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3907617) q[0];
sx q[0];
rz(-1.9347235) q[0];
sx q[0];
rz(-2.3857351) q[0];
rz(2.4769056) q[1];
sx q[1];
rz(-1.942778) q[1];
sx q[1];
rz(1.8519148) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3493746) q[0];
sx q[0];
rz(-1.5070033) q[0];
sx q[0];
rz(-2.347058) q[0];
rz(1.2694744) q[2];
sx q[2];
rz(-1.0507283) q[2];
sx q[2];
rz(-3.0406496) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3933223) q[1];
sx q[1];
rz(-0.38667187) q[1];
sx q[1];
rz(2.3063117) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4699072) q[3];
sx q[3];
rz(-1.2630188) q[3];
sx q[3];
rz(0.65478126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0083302) q[2];
sx q[2];
rz(-1.8066758) q[2];
sx q[2];
rz(-1.5745715) q[2];
rz(-1.745892) q[3];
sx q[3];
rz(-1.402366) q[3];
sx q[3];
rz(0.76656669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.72346) q[0];
sx q[0];
rz(-0.27579951) q[0];
sx q[0];
rz(-2.7511399) q[0];
rz(-2.0402724) q[1];
sx q[1];
rz(-1.3464255) q[1];
sx q[1];
rz(2.2094545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.348968) q[0];
sx q[0];
rz(-1.5431964) q[0];
sx q[0];
rz(-1.6076002) q[0];
rz(-2.6022692) q[2];
sx q[2];
rz(-1.5289837) q[2];
sx q[2];
rz(1.5635179) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.336672) q[1];
sx q[1];
rz(-1.3778566) q[1];
sx q[1];
rz(2.8590917) q[1];
rz(-0.45472179) q[3];
sx q[3];
rz(-2.623436) q[3];
sx q[3];
rz(1.3738812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0497047) q[2];
sx q[2];
rz(-2.4447417) q[2];
sx q[2];
rz(1.4616802) q[2];
rz(1.0852496) q[3];
sx q[3];
rz(-2.0094252) q[3];
sx q[3];
rz(2.5011044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.28814462) q[0];
sx q[0];
rz(-2.4818821) q[0];
sx q[0];
rz(1.1261384) q[0];
rz(2.0938734) q[1];
sx q[1];
rz(-2.5151217) q[1];
sx q[1];
rz(1.610021) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0275912) q[0];
sx q[0];
rz(-0.023311255) q[0];
sx q[0];
rz(3.0816881) q[0];
rz(-pi) q[1];
rz(0.44785277) q[2];
sx q[2];
rz(-0.50074762) q[2];
sx q[2];
rz(-0.66381493) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6781552) q[1];
sx q[1];
rz(-3.040294) q[1];
sx q[1];
rz(-0.8763635) q[1];
x q[2];
rz(0.13981847) q[3];
sx q[3];
rz(-0.56666683) q[3];
sx q[3];
rz(1.293928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62022007) q[2];
sx q[2];
rz(-2.3282101) q[2];
sx q[2];
rz(-0.20624557) q[2];
rz(-0.5591875) q[3];
sx q[3];
rz(-1.6181479) q[3];
sx q[3];
rz(2.0172113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6013721) q[0];
sx q[0];
rz(-1.2982397) q[0];
sx q[0];
rz(0.44197767) q[0];
rz(1.9961009) q[1];
sx q[1];
rz(-1.835123) q[1];
sx q[1];
rz(1.8225972) q[1];
rz(-0.27362846) q[2];
sx q[2];
rz(-1.3009334) q[2];
sx q[2];
rz(1.0197659) q[2];
rz(1.3973936) q[3];
sx q[3];
rz(-2.2874293) q[3];
sx q[3];
rz(2.3746272) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
