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
rz(-2.0353844) q[0];
sx q[0];
rz(-0.68618965) q[0];
sx q[0];
rz(1.9880779) q[0];
rz(-1.8312307) q[1];
sx q[1];
rz(-2.6481833) q[1];
sx q[1];
rz(-2.6931813) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7767186) q[0];
sx q[0];
rz(-1.5176511) q[0];
sx q[0];
rz(-2.4928635) q[0];
rz(-pi) q[1];
rz(-0.37031074) q[2];
sx q[2];
rz(-2.8680621) q[2];
sx q[2];
rz(2.6639378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7192711) q[1];
sx q[1];
rz(-2.0930674) q[1];
sx q[1];
rz(-2.2838255) q[1];
x q[2];
rz(-0.5011933) q[3];
sx q[3];
rz(-1.5287011) q[3];
sx q[3];
rz(0.30748366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.51556921) q[2];
sx q[2];
rz(-1.7531351) q[2];
sx q[2];
rz(2.5173729) q[2];
rz(-0.37912399) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(-2.8930801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9666331) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(-2.1061184) q[0];
rz(1.8677208) q[1];
sx q[1];
rz(-0.73718166) q[1];
sx q[1];
rz(-2.6587291) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.14014) q[0];
sx q[0];
rz(-1.4743544) q[0];
sx q[0];
rz(-1.6167377) q[0];
x q[1];
rz(2.5245978) q[2];
sx q[2];
rz(-2.9386387) q[2];
sx q[2];
rz(2.647612) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.711593) q[1];
sx q[1];
rz(-1.5801799) q[1];
sx q[1];
rz(1.1389772) q[1];
x q[2];
rz(-1.6522406) q[3];
sx q[3];
rz(-1.2508498) q[3];
sx q[3];
rz(-2.4575352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1172993) q[2];
sx q[2];
rz(-1.5017941) q[2];
sx q[2];
rz(-1.2705605) q[2];
rz(-0.67000669) q[3];
sx q[3];
rz(-2.5195401) q[3];
sx q[3];
rz(-0.89699927) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4033177) q[0];
sx q[0];
rz(-0.44602317) q[0];
sx q[0];
rz(-2.4025412) q[0];
rz(2.1801379) q[1];
sx q[1];
rz(-2.8196204) q[1];
sx q[1];
rz(-2.3846073) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8563499) q[0];
sx q[0];
rz(-2.3654571) q[0];
sx q[0];
rz(-1.9479284) q[0];
rz(-1.1951094) q[2];
sx q[2];
rz(-2.2268112) q[2];
sx q[2];
rz(-0.13740787) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0464335) q[1];
sx q[1];
rz(-2.1951402) q[1];
sx q[1];
rz(-1.4051564) q[1];
rz(-1.8784896) q[3];
sx q[3];
rz(-1.6374) q[3];
sx q[3];
rz(2.8860725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6387393) q[2];
sx q[2];
rz(-2.2245378) q[2];
sx q[2];
rz(0.70510954) q[2];
rz(-2.8201568) q[3];
sx q[3];
rz(-1.0175846) q[3];
sx q[3];
rz(-0.41079918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0099156378) q[0];
sx q[0];
rz(-0.83808815) q[0];
sx q[0];
rz(1.2257082) q[0];
rz(-1.907584) q[1];
sx q[1];
rz(-2.1261413) q[1];
sx q[1];
rz(0.066224901) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8787251) q[0];
sx q[0];
rz(-1.4811133) q[0];
sx q[0];
rz(-1.8171726) q[0];
x q[1];
rz(-0.29860626) q[2];
sx q[2];
rz(-2.3176607) q[2];
sx q[2];
rz(0.9048942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45715047) q[1];
sx q[1];
rz(-1.8891786) q[1];
sx q[1];
rz(-1.063698) q[1];
rz(2.7755967) q[3];
sx q[3];
rz(-1.4533236) q[3];
sx q[3];
rz(2.7155196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0464728) q[2];
sx q[2];
rz(-2.1472223) q[2];
sx q[2];
rz(2.7009916) q[2];
rz(-2.1353841) q[3];
sx q[3];
rz(-1.3738084) q[3];
sx q[3];
rz(0.83427507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7202268) q[0];
sx q[0];
rz(-0.81291968) q[0];
sx q[0];
rz(-1.6356069) q[0];
rz(1.6141363) q[1];
sx q[1];
rz(-1.9215877) q[1];
sx q[1];
rz(1.6857326) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2746873) q[0];
sx q[0];
rz(-0.84466776) q[0];
sx q[0];
rz(0.19498904) q[0];
x q[1];
rz(1.6873932) q[2];
sx q[2];
rz(-1.3449838) q[2];
sx q[2];
rz(1.9436398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2537186) q[1];
sx q[1];
rz(-1.5504596) q[1];
sx q[1];
rz(0.65171297) q[1];
rz(2.9221228) q[3];
sx q[3];
rz(-2.5220036) q[3];
sx q[3];
rz(-1.0426903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.27318925) q[2];
sx q[2];
rz(-1.578178) q[2];
sx q[2];
rz(-2.2301162) q[2];
rz(2.2677926) q[3];
sx q[3];
rz(-0.74207145) q[3];
sx q[3];
rz(-0.0066643683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28354302) q[0];
sx q[0];
rz(-0.25074211) q[0];
sx q[0];
rz(0.13033303) q[0];
rz(2.6538972) q[1];
sx q[1];
rz(-0.54134381) q[1];
sx q[1];
rz(0.74388751) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.795563) q[0];
sx q[0];
rz(-1.6206121) q[0];
sx q[0];
rz(1.4883785) q[0];
x q[1];
rz(-3.0024372) q[2];
sx q[2];
rz(-2.1916813) q[2];
sx q[2];
rz(-2.5992952) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86539093) q[1];
sx q[1];
rz(-1.2873833) q[1];
sx q[1];
rz(2.1332425) q[1];
rz(-0.48527511) q[3];
sx q[3];
rz(-1.4111184) q[3];
sx q[3];
rz(-1.6078469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.1193715) q[2];
sx q[2];
rz(-2.2890942) q[2];
sx q[2];
rz(-0.96088299) q[2];
rz(-0.39572257) q[3];
sx q[3];
rz(-1.8053677) q[3];
sx q[3];
rz(0.1161639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.457086) q[0];
sx q[0];
rz(-2.0261903) q[0];
sx q[0];
rz(2.0945666) q[0];
rz(-0.60415769) q[1];
sx q[1];
rz(-1.2774757) q[1];
sx q[1];
rz(-2.7488757) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1395124) q[0];
sx q[0];
rz(-0.89186984) q[0];
sx q[0];
rz(2.9154791) q[0];
rz(-1.4899859) q[2];
sx q[2];
rz(-2.5210292) q[2];
sx q[2];
rz(2.0106237) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3369357) q[1];
sx q[1];
rz(-1.5439022) q[1];
sx q[1];
rz(-2.8662445) q[1];
rz(-pi) q[2];
rz(2.8128871) q[3];
sx q[3];
rz(-1.8227326) q[3];
sx q[3];
rz(0.011289277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.087223209) q[2];
sx q[2];
rz(-1.9471709) q[2];
sx q[2];
rz(1.3090022) q[2];
rz(1.7878112) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(-0.27031171) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29276174) q[0];
sx q[0];
rz(-1.6852385) q[0];
sx q[0];
rz(0.30717474) q[0];
rz(-1.3245026) q[1];
sx q[1];
rz(-2.7565286) q[1];
sx q[1];
rz(-0.90075341) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9884315) q[0];
sx q[0];
rz(-1.5539196) q[0];
sx q[0];
rz(-0.26217006) q[0];
rz(2.6866954) q[2];
sx q[2];
rz(-0.96856801) q[2];
sx q[2];
rz(1.0554316) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.72462725) q[1];
sx q[1];
rz(-0.46675439) q[1];
sx q[1];
rz(-0.80893597) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6314478) q[3];
sx q[3];
rz(-1.4561781) q[3];
sx q[3];
rz(1.7607911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33425346) q[2];
sx q[2];
rz(-0.05143493) q[2];
sx q[2];
rz(-0.61484289) q[2];
rz(1.0963415) q[3];
sx q[3];
rz(-2.5494826) q[3];
sx q[3];
rz(-0.69826564) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39356247) q[0];
sx q[0];
rz(-0.63069558) q[0];
sx q[0];
rz(0.50931859) q[0];
rz(2.9604984) q[1];
sx q[1];
rz(-1.7362005) q[1];
sx q[1];
rz(-2.1720355) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80705331) q[0];
sx q[0];
rz(-1.4167804) q[0];
sx q[0];
rz(-0.75169433) q[0];
rz(0.050886919) q[2];
sx q[2];
rz(-2.0774088) q[2];
sx q[2];
rz(-1.8614872) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0694094) q[1];
sx q[1];
rz(-1.9095712) q[1];
sx q[1];
rz(1.4563926) q[1];
rz(-pi) q[2];
rz(-2.2990555) q[3];
sx q[3];
rz(-1.304783) q[3];
sx q[3];
rz(-1.0501978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8934882) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(-1.0815557) q[2];
rz(-0.23165101) q[3];
sx q[3];
rz(-1.3737498) q[3];
sx q[3];
rz(2.933568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3807826) q[0];
sx q[0];
rz(-2.91687) q[0];
sx q[0];
rz(0.64714062) q[0];
rz(0.45516792) q[1];
sx q[1];
rz(-1.4249233) q[1];
sx q[1];
rz(2.3796577) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3477701) q[0];
sx q[0];
rz(-0.88632727) q[0];
sx q[0];
rz(1.9467627) q[0];
rz(-pi) q[1];
rz(0.64024957) q[2];
sx q[2];
rz(-2.2138322) q[2];
sx q[2];
rz(-1.6563479) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6651407) q[1];
sx q[1];
rz(-2.2287205) q[1];
sx q[1];
rz(1.614078) q[1];
rz(-1.5555218) q[3];
sx q[3];
rz(-0.46746436) q[3];
sx q[3];
rz(1.8031507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0882988) q[2];
sx q[2];
rz(-0.21685728) q[2];
sx q[2];
rz(-2.2376412) q[2];
rz(-1.5229185) q[3];
sx q[3];
rz(-2.121033) q[3];
sx q[3];
rz(1.9797549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91697964) q[0];
sx q[0];
rz(-1.9539178) q[0];
sx q[0];
rz(1.7896347) q[0];
rz(0.12679535) q[1];
sx q[1];
rz(-2.3190111) q[1];
sx q[1];
rz(0.93228985) q[1];
rz(-1.1197208) q[2];
sx q[2];
rz(-1.3430165) q[2];
sx q[2];
rz(-2.7271885) q[2];
rz(-0.0062777304) q[3];
sx q[3];
rz(-2.7115887) q[3];
sx q[3];
rz(-0.030621519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
