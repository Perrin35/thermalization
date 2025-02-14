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
rz(1.1062082) q[0];
sx q[0];
rz(-2.455403) q[0];
sx q[0];
rz(-1.9880779) q[0];
rz(-1.8312307) q[1];
sx q[1];
rz(-2.6481833) q[1];
sx q[1];
rz(0.44841132) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9759244) q[0];
sx q[0];
rz(-0.92313719) q[0];
sx q[0];
rz(-1.5041385) q[0];
x q[1];
rz(2.7712819) q[2];
sx q[2];
rz(-0.2735306) q[2];
sx q[2];
rz(0.47765484) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37473512) q[1];
sx q[1];
rz(-2.285706) q[1];
sx q[1];
rz(0.84918569) q[1];
x q[2];
rz(3.0541522) q[3];
sx q[3];
rz(-2.6387847) q[3];
sx q[3];
rz(1.3399762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51556921) q[2];
sx q[2];
rz(-1.3884576) q[2];
sx q[2];
rz(-2.5173729) q[2];
rz(0.37912399) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(2.8930801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9666331) q[0];
sx q[0];
rz(-2.3269854) q[0];
sx q[0];
rz(1.0354743) q[0];
rz(1.8677208) q[1];
sx q[1];
rz(-0.73718166) q[1];
sx q[1];
rz(0.48286352) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55601701) q[0];
sx q[0];
rz(-0.10679467) q[0];
sx q[0];
rz(0.44321816) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5245978) q[2];
sx q[2];
rz(-0.20295396) q[2];
sx q[2];
rz(-2.647612) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9804364) q[1];
sx q[1];
rz(-0.43191467) q[1];
sx q[1];
rz(-1.5932139) q[1];
rz(-pi) q[2];
x q[2];
rz(1.489352) q[3];
sx q[3];
rz(-1.8907428) q[3];
sx q[3];
rz(2.4575352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.024293385) q[2];
sx q[2];
rz(-1.5017941) q[2];
sx q[2];
rz(1.8710322) q[2];
rz(-2.471586) q[3];
sx q[3];
rz(-0.62205258) q[3];
sx q[3];
rz(2.2445934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73827493) q[0];
sx q[0];
rz(-0.44602317) q[0];
sx q[0];
rz(-0.73905149) q[0];
rz(-0.96145472) q[1];
sx q[1];
rz(-2.8196204) q[1];
sx q[1];
rz(-2.3846073) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2852427) q[0];
sx q[0];
rz(-2.3654571) q[0];
sx q[0];
rz(1.9479284) q[0];
x q[1];
rz(-0.44481014) q[2];
sx q[2];
rz(-0.74193566) q[2];
sx q[2];
rz(2.7049989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0464335) q[1];
sx q[1];
rz(-2.1951402) q[1];
sx q[1];
rz(1.7364362) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0717172) q[3];
sx q[3];
rz(-1.2638076) q[3];
sx q[3];
rz(1.8474634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6387393) q[2];
sx q[2];
rz(-2.2245378) q[2];
sx q[2];
rz(2.4364831) q[2];
rz(2.8201568) q[3];
sx q[3];
rz(-2.1240081) q[3];
sx q[3];
rz(-0.41079918) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.131677) q[0];
sx q[0];
rz(-2.3035045) q[0];
sx q[0];
rz(-1.2257082) q[0];
rz(-1.2340087) q[1];
sx q[1];
rz(-2.1261413) q[1];
sx q[1];
rz(-0.066224901) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8561838) q[0];
sx q[0];
rz(-1.8161621) q[0];
sx q[0];
rz(3.0491327) q[0];
rz(-pi) q[1];
rz(1.8784857) q[2];
sx q[2];
rz(-0.79350797) q[2];
sx q[2];
rz(-0.47946489) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.855763) q[1];
sx q[1];
rz(-2.0501839) q[1];
sx q[1];
rz(-2.781032) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31853557) q[3];
sx q[3];
rz(-2.7580166) q[3];
sx q[3];
rz(-0.84795241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0464728) q[2];
sx q[2];
rz(-0.9943704) q[2];
sx q[2];
rz(2.7009916) q[2];
rz(2.1353841) q[3];
sx q[3];
rz(-1.7677842) q[3];
sx q[3];
rz(-2.3073176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7202268) q[0];
sx q[0];
rz(-2.328673) q[0];
sx q[0];
rz(-1.6356069) q[0];
rz(1.5274564) q[1];
sx q[1];
rz(-1.2200049) q[1];
sx q[1];
rz(1.6857326) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86690534) q[0];
sx q[0];
rz(-0.84466776) q[0];
sx q[0];
rz(0.19498904) q[0];
rz(-pi) q[1];
rz(-2.6728476) q[2];
sx q[2];
rz(-0.25368099) q[2];
sx q[2];
rz(1.6799334) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8511635) q[1];
sx q[1];
rz(-2.4896087) q[1];
sx q[1];
rz(-0.033520582) q[1];
rz(-pi) q[2];
rz(1.7248575) q[3];
sx q[3];
rz(-0.96821456) q[3];
sx q[3];
rz(1.3101206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8684034) q[2];
sx q[2];
rz(-1.578178) q[2];
sx q[2];
rz(0.91147649) q[2];
rz(0.8738001) q[3];
sx q[3];
rz(-0.74207145) q[3];
sx q[3];
rz(-3.1349283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28354302) q[0];
sx q[0];
rz(-2.8908505) q[0];
sx q[0];
rz(3.0112596) q[0];
rz(2.6538972) q[1];
sx q[1];
rz(-0.54134381) q[1];
sx q[1];
rz(-2.3977051) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.795563) q[0];
sx q[0];
rz(-1.6206121) q[0];
sx q[0];
rz(1.6532142) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3792453) q[2];
sx q[2];
rz(-2.5073176) q[2];
sx q[2];
rz(-0.77855643) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53092693) q[1];
sx q[1];
rz(-2.1082889) q[1];
sx q[1];
rz(2.8100138) q[1];
rz(-pi) q[2];
rz(2.6563175) q[3];
sx q[3];
rz(-1.4111184) q[3];
sx q[3];
rz(-1.6078469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0222212) q[2];
sx q[2];
rz(-0.85249844) q[2];
sx q[2];
rz(-0.96088299) q[2];
rz(2.7458701) q[3];
sx q[3];
rz(-1.8053677) q[3];
sx q[3];
rz(-3.0254288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6845067) q[0];
sx q[0];
rz(-1.1154024) q[0];
sx q[0];
rz(-2.0945666) q[0];
rz(2.537435) q[1];
sx q[1];
rz(-1.2774757) q[1];
sx q[1];
rz(0.39271694) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1395124) q[0];
sx q[0];
rz(-2.2497228) q[0];
sx q[0];
rz(2.9154791) q[0];
x q[1];
rz(-1.6516067) q[2];
sx q[2];
rz(-0.62056345) q[2];
sx q[2];
rz(-1.1309689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8046569) q[1];
sx q[1];
rz(-1.5976904) q[1];
sx q[1];
rz(-2.8662445) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.673135) q[3];
sx q[3];
rz(-2.7302448) q[3];
sx q[3];
rz(2.213221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0543694) q[2];
sx q[2];
rz(-1.1944218) q[2];
sx q[2];
rz(-1.8325904) q[2];
rz(-1.3537815) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(-0.27031171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.8488309) q[0];
sx q[0];
rz(-1.6852385) q[0];
sx q[0];
rz(-2.8344179) q[0];
rz(1.3245026) q[1];
sx q[1];
rz(-0.38506404) q[1];
sx q[1];
rz(-0.90075341) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4221638) q[0];
sx q[0];
rz(-1.3086645) q[0];
sx q[0];
rz(1.5882701) q[0];
x q[1];
rz(-1.0020761) q[2];
sx q[2];
rz(-2.4042712) q[2];
sx q[2];
rz(-1.7676644) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5985721) q[1];
sx q[1];
rz(-1.2391587) q[1];
sx q[1];
rz(2.806862) q[1];
rz(-pi) q[2];
rz(-0.19260223) q[3];
sx q[3];
rz(-2.5012272) q[3];
sx q[3];
rz(0.34512025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.33425346) q[2];
sx q[2];
rz(-3.0901577) q[2];
sx q[2];
rz(-2.5267498) q[2];
rz(-2.0452512) q[3];
sx q[3];
rz(-2.5494826) q[3];
sx q[3];
rz(-0.69826564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39356247) q[0];
sx q[0];
rz(-0.63069558) q[0];
sx q[0];
rz(-2.6322741) q[0];
rz(-2.9604984) q[1];
sx q[1];
rz(-1.4053922) q[1];
sx q[1];
rz(0.96955713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60107952) q[0];
sx q[0];
rz(-0.76427312) q[0];
sx q[0];
rz(2.918052) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0779586) q[2];
sx q[2];
rz(-1.5263057) q[2];
sx q[2];
rz(-0.26598334) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6020376) q[1];
sx q[1];
rz(-1.6786715) q[1];
sx q[1];
rz(-0.34083582) q[1];
rz(-pi) q[2];
rz(2.7915482) q[3];
sx q[3];
rz(-2.2681103) q[3];
sx q[3];
rz(0.75087912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8934882) q[2];
sx q[2];
rz(-2.9073145) q[2];
sx q[2];
rz(-2.060037) q[2];
rz(-0.23165101) q[3];
sx q[3];
rz(-1.3737498) q[3];
sx q[3];
rz(-0.20802465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.76081) q[0];
sx q[0];
rz(-2.91687) q[0];
sx q[0];
rz(-0.64714062) q[0];
rz(0.45516792) q[1];
sx q[1];
rz(-1.7166694) q[1];
sx q[1];
rz(0.761935) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9058911) q[0];
sx q[0];
rz(-0.76607031) q[0];
sx q[0];
rz(-0.42278843) q[0];
rz(-pi) q[1];
rz(2.5013431) q[2];
sx q[2];
rz(-2.2138322) q[2];
sx q[2];
rz(1.6563479) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0207723) q[1];
sx q[1];
rz(-1.6050395) q[1];
sx q[1];
rz(2.4832151) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1033789) q[3];
sx q[3];
rz(-1.5776792) q[3];
sx q[3];
rz(0.21871834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.053293856) q[2];
sx q[2];
rz(-2.9247354) q[2];
sx q[2];
rz(2.2376412) q[2];
rz(1.6186742) q[3];
sx q[3];
rz(-2.121033) q[3];
sx q[3];
rz(-1.1618377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91697964) q[0];
sx q[0];
rz(-1.9539178) q[0];
sx q[0];
rz(1.7896347) q[0];
rz(3.0147973) q[1];
sx q[1];
rz(-0.8225816) q[1];
sx q[1];
rz(-2.2093028) q[1];
rz(-2.0595111) q[2];
sx q[2];
rz(-2.6398224) q[2];
sx q[2];
rz(2.421438) q[2];
rz(-1.5736754) q[3];
sx q[3];
rz(-1.1408014) q[3];
sx q[3];
rz(-0.037527966) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
