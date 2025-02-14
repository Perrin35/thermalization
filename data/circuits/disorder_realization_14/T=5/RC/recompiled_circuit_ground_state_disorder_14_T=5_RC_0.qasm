OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2603962) q[0];
sx q[0];
rz(-1.2328147) q[0];
sx q[0];
rz(0.32786274) q[0];
rz(-0.39363632) q[1];
sx q[1];
rz(4.7596158) q[1];
sx q[1];
rz(14.16852) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6820591) q[0];
sx q[0];
rz(-1.5757808) q[0];
sx q[0];
rz(-0.2010767) q[0];
rz(-1.0225564) q[2];
sx q[2];
rz(-1.3371144) q[2];
sx q[2];
rz(1.7811687) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1069168) q[1];
sx q[1];
rz(-0.44245103) q[1];
sx q[1];
rz(-0.4974858) q[1];
rz(-pi) q[2];
rz(2.7490535) q[3];
sx q[3];
rz(-1.2219567) q[3];
sx q[3];
rz(-0.67482812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9615122) q[2];
sx q[2];
rz(-0.95246685) q[2];
sx q[2];
rz(2.3251779) q[2];
rz(2.4444729) q[3];
sx q[3];
rz(-0.37728798) q[3];
sx q[3];
rz(-0.94059801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3546251) q[0];
sx q[0];
rz(-2.0836199) q[0];
sx q[0];
rz(2.4355198) q[0];
rz(-2.970703) q[1];
sx q[1];
rz(-2.6306174) q[1];
sx q[1];
rz(-0.99542803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1814897) q[0];
sx q[0];
rz(-1.8892291) q[0];
sx q[0];
rz(0.35950999) q[0];
rz(-pi) q[1];
rz(-2.273492) q[2];
sx q[2];
rz(-1.1836402) q[2];
sx q[2];
rz(1.1328982) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5157516) q[1];
sx q[1];
rz(-1.7236354) q[1];
sx q[1];
rz(1.7865074) q[1];
x q[2];
rz(-1.6042406) q[3];
sx q[3];
rz(-2.20801) q[3];
sx q[3];
rz(-0.49195615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5848315) q[2];
sx q[2];
rz(-2.2274667) q[2];
sx q[2];
rz(0.27908984) q[2];
rz(-0.68785214) q[3];
sx q[3];
rz(-2.9180241) q[3];
sx q[3];
rz(-2.1384625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2897709) q[0];
sx q[0];
rz(-2.2230447) q[0];
sx q[0];
rz(-0.37891349) q[0];
rz(-1.1907578) q[1];
sx q[1];
rz(-0.36668101) q[1];
sx q[1];
rz(-2.3653638) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8541753) q[0];
sx q[0];
rz(-0.45817927) q[0];
sx q[0];
rz(-1.929438) q[0];
rz(1.7227145) q[2];
sx q[2];
rz(-1.7524652) q[2];
sx q[2];
rz(2.7318397) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32404985) q[1];
sx q[1];
rz(-2.2607798) q[1];
sx q[1];
rz(1.2381366) q[1];
rz(-2.0950731) q[3];
sx q[3];
rz(-1.610743) q[3];
sx q[3];
rz(1.0429045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0188521) q[2];
sx q[2];
rz(-1.9032225) q[2];
sx q[2];
rz(-2.5416601) q[2];
rz(-1.2395202) q[3];
sx q[3];
rz(-1.7299088) q[3];
sx q[3];
rz(2.5304573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3032853) q[0];
sx q[0];
rz(-0.67876434) q[0];
sx q[0];
rz(1.432206) q[0];
rz(-0.86668658) q[1];
sx q[1];
rz(-1.7984093) q[1];
sx q[1];
rz(2.0492882) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1537395) q[0];
sx q[0];
rz(-0.9747552) q[0];
sx q[0];
rz(-0.88576646) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.230365) q[2];
sx q[2];
rz(-1.22515) q[2];
sx q[2];
rz(-2.4663053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0580313) q[1];
sx q[1];
rz(-1.5688854) q[1];
sx q[1];
rz(-0.46431904) q[1];
rz(-pi) q[2];
rz(-0.82073145) q[3];
sx q[3];
rz(-0.98086517) q[3];
sx q[3];
rz(1.8128519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0852574) q[2];
sx q[2];
rz(-1.9838355) q[2];
sx q[2];
rz(1.2752656) q[2];
rz(1.3750252) q[3];
sx q[3];
rz(-1.2396038) q[3];
sx q[3];
rz(-2.1084771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.035868693) q[0];
sx q[0];
rz(-2.2585456) q[0];
sx q[0];
rz(1.7778273) q[0];
rz(2.5719602) q[1];
sx q[1];
rz(-0.49915794) q[1];
sx q[1];
rz(2.6380576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.106905) q[0];
sx q[0];
rz(-1.4464753) q[0];
sx q[0];
rz(0.19363815) q[0];
rz(-pi) q[1];
rz(0.72959186) q[2];
sx q[2];
rz(-2.1520832) q[2];
sx q[2];
rz(0.94129291) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.40696661) q[1];
sx q[1];
rz(-1.5583724) q[1];
sx q[1];
rz(1.4895275) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3713417) q[3];
sx q[3];
rz(-1.3512843) q[3];
sx q[3];
rz(1.7444004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2786431) q[2];
sx q[2];
rz(-2.7549665) q[2];
sx q[2];
rz(-1.6031727) q[2];
rz(-0.1362416) q[3];
sx q[3];
rz(-1.4038266) q[3];
sx q[3];
rz(-1.2374102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513026) q[0];
sx q[0];
rz(-2.5232115) q[0];
sx q[0];
rz(2.8001617) q[0];
rz(0.731172) q[1];
sx q[1];
rz(-0.57594222) q[1];
sx q[1];
rz(2.0700571) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7853189) q[0];
sx q[0];
rz(-0.59138227) q[0];
sx q[0];
rz(-0.7749556) q[0];
rz(3.0546636) q[2];
sx q[2];
rz(-2.5247716) q[2];
sx q[2];
rz(1.7380184) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2372969) q[1];
sx q[1];
rz(-1.132442) q[1];
sx q[1];
rz(2.3890952) q[1];
rz(-pi) q[2];
rz(-3.1370779) q[3];
sx q[3];
rz(-2.0751953) q[3];
sx q[3];
rz(2.2405091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7858481) q[2];
sx q[2];
rz(-1.5037437) q[2];
sx q[2];
rz(-0.49622932) q[2];
rz(-1.548454) q[3];
sx q[3];
rz(-1.691247) q[3];
sx q[3];
rz(2.1350071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9736495) q[0];
sx q[0];
rz(-1.4730467) q[0];
sx q[0];
rz(-3.1267401) q[0];
rz(1.7282093) q[1];
sx q[1];
rz(-1.1057066) q[1];
sx q[1];
rz(2.3766439) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0524372) q[0];
sx q[0];
rz(-0.90298072) q[0];
sx q[0];
rz(-2.3994156) q[0];
rz(-1.9221481) q[2];
sx q[2];
rz(-2.2235011) q[2];
sx q[2];
rz(1.2437133) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97100869) q[1];
sx q[1];
rz(-0.39675823) q[1];
sx q[1];
rz(0.011406974) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21723218) q[3];
sx q[3];
rz(-1.1558371) q[3];
sx q[3];
rz(0.46997786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0726274) q[2];
sx q[2];
rz(-1.5494538) q[2];
sx q[2];
rz(2.3907982) q[2];
rz(-2.1454861) q[3];
sx q[3];
rz(-0.80544296) q[3];
sx q[3];
rz(2.6944845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.65425101) q[0];
sx q[0];
rz(-0.60289201) q[0];
sx q[0];
rz(2.531429) q[0];
rz(-0.24972406) q[1];
sx q[1];
rz(-1.4796673) q[1];
sx q[1];
rz(3.0009559) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6941255) q[0];
sx q[0];
rz(-0.84703723) q[0];
sx q[0];
rz(-2.1985017) q[0];
rz(-pi) q[1];
rz(-2.6360077) q[2];
sx q[2];
rz(-1.4806432) q[2];
sx q[2];
rz(-0.32109993) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.096123213) q[1];
sx q[1];
rz(-0.28633196) q[1];
sx q[1];
rz(-2.8827122) q[1];
x q[2];
rz(-0.76042995) q[3];
sx q[3];
rz(-1.7685585) q[3];
sx q[3];
rz(0.20748479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3934624) q[2];
sx q[2];
rz(-0.40337864) q[2];
sx q[2];
rz(2.051579) q[2];
rz(1.9213093) q[3];
sx q[3];
rz(-2.0956495) q[3];
sx q[3];
rz(2.2530341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1423993) q[0];
sx q[0];
rz(-1.1141454) q[0];
sx q[0];
rz(-0.92495579) q[0];
rz(-1.6389182) q[1];
sx q[1];
rz(-1.0097367) q[1];
sx q[1];
rz(-0.22770539) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0696237) q[0];
sx q[0];
rz(-1.6771131) q[0];
sx q[0];
rz(-1.0384472) q[0];
rz(-pi) q[1];
rz(-0.41267379) q[2];
sx q[2];
rz(-2.5209171) q[2];
sx q[2];
rz(0.35457573) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5081577) q[1];
sx q[1];
rz(-0.40649182) q[1];
sx q[1];
rz(-2.6395931) q[1];
rz(-0.27276178) q[3];
sx q[3];
rz(-2.5099059) q[3];
sx q[3];
rz(-0.28608337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0324675) q[2];
sx q[2];
rz(-0.6296857) q[2];
sx q[2];
rz(-2.6510748) q[2];
rz(1.0444752) q[3];
sx q[3];
rz(-0.46382469) q[3];
sx q[3];
rz(-3.131955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5661285) q[0];
sx q[0];
rz(-0.45970356) q[0];
sx q[0];
rz(1.2658966) q[0];
rz(-1.6195541) q[1];
sx q[1];
rz(-1.8892037) q[1];
sx q[1];
rz(0.44580805) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4067526) q[0];
sx q[0];
rz(-1.8963843) q[0];
sx q[0];
rz(2.1181137) q[0];
x q[1];
rz(-1.61779) q[2];
sx q[2];
rz(-1.0596152) q[2];
sx q[2];
rz(1.5906478) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.5490732) q[1];
sx q[1];
rz(-2.5142094) q[1];
sx q[1];
rz(-2.907014) q[1];
rz(-pi) q[2];
rz(0.021546797) q[3];
sx q[3];
rz(-1.2974707) q[3];
sx q[3];
rz(-2.1692028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1294935) q[2];
sx q[2];
rz(-0.98126498) q[2];
sx q[2];
rz(2.1141466) q[2];
rz(-1.8217314) q[3];
sx q[3];
rz(-0.3052932) q[3];
sx q[3];
rz(-2.6125438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8158067) q[0];
sx q[0];
rz(-1.7740213) q[0];
sx q[0];
rz(-1.0446145) q[0];
rz(-2.2141937) q[1];
sx q[1];
rz(-1.2694042) q[1];
sx q[1];
rz(-0.70613695) q[1];
rz(-1.883222) q[2];
sx q[2];
rz(-1.2580522) q[2];
sx q[2];
rz(-3.0589041) q[2];
rz(1.2835684) q[3];
sx q[3];
rz(-0.97097048) q[3];
sx q[3];
rz(0.91965796) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
