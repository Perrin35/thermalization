OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3988848) q[0];
sx q[0];
rz(-2.3595915) q[0];
sx q[0];
rz(-1.8703823) q[0];
rz(0.27702364) q[1];
sx q[1];
rz(-0.47200534) q[1];
sx q[1];
rz(-3.1402052) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1479552) q[0];
sx q[0];
rz(-1.4076828) q[0];
sx q[0];
rz(1.9440584) q[0];
x q[1];
rz(0.44714655) q[2];
sx q[2];
rz(-1.6086173) q[2];
sx q[2];
rz(-2.56074) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5056155) q[1];
sx q[1];
rz(-0.50230366) q[1];
sx q[1];
rz(0.14640267) q[1];
rz(-pi) q[2];
rz(1.1997585) q[3];
sx q[3];
rz(-1.3317809) q[3];
sx q[3];
rz(1.0306851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9871621) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(-0.74938613) q[2];
rz(-2.1253712) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(0.40482503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4352903) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-0.97066561) q[0];
rz(-2.1043815) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(-0.81545365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6204651) q[0];
sx q[0];
rz(-1.5107811) q[0];
sx q[0];
rz(-0.64882664) q[0];
rz(-pi) q[1];
rz(2.4120861) q[2];
sx q[2];
rz(-1.365005) q[2];
sx q[2];
rz(1.4480928) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26317877) q[1];
sx q[1];
rz(-2.8385332) q[1];
sx q[1];
rz(3.0221699) q[1];
rz(-1.3482413) q[3];
sx q[3];
rz(-2.2046304) q[3];
sx q[3];
rz(2.8702877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4619535) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(-0.63278502) q[2];
rz(-1.9880382) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(-2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(-0.87483037) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(-0.99951807) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8529352) q[0];
sx q[0];
rz(-1.7203727) q[0];
sx q[0];
rz(1.268671) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61727662) q[2];
sx q[2];
rz(-1.0059788) q[2];
sx q[2];
rz(1.9895983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1516583) q[1];
sx q[1];
rz(-2.5807568) q[1];
sx q[1];
rz(2.6480617) q[1];
rz(-3.0619377) q[3];
sx q[3];
rz(-1.0632535) q[3];
sx q[3];
rz(1.5505276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6040566) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(-0.58004722) q[2];
rz(2.3245658) q[3];
sx q[3];
rz(-1.3823119) q[3];
sx q[3];
rz(-1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(2.2312009) q[0];
rz(-2.6903649) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(0.27483637) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9592322) q[0];
sx q[0];
rz(-1.649547) q[0];
sx q[0];
rz(-1.4819281) q[0];
rz(-pi) q[1];
rz(-1.5675797) q[2];
sx q[2];
rz(-2.6811757) q[2];
sx q[2];
rz(-2.761063) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7819314) q[1];
sx q[1];
rz(-1.1210103) q[1];
sx q[1];
rz(0.91314544) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0148388) q[3];
sx q[3];
rz(-2.2634014) q[3];
sx q[3];
rz(2.9664489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.794902) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(-0.16170734) q[3];
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
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4836924) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(2.143798) q[0];
rz(0.18355852) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(-1.6246187) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51380101) q[0];
sx q[0];
rz(-2.1876946) q[0];
sx q[0];
rz(-0.4517201) q[0];
rz(-pi) q[1];
rz(2.5324608) q[2];
sx q[2];
rz(-2.5172148) q[2];
sx q[2];
rz(0.49027157) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9285674) q[1];
sx q[1];
rz(-1.9145609) q[1];
sx q[1];
rz(-0.43069559) q[1];
x q[2];
rz(-2.3715641) q[3];
sx q[3];
rz(-0.75776811) q[3];
sx q[3];
rz(-3.1185574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(0.3240164) q[2];
rz(1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(-1.2639686) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(2.8009159) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9522889) q[0];
sx q[0];
rz(-1.6366819) q[0];
sx q[0];
rz(-0.031866372) q[0];
rz(0.31008115) q[2];
sx q[2];
rz(-2.3611464) q[2];
sx q[2];
rz(-2.9030637) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1250455) q[1];
sx q[1];
rz(-0.53495896) q[1];
sx q[1];
rz(-0.46334456) q[1];
x q[2];
rz(-0.42231456) q[3];
sx q[3];
rz(-2.0229682) q[3];
sx q[3];
rz(0.57002588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.95057758) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(2.3699956) q[2];
rz(-2.5937882) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9724378) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(-0.34564885) q[0];
rz(-0.06282839) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(-2.6766434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15360399) q[0];
sx q[0];
rz(-2.0360332) q[0];
sx q[0];
rz(2.9265755) q[0];
x q[1];
rz(-0.99636997) q[2];
sx q[2];
rz(-1.9949706) q[2];
sx q[2];
rz(3.0976354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.834224) q[1];
sx q[1];
rz(-2.8982179) q[1];
sx q[1];
rz(2.6878396) q[1];
rz(3.1157324) q[3];
sx q[3];
rz(-1.5246632) q[3];
sx q[3];
rz(-2.9010454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0039625) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(2.8239992) q[2];
rz(-0.57146227) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(-2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5193609) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(2.8572594) q[0];
rz(2.590086) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(-0.078358738) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1348304) q[0];
sx q[0];
rz(-1.6400596) q[0];
sx q[0];
rz(2.3368821) q[0];
rz(-pi) q[1];
rz(2.9044754) q[2];
sx q[2];
rz(-1.8097005) q[2];
sx q[2];
rz(-0.6616) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2625418) q[1];
sx q[1];
rz(-1.3971551) q[1];
sx q[1];
rz(2.2294728) q[1];
rz(-pi) q[2];
rz(3.1073242) q[3];
sx q[3];
rz(-1.4141603) q[3];
sx q[3];
rz(-0.53965118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4006965) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(0.25137869) q[2];
rz(-0.58319432) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(1.4096227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79779977) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(0.051368512) q[0];
rz(-0.92357606) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-0.87402469) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6286205) q[0];
sx q[0];
rz(-2.3734833) q[0];
sx q[0];
rz(0.012499768) q[0];
rz(-0.36058493) q[2];
sx q[2];
rz(-1.1869831) q[2];
sx q[2];
rz(-1.4521445) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8928788) q[1];
sx q[1];
rz(-1.6422179) q[1];
sx q[1];
rz(2.0715269) q[1];
x q[2];
rz(-0.60871082) q[3];
sx q[3];
rz(-1.8551644) q[3];
sx q[3];
rz(-1.9539208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.727227) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(0.61974636) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(1.363389) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0062362) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(-0.18877098) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.1788517) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2154685) q[0];
sx q[0];
rz(-1.3660396) q[0];
sx q[0];
rz(-1.2925076) q[0];
x q[1];
rz(-0.63209052) q[2];
sx q[2];
rz(-2.9846016) q[2];
sx q[2];
rz(-1.1319515) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5621592) q[1];
sx q[1];
rz(-2.4907618) q[1];
sx q[1];
rz(1.8821554) q[1];
x q[2];
rz(0.3801109) q[3];
sx q[3];
rz(-2.5862525) q[3];
sx q[3];
rz(-1.4243319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0835691) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(0.87456885) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(-2.3836366) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(1.6451251) q[2];
sx q[2];
rz(-2.7004514) q[2];
sx q[2];
rz(-2.2218291) q[2];
rz(-1.422613) q[3];
sx q[3];
rz(-1.8437181) q[3];
sx q[3];
rz(-3.0317882) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
