OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0410864) q[0];
sx q[0];
rz(-0.20809986) q[0];
sx q[0];
rz(-0.40325525) q[0];
rz(2.9085605) q[1];
sx q[1];
rz(-1.7014528) q[1];
sx q[1];
rz(-0.22388248) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1772645) q[0];
sx q[0];
rz(-1.3601662) q[0];
sx q[0];
rz(-0.62646477) q[0];
rz(-2.2751132) q[2];
sx q[2];
rz(-1.472855) q[2];
sx q[2];
rz(0.78664727) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63131166) q[1];
sx q[1];
rz(-1.3461539) q[1];
sx q[1];
rz(-0.90490492) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81108629) q[3];
sx q[3];
rz(-1.038365) q[3];
sx q[3];
rz(0.48871751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.69748059) q[2];
sx q[2];
rz(-2.8287502) q[2];
sx q[2];
rz(0.035813896) q[2];
rz(2.4114285) q[3];
sx q[3];
rz(-1.9883479) q[3];
sx q[3];
rz(-0.20619503) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8908454) q[0];
sx q[0];
rz(-0.99531168) q[0];
sx q[0];
rz(0.19749755) q[0];
rz(-2.6644871) q[1];
sx q[1];
rz(-2.4767866) q[1];
sx q[1];
rz(-0.8055996) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44551906) q[0];
sx q[0];
rz(-1.7564) q[0];
sx q[0];
rz(-2.473102) q[0];
rz(2.9764011) q[2];
sx q[2];
rz(-0.38536638) q[2];
sx q[2];
rz(-1.2549653) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9220613) q[1];
sx q[1];
rz(-2.1967255) q[1];
sx q[1];
rz(2.3971167) q[1];
rz(-pi) q[2];
rz(-2.9904705) q[3];
sx q[3];
rz(-1.7100705) q[3];
sx q[3];
rz(0.81859156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3104559) q[2];
sx q[2];
rz(-1.4996935) q[2];
sx q[2];
rz(-1.7384701) q[2];
rz(0.63203114) q[3];
sx q[3];
rz(-2.8681614) q[3];
sx q[3];
rz(0.12701756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(1.8359351) q[0];
sx q[0];
rz(-2.5129565) q[0];
sx q[0];
rz(2.054731) q[0];
rz(-0.19733363) q[1];
sx q[1];
rz(-1.5060164) q[1];
sx q[1];
rz(2.8089583) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4950748) q[0];
sx q[0];
rz(-1.7103638) q[0];
sx q[0];
rz(-2.0049176) q[0];
rz(-pi) q[1];
rz(1.640075) q[2];
sx q[2];
rz(-2.4054619) q[2];
sx q[2];
rz(-2.2683805) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8197934) q[1];
sx q[1];
rz(-1.7463356) q[1];
sx q[1];
rz(-1.790253) q[1];
rz(0.25742297) q[3];
sx q[3];
rz(-0.74353131) q[3];
sx q[3];
rz(-0.0071255077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3745594) q[2];
sx q[2];
rz(-1.550753) q[2];
sx q[2];
rz(-1.4683051) q[2];
rz(-0.0091920216) q[3];
sx q[3];
rz(-2.6720948) q[3];
sx q[3];
rz(2.3495242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62995768) q[0];
sx q[0];
rz(-0.63794962) q[0];
sx q[0];
rz(1.6988423) q[0];
rz(1.0244145) q[1];
sx q[1];
rz(-1.9099648) q[1];
sx q[1];
rz(-2.3146497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8222745) q[0];
sx q[0];
rz(-1.3869023) q[0];
sx q[0];
rz(3.0660423) q[0];
rz(-2.9873136) q[2];
sx q[2];
rz(-0.88866975) q[2];
sx q[2];
rz(-0.23274225) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6232439) q[1];
sx q[1];
rz(-2.001013) q[1];
sx q[1];
rz(-1.3395549) q[1];
rz(-pi) q[2];
rz(2.2531281) q[3];
sx q[3];
rz(-1.4766221) q[3];
sx q[3];
rz(1.3741796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3247165) q[2];
sx q[2];
rz(-2.3574895) q[2];
sx q[2];
rz(1.8640222) q[2];
rz(0.68814021) q[3];
sx q[3];
rz(-1.0253996) q[3];
sx q[3];
rz(-2.1791606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14030305) q[0];
sx q[0];
rz(-1.7004509) q[0];
sx q[0];
rz(2.9241614) q[0];
rz(-2.8561719) q[1];
sx q[1];
rz(-2.5358584) q[1];
sx q[1];
rz(2.6434456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0131684) q[0];
sx q[0];
rz(-1.8966872) q[0];
sx q[0];
rz(-0.89114916) q[0];
rz(-pi) q[1];
rz(-2.010538) q[2];
sx q[2];
rz(-0.56098191) q[2];
sx q[2];
rz(-1.2351241) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5707055) q[1];
sx q[1];
rz(-2.0556309) q[1];
sx q[1];
rz(-0.2356694) q[1];
rz(-1.3181456) q[3];
sx q[3];
rz(-1.9985191) q[3];
sx q[3];
rz(-1.778217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.51776) q[2];
sx q[2];
rz(-2.0443003) q[2];
sx q[2];
rz(2.9916054) q[2];
rz(0.013966694) q[3];
sx q[3];
rz(-1.59168) q[3];
sx q[3];
rz(-0.51378957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49195313) q[0];
sx q[0];
rz(-2.5998901) q[0];
sx q[0];
rz(2.3288222) q[0];
rz(-1.7236408) q[1];
sx q[1];
rz(-1.6287454) q[1];
sx q[1];
rz(2.0779804) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9169711) q[0];
sx q[0];
rz(-1.2422891) q[0];
sx q[0];
rz(-2.5909831) q[0];
x q[1];
rz(2.9779149) q[2];
sx q[2];
rz(-1.5874169) q[2];
sx q[2];
rz(1.4281685) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.221595) q[1];
sx q[1];
rz(-2.1694899) q[1];
sx q[1];
rz(-1.5845951) q[1];
rz(1.316458) q[3];
sx q[3];
rz(-1.7539548) q[3];
sx q[3];
rz(-1.1381799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.873988) q[2];
sx q[2];
rz(-2.5817817) q[2];
sx q[2];
rz(1.7248636) q[2];
rz(-3.0806165) q[3];
sx q[3];
rz(-2.2305326) q[3];
sx q[3];
rz(-1.6772259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7199719) q[0];
sx q[0];
rz(-0.82141972) q[0];
sx q[0];
rz(-0.11909568) q[0];
rz(0.47856092) q[1];
sx q[1];
rz(-0.81399545) q[1];
sx q[1];
rz(-0.68971577) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0937492) q[0];
sx q[0];
rz(-0.99573538) q[0];
sx q[0];
rz(2.439365) q[0];
rz(-0.63197559) q[2];
sx q[2];
rz(-2.5594829) q[2];
sx q[2];
rz(1.4402953) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8552095) q[1];
sx q[1];
rz(-0.70806009) q[1];
sx q[1];
rz(-1.5016012) q[1];
x q[2];
rz(-0.85817091) q[3];
sx q[3];
rz(-1.6827876) q[3];
sx q[3];
rz(1.6046765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0488284) q[2];
sx q[2];
rz(-3.0807639) q[2];
sx q[2];
rz(2.0737341) q[2];
rz(2.974406) q[3];
sx q[3];
rz(-1.7696295) q[3];
sx q[3];
rz(-0.13944496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-1.6086513) q[0];
sx q[0];
rz(-2.3280188) q[0];
sx q[0];
rz(2.2331878) q[0];
rz(0.99336973) q[1];
sx q[1];
rz(-0.16255957) q[1];
sx q[1];
rz(2.2672674) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4607166) q[0];
sx q[0];
rz(-1.5534288) q[0];
sx q[0];
rz(3.0094023) q[0];
rz(-pi) q[1];
x q[1];
rz(0.019432391) q[2];
sx q[2];
rz(-1.3131058) q[2];
sx q[2];
rz(-1.9382221) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.43046865) q[1];
sx q[1];
rz(-1.5114771) q[1];
sx q[1];
rz(1.7064246) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7818412) q[3];
sx q[3];
rz(-1.0488516) q[3];
sx q[3];
rz(-0.2019384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5270762) q[2];
sx q[2];
rz(-0.15915844) q[2];
sx q[2];
rz(2.2873774) q[2];
rz(2.6863875) q[3];
sx q[3];
rz(-2.5762317) q[3];
sx q[3];
rz(0.2555041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47827569) q[0];
sx q[0];
rz(-1.9238967) q[0];
sx q[0];
rz(-1.4772557) q[0];
rz(-1.5669589) q[1];
sx q[1];
rz(-0.68395558) q[1];
sx q[1];
rz(2.4050567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7304717) q[0];
sx q[0];
rz(-0.53148848) q[0];
sx q[0];
rz(-1.0941675) q[0];
rz(-pi) q[1];
rz(-1.4669424) q[2];
sx q[2];
rz(-1.1919293) q[2];
sx q[2];
rz(1.1844289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.354035) q[1];
sx q[1];
rz(-1.8144467) q[1];
sx q[1];
rz(2.8188348) q[1];
x q[2];
rz(-3.1062713) q[3];
sx q[3];
rz(-0.97408726) q[3];
sx q[3];
rz(-2.2281856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3972724) q[2];
sx q[2];
rz(-2.263415) q[2];
sx q[2];
rz(1.6142023) q[2];
rz(-0.39198908) q[3];
sx q[3];
rz(-1.1992998) q[3];
sx q[3];
rz(3.0919302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.117332) q[0];
sx q[0];
rz(-1.4735104) q[0];
sx q[0];
rz(0.20183739) q[0];
rz(-2.9182538) q[1];
sx q[1];
rz(-2.6170862) q[1];
sx q[1];
rz(2.7490659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46014547) q[0];
sx q[0];
rz(-1.1090312) q[0];
sx q[0];
rz(0.11790922) q[0];
rz(-pi) q[1];
rz(1.5129651) q[2];
sx q[2];
rz(-2.4041345) q[2];
sx q[2];
rz(-2.2493169) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0328249) q[1];
sx q[1];
rz(-2.4420173) q[1];
sx q[1];
rz(-0.8731858) q[1];
rz(-0.48247561) q[3];
sx q[3];
rz(-1.8419208) q[3];
sx q[3];
rz(1.9737873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4203804) q[2];
sx q[2];
rz(-0.98782295) q[2];
sx q[2];
rz(2.5610899) q[2];
rz(2.765559) q[3];
sx q[3];
rz(-2.1707363) q[3];
sx q[3];
rz(1.6002801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9911983) q[0];
sx q[0];
rz(-1.5924441) q[0];
sx q[0];
rz(1.7572255) q[0];
rz(-1.6509542) q[1];
sx q[1];
rz(-1.2532267) q[1];
sx q[1];
rz(0.86029235) q[1];
rz(-2.2223086) q[2];
sx q[2];
rz(-2.2625661) q[2];
sx q[2];
rz(1.7023466) q[2];
rz(0.31181351) q[3];
sx q[3];
rz(-1.5494294) q[3];
sx q[3];
rz(2.6904306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
