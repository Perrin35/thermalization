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
rz(1.4229245) q[0];
sx q[0];
rz(-2.0473502) q[0];
sx q[0];
rz(-0.25804582) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(-1.300783) q[1];
sx q[1];
rz(0.25564495) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6545171) q[0];
sx q[0];
rz(-2.1571113) q[0];
sx q[0];
rz(-0.70588995) q[0];
rz(-pi) q[1];
rz(-1.5077816) q[2];
sx q[2];
rz(-1.1828701) q[2];
sx q[2];
rz(0.20252075) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2985743) q[1];
sx q[1];
rz(-2.9461423) q[1];
sx q[1];
rz(-1.9324383) q[1];
rz(-pi) q[2];
rz(-1.3575451) q[3];
sx q[3];
rz(-2.0590326) q[3];
sx q[3];
rz(-2.8247716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51097441) q[2];
sx q[2];
rz(-1.3773842) q[2];
sx q[2];
rz(-2.0261436) q[2];
rz(0.014178064) q[3];
sx q[3];
rz(-1.3445798) q[3];
sx q[3];
rz(0.43629638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2071335) q[0];
sx q[0];
rz(-2.5196228) q[0];
sx q[0];
rz(1.2472664) q[0];
rz(-0.44218749) q[1];
sx q[1];
rz(-1.7555883) q[1];
sx q[1];
rz(-0.8173379) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7572875) q[0];
sx q[0];
rz(-1.7624859) q[0];
sx q[0];
rz(-1.0377645) q[0];
x q[1];
rz(-1.1306612) q[2];
sx q[2];
rz(-1.3190184) q[2];
sx q[2];
rz(1.8735069) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8838447) q[1];
sx q[1];
rz(-2.2040329) q[1];
sx q[1];
rz(-1.3650989) q[1];
x q[2];
rz(-0.76089184) q[3];
sx q[3];
rz(-2.0756654) q[3];
sx q[3];
rz(0.99772108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77357972) q[2];
sx q[2];
rz(-2.7648338) q[2];
sx q[2];
rz(2.9212941) q[2];
rz(1.9593272) q[3];
sx q[3];
rz(-1.9672829) q[3];
sx q[3];
rz(-0.063145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.050506266) q[0];
sx q[0];
rz(-1.3171221) q[0];
sx q[0];
rz(-1.1385981) q[0];
rz(2.7298722) q[1];
sx q[1];
rz(-2.3717334) q[1];
sx q[1];
rz(-1.0134816) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0419755) q[0];
sx q[0];
rz(-1.3177773) q[0];
sx q[0];
rz(1.5554886) q[0];
rz(0.13229741) q[2];
sx q[2];
rz(-1.1331285) q[2];
sx q[2];
rz(0.79298151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.93597465) q[1];
sx q[1];
rz(-2.0009319) q[1];
sx q[1];
rz(0.79367743) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0115836) q[3];
sx q[3];
rz(-0.84835669) q[3];
sx q[3];
rz(0.60628451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98075214) q[2];
sx q[2];
rz(-1.1557121) q[2];
sx q[2];
rz(-1.2551003) q[2];
rz(-0.27215019) q[3];
sx q[3];
rz(-2.3641219) q[3];
sx q[3];
rz(0.14061558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.960152) q[0];
sx q[0];
rz(-0.93638268) q[0];
sx q[0];
rz(2.5312359) q[0];
rz(2.4329674) q[1];
sx q[1];
rz(-2.1253864) q[1];
sx q[1];
rz(1.5707387) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33922903) q[0];
sx q[0];
rz(-1.4177563) q[0];
sx q[0];
rz(1.7101076) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6897292) q[2];
sx q[2];
rz(-1.2858675) q[2];
sx q[2];
rz(-2.1708084) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6923234) q[1];
sx q[1];
rz(-0.27579112) q[1];
sx q[1];
rz(2.9640366) q[1];
x q[2];
rz(-1.1669772) q[3];
sx q[3];
rz(-2.32956) q[3];
sx q[3];
rz(2.5126575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2107971) q[2];
sx q[2];
rz(-2.1288629) q[2];
sx q[2];
rz(-2.5861758) q[2];
rz(-0.72426116) q[3];
sx q[3];
rz(-1.9925502) q[3];
sx q[3];
rz(-0.65653062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50437462) q[0];
sx q[0];
rz(-1.1906304) q[0];
sx q[0];
rz(2.7727238) q[0];
rz(0.24636191) q[1];
sx q[1];
rz(-1.8135704) q[1];
sx q[1];
rz(1.7074283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95043889) q[0];
sx q[0];
rz(-2.3847347) q[0];
sx q[0];
rz(-3.1340772) q[0];
x q[1];
rz(-0.35583115) q[2];
sx q[2];
rz(-0.47105481) q[2];
sx q[2];
rz(1.6182773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37232698) q[1];
sx q[1];
rz(-2.7397836) q[1];
sx q[1];
rz(0.010784464) q[1];
x q[2];
rz(-0.09999545) q[3];
sx q[3];
rz(-1.3444364) q[3];
sx q[3];
rz(-2.5842427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4904334) q[2];
sx q[2];
rz(-1.3110524) q[2];
sx q[2];
rz(2.0595713) q[2];
rz(-0.79484445) q[3];
sx q[3];
rz(-1.4719897) q[3];
sx q[3];
rz(1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.865888) q[0];
sx q[0];
rz(-3.0401433) q[0];
sx q[0];
rz(-0.24359447) q[0];
rz(-2.1547735) q[1];
sx q[1];
rz(-2.3934264) q[1];
sx q[1];
rz(2.2056244) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2631131) q[0];
sx q[0];
rz(-1.8917483) q[0];
sx q[0];
rz(-1.9892938) q[0];
rz(-pi) q[1];
rz(-0.34910874) q[2];
sx q[2];
rz(-0.85452467) q[2];
sx q[2];
rz(-0.30308613) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82560357) q[1];
sx q[1];
rz(-2.2554419) q[1];
sx q[1];
rz(2.0988093) q[1];
rz(-3.1201911) q[3];
sx q[3];
rz(-0.97107065) q[3];
sx q[3];
rz(0.76364005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.662107) q[2];
sx q[2];
rz(-2.0027436) q[2];
sx q[2];
rz(2.7916419) q[2];
rz(-2.5838666) q[3];
sx q[3];
rz(-1.9316659) q[3];
sx q[3];
rz(0.85132712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4418942) q[0];
sx q[0];
rz(-0.63755578) q[0];
sx q[0];
rz(-0.30174524) q[0];
rz(0.31983495) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(-2.005827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58590305) q[0];
sx q[0];
rz(-1.6246038) q[0];
sx q[0];
rz(-0.64482989) q[0];
rz(-pi) q[1];
rz(0.91355027) q[2];
sx q[2];
rz(-0.74591178) q[2];
sx q[2];
rz(-0.28829703) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7702076) q[1];
sx q[1];
rz(-1.4911629) q[1];
sx q[1];
rz(2.1712028) q[1];
rz(-pi) q[2];
rz(1.6834176) q[3];
sx q[3];
rz(-0.78312342) q[3];
sx q[3];
rz(2.6118731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7319506) q[2];
sx q[2];
rz(-2.4739517) q[2];
sx q[2];
rz(-2.9988334) q[2];
rz(-1.7536633) q[3];
sx q[3];
rz(-1.1889435) q[3];
sx q[3];
rz(2.4587542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9614354) q[0];
sx q[0];
rz(-3.1249983) q[0];
sx q[0];
rz(2.5855682) q[0];
rz(-0.22932886) q[1];
sx q[1];
rz(-1.2622958) q[1];
sx q[1];
rz(1.75846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.589298) q[0];
sx q[0];
rz(-2.4032666) q[0];
sx q[0];
rz(3.0895322) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6830693) q[2];
sx q[2];
rz(-1.5178198) q[2];
sx q[2];
rz(0.76901877) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.30547678) q[1];
sx q[1];
rz(-0.51199847) q[1];
sx q[1];
rz(2.192537) q[1];
rz(-pi) q[2];
rz(-1.1874299) q[3];
sx q[3];
rz(-1.995242) q[3];
sx q[3];
rz(-0.92680537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51656276) q[2];
sx q[2];
rz(-2.8690858) q[2];
sx q[2];
rz(-1.2410835) q[2];
rz(0.062601335) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(0.82714287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1271707) q[0];
sx q[0];
rz(-1.4143455) q[0];
sx q[0];
rz(-2.993809) q[0];
rz(-2.6328909) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(-2.7630189) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0914144) q[0];
sx q[0];
rz(-1.1568406) q[0];
sx q[0];
rz(1.3148091) q[0];
rz(-pi) q[1];
rz(-0.97855391) q[2];
sx q[2];
rz(-1.2662953) q[2];
sx q[2];
rz(-2.1155807) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55146688) q[1];
sx q[1];
rz(-2.6762137) q[1];
sx q[1];
rz(-0.82222934) q[1];
x q[2];
rz(0.76102961) q[3];
sx q[3];
rz(-1.6609523) q[3];
sx q[3];
rz(-2.6554299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1754237) q[2];
sx q[2];
rz(-0.95169008) q[2];
sx q[2];
rz(1.2790722) q[2];
rz(-1.4281979) q[3];
sx q[3];
rz(-1.0019852) q[3];
sx q[3];
rz(0.14555791) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3928423) q[0];
sx q[0];
rz(-0.57806438) q[0];
sx q[0];
rz(0.064099126) q[0];
rz(-0.52866689) q[1];
sx q[1];
rz(-1.8730947) q[1];
sx q[1];
rz(0.97506964) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8846183) q[0];
sx q[0];
rz(-1.8351089) q[0];
sx q[0];
rz(-2.6340911) q[0];
rz(2.3401572) q[2];
sx q[2];
rz(-0.74432997) q[2];
sx q[2];
rz(-2.3542885) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9916315) q[1];
sx q[1];
rz(-1.35222) q[1];
sx q[1];
rz(0.44709713) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24419489) q[3];
sx q[3];
rz(-1.7932442) q[3];
sx q[3];
rz(-2.4110068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9958682) q[2];
sx q[2];
rz(-0.66103649) q[2];
sx q[2];
rz(-2.5346942) q[2];
rz(2.8912344) q[3];
sx q[3];
rz(-1.4162049) q[3];
sx q[3];
rz(2.6355766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.09457) q[0];
sx q[0];
rz(-1.2170412) q[0];
sx q[0];
rz(1.8893597) q[0];
rz(0.49867123) q[1];
sx q[1];
rz(-1.55232) q[1];
sx q[1];
rz(1.3943863) q[1];
rz(2.5468536) q[2];
sx q[2];
rz(-2.1872216) q[2];
sx q[2];
rz(-2.8333153) q[2];
rz(2.102643) q[3];
sx q[3];
rz(-2.9038351) q[3];
sx q[3];
rz(-2.6701715) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
