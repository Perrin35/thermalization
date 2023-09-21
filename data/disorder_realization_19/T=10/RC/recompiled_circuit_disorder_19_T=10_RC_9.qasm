OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0269545) q[0];
sx q[0];
rz(-1.6898328) q[0];
sx q[0];
rz(0.64557689) q[0];
rz(0.37880701) q[1];
sx q[1];
rz(4.9103476) q[1];
sx q[1];
rz(11.068439) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49255532) q[0];
sx q[0];
rz(-1.7151924) q[0];
sx q[0];
rz(-1.210968) q[0];
rz(-pi) q[1];
x q[1];
rz(1.475392) q[2];
sx q[2];
rz(-2.5089052) q[2];
sx q[2];
rz(2.9205517) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5421996) q[1];
sx q[1];
rz(-1.3560055) q[1];
sx q[1];
rz(2.3989124) q[1];
rz(-2.492339) q[3];
sx q[3];
rz(-1.89562) q[3];
sx q[3];
rz(2.0089825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2215185) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(-2.8011838) q[2];
rz(0.83299625) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67265636) q[0];
sx q[0];
rz(-2.7514973) q[0];
sx q[0];
rz(2.5966068) q[0];
rz(-2.2333721) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(2.3166336) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93281125) q[0];
sx q[0];
rz(-2.6926846) q[0];
sx q[0];
rz(-2.0138028) q[0];
rz(-1.5603793) q[2];
sx q[2];
rz(-1.1679107) q[2];
sx q[2];
rz(-1.252623) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1383914) q[1];
sx q[1];
rz(-1.0111286) q[1];
sx q[1];
rz(-1.4336587) q[1];
rz(-pi) q[2];
rz(2.1373848) q[3];
sx q[3];
rz(-1.7363747) q[3];
sx q[3];
rz(0.86629888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5543582) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(2.9651802) q[2];
rz(0.13088626) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(0.46974716) q[0];
rz(-1.6167971) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(0.038539561) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8174724) q[0];
sx q[0];
rz(-1.5714374) q[0];
sx q[0];
rz(-1.5630432) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23080319) q[2];
sx q[2];
rz(-1.7204086) q[2];
sx q[2];
rz(-2.9485867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5765499) q[1];
sx q[1];
rz(-1.0817263) q[1];
sx q[1];
rz(3.0719041) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3452957) q[3];
sx q[3];
rz(-2.3781804) q[3];
sx q[3];
rz(0.89087668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5588351) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(-2.2290686) q[2];
rz(-1.8330666) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38347605) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(0.34969774) q[0];
rz(-1.8967459) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(-2.9464338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0168266) q[0];
sx q[0];
rz(-0.15206465) q[0];
sx q[0];
rz(2.2269339) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48093502) q[2];
sx q[2];
rz(-1.0785042) q[2];
sx q[2];
rz(0.42602691) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2745167) q[1];
sx q[1];
rz(-1.3356326) q[1];
sx q[1];
rz(-2.0348674) q[1];
rz(-pi) q[2];
rz(1.5870729) q[3];
sx q[3];
rz(-2.2098594) q[3];
sx q[3];
rz(-1.0234327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23094709) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(-1.8738497) q[2];
rz(-1.0686482) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0347663) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(0.1396133) q[0];
rz(-1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(2.7635014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75487126) q[0];
sx q[0];
rz(-2.5140962) q[0];
sx q[0];
rz(1.0794712) q[0];
x q[1];
rz(0.81850448) q[2];
sx q[2];
rz(-2.2193925) q[2];
sx q[2];
rz(-1.6550145) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0542478) q[1];
sx q[1];
rz(-1.9964881) q[1];
sx q[1];
rz(2.737983) q[1];
rz(-pi) q[2];
rz(-0.46697163) q[3];
sx q[3];
rz(-1.3936371) q[3];
sx q[3];
rz(2.5339479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87171626) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(-2.2536229) q[2];
rz(-2.1652083) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(-2.1358657) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3865005) q[0];
sx q[0];
rz(-0.077682406) q[0];
sx q[0];
rz(-2.7639672) q[0];
rz(0.32304421) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(-2.2264218) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8422896) q[0];
sx q[0];
rz(-1.8931086) q[0];
sx q[0];
rz(1.4847941) q[0];
x q[1];
rz(0.50778163) q[2];
sx q[2];
rz(-2.0785329) q[2];
sx q[2];
rz(2.7024262) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8296216) q[1];
sx q[1];
rz(-1.8835856) q[1];
sx q[1];
rz(-1.642986) q[1];
rz(1.4573426) q[3];
sx q[3];
rz(-0.87943422) q[3];
sx q[3];
rz(-1.7951579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6043828) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(2.3510695) q[2];
rz(0.28997713) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.73615605) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(-1.908318) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(1.4253915) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7084301) q[0];
sx q[0];
rz(-1.7571974) q[0];
sx q[0];
rz(2.7011046) q[0];
rz(-0.44600365) q[2];
sx q[2];
rz(-1.6288695) q[2];
sx q[2];
rz(0.29689483) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.094147625) q[1];
sx q[1];
rz(-0.74456753) q[1];
sx q[1];
rz(-2.5658957) q[1];
rz(-pi) q[2];
rz(-0.013399259) q[3];
sx q[3];
rz(-2.074476) q[3];
sx q[3];
rz(-2.3434533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9795064) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(1.210775) q[2];
rz(-0.0018421729) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2974671) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-3.0650744) q[0];
rz(-0.67529768) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(-0.75235596) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233326) q[0];
sx q[0];
rz(-1.8857297) q[0];
sx q[0];
rz(0.53091913) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2194527) q[2];
sx q[2];
rz(-1.4553918) q[2];
sx q[2];
rz(2.5285072) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89531089) q[1];
sx q[1];
rz(-1.4824251) q[1];
sx q[1];
rz(-1.6842711) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5520677) q[3];
sx q[3];
rz(-1.7426963) q[3];
sx q[3];
rz(-0.37026438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.9383119) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(3.1414462) q[2];
rz(-1.1095307) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9922441) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(2.1355656) q[0];
rz(0.26793119) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(-1.3148274) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49004236) q[0];
sx q[0];
rz(-2.393912) q[0];
sx q[0];
rz(-0.3726532) q[0];
rz(-pi) q[1];
rz(0.36395145) q[2];
sx q[2];
rz(-1.4462573) q[2];
sx q[2];
rz(-3.0968015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.41373738) q[1];
sx q[1];
rz(-1.1579885) q[1];
sx q[1];
rz(2.2933943) q[1];
rz(-pi) q[2];
rz(0.0074578961) q[3];
sx q[3];
rz(-2.3141626) q[3];
sx q[3];
rz(-2.5936562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7018147) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(2.8961199) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4196639) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(-0.87896705) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(0.39189664) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.648801) q[0];
sx q[0];
rz(-1.0133044) q[0];
sx q[0];
rz(3.1317943) q[0];
rz(-pi) q[1];
x q[1];
rz(1.395412) q[2];
sx q[2];
rz(-2.4782964) q[2];
sx q[2];
rz(-1.7593918) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7614903) q[1];
sx q[1];
rz(-0.83161608) q[1];
sx q[1];
rz(-2.889061) q[1];
rz(-0.55962555) q[3];
sx q[3];
rz(-1.8436173) q[3];
sx q[3];
rz(0.19459693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5320756) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(-1.1575451) q[2];
rz(-0.55084294) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3863603) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(-1.339284) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(-0.77970589) q[2];
sx q[2];
rz(-1.69366) q[2];
sx q[2];
rz(2.2488307) q[2];
rz(2.1035351) q[3];
sx q[3];
rz(-2.5235015) q[3];
sx q[3];
rz(-1.5550767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
