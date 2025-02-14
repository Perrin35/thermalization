OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7336361) q[0];
sx q[0];
rz(-0.1853369) q[0];
sx q[0];
rz(1.4987401) q[0];
rz(2.9474131) q[1];
sx q[1];
rz(-2.7979538) q[1];
sx q[1];
rz(-2.763881) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.728017) q[0];
sx q[0];
rz(-2.8110162) q[0];
sx q[0];
rz(0.74329336) q[0];
rz(-pi) q[1];
x q[1];
rz(2.333856) q[2];
sx q[2];
rz(-0.57538549) q[2];
sx q[2];
rz(-2.6859716) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9895674) q[1];
sx q[1];
rz(-1.6049687) q[1];
sx q[1];
rz(-0.65712838) q[1];
rz(3.1123766) q[3];
sx q[3];
rz(-1.3061285) q[3];
sx q[3];
rz(-0.1614557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7517884) q[2];
sx q[2];
rz(-1.7916388) q[2];
sx q[2];
rz(-2.8765615) q[2];
rz(0.42936471) q[3];
sx q[3];
rz(-1.2484442) q[3];
sx q[3];
rz(-3.140669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65381831) q[0];
sx q[0];
rz(-0.14350292) q[0];
sx q[0];
rz(1.8408467) q[0];
rz(0.067151345) q[1];
sx q[1];
rz(-0.90922272) q[1];
sx q[1];
rz(0.86004177) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3305682) q[0];
sx q[0];
rz(-2.1386792) q[0];
sx q[0];
rz(-0.83356838) q[0];
x q[1];
rz(0.70921398) q[2];
sx q[2];
rz(-1.0150036) q[2];
sx q[2];
rz(1.1467939) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7886161) q[1];
sx q[1];
rz(-1.9771655) q[1];
sx q[1];
rz(-0.95384903) q[1];
rz(1.8710526) q[3];
sx q[3];
rz(-1.7466063) q[3];
sx q[3];
rz(2.6838357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.63090515) q[2];
sx q[2];
rz(-2.5248933) q[2];
sx q[2];
rz(2.5751233) q[2];
rz(-2.412292) q[3];
sx q[3];
rz(-2.2972378) q[3];
sx q[3];
rz(-2.9384889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1931964) q[0];
sx q[0];
rz(-1.0862792) q[0];
sx q[0];
rz(-2.7753944) q[0];
rz(-1.5858448) q[1];
sx q[1];
rz(-2.0987174) q[1];
sx q[1];
rz(-3.0911456) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9589466) q[0];
sx q[0];
rz(-0.0894657) q[0];
sx q[0];
rz(-2.0561809) q[0];
x q[1];
rz(2.5620805) q[2];
sx q[2];
rz(-0.77634927) q[2];
sx q[2];
rz(-3.1402328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37990268) q[1];
sx q[1];
rz(-2.8643048) q[1];
sx q[1];
rz(-0.65245858) q[1];
rz(-pi) q[2];
rz(-0.14536038) q[3];
sx q[3];
rz(-1.1499377) q[3];
sx q[3];
rz(-3.1042254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0744276) q[2];
sx q[2];
rz(-1.0600435) q[2];
sx q[2];
rz(1.3107497) q[2];
rz(1.358076) q[3];
sx q[3];
rz(-0.57099968) q[3];
sx q[3];
rz(-1.8324435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27291372) q[0];
sx q[0];
rz(-2.9153115) q[0];
sx q[0];
rz(1.7568461) q[0];
rz(1.9028496) q[1];
sx q[1];
rz(-1.9107995) q[1];
sx q[1];
rz(-1.967427) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5730302) q[0];
sx q[0];
rz(-1.996604) q[0];
sx q[0];
rz(-2.9338475) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8644845) q[2];
sx q[2];
rz(-1.8113675) q[2];
sx q[2];
rz(1.5391369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.13951835) q[1];
sx q[1];
rz(-1.0660831) q[1];
sx q[1];
rz(0.37481777) q[1];
rz(-pi) q[2];
rz(-0.87808319) q[3];
sx q[3];
rz(-1.2125848) q[3];
sx q[3];
rz(-0.75853759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43526402) q[2];
sx q[2];
rz(-1.3202983) q[2];
sx q[2];
rz(3.1026133) q[2];
rz(-0.6428166) q[3];
sx q[3];
rz(-1.0182074) q[3];
sx q[3];
rz(2.5207998) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3103127) q[0];
sx q[0];
rz(-1.0142925) q[0];
sx q[0];
rz(-3.1211299) q[0];
rz(-2.9488355) q[1];
sx q[1];
rz(-0.61734504) q[1];
sx q[1];
rz(1.9761168) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3912688) q[0];
sx q[0];
rz(-1.170715) q[0];
sx q[0];
rz(-2.4720936) q[0];
rz(-0.67365174) q[2];
sx q[2];
rz(-1.9314226) q[2];
sx q[2];
rz(1.1471105) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1736892) q[1];
sx q[1];
rz(-0.5040796) q[1];
sx q[1];
rz(-1.6126812) q[1];
rz(1.194077) q[3];
sx q[3];
rz(-2.93749) q[3];
sx q[3];
rz(-0.78777504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3141025) q[2];
sx q[2];
rz(-0.9919439) q[2];
sx q[2];
rz(2.6143383) q[2];
rz(2.5984247) q[3];
sx q[3];
rz(-0.68734622) q[3];
sx q[3];
rz(2.7988722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0813893) q[0];
sx q[0];
rz(-2.9855766) q[0];
sx q[0];
rz(-2.7348837) q[0];
rz(-1.1609062) q[1];
sx q[1];
rz(-0.91861594) q[1];
sx q[1];
rz(2.1312174) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.564931) q[0];
sx q[0];
rz(-1.0258249) q[0];
sx q[0];
rz(-0.14133639) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.819028) q[2];
sx q[2];
rz(-1.7062213) q[2];
sx q[2];
rz(2.5544006) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53394371) q[1];
sx q[1];
rz(-1.2022809) q[1];
sx q[1];
rz(-1.7979969) q[1];
rz(-1.3090735) q[3];
sx q[3];
rz(-0.78262586) q[3];
sx q[3];
rz(1.7109551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0747718) q[2];
sx q[2];
rz(-1.2819042) q[2];
sx q[2];
rz(2.1066966) q[2];
rz(-1.7847938) q[3];
sx q[3];
rz(-2.5467338) q[3];
sx q[3];
rz(2.518173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9676301) q[0];
sx q[0];
rz(-1.5246464) q[0];
sx q[0];
rz(-2.8136643) q[0];
rz(0.11218849) q[1];
sx q[1];
rz(-1.9393549) q[1];
sx q[1];
rz(2.1690538) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9754494) q[0];
sx q[0];
rz(-1.219141) q[0];
sx q[0];
rz(-0.97197284) q[0];
rz(-pi) q[1];
rz(-2.8424524) q[2];
sx q[2];
rz(-1.1465596) q[2];
sx q[2];
rz(-1.5628634) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9303846) q[1];
sx q[1];
rz(-0.58961419) q[1];
sx q[1];
rz(-3.1077887) q[1];
rz(-pi) q[2];
rz(2.6554606) q[3];
sx q[3];
rz(-2.4555169) q[3];
sx q[3];
rz(-3.0239575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5358676) q[2];
sx q[2];
rz(-0.87656993) q[2];
sx q[2];
rz(-2.4884339) q[2];
rz(0.43290916) q[3];
sx q[3];
rz(-2.631729) q[3];
sx q[3];
rz(0.21231095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9369649) q[0];
sx q[0];
rz(-2.300394) q[0];
sx q[0];
rz(2.9168108) q[0];
rz(-0.79554355) q[1];
sx q[1];
rz(-1.8276151) q[1];
sx q[1];
rz(2.0733817) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8248937) q[0];
sx q[0];
rz(-0.84352701) q[0];
sx q[0];
rz(-2.4466912) q[0];
x q[1];
rz(0.62309391) q[2];
sx q[2];
rz(-0.87298191) q[2];
sx q[2];
rz(2.0562003) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.28212122) q[1];
sx q[1];
rz(-0.88338125) q[1];
sx q[1];
rz(-1.6038546) q[1];
x q[2];
rz(-3.0826236) q[3];
sx q[3];
rz(-1.7921721) q[3];
sx q[3];
rz(-0.27316545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75605679) q[2];
sx q[2];
rz(-1.7588561) q[2];
sx q[2];
rz(-2.5059911) q[2];
rz(0.33291891) q[3];
sx q[3];
rz(-1.0115441) q[3];
sx q[3];
rz(-2.6788768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(1.8024837) q[0];
sx q[0];
rz(-3.1153296) q[0];
sx q[0];
rz(3.0185757) q[0];
rz(-1.7440375) q[1];
sx q[1];
rz(-1.7536283) q[1];
sx q[1];
rz(-0.77973286) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62209475) q[0];
sx q[0];
rz(-0.72185282) q[0];
sx q[0];
rz(-1.1546385) q[0];
x q[1];
rz(0.17895584) q[2];
sx q[2];
rz(-2.4479986) q[2];
sx q[2];
rz(-0.041698448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18786622) q[1];
sx q[1];
rz(-0.426891) q[1];
sx q[1];
rz(-2.0679451) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.443583) q[3];
sx q[3];
rz(-1.4641342) q[3];
sx q[3];
rz(-0.17388177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66703779) q[2];
sx q[2];
rz(-1.2660618) q[2];
sx q[2];
rz(0.090465821) q[2];
rz(2.7247834) q[3];
sx q[3];
rz(-0.79206812) q[3];
sx q[3];
rz(2.1478103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4823293) q[0];
sx q[0];
rz(-1.3220795) q[0];
sx q[0];
rz(-3.0521159) q[0];
rz(0.80884519) q[1];
sx q[1];
rz(-0.50061148) q[1];
sx q[1];
rz(-1.0577177) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9442975) q[0];
sx q[0];
rz(-0.57387251) q[0];
sx q[0];
rz(2.3969335) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5799149) q[2];
sx q[2];
rz(-2.5895725) q[2];
sx q[2];
rz(0.43285757) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6597418) q[1];
sx q[1];
rz(-2.5999477) q[1];
sx q[1];
rz(-0.84195824) q[1];
x q[2];
rz(-1.5461033) q[3];
sx q[3];
rz(-2.3273483) q[3];
sx q[3];
rz(0.3972185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7384501) q[2];
sx q[2];
rz(-0.53692997) q[2];
sx q[2];
rz(-2.7682847) q[2];
rz(-2.8849854) q[3];
sx q[3];
rz(-1.720287) q[3];
sx q[3];
rz(-3.0671425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8858717) q[0];
sx q[0];
rz(-1.5626361) q[0];
sx q[0];
rz(1.7241021) q[0];
rz(1.7120842) q[1];
sx q[1];
rz(-2.7919339) q[1];
sx q[1];
rz(-2.3333593) q[1];
rz(-3.1222432) q[2];
sx q[2];
rz(-0.45889284) q[2];
sx q[2];
rz(1.5677551) q[2];
rz(-0.96961602) q[3];
sx q[3];
rz(-1.5388699) q[3];
sx q[3];
rz(-0.13997302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
