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
rz(2.7569438) q[0];
sx q[0];
rz(-0.83847133) q[0];
sx q[0];
rz(-0.60929259) q[0];
rz(0.72120178) q[1];
sx q[1];
rz(-0.92858044) q[1];
sx q[1];
rz(-2.0452926) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1288777) q[0];
sx q[0];
rz(-2.4372074) q[0];
sx q[0];
rz(-1.7449858) q[0];
x q[1];
rz(2.6651971) q[2];
sx q[2];
rz(-2.1557249) q[2];
sx q[2];
rz(-2.3830519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3071202) q[1];
sx q[1];
rz(-2.5398919) q[1];
sx q[1];
rz(-0.51096075) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9864086) q[3];
sx q[3];
rz(-1.461457) q[3];
sx q[3];
rz(-0.21462277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.91813749) q[2];
sx q[2];
rz(-1.3104985) q[2];
sx q[2];
rz(1.9523331) q[2];
rz(0.39092815) q[3];
sx q[3];
rz(-1.1632185) q[3];
sx q[3];
rz(-1.1322359) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44553462) q[0];
sx q[0];
rz(-1.3548509) q[0];
sx q[0];
rz(1.6433486) q[0];
rz(-0.6779201) q[1];
sx q[1];
rz(-1.3005715) q[1];
sx q[1];
rz(-0.71151412) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16879119) q[0];
sx q[0];
rz(-1.1074705) q[0];
sx q[0];
rz(0.40131779) q[0];
rz(-pi) q[1];
rz(1.3337801) q[2];
sx q[2];
rz(-2.0227602) q[2];
sx q[2];
rz(2.1353561) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9715648) q[1];
sx q[1];
rz(-1.6508647) q[1];
sx q[1];
rz(0.34645924) q[1];
rz(-pi) q[2];
x q[2];
rz(0.031458843) q[3];
sx q[3];
rz(-0.94740564) q[3];
sx q[3];
rz(-1.7707847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.70637643) q[2];
sx q[2];
rz(-1.3560359) q[2];
sx q[2];
rz(1.0630652) q[2];
rz(1.4937909) q[3];
sx q[3];
rz(-2.6726275) q[3];
sx q[3];
rz(-0.74023214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4110334) q[0];
sx q[0];
rz(-0.39091245) q[0];
sx q[0];
rz(0.60182369) q[0];
rz(2.9959294) q[1];
sx q[1];
rz(-1.0258976) q[1];
sx q[1];
rz(-2.513733) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8191166) q[0];
sx q[0];
rz(-1.5014582) q[0];
sx q[0];
rz(3.022911) q[0];
rz(1.026903) q[2];
sx q[2];
rz(-2.5182407) q[2];
sx q[2];
rz(2.198213) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16937373) q[1];
sx q[1];
rz(-1.7986606) q[1];
sx q[1];
rz(1.9337173) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95730036) q[3];
sx q[3];
rz(-1.9578551) q[3];
sx q[3];
rz(1.8107896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7566028) q[2];
sx q[2];
rz(-2.0774697) q[2];
sx q[2];
rz(0.21558726) q[2];
rz(-2.0032517) q[3];
sx q[3];
rz(-1.3201069) q[3];
sx q[3];
rz(2.5707572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0067921) q[0];
sx q[0];
rz(-1.726806) q[0];
sx q[0];
rz(0.11882812) q[0];
rz(1.6670594) q[1];
sx q[1];
rz(-2.2524565) q[1];
sx q[1];
rz(0.80783358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87619239) q[0];
sx q[0];
rz(-1.0256488) q[0];
sx q[0];
rz(2.4643174) q[0];
x q[1];
rz(1.7052682) q[2];
sx q[2];
rz(-2.8155012) q[2];
sx q[2];
rz(1.6539314) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1203907) q[1];
sx q[1];
rz(-2.5087416) q[1];
sx q[1];
rz(-1.2854043) q[1];
rz(-pi) q[2];
rz(-2.4633154) q[3];
sx q[3];
rz(-1.5587574) q[3];
sx q[3];
rz(1.2348242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5011751) q[2];
sx q[2];
rz(-1.4150323) q[2];
sx q[2];
rz(0.51587063) q[2];
rz(-3.0473895) q[3];
sx q[3];
rz(-0.7754063) q[3];
sx q[3];
rz(1.5883821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3696988) q[0];
sx q[0];
rz(-1.1006681) q[0];
sx q[0];
rz(1.1871673) q[0];
rz(2.5502603) q[1];
sx q[1];
rz(-1.2966803) q[1];
sx q[1];
rz(-2.3147413) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39385228) q[0];
sx q[0];
rz(-1.6918139) q[0];
sx q[0];
rz(2.7583102) q[0];
rz(1.44136) q[2];
sx q[2];
rz(-2.5797306) q[2];
sx q[2];
rz(-0.33298408) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0650658) q[1];
sx q[1];
rz(-0.69628104) q[1];
sx q[1];
rz(-1.7743006) q[1];
x q[2];
rz(2.1416792) q[3];
sx q[3];
rz(-1.1383082) q[3];
sx q[3];
rz(-1.692358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7190711) q[2];
sx q[2];
rz(-2.8870388) q[2];
sx q[2];
rz(-2.1264326) q[2];
rz(-0.90967956) q[3];
sx q[3];
rz(-1.2404975) q[3];
sx q[3];
rz(-0.66143405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5231617) q[0];
sx q[0];
rz(-1.2105415) q[0];
sx q[0];
rz(-0.30995187) q[0];
rz(3.1228206) q[1];
sx q[1];
rz(-1.680178) q[1];
sx q[1];
rz(3.054256) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76742889) q[0];
sx q[0];
rz(-1.5383676) q[0];
sx q[0];
rz(-0.57291605) q[0];
rz(-pi) q[1];
x q[1];
rz(2.491966) q[2];
sx q[2];
rz(-0.98043331) q[2];
sx q[2];
rz(0.43250674) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8467056) q[1];
sx q[1];
rz(-2.069983) q[1];
sx q[1];
rz(0.45392068) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15743001) q[3];
sx q[3];
rz(-2.0937908) q[3];
sx q[3];
rz(-0.98539814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7901018) q[2];
sx q[2];
rz(-0.47372207) q[2];
sx q[2];
rz(-0.52004415) q[2];
rz(-0.13488787) q[3];
sx q[3];
rz(-1.3210195) q[3];
sx q[3];
rz(1.7322056) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0670369) q[0];
sx q[0];
rz(-1.073607) q[0];
sx q[0];
rz(0.38597646) q[0];
rz(-2.0416253) q[1];
sx q[1];
rz(-2.4620582) q[1];
sx q[1];
rz(-2.3540672) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58819492) q[0];
sx q[0];
rz(-2.6752691) q[0];
sx q[0];
rz(-2.2539) q[0];
x q[1];
rz(2.6359796) q[2];
sx q[2];
rz(-2.4596301) q[2];
sx q[2];
rz(-1.4346892) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13788504) q[1];
sx q[1];
rz(-2.2785997) q[1];
sx q[1];
rz(-0.46143175) q[1];
rz(-pi) q[2];
rz(-2.1925364) q[3];
sx q[3];
rz(-0.28655401) q[3];
sx q[3];
rz(1.3499945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0218574) q[2];
sx q[2];
rz(-1.3081552) q[2];
sx q[2];
rz(-1.5137399) q[2];
rz(-1.9205903) q[3];
sx q[3];
rz(-1.9767714) q[3];
sx q[3];
rz(0.47202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30962238) q[0];
sx q[0];
rz(-1.4291052) q[0];
sx q[0];
rz(0.23304644) q[0];
rz(-1.4876935) q[1];
sx q[1];
rz(-2.1079886) q[1];
sx q[1];
rz(2.1344562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45389913) q[0];
sx q[0];
rz(-1.0570044) q[0];
sx q[0];
rz(-1.9625241) q[0];
rz(1.8677668) q[2];
sx q[2];
rz(-2.179702) q[2];
sx q[2];
rz(-1.1394274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1888652) q[1];
sx q[1];
rz(-0.42546526) q[1];
sx q[1];
rz(1.2098321) q[1];
rz(-1.5910025) q[3];
sx q[3];
rz(-1.8148225) q[3];
sx q[3];
rz(-0.013210162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86203376) q[2];
sx q[2];
rz(-1.1125914) q[2];
sx q[2];
rz(2.2588008) q[2];
rz(-0.27859303) q[3];
sx q[3];
rz(-1.168058) q[3];
sx q[3];
rz(-2.6826503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19926628) q[0];
sx q[0];
rz(-2.6658604) q[0];
sx q[0];
rz(1.5698154) q[0];
rz(0.78872952) q[1];
sx q[1];
rz(-0.9659583) q[1];
sx q[1];
rz(-2.4615361) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4231789) q[0];
sx q[0];
rz(-2.3299651) q[0];
sx q[0];
rz(1.4885097) q[0];
rz(-pi) q[1];
rz(-3.0228258) q[2];
sx q[2];
rz(-1.4450018) q[2];
sx q[2];
rz(-2.072352) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9763261) q[1];
sx q[1];
rz(-1.191947) q[1];
sx q[1];
rz(2.5514609) q[1];
rz(-pi) q[2];
rz(-0.84362883) q[3];
sx q[3];
rz(-1.9988201) q[3];
sx q[3];
rz(-2.6035863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1532229) q[2];
sx q[2];
rz(-0.92066568) q[2];
sx q[2];
rz(1.9752768) q[2];
rz(1.9370646) q[3];
sx q[3];
rz(-1.322999) q[3];
sx q[3];
rz(-2.5148463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.523943) q[0];
sx q[0];
rz(-0.88215041) q[0];
sx q[0];
rz(-0.43706617) q[0];
rz(-2.3433459) q[1];
sx q[1];
rz(-2.6415446) q[1];
sx q[1];
rz(0.22013586) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6740538) q[0];
sx q[0];
rz(-2.3706145) q[0];
sx q[0];
rz(-1.5037554) q[0];
rz(-1.885845) q[2];
sx q[2];
rz(-1.7654668) q[2];
sx q[2];
rz(2.6597629) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5497351) q[1];
sx q[1];
rz(-2.3384574) q[1];
sx q[1];
rz(-1.1863043) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6159358) q[3];
sx q[3];
rz(-2.3718826) q[3];
sx q[3];
rz(1.5509645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2685214) q[2];
sx q[2];
rz(-1.2348509) q[2];
sx q[2];
rz(0.29042563) q[2];
rz(2.5051266) q[3];
sx q[3];
rz(-1.7593743) q[3];
sx q[3];
rz(1.9071473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81454043) q[0];
sx q[0];
rz(-2.4335813) q[0];
sx q[0];
rz(-2.0306564) q[0];
rz(-3.0771599) q[1];
sx q[1];
rz(-1.4705407) q[1];
sx q[1];
rz(2.4992117) q[1];
rz(1.5242794) q[2];
sx q[2];
rz(-0.59534335) q[2];
sx q[2];
rz(-2.6186242) q[2];
rz(-2.6024184) q[3];
sx q[3];
rz(-1.3414345) q[3];
sx q[3];
rz(0.10673005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
