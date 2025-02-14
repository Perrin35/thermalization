OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.34831369) q[0];
sx q[0];
rz(8.0170595) q[0];
sx q[0];
rz(9.4693139) q[0];
rz(-1.1107923) q[1];
sx q[1];
rz(5.0587237) q[1];
sx q[1];
rz(8.6375477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6605741) q[0];
sx q[0];
rz(-2.3957167) q[0];
sx q[0];
rz(-0.77156271) q[0];
rz(-pi) q[1];
rz(-0.039723176) q[2];
sx q[2];
rz(-1.2457704) q[2];
sx q[2];
rz(2.9892493) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8742919) q[1];
sx q[1];
rz(-2.1484765) q[1];
sx q[1];
rz(3.0588052) q[1];
rz(-0.23870339) q[3];
sx q[3];
rz(-1.9145962) q[3];
sx q[3];
rz(1.3581004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.74701509) q[2];
sx q[2];
rz(-1.0801103) q[2];
sx q[2];
rz(0.75418312) q[2];
rz(-1.7265823) q[3];
sx q[3];
rz(-0.49608803) q[3];
sx q[3];
rz(0.32895857) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0542145) q[0];
sx q[0];
rz(-2.7890451) q[0];
sx q[0];
rz(1.8074328) q[0];
rz(-1.2913903) q[1];
sx q[1];
rz(-1.038237) q[1];
sx q[1];
rz(2.2659567) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74036769) q[0];
sx q[0];
rz(-0.90310192) q[0];
sx q[0];
rz(0.69913787) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0793117) q[2];
sx q[2];
rz(-1.0750293) q[2];
sx q[2];
rz(-2.371126) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39159472) q[1];
sx q[1];
rz(-2.007496) q[1];
sx q[1];
rz(-0.32086793) q[1];
rz(1.3631212) q[3];
sx q[3];
rz(-1.2856348) q[3];
sx q[3];
rz(2.7792197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.69677532) q[2];
sx q[2];
rz(-2.4397218) q[2];
sx q[2];
rz(2.2689421) q[2];
rz(-0.78898346) q[3];
sx q[3];
rz(-0.9011457) q[3];
sx q[3];
rz(-0.22855973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.2149684) q[0];
sx q[0];
rz(-1.3525532) q[0];
sx q[0];
rz(3.1335926) q[0];
rz(0.81832111) q[1];
sx q[1];
rz(-0.74940959) q[1];
sx q[1];
rz(1.2368894) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0950355) q[0];
sx q[0];
rz(-0.76991428) q[0];
sx q[0];
rz(2.0410935) q[0];
rz(-pi) q[1];
rz(-2.0047999) q[2];
sx q[2];
rz(-1.889359) q[2];
sx q[2];
rz(-0.54036507) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.19265402) q[1];
sx q[1];
rz(-2.1100419) q[1];
sx q[1];
rz(1.9303028) q[1];
x q[2];
rz(1.453492) q[3];
sx q[3];
rz(-1.8059732) q[3];
sx q[3];
rz(-1.4405574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7311953) q[2];
sx q[2];
rz(-0.13430139) q[2];
sx q[2];
rz(-2.2006939) q[2];
rz(-1.1040374) q[3];
sx q[3];
rz(-1.7887812) q[3];
sx q[3];
rz(1.997939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10821548) q[0];
sx q[0];
rz(-2.5984851) q[0];
sx q[0];
rz(0.23750842) q[0];
rz(0.6595276) q[1];
sx q[1];
rz(-0.57412761) q[1];
sx q[1];
rz(-3.1245756) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1579817) q[0];
sx q[0];
rz(-1.6189623) q[0];
sx q[0];
rz(1.7397797) q[0];
rz(-pi) q[1];
rz(-1.7791564) q[2];
sx q[2];
rz(-2.0063498) q[2];
sx q[2];
rz(2.5591133) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.915565) q[1];
sx q[1];
rz(-1.6751958) q[1];
sx q[1];
rz(-0.56609367) q[1];
rz(2.4260826) q[3];
sx q[3];
rz(-0.9431211) q[3];
sx q[3];
rz(-1.1174517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90618769) q[2];
sx q[2];
rz(-1.2013288) q[2];
sx q[2];
rz(-2.2271633) q[2];
rz(-2.7794481) q[3];
sx q[3];
rz(-0.80949628) q[3];
sx q[3];
rz(2.5563498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(1.2528766) q[0];
sx q[0];
rz(-1.7676366) q[0];
sx q[0];
rz(0.1121029) q[0];
rz(2.0762699) q[1];
sx q[1];
rz(-1.0708829) q[1];
sx q[1];
rz(-1.0692474) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4010942) q[0];
sx q[0];
rz(-1.8938657) q[0];
sx q[0];
rz(-1.857398) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5324131) q[2];
sx q[2];
rz(-1.1523917) q[2];
sx q[2];
rz(2.6765228) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4878975) q[1];
sx q[1];
rz(-1.2733165) q[1];
sx q[1];
rz(-1.9031699) q[1];
x q[2];
rz(1.7642683) q[3];
sx q[3];
rz(-2.5335823) q[3];
sx q[3];
rz(1.6093169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7142882) q[2];
sx q[2];
rz(-0.71791831) q[2];
sx q[2];
rz(2.60738) q[2];
rz(0.30321768) q[3];
sx q[3];
rz(-1.5568308) q[3];
sx q[3];
rz(1.5155972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41820207) q[0];
sx q[0];
rz(-1.8573107) q[0];
sx q[0];
rz(0.89163017) q[0];
rz(1.760969) q[1];
sx q[1];
rz(-1.2029519) q[1];
sx q[1];
rz(-0.92323971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1311125) q[0];
sx q[0];
rz(-1.4725871) q[0];
sx q[0];
rz(1.5183856) q[0];
x q[1];
rz(2.7545287) q[2];
sx q[2];
rz(-1.3157857) q[2];
sx q[2];
rz(2.2789291) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6369317) q[1];
sx q[1];
rz(-1.8430222) q[1];
sx q[1];
rz(-2.9039761) q[1];
x q[2];
rz(0.51839101) q[3];
sx q[3];
rz(-0.71037358) q[3];
sx q[3];
rz(0.47391665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0455857) q[2];
sx q[2];
rz(-0.43232375) q[2];
sx q[2];
rz(1.879479) q[2];
rz(-1.9803842) q[3];
sx q[3];
rz(-1.592417) q[3];
sx q[3];
rz(-1.8956634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4534509) q[0];
sx q[0];
rz(-0.83204404) q[0];
sx q[0];
rz(1.6081109) q[0];
rz(-2.7104132) q[1];
sx q[1];
rz(-0.9175514) q[1];
sx q[1];
rz(1.7327488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4154177) q[0];
sx q[0];
rz(-1.2841932) q[0];
sx q[0];
rz(0.64985259) q[0];
x q[1];
rz(1.1860851) q[2];
sx q[2];
rz(-1.4885579) q[2];
sx q[2];
rz(-0.68160328) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3235493) q[1];
sx q[1];
rz(-2.9406639) q[1];
sx q[1];
rz(2.4646321) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.124229) q[3];
sx q[3];
rz(-2.0696739) q[3];
sx q[3];
rz(-1.4829829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2048753) q[2];
sx q[2];
rz(-1.1254346) q[2];
sx q[2];
rz(1.3196866) q[2];
rz(3.1070869) q[3];
sx q[3];
rz(-1.6104108) q[3];
sx q[3];
rz(1.7721133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40783229) q[0];
sx q[0];
rz(-0.5540846) q[0];
sx q[0];
rz(1.2169417) q[0];
rz(0.92562428) q[1];
sx q[1];
rz(-1.7763014) q[1];
sx q[1];
rz(-1.2006753) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72382765) q[0];
sx q[0];
rz(-0.96207959) q[0];
sx q[0];
rz(-0.62792741) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.578892) q[2];
sx q[2];
rz(-2.4614053) q[2];
sx q[2];
rz(-2.0923751) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0266704) q[1];
sx q[1];
rz(-1.5296796) q[1];
sx q[1];
rz(-2.369472) q[1];
rz(-pi) q[2];
rz(-1.1481836) q[3];
sx q[3];
rz(-1.234493) q[3];
sx q[3];
rz(0.5748394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.919148) q[2];
sx q[2];
rz(-1.5347975) q[2];
sx q[2];
rz(-2.0929125) q[2];
rz(-0.18516304) q[3];
sx q[3];
rz(-1.984963) q[3];
sx q[3];
rz(-3.0345501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1222526) q[0];
sx q[0];
rz(-0.66146079) q[0];
sx q[0];
rz(-0.81073236) q[0];
rz(-0.73075378) q[1];
sx q[1];
rz(-2.0850875) q[1];
sx q[1];
rz(-0.96053851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86772462) q[0];
sx q[0];
rz(-1.5698213) q[0];
sx q[0];
rz(-0.0033897059) q[0];
x q[1];
rz(2.2645017) q[2];
sx q[2];
rz(-0.6232647) q[2];
sx q[2];
rz(-2.7359112) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3523622) q[1];
sx q[1];
rz(-1.8327729) q[1];
sx q[1];
rz(2.1195223) q[1];
rz(-1.9979565) q[3];
sx q[3];
rz(-1.8064326) q[3];
sx q[3];
rz(0.7251265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9591799) q[2];
sx q[2];
rz(-1.7824087) q[2];
sx q[2];
rz(-1.5394999) q[2];
rz(2.0942073) q[3];
sx q[3];
rz(-1.7644707) q[3];
sx q[3];
rz(1.876095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3613116) q[0];
sx q[0];
rz(-2.2762716) q[0];
sx q[0];
rz(0.58746946) q[0];
rz(2.7777708) q[1];
sx q[1];
rz(-1.1452585) q[1];
sx q[1];
rz(2.2344373) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.912187) q[0];
sx q[0];
rz(-2.7813781) q[0];
sx q[0];
rz(0.2304669) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67423363) q[2];
sx q[2];
rz(-2.5620954) q[2];
sx q[2];
rz(-0.63636875) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.75395155) q[1];
sx q[1];
rz(-2.6389737) q[1];
sx q[1];
rz(0.44746621) q[1];
rz(0.085368319) q[3];
sx q[3];
rz(-1.2078346) q[3];
sx q[3];
rz(-0.1051108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5017447) q[2];
sx q[2];
rz(-0.088968337) q[2];
sx q[2];
rz(-2.7196344) q[2];
rz(0.0095327775) q[3];
sx q[3];
rz(-1.5671174) q[3];
sx q[3];
rz(-0.52614051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6137153) q[0];
sx q[0];
rz(-1.2953225) q[0];
sx q[0];
rz(-3.0184826) q[0];
rz(-2.4404424) q[1];
sx q[1];
rz(-1.9155365) q[1];
sx q[1];
rz(1.6446) q[1];
rz(-1.9385626) q[2];
sx q[2];
rz(-2.2333973) q[2];
sx q[2];
rz(-1.0495627) q[2];
rz(1.9811859) q[3];
sx q[3];
rz(-1.3153362) q[3];
sx q[3];
rz(1.5462331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
