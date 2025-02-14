OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.793279) q[0];
sx q[0];
rz(-1.7338742) q[0];
sx q[0];
rz(-3.0970567) q[0];
rz(-1.1107923) q[1];
sx q[1];
rz(5.0587237) q[1];
sx q[1];
rz(8.6375477) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6985748) q[0];
sx q[0];
rz(-2.0788143) q[0];
sx q[0];
rz(-0.99850151) q[0];
x q[1];
rz(1.8960612) q[2];
sx q[2];
rz(-1.533154) q[2];
sx q[2];
rz(1.4311439) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2582142) q[1];
sx q[1];
rz(-1.5014663) q[1];
sx q[1];
rz(2.1500471) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23870339) q[3];
sx q[3];
rz(-1.2269964) q[3];
sx q[3];
rz(-1.7834922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3945776) q[2];
sx q[2];
rz(-1.0801103) q[2];
sx q[2];
rz(-0.75418312) q[2];
rz(1.4150103) q[3];
sx q[3];
rz(-2.6455046) q[3];
sx q[3];
rz(2.8126341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0542145) q[0];
sx q[0];
rz(-2.7890451) q[0];
sx q[0];
rz(-1.8074328) q[0];
rz(1.8502024) q[1];
sx q[1];
rz(-1.038237) q[1];
sx q[1];
rz(-0.87563595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74036769) q[0];
sx q[0];
rz(-0.90310192) q[0];
sx q[0];
rz(-0.69913787) q[0];
x q[1];
rz(3.0793117) q[2];
sx q[2];
rz(-2.0665633) q[2];
sx q[2];
rz(2.371126) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7499979) q[1];
sx q[1];
rz(-1.1340967) q[1];
sx q[1];
rz(-0.32086793) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61296566) q[3];
sx q[3];
rz(-0.35109441) q[3];
sx q[3];
rz(-2.8610817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.69677532) q[2];
sx q[2];
rz(-2.4397218) q[2];
sx q[2];
rz(-0.87265054) q[2];
rz(-0.78898346) q[3];
sx q[3];
rz(-2.240447) q[3];
sx q[3];
rz(0.22855973) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2149684) q[0];
sx q[0];
rz(-1.7890395) q[0];
sx q[0];
rz(0.0080000814) q[0];
rz(2.3232715) q[1];
sx q[1];
rz(-2.3921831) q[1];
sx q[1];
rz(-1.9047033) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2674026) q[0];
sx q[0];
rz(-1.2498901) q[0];
sx q[0];
rz(0.85808922) q[0];
x q[1];
rz(-2.7929467) q[2];
sx q[2];
rz(-1.1599891) q[2];
sx q[2];
rz(2.2553159) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3171009) q[1];
sx q[1];
rz(-2.503509) q[1];
sx q[1];
rz(2.610105) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6872356) q[3];
sx q[3];
rz(-0.26232346) q[3];
sx q[3];
rz(1.9087877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7311953) q[2];
sx q[2];
rz(-3.0072913) q[2];
sx q[2];
rz(-0.94089874) q[2];
rz(-2.0375552) q[3];
sx q[3];
rz(-1.7887812) q[3];
sx q[3];
rz(-1.997939) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0333772) q[0];
sx q[0];
rz(-0.5431076) q[0];
sx q[0];
rz(-0.23750842) q[0];
rz(-2.4820651) q[1];
sx q[1];
rz(-0.57412761) q[1];
sx q[1];
rz(0.017017078) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7369923) q[0];
sx q[0];
rz(-1.4020108) q[0];
sx q[0];
rz(-0.048860839) q[0];
x q[1];
rz(2.6976162) q[2];
sx q[2];
rz(-1.7594565) q[2];
sx q[2];
rz(-0.89934811) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.81777311) q[1];
sx q[1];
rz(-2.5669879) q[1];
sx q[1];
rz(2.9486604) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83570273) q[3];
sx q[3];
rz(-0.91360211) q[3];
sx q[3];
rz(3.0007255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.235405) q[2];
sx q[2];
rz(-1.9402639) q[2];
sx q[2];
rz(-0.91442937) q[2];
rz(0.36214456) q[3];
sx q[3];
rz(-2.3320964) q[3];
sx q[3];
rz(-2.5563498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-1.2528766) q[0];
sx q[0];
rz(-1.7676366) q[0];
sx q[0];
rz(3.0294898) q[0];
rz(-1.0653227) q[1];
sx q[1];
rz(-1.0708829) q[1];
sx q[1];
rz(-1.0692474) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0651848) q[0];
sx q[0];
rz(-1.2994081) q[0];
sx q[0];
rz(0.33581622) q[0];
rz(-1.094355) q[2];
sx q[2];
rz(-2.0531056) q[2];
sx q[2];
rz(-2.2708238) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4878975) q[1];
sx q[1];
rz(-1.2733165) q[1];
sx q[1];
rz(-1.2384227) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97158708) q[3];
sx q[3];
rz(-1.4607444) q[3];
sx q[3];
rz(-0.19796619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4273044) q[2];
sx q[2];
rz(-0.71791831) q[2];
sx q[2];
rz(-2.60738) q[2];
rz(2.838375) q[3];
sx q[3];
rz(-1.5847619) q[3];
sx q[3];
rz(1.5155972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41820207) q[0];
sx q[0];
rz(-1.8573107) q[0];
sx q[0];
rz(2.2499625) q[0];
rz(-1.760969) q[1];
sx q[1];
rz(-1.2029519) q[1];
sx q[1];
rz(-2.2183529) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5864201) q[0];
sx q[0];
rz(-1.6229543) q[0];
sx q[0];
rz(3.0432493) q[0];
x q[1];
rz(1.8452066) q[2];
sx q[2];
rz(-1.9447127) q[2];
sx q[2];
rz(-0.81059882) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1404851) q[1];
sx q[1];
rz(-1.342089) q[1];
sx q[1];
rz(1.8504935) q[1];
rz(-2.6232016) q[3];
sx q[3];
rz(-2.4312191) q[3];
sx q[3];
rz(2.667676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.096006958) q[2];
sx q[2];
rz(-2.7092689) q[2];
sx q[2];
rz(-1.2621137) q[2];
rz(-1.1612085) q[3];
sx q[3];
rz(-1.592417) q[3];
sx q[3];
rz(1.8956634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4534509) q[0];
sx q[0];
rz(-0.83204404) q[0];
sx q[0];
rz(-1.6081109) q[0];
rz(-2.7104132) q[1];
sx q[1];
rz(-0.9175514) q[1];
sx q[1];
rz(-1.4088438) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20074318) q[0];
sx q[0];
rz(-2.4398167) q[0];
sx q[0];
rz(0.45325847) q[0];
rz(0.088690745) q[2];
sx q[2];
rz(-1.1874532) q[2];
sx q[2];
rz(-2.2191522) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.085693749) q[1];
sx q[1];
rz(-1.4454464) q[1];
sx q[1];
rz(-0.15744822) q[1];
rz(-1.3471782) q[3];
sx q[3];
rz(-0.51284344) q[3];
sx q[3];
rz(1.7382857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93671736) q[2];
sx q[2];
rz(-2.016158) q[2];
sx q[2];
rz(1.8219061) q[2];
rz(-3.1070869) q[3];
sx q[3];
rz(-1.6104108) q[3];
sx q[3];
rz(-1.7721133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7337604) q[0];
sx q[0];
rz(-0.5540846) q[0];
sx q[0];
rz(-1.2169417) q[0];
rz(-0.92562428) q[1];
sx q[1];
rz(-1.3652912) q[1];
sx q[1];
rz(1.9409174) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.417765) q[0];
sx q[0];
rz(-2.1795131) q[0];
sx q[0];
rz(2.5136652) q[0];
rz(-pi) q[1];
rz(-1.1633792) q[2];
sx q[2];
rz(-1.0098741) q[2];
sx q[2];
rz(1.7307868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7279049) q[1];
sx q[1];
rz(-2.3686045) q[1];
sx q[1];
rz(-3.0826921) q[1];
rz(0.36603277) q[3];
sx q[3];
rz(-1.1732374) q[3];
sx q[3];
rz(-1.9983069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.22244464) q[2];
sx q[2];
rz(-1.6067952) q[2];
sx q[2];
rz(-2.0929125) q[2];
rz(0.18516304) q[3];
sx q[3];
rz(-1.984963) q[3];
sx q[3];
rz(-0.10704253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1222526) q[0];
sx q[0];
rz(-2.4801319) q[0];
sx q[0];
rz(-2.3308603) q[0];
rz(0.73075378) q[1];
sx q[1];
rz(-2.0850875) q[1];
sx q[1];
rz(0.96053851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.273868) q[0];
sx q[0];
rz(-1.5698213) q[0];
sx q[0];
rz(-0.0033897059) q[0];
rz(-pi) q[1];
rz(-0.4308295) q[2];
sx q[2];
rz(-2.0362051) q[2];
sx q[2];
rz(-2.7499466) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.78923046) q[1];
sx q[1];
rz(-1.3088198) q[1];
sx q[1];
rz(2.1195223) q[1];
rz(-pi) q[2];
rz(2.8836684) q[3];
sx q[3];
rz(-1.9854133) q[3];
sx q[3];
rz(0.95154155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9591799) q[2];
sx q[2];
rz(-1.3591839) q[2];
sx q[2];
rz(-1.6020927) q[2];
rz(-2.0942073) q[3];
sx q[3];
rz(-1.3771219) q[3];
sx q[3];
rz(1.876095) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3613116) q[0];
sx q[0];
rz(-2.2762716) q[0];
sx q[0];
rz(-2.5541232) q[0];
rz(-0.36382183) q[1];
sx q[1];
rz(-1.1452585) q[1];
sx q[1];
rz(2.2344373) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5575378) q[0];
sx q[0];
rz(-1.6514002) q[0];
sx q[0];
rz(0.35146468) q[0];
rz(0.67423363) q[2];
sx q[2];
rz(-0.57949726) q[2];
sx q[2];
rz(-0.63636875) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.75395155) q[1];
sx q[1];
rz(-0.50261897) q[1];
sx q[1];
rz(0.44746621) q[1];
rz(-1.3499522) q[3];
sx q[3];
rz(-0.37243249) q[3];
sx q[3];
rz(0.34162921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63984799) q[2];
sx q[2];
rz(-3.0526243) q[2];
sx q[2];
rz(0.42195827) q[2];
rz(-0.0095327775) q[3];
sx q[3];
rz(-1.5671174) q[3];
sx q[3];
rz(-2.6154521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.20303) q[2];
sx q[2];
rz(-0.90819539) q[2];
sx q[2];
rz(2.09203) q[2];
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
