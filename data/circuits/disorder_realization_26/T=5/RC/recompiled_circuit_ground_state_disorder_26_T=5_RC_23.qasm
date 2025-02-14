OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.77312624) q[0];
sx q[0];
rz(2.3946895) q[0];
sx q[0];
rz(11.725732) q[0];
rz(3.2631915) q[1];
sx q[1];
rz(-1.8687948) q[1];
sx q[1];
rz(12.267332) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5344319) q[0];
sx q[0];
rz(-0.91437712) q[0];
sx q[0];
rz(-0.54404152) q[0];
rz(1.4242572) q[2];
sx q[2];
rz(-2.6508109) q[2];
sx q[2];
rz(-2.7637568) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12451367) q[1];
sx q[1];
rz(-2.8610367) q[1];
sx q[1];
rz(-2.6002167) q[1];
rz(-pi) q[2];
rz(1.1160581) q[3];
sx q[3];
rz(-1.4423496) q[3];
sx q[3];
rz(-1.0806395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76089871) q[2];
sx q[2];
rz(-2.0913405) q[2];
sx q[2];
rz(2.8139581) q[2];
rz(1.3753752) q[3];
sx q[3];
rz(-1.7211434) q[3];
sx q[3];
rz(0.49427858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7647917) q[0];
sx q[0];
rz(-1.5400274) q[0];
sx q[0];
rz(-2.2862527) q[0];
rz(3.1335462) q[1];
sx q[1];
rz(-1.2866373) q[1];
sx q[1];
rz(-0.45113742) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83734918) q[0];
sx q[0];
rz(-2.3991802) q[0];
sx q[0];
rz(-1.3377331) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9443953) q[2];
sx q[2];
rz(-2.3596546) q[2];
sx q[2];
rz(0.3202084) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0939744) q[1];
sx q[1];
rz(-2.6674358) q[1];
sx q[1];
rz(-0.38800254) q[1];
rz(-pi) q[2];
rz(0.47685949) q[3];
sx q[3];
rz(-1.9369159) q[3];
sx q[3];
rz(1.2211293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.96192876) q[2];
sx q[2];
rz(-2.0018061) q[2];
sx q[2];
rz(3.0120604) q[2];
rz(2.9642963) q[3];
sx q[3];
rz(-0.57759053) q[3];
sx q[3];
rz(3.0887443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7496846) q[0];
sx q[0];
rz(-2.7295697) q[0];
sx q[0];
rz(-2.5823197) q[0];
rz(-0.094712146) q[1];
sx q[1];
rz(-1.5385224) q[1];
sx q[1];
rz(-0.47725484) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6420631) q[0];
sx q[0];
rz(-1.9669878) q[0];
sx q[0];
rz(-0.27735932) q[0];
x q[1];
rz(-1.7202708) q[2];
sx q[2];
rz(-1.9286619) q[2];
sx q[2];
rz(2.0400816) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2169358) q[1];
sx q[1];
rz(-0.81243304) q[1];
sx q[1];
rz(1.8820018) q[1];
rz(-pi) q[2];
rz(2.151793) q[3];
sx q[3];
rz(-0.81077164) q[3];
sx q[3];
rz(-0.9188942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7462848) q[2];
sx q[2];
rz(-1.2664653) q[2];
sx q[2];
rz(0.9300119) q[2];
rz(-1.8837455) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(3.0959082) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3607218) q[0];
sx q[0];
rz(-1.5683132) q[0];
sx q[0];
rz(-2.1029396) q[0];
rz(-0.61141283) q[1];
sx q[1];
rz(-2.2544506) q[1];
sx q[1];
rz(2.2183653) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1194342) q[0];
sx q[0];
rz(-0.6878449) q[0];
sx q[0];
rz(1.0115959) q[0];
x q[1];
rz(1.153062) q[2];
sx q[2];
rz(-1.5600403) q[2];
sx q[2];
rz(-2.4962184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.23455305) q[1];
sx q[1];
rz(-0.72826339) q[1];
sx q[1];
rz(1.4437463) q[1];
rz(-pi) q[2];
rz(0.65030469) q[3];
sx q[3];
rz(-1.7803444) q[3];
sx q[3];
rz(2.9858231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5103147) q[2];
sx q[2];
rz(-2.6738561) q[2];
sx q[2];
rz(-0.54747096) q[2];
rz(-0.064420961) q[3];
sx q[3];
rz(-1.3037325) q[3];
sx q[3];
rz(0.87375435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15544686) q[0];
sx q[0];
rz(-2.2387945) q[0];
sx q[0];
rz(1.4601532) q[0];
rz(2.1414781) q[1];
sx q[1];
rz(-2.2747048) q[1];
sx q[1];
rz(-0.11988457) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016235624) q[0];
sx q[0];
rz(-0.88847697) q[0];
sx q[0];
rz(-2.8518139) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2861459) q[2];
sx q[2];
rz(-2.3513633) q[2];
sx q[2];
rz(1.2145192) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8441808) q[1];
sx q[1];
rz(-1.4919123) q[1];
sx q[1];
rz(-0.21551883) q[1];
x q[2];
rz(0.67976953) q[3];
sx q[3];
rz(-0.7228557) q[3];
sx q[3];
rz(-2.2732609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1025866) q[2];
sx q[2];
rz(-0.90709364) q[2];
sx q[2];
rz(2.8809123) q[2];
rz(0.87812224) q[3];
sx q[3];
rz(-1.9120646) q[3];
sx q[3];
rz(2.3584283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53073019) q[0];
sx q[0];
rz(-0.090228883) q[0];
sx q[0];
rz(1.0020142) q[0];
rz(-0.37725457) q[1];
sx q[1];
rz(-2.1899624) q[1];
sx q[1];
rz(0.9446876) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1169918) q[0];
sx q[0];
rz(-2.5178066) q[0];
sx q[0];
rz(2.7626415) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89904191) q[2];
sx q[2];
rz(-2.0625249) q[2];
sx q[2];
rz(1.1373625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3516253) q[1];
sx q[1];
rz(-1.891544) q[1];
sx q[1];
rz(-1.7299537) q[1];
rz(-1.0287337) q[3];
sx q[3];
rz(-1.4698062) q[3];
sx q[3];
rz(0.93588692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3433596) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(3.0976683) q[2];
rz(-1.6262866) q[3];
sx q[3];
rz(-1.1224727) q[3];
sx q[3];
rz(-1.6759492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8632904) q[0];
sx q[0];
rz(-0.55178061) q[0];
sx q[0];
rz(-0.58468753) q[0];
rz(-2.0901285) q[1];
sx q[1];
rz(-2.3221071) q[1];
sx q[1];
rz(-0.47031602) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6795652) q[0];
sx q[0];
rz(-0.60067486) q[0];
sx q[0];
rz(-2.8810049) q[0];
x q[1];
rz(2.3114071) q[2];
sx q[2];
rz(-0.52415327) q[2];
sx q[2];
rz(-0.94607991) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2863949) q[1];
sx q[1];
rz(-0.87955399) q[1];
sx q[1];
rz(-2.046584) q[1];
rz(0.11300762) q[3];
sx q[3];
rz(-1.5829493) q[3];
sx q[3];
rz(3.0719824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8768002) q[2];
sx q[2];
rz(-2.6046643) q[2];
sx q[2];
rz(0.31141591) q[2];
rz(0.16820678) q[3];
sx q[3];
rz(-1.5663389) q[3];
sx q[3];
rz(0.28347191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6300221) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(-2.2054963) q[0];
rz(-0.49631897) q[1];
sx q[1];
rz(-2.4744108) q[1];
sx q[1];
rz(0.47028968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5514126) q[0];
sx q[0];
rz(-2.0594993) q[0];
sx q[0];
rz(-1.2312698) q[0];
x q[1];
rz(-1.5835205) q[2];
sx q[2];
rz(-0.1977405) q[2];
sx q[2];
rz(1.5997353) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6185604) q[1];
sx q[1];
rz(-1.9123239) q[1];
sx q[1];
rz(3.0681472) q[1];
rz(-2.8652898) q[3];
sx q[3];
rz(-1.4014162) q[3];
sx q[3];
rz(2.9864431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63142598) q[2];
sx q[2];
rz(-1.3672071) q[2];
sx q[2];
rz(-1.4754971) q[2];
rz(-2.3462319) q[3];
sx q[3];
rz(-0.16203351) q[3];
sx q[3];
rz(1.8719261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0610166) q[0];
sx q[0];
rz(-0.59122714) q[0];
sx q[0];
rz(-3.0812145) q[0];
rz(-0.16054842) q[1];
sx q[1];
rz(-1.6146654) q[1];
sx q[1];
rz(-2.1626332) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5125324) q[0];
sx q[0];
rz(-2.1950025) q[0];
sx q[0];
rz(1.7140304) q[0];
rz(1.9039745) q[2];
sx q[2];
rz(-0.71014437) q[2];
sx q[2];
rz(1.6802481) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97040788) q[1];
sx q[1];
rz(-2.0386341) q[1];
sx q[1];
rz(-0.68433185) q[1];
rz(-pi) q[2];
rz(0.99154226) q[3];
sx q[3];
rz(-1.5421621) q[3];
sx q[3];
rz(2.7987409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0673361) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(-3.0687029) q[2];
rz(0.59761754) q[3];
sx q[3];
rz(-1.7643192) q[3];
sx q[3];
rz(0.0028751956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6947967) q[0];
sx q[0];
rz(-2.1294761) q[0];
sx q[0];
rz(-0.52892518) q[0];
rz(-2.9534598) q[1];
sx q[1];
rz(-2.4362322) q[1];
sx q[1];
rz(-0.58427748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0952632) q[0];
sx q[0];
rz(-1.3130616) q[0];
sx q[0];
rz(-1.3168174) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1816013) q[2];
sx q[2];
rz(-1.1616594) q[2];
sx q[2];
rz(-1.429806) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0127231) q[1];
sx q[1];
rz(-1.979961) q[1];
sx q[1];
rz(1.6175458) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9464302) q[3];
sx q[3];
rz(-2.4886819) q[3];
sx q[3];
rz(-0.87797166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3820485) q[2];
sx q[2];
rz(-1.0039165) q[2];
sx q[2];
rz(-3.0697401) q[2];
rz(-1.0233277) q[3];
sx q[3];
rz(-1.5691248) q[3];
sx q[3];
rz(0.48880997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95579424) q[0];
sx q[0];
rz(-0.6577984) q[0];
sx q[0];
rz(-0.26554769) q[0];
rz(2.4664948) q[1];
sx q[1];
rz(-1.5470807) q[1];
sx q[1];
rz(3.0463228) q[1];
rz(-3.0037389) q[2];
sx q[2];
rz(-0.51904917) q[2];
sx q[2];
rz(2.9067155) q[2];
rz(-2.3768589) q[3];
sx q[3];
rz(-2.911829) q[3];
sx q[3];
rz(-3.046934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
