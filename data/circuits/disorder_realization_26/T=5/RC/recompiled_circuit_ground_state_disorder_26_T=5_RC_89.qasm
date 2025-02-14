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
rz(-0.74690312) q[0];
sx q[0];
rz(-2.3009543) q[0];
rz(0.12159881) q[1];
sx q[1];
rz(-1.2727979) q[1];
sx q[1];
rz(-2.8425541) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6071607) q[0];
sx q[0];
rz(-2.2272155) q[0];
sx q[0];
rz(0.54404152) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7173355) q[2];
sx q[2];
rz(-2.6508109) q[2];
sx q[2];
rz(2.7637568) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92235293) q[1];
sx q[1];
rz(-1.7139707) q[1];
sx q[1];
rz(2.8994865) q[1];
rz(-pi) q[2];
rz(2.998803) q[3];
sx q[3];
rz(-2.0215109) q[3];
sx q[3];
rz(-2.5888909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76089871) q[2];
sx q[2];
rz(-2.0913405) q[2];
sx q[2];
rz(-2.8139581) q[2];
rz(1.3753752) q[3];
sx q[3];
rz(-1.4204493) q[3];
sx q[3];
rz(2.6473141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7647917) q[0];
sx q[0];
rz(-1.6015653) q[0];
sx q[0];
rz(2.2862527) q[0];
rz(3.1335462) q[1];
sx q[1];
rz(-1.2866373) q[1];
sx q[1];
rz(2.6904552) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52569235) q[0];
sx q[0];
rz(-2.2886758) q[0];
sx q[0];
rz(-2.9327716) q[0];
rz(-pi) q[1];
rz(-1.3786267) q[2];
sx q[2];
rz(-0.8079257) q[2];
sx q[2];
rz(0.59468036) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0939744) q[1];
sx q[1];
rz(-2.6674358) q[1];
sx q[1];
rz(0.38800254) q[1];
rz(-pi) q[2];
rz(1.1633918) q[3];
sx q[3];
rz(-2.0137069) q[3];
sx q[3];
rz(-0.16678424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.96192876) q[2];
sx q[2];
rz(-1.1397866) q[2];
sx q[2];
rz(0.12953225) q[2];
rz(-0.1772964) q[3];
sx q[3];
rz(-2.5640021) q[3];
sx q[3];
rz(0.052848335) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.391908) q[0];
sx q[0];
rz(-0.41202298) q[0];
sx q[0];
rz(-2.5823197) q[0];
rz(0.094712146) q[1];
sx q[1];
rz(-1.5385224) q[1];
sx q[1];
rz(0.47725484) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9608979) q[0];
sx q[0];
rz(-1.8261693) q[0];
sx q[0];
rz(-1.9810505) q[0];
x q[1];
rz(-2.7800326) q[2];
sx q[2];
rz(-1.430856) q[2];
sx q[2];
rz(-0.41658336) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2169358) q[1];
sx q[1];
rz(-0.81243304) q[1];
sx q[1];
rz(1.8820018) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52366728) q[3];
sx q[3];
rz(-0.91991495) q[3];
sx q[3];
rz(-1.6802579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39530784) q[2];
sx q[2];
rz(-1.8751273) q[2];
sx q[2];
rz(2.2115808) q[2];
rz(1.2578472) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(-0.045684489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3607218) q[0];
sx q[0];
rz(-1.5683132) q[0];
sx q[0];
rz(2.1029396) q[0];
rz(0.61141283) q[1];
sx q[1];
rz(-2.2544506) q[1];
sx q[1];
rz(0.92322737) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1425689) q[0];
sx q[0];
rz(-1.9143189) q[0];
sx q[0];
rz(-0.96238636) q[0];
rz(-pi) q[1];
rz(-3.1298248) q[2];
sx q[2];
rz(-1.9885049) q[2];
sx q[2];
rz(0.92064806) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7103017) q[1];
sx q[1];
rz(-1.4863621) q[1];
sx q[1];
rz(-0.84655098) q[1];
rz(-pi) q[2];
x q[2];
rz(2.491288) q[3];
sx q[3];
rz(-1.3612483) q[3];
sx q[3];
rz(2.9858231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63127798) q[2];
sx q[2];
rz(-0.46773657) q[2];
sx q[2];
rz(0.54747096) q[2];
rz(-0.064420961) q[3];
sx q[3];
rz(-1.8378601) q[3];
sx q[3];
rz(2.2678383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.6814394) q[0];
rz(2.1414781) q[1];
sx q[1];
rz(-0.8668879) q[1];
sx q[1];
rz(0.11988457) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7728876) q[0];
sx q[0];
rz(-1.3471706) q[0];
sx q[0];
rz(2.2740433) q[0];
rz(-2.5566543) q[2];
sx q[2];
rz(-1.0046994) q[2];
sx q[2];
rz(-1.0370129) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29741187) q[1];
sx q[1];
rz(-1.6496804) q[1];
sx q[1];
rz(0.21551883) q[1];
rz(-pi) q[2];
rz(-0.67976953) q[3];
sx q[3];
rz(-0.7228557) q[3];
sx q[3];
rz(2.2732609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.039006058) q[2];
sx q[2];
rz(-2.234499) q[2];
sx q[2];
rz(-2.8809123) q[2];
rz(-0.87812224) q[3];
sx q[3];
rz(-1.9120646) q[3];
sx q[3];
rz(0.78316435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53073019) q[0];
sx q[0];
rz(-3.0513638) q[0];
sx q[0];
rz(-1.0020142) q[0];
rz(-0.37725457) q[1];
sx q[1];
rz(-0.95163029) q[1];
sx q[1];
rz(-0.9446876) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1169918) q[0];
sx q[0];
rz(-2.5178066) q[0];
sx q[0];
rz(-0.37895112) q[0];
rz(-2.5414921) q[2];
sx q[2];
rz(-2.1515111) q[2];
sx q[2];
rz(2.3490259) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2608883) q[1];
sx q[1];
rz(-0.35683888) q[1];
sx q[1];
rz(-0.44512213) q[1];
rz(-1.3768436) q[3];
sx q[3];
rz(-0.55046457) q[3];
sx q[3];
rz(2.6725519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79823309) q[2];
sx q[2];
rz(-1.6075906) q[2];
sx q[2];
rz(0.043924335) q[2];
rz(-1.6262866) q[3];
sx q[3];
rz(-1.1224727) q[3];
sx q[3];
rz(1.4656434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2783022) q[0];
sx q[0];
rz(-0.55178061) q[0];
sx q[0];
rz(-0.58468753) q[0];
rz(2.0901285) q[1];
sx q[1];
rz(-2.3221071) q[1];
sx q[1];
rz(-2.6712766) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7746587) q[0];
sx q[0];
rz(-2.1485188) q[0];
sx q[0];
rz(1.3960741) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1675179) q[2];
sx q[2];
rz(-1.2263311) q[2];
sx q[2];
rz(-0.044805077) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17543381) q[1];
sx q[1];
rz(-0.81635469) q[1];
sx q[1];
rz(-0.50558009) q[1];
rz(1.5830273) q[3];
sx q[3];
rz(-1.6837956) q[3];
sx q[3];
rz(1.6390273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26479244) q[2];
sx q[2];
rz(-2.6046643) q[2];
sx q[2];
rz(2.8301767) q[2];
rz(-0.16820678) q[3];
sx q[3];
rz(-1.5663389) q[3];
sx q[3];
rz(-0.28347191) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6300221) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(-0.93609634) q[0];
rz(-0.49631897) q[1];
sx q[1];
rz(-0.66718188) q[1];
sx q[1];
rz(-0.47028968) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5514126) q[0];
sx q[0];
rz(-2.0594993) q[0];
sx q[0];
rz(1.9103229) q[0];
rz(-pi) q[1];
rz(3.1390433) q[2];
sx q[2];
rz(-1.7685206) q[2];
sx q[2];
rz(-1.5288803) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3067842) q[1];
sx q[1];
rz(-0.34903279) q[1];
sx q[1];
rz(-1.7743737) q[1];
x q[2];
rz(1.3948729) q[3];
sx q[3];
rz(-1.8430437) q[3];
sx q[3];
rz(-1.3678838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.63142598) q[2];
sx q[2];
rz(-1.3672071) q[2];
sx q[2];
rz(1.6660956) q[2];
rz(-0.79536074) q[3];
sx q[3];
rz(-0.16203351) q[3];
sx q[3];
rz(-1.8719261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0610166) q[0];
sx q[0];
rz(-0.59122714) q[0];
sx q[0];
rz(3.0812145) q[0];
rz(-2.9810442) q[1];
sx q[1];
rz(-1.5269273) q[1];
sx q[1];
rz(0.9789595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2706021) q[0];
sx q[0];
rz(-0.63828642) q[0];
sx q[0];
rz(2.9459475) q[0];
rz(-pi) q[1];
rz(-1.2376182) q[2];
sx q[2];
rz(-0.71014437) q[2];
sx q[2];
rz(1.6802481) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1886602) q[1];
sx q[1];
rz(-0.97129909) q[1];
sx q[1];
rz(-2.1486234) q[1];
rz(-1.5185202) q[3];
sx q[3];
rz(-2.5617122) q[3];
sx q[3];
rz(-1.9573905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0673361) q[2];
sx q[2];
rz(-0.29198519) q[2];
sx q[2];
rz(3.0687029) q[2];
rz(0.59761754) q[3];
sx q[3];
rz(-1.3772734) q[3];
sx q[3];
rz(-0.0028751956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6947967) q[0];
sx q[0];
rz(-2.1294761) q[0];
sx q[0];
rz(2.6126675) q[0];
rz(2.9534598) q[1];
sx q[1];
rz(-0.70536047) q[1];
sx q[1];
rz(2.5573152) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3007418) q[0];
sx q[0];
rz(-0.35984958) q[0];
sx q[0];
rz(-0.76143439) q[0];
rz(-pi) q[1];
rz(2.4227117) q[2];
sx q[2];
rz(-0.55698697) q[2];
sx q[2];
rz(-2.5124036) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0127231) q[1];
sx q[1];
rz(-1.1616316) q[1];
sx q[1];
rz(-1.6175458) q[1];
rz(-pi) q[2];
rz(0.27354555) q[3];
sx q[3];
rz(-0.97018948) q[3];
sx q[3];
rz(2.724444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.75954413) q[2];
sx q[2];
rz(-2.1376762) q[2];
sx q[2];
rz(-0.071852597) q[2];
rz(-1.0233277) q[3];
sx q[3];
rz(-1.5724678) q[3];
sx q[3];
rz(2.6527827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95579424) q[0];
sx q[0];
rz(-2.4837942) q[0];
sx q[0];
rz(2.876045) q[0];
rz(-2.4664948) q[1];
sx q[1];
rz(-1.594512) q[1];
sx q[1];
rz(-0.095269861) q[1];
rz(0.51495348) q[2];
sx q[2];
rz(-1.5025768) q[2];
sx q[2];
rz(1.2160355) q[2];
rz(-0.76473372) q[3];
sx q[3];
rz(-0.22976362) q[3];
sx q[3];
rz(0.094658628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
