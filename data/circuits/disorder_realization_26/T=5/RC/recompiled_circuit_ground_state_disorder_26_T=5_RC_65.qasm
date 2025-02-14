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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75343412) q[0];
sx q[0];
rz(-2.3155766) q[0];
sx q[0];
rz(-2.1623934) q[0];
x q[1];
rz(0.077871696) q[2];
sx q[2];
rz(-2.0558594) q[2];
sx q[2];
rz(-2.9295706) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68366508) q[1];
sx q[1];
rz(-1.8103765) q[1];
sx q[1];
rz(-1.423382) q[1];
x q[2];
rz(-1.2848008) q[3];
sx q[3];
rz(-0.47130775) q[3];
sx q[3];
rz(-0.23391868) q[3];
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
rz(-0.32763457) q[2];
rz(-1.3753752) q[3];
sx q[3];
rz(-1.4204493) q[3];
sx q[3];
rz(-2.6473141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.37680092) q[0];
sx q[0];
rz(-1.6015653) q[0];
sx q[0];
rz(0.85533992) q[0];
rz(3.1335462) q[1];
sx q[1];
rz(-1.8549553) q[1];
sx q[1];
rz(0.45113742) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52569235) q[0];
sx q[0];
rz(-2.2886758) q[0];
sx q[0];
rz(0.20882102) q[0];
rz(-pi) q[1];
rz(2.3694384) q[2];
sx q[2];
rz(-1.4322965) q[2];
sx q[2];
rz(-1.1097627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0939744) q[1];
sx q[1];
rz(-2.6674358) q[1];
sx q[1];
rz(-0.38800254) q[1];
x q[2];
rz(-2.4456817) q[3];
sx q[3];
rz(-0.59248052) q[3];
sx q[3];
rz(0.95564465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1796639) q[2];
sx q[2];
rz(-2.0018061) q[2];
sx q[2];
rz(-0.12953225) q[2];
rz(0.1772964) q[3];
sx q[3];
rz(-2.5640021) q[3];
sx q[3];
rz(-0.052848335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7496846) q[0];
sx q[0];
rz(-0.41202298) q[0];
sx q[0];
rz(-2.5823197) q[0];
rz(0.094712146) q[1];
sx q[1];
rz(-1.5385224) q[1];
sx q[1];
rz(0.47725484) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4995296) q[0];
sx q[0];
rz(-1.1746049) q[0];
sx q[0];
rz(2.8642333) q[0];
rz(-pi) q[1];
rz(2.7800326) q[2];
sx q[2];
rz(-1.7107367) q[2];
sx q[2];
rz(-0.41658336) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.4871769) q[1];
sx q[1];
rz(-0.80779874) q[1];
sx q[1];
rz(0.31262763) q[1];
rz(-pi) q[2];
rz(-0.52366728) q[3];
sx q[3];
rz(-0.91991495) q[3];
sx q[3];
rz(1.4613348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39530784) q[2];
sx q[2];
rz(-1.2664653) q[2];
sx q[2];
rz(2.2115808) q[2];
rz(-1.8837455) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(3.0959082) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3607218) q[0];
sx q[0];
rz(-1.5732795) q[0];
sx q[0];
rz(2.1029396) q[0];
rz(0.61141283) q[1];
sx q[1];
rz(-0.88714209) q[1];
sx q[1];
rz(2.2183653) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1425689) q[0];
sx q[0];
rz(-1.9143189) q[0];
sx q[0];
rz(0.96238636) q[0];
rz(-pi) q[1];
rz(1.5442886) q[2];
sx q[2];
rz(-0.4178646) q[2];
sx q[2];
rz(-0.94965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.065040914) q[1];
sx q[1];
rz(-0.84969798) q[1];
sx q[1];
rz(-0.11252071) q[1];
rz(-0.33780725) q[3];
sx q[3];
rz(-2.4630513) q[3];
sx q[3];
rz(-1.1480918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5103147) q[2];
sx q[2];
rz(-0.46773657) q[2];
sx q[2];
rz(-0.54747096) q[2];
rz(-3.0771717) q[3];
sx q[3];
rz(-1.8378601) q[3];
sx q[3];
rz(-2.2678383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-0.8668879) q[1];
sx q[1];
rz(-3.0217081) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016235624) q[0];
sx q[0];
rz(-2.2531157) q[0];
sx q[0];
rz(0.28977878) q[0];
rz(-2.5566543) q[2];
sx q[2];
rz(-2.1368933) q[2];
sx q[2];
rz(1.0370129) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8854578) q[1];
sx q[1];
rz(-1.3559582) q[1];
sx q[1];
rz(1.4900521) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67976953) q[3];
sx q[3];
rz(-2.418737) q[3];
sx q[3];
rz(0.86833176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1025866) q[2];
sx q[2];
rz(-2.234499) q[2];
sx q[2];
rz(-0.26068035) q[2];
rz(-2.2634704) q[3];
sx q[3];
rz(-1.2295281) q[3];
sx q[3];
rz(0.78316435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6108625) q[0];
sx q[0];
rz(-0.090228883) q[0];
sx q[0];
rz(1.0020142) q[0];
rz(2.7643381) q[1];
sx q[1];
rz(-0.95163029) q[1];
sx q[1];
rz(2.196905) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8588036) q[0];
sx q[0];
rz(-1.7886046) q[0];
sx q[0];
rz(0.58932886) q[0];
rz(-pi) q[1];
rz(2.2425507) q[2];
sx q[2];
rz(-2.0625249) q[2];
sx q[2];
rz(-1.1373625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73026555) q[1];
sx q[1];
rz(-1.7217727) q[1];
sx q[1];
rz(-0.32457268) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0287337) q[3];
sx q[3];
rz(-1.6717864) q[3];
sx q[3];
rz(0.93588692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.79823309) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(3.0976683) q[2];
rz(-1.6262866) q[3];
sx q[3];
rz(-1.1224727) q[3];
sx q[3];
rz(1.4656434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8632904) q[0];
sx q[0];
rz(-2.589812) q[0];
sx q[0];
rz(2.5569051) q[0];
rz(-2.0901285) q[1];
sx q[1];
rz(-2.3221071) q[1];
sx q[1];
rz(-0.47031602) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6795652) q[0];
sx q[0];
rz(-0.60067486) q[0];
sx q[0];
rz(-2.8810049) q[0];
x q[1];
rz(-1.1675179) q[2];
sx q[2];
rz(-1.2263311) q[2];
sx q[2];
rz(3.0967876) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0329758) q[1];
sx q[1];
rz(-1.931463) q[1];
sx q[1];
rz(2.3920139) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11300762) q[3];
sx q[3];
rz(-1.5829493) q[3];
sx q[3];
rz(-0.069610217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.26479244) q[2];
sx q[2];
rz(-0.53692836) q[2];
sx q[2];
rz(0.31141591) q[2];
rz(-2.9733859) q[3];
sx q[3];
rz(-1.5752537) q[3];
sx q[3];
rz(-0.28347191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51157057) q[0];
sx q[0];
rz(-0.011307414) q[0];
sx q[0];
rz(-0.93609634) q[0];
rz(-0.49631897) q[1];
sx q[1];
rz(-0.66718188) q[1];
sx q[1];
rz(-0.47028968) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9966489) q[0];
sx q[0];
rz(-1.2723288) q[0];
sx q[0];
rz(0.51346438) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3730714) q[2];
sx q[2];
rz(-1.573296) q[2];
sx q[2];
rz(-3.1001774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6185604) q[1];
sx q[1];
rz(-1.2292687) q[1];
sx q[1];
rz(3.0681472) q[1];
rz(-2.5816282) q[3];
sx q[3];
rz(-0.32295152) q[3];
sx q[3];
rz(-1.1895869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5101667) q[2];
sx q[2];
rz(-1.3672071) q[2];
sx q[2];
rz(1.6660956) q[2];
rz(-2.3462319) q[3];
sx q[3];
rz(-2.9795591) q[3];
sx q[3];
rz(-1.8719261) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0610166) q[0];
sx q[0];
rz(-0.59122714) q[0];
sx q[0];
rz(3.0812145) q[0];
rz(-0.16054842) q[1];
sx q[1];
rz(-1.5269273) q[1];
sx q[1];
rz(2.1626332) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5125324) q[0];
sx q[0];
rz(-0.94659014) q[0];
sx q[0];
rz(1.4275622) q[0];
x q[1];
rz(0.8884807) q[2];
sx q[2];
rz(-1.7856626) q[2];
sx q[2];
rz(-0.14718283) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95293249) q[1];
sx q[1];
rz(-0.97129909) q[1];
sx q[1];
rz(0.99296928) q[1];
rz(-0.99154226) q[3];
sx q[3];
rz(-1.5994306) q[3];
sx q[3];
rz(-0.34285173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0742566) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(-3.0687029) q[2];
rz(-2.5439751) q[3];
sx q[3];
rz(-1.7643192) q[3];
sx q[3];
rz(-3.1387175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.446796) q[0];
sx q[0];
rz(-2.1294761) q[0];
sx q[0];
rz(-2.6126675) q[0];
rz(0.18813285) q[1];
sx q[1];
rz(-0.70536047) q[1];
sx q[1];
rz(0.58427748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8408509) q[0];
sx q[0];
rz(-0.35984958) q[0];
sx q[0];
rz(2.3801583) q[0];
rz(-1.9599914) q[2];
sx q[2];
rz(-1.9799332) q[2];
sx q[2];
rz(1.429806) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6021298) q[1];
sx q[1];
rz(-1.6136843) q[1];
sx q[1];
rz(-2.7320288) q[1];
x q[2];
rz(-1.9464302) q[3];
sx q[3];
rz(-2.4886819) q[3];
sx q[3];
rz(0.87797166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.75954413) q[2];
sx q[2];
rz(-2.1376762) q[2];
sx q[2];
rz(-0.071852597) q[2];
rz(1.0233277) q[3];
sx q[3];
rz(-1.5691248) q[3];
sx q[3];
rz(-0.48880997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.95579424) q[0];
sx q[0];
rz(-0.6577984) q[0];
sx q[0];
rz(-0.26554769) q[0];
rz(0.67509782) q[1];
sx q[1];
rz(-1.594512) q[1];
sx q[1];
rz(-0.095269861) q[1];
rz(-2.6266392) q[2];
sx q[2];
rz(-1.5025768) q[2];
sx q[2];
rz(1.2160355) q[2];
rz(-1.731338) q[3];
sx q[3];
rz(-1.7358801) q[3];
sx q[3];
rz(0.87270234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
