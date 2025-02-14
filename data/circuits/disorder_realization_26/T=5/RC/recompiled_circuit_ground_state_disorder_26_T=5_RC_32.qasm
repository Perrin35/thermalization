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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75343412) q[0];
sx q[0];
rz(-2.3155766) q[0];
sx q[0];
rz(0.97919925) q[0];
rz(3.063721) q[2];
sx q[2];
rz(-1.0857333) q[2];
sx q[2];
rz(0.21202206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12451367) q[1];
sx q[1];
rz(-2.8610367) q[1];
sx q[1];
rz(0.54137595) q[1];
rz(-pi) q[2];
rz(-0.14278966) q[3];
sx q[3];
rz(-2.0215109) q[3];
sx q[3];
rz(0.55270178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76089871) q[2];
sx q[2];
rz(-1.0502522) q[2];
sx q[2];
rz(0.32763457) q[2];
rz(-1.7662175) q[3];
sx q[3];
rz(-1.4204493) q[3];
sx q[3];
rz(-0.49427858) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37680092) q[0];
sx q[0];
rz(-1.5400274) q[0];
sx q[0];
rz(0.85533992) q[0];
rz(-3.1335462) q[1];
sx q[1];
rz(-1.8549553) q[1];
sx q[1];
rz(-0.45113742) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2349885) q[0];
sx q[0];
rz(-1.7275817) q[0];
sx q[0];
rz(-0.84201987) q[0];
x q[1];
rz(1.3786267) q[2];
sx q[2];
rz(-2.333667) q[2];
sx q[2];
rz(-2.5469123) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6169713) q[1];
sx q[1];
rz(-1.1344304) q[1];
sx q[1];
rz(1.379016) q[1];
rz(-pi) q[2];
rz(1.1633918) q[3];
sx q[3];
rz(-1.1278858) q[3];
sx q[3];
rz(-2.9748084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.96192876) q[2];
sx q[2];
rz(-2.0018061) q[2];
sx q[2];
rz(0.12953225) q[2];
rz(-2.9642963) q[3];
sx q[3];
rz(-0.57759053) q[3];
sx q[3];
rz(0.052848335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.7496846) q[0];
sx q[0];
rz(-2.7295697) q[0];
sx q[0];
rz(0.55927292) q[0];
rz(3.0468805) q[1];
sx q[1];
rz(-1.6030703) q[1];
sx q[1];
rz(-2.6643378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1806948) q[0];
sx q[0];
rz(-1.3154234) q[0];
sx q[0];
rz(-1.1605422) q[0];
x q[1];
rz(1.7202708) q[2];
sx q[2];
rz(-1.2129307) q[2];
sx q[2];
rz(-1.1015111) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86384799) q[1];
sx q[1];
rz(-1.3466292) q[1];
sx q[1];
rz(0.78296354) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84945143) q[3];
sx q[3];
rz(-1.1616544) q[3];
sx q[3];
rz(0.22709286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39530784) q[2];
sx q[2];
rz(-1.2664653) q[2];
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
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.88714209) q[1];
sx q[1];
rz(0.92322737) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99902376) q[0];
sx q[0];
rz(-1.2272738) q[0];
sx q[0];
rz(2.1792063) q[0];
rz(1.153062) q[2];
sx q[2];
rz(-1.5815524) q[2];
sx q[2];
rz(-0.64537424) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9070396) q[1];
sx q[1];
rz(-2.4133293) q[1];
sx q[1];
rz(-1.6978463) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33780725) q[3];
sx q[3];
rz(-0.67854133) q[3];
sx q[3];
rz(-1.9935009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63127798) q[2];
sx q[2];
rz(-0.46773657) q[2];
sx q[2];
rz(0.54747096) q[2];
rz(3.0771717) q[3];
sx q[3];
rz(-1.8378601) q[3];
sx q[3];
rz(2.2678383) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15544686) q[0];
sx q[0];
rz(-0.90279818) q[0];
sx q[0];
rz(1.4601532) q[0];
rz(-2.1414781) q[1];
sx q[1];
rz(-2.2747048) q[1];
sx q[1];
rz(0.11988457) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45792199) q[0];
sx q[0];
rz(-2.409465) q[0];
sx q[0];
rz(1.2326272) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91953711) q[2];
sx q[2];
rz(-2.0555758) q[2];
sx q[2];
rz(-2.9491021) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.29741187) q[1];
sx q[1];
rz(-1.4919123) q[1];
sx q[1];
rz(-2.9260738) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0644887) q[3];
sx q[3];
rz(-2.1112006) q[3];
sx q[3];
rz(-3.0960954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.039006058) q[2];
sx q[2];
rz(-0.90709364) q[2];
sx q[2];
rz(-0.26068035) q[2];
rz(0.87812224) q[3];
sx q[3];
rz(-1.2295281) q[3];
sx q[3];
rz(0.78316435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53073019) q[0];
sx q[0];
rz(-3.0513638) q[0];
sx q[0];
rz(-2.1395785) q[0];
rz(0.37725457) q[1];
sx q[1];
rz(-2.1899624) q[1];
sx q[1];
rz(2.196905) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2827891) q[0];
sx q[0];
rz(-1.3529881) q[0];
sx q[0];
rz(2.5522638) q[0];
rz(2.5414921) q[2];
sx q[2];
rz(-0.99008152) q[2];
sx q[2];
rz(-0.7925668) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4113271) q[1];
sx q[1];
rz(-1.7217727) q[1];
sx q[1];
rz(-2.81702) q[1];
rz(-pi) q[2];
rz(-0.11774534) q[3];
sx q[3];
rz(-1.0317993) q[3];
sx q[3];
rz(-0.6955516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79823309) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(-0.043924335) q[2];
rz(1.515306) q[3];
sx q[3];
rz(-1.1224727) q[3];
sx q[3];
rz(1.4656434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8632904) q[0];
sx q[0];
rz(-2.589812) q[0];
sx q[0];
rz(-0.58468753) q[0];
rz(-2.0901285) q[1];
sx q[1];
rz(-0.81948558) q[1];
sx q[1];
rz(-2.6712766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4620275) q[0];
sx q[0];
rz(-2.5409178) q[0];
sx q[0];
rz(-2.8810049) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1675179) q[2];
sx q[2];
rz(-1.9152616) q[2];
sx q[2];
rz(-0.044805077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1086169) q[1];
sx q[1];
rz(-1.2101296) q[1];
sx q[1];
rz(2.3920139) q[1];
rz(-pi) q[2];
rz(-1.5585654) q[3];
sx q[3];
rz(-1.6837956) q[3];
sx q[3];
rz(1.6390273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6300221) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(0.93609634) q[0];
rz(2.6452737) q[1];
sx q[1];
rz(-0.66718188) q[1];
sx q[1];
rz(-0.47028968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1449438) q[0];
sx q[0];
rz(-1.2723288) q[0];
sx q[0];
rz(2.6281283) q[0];
rz(-1.5580721) q[2];
sx q[2];
rz(-0.1977405) q[2];
sx q[2];
rz(1.5418574) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.072402231) q[1];
sx q[1];
rz(-1.5015998) q[1];
sx q[1];
rz(1.9131768) q[1];
rz(-pi) q[2];
rz(1.7467198) q[3];
sx q[3];
rz(-1.8430437) q[3];
sx q[3];
rz(1.3678838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5101667) q[2];
sx q[2];
rz(-1.3672071) q[2];
sx q[2];
rz(1.4754971) q[2];
rz(-2.3462319) q[3];
sx q[3];
rz(-2.9795591) q[3];
sx q[3];
rz(-1.8719261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.0805761) q[0];
sx q[0];
rz(-0.59122714) q[0];
sx q[0];
rz(-0.06037816) q[0];
rz(2.9810442) q[1];
sx q[1];
rz(-1.6146654) q[1];
sx q[1];
rz(0.9789595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5125324) q[0];
sx q[0];
rz(-0.94659014) q[0];
sx q[0];
rz(-1.4275622) q[0];
rz(-pi) q[1];
rz(2.8674815) q[2];
sx q[2];
rz(-0.90702552) q[2];
sx q[2];
rz(1.2520777) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0461523) q[1];
sx q[1];
rz(-2.3344724) q[1];
sx q[1];
rz(2.4673106) q[1];
rz(-pi) q[2];
rz(3.1073808) q[3];
sx q[3];
rz(-0.99181038) q[3];
sx q[3];
rz(1.8949231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0673361) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(3.0687029) q[2];
rz(-0.59761754) q[3];
sx q[3];
rz(-1.3772734) q[3];
sx q[3];
rz(0.0028751956) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.446796) q[0];
sx q[0];
rz(-1.0121166) q[0];
sx q[0];
rz(-0.52892518) q[0];
rz(0.18813285) q[1];
sx q[1];
rz(-2.4362322) q[1];
sx q[1];
rz(-0.58427748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3007418) q[0];
sx q[0];
rz(-0.35984958) q[0];
sx q[0];
rz(2.3801583) q[0];
rz(-2.7033349) q[2];
sx q[2];
rz(-1.9263679) q[2];
sx q[2];
rz(-0.3027161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0127231) q[1];
sx q[1];
rz(-1.979961) q[1];
sx q[1];
rz(1.6175458) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95240611) q[3];
sx q[3];
rz(-1.795553) q[3];
sx q[3];
rz(0.99638961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3820485) q[2];
sx q[2];
rz(-2.1376762) q[2];
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
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1857984) q[0];
sx q[0];
rz(-2.4837942) q[0];
sx q[0];
rz(2.876045) q[0];
rz(-2.4664948) q[1];
sx q[1];
rz(-1.594512) q[1];
sx q[1];
rz(-0.095269861) q[1];
rz(3.0037389) q[2];
sx q[2];
rz(-2.6225435) q[2];
sx q[2];
rz(-0.23487716) q[2];
rz(-0.16719462) q[3];
sx q[3];
rz(-1.4124558) q[3];
sx q[3];
rz(2.4168933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
