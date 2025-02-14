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
rz(3.2631915) q[1];
sx q[1];
rz(-1.8687948) q[1];
sx q[1];
rz(12.267332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6071607) q[0];
sx q[0];
rz(-0.91437712) q[0];
sx q[0];
rz(-2.5975511) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0844803) q[2];
sx q[2];
rz(-1.5019226) q[2];
sx q[2];
rz(-1.322408) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.017079) q[1];
sx q[1];
rz(-0.28055596) q[1];
sx q[1];
rz(-2.6002167) q[1];
rz(1.2848008) q[3];
sx q[3];
rz(-0.47130775) q[3];
sx q[3];
rz(-2.907674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.76089871) q[2];
sx q[2];
rz(-2.0913405) q[2];
sx q[2];
rz(-0.32763457) q[2];
rz(1.3753752) q[3];
sx q[3];
rz(-1.7211434) q[3];
sx q[3];
rz(-2.6473141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7647917) q[0];
sx q[0];
rz(-1.5400274) q[0];
sx q[0];
rz(-2.2862527) q[0];
rz(0.0080464706) q[1];
sx q[1];
rz(-1.2866373) q[1];
sx q[1];
rz(0.45113742) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3042435) q[0];
sx q[0];
rz(-2.3991802) q[0];
sx q[0];
rz(1.3377331) q[0];
rz(0.19719736) q[2];
sx q[2];
rz(-0.78193808) q[2];
sx q[2];
rz(-0.3202084) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0135368) q[1];
sx q[1];
rz(-1.3971796) q[1];
sx q[1];
rz(-0.44349576) q[1];
x q[2];
rz(-2.4456817) q[3];
sx q[3];
rz(-2.5491121) q[3];
sx q[3];
rz(-0.95564465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96192876) q[2];
sx q[2];
rz(-1.1397866) q[2];
sx q[2];
rz(-3.0120604) q[2];
rz(2.9642963) q[3];
sx q[3];
rz(-2.5640021) q[3];
sx q[3];
rz(-3.0887443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7496846) q[0];
sx q[0];
rz(-2.7295697) q[0];
sx q[0];
rz(2.5823197) q[0];
rz(-0.094712146) q[1];
sx q[1];
rz(-1.5385224) q[1];
sx q[1];
rz(2.6643378) q[1];
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
rz(-pi) q[1];
x q[1];
rz(1.4213218) q[2];
sx q[2];
rz(-1.2129307) q[2];
sx q[2];
rz(1.1015111) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86384799) q[1];
sx q[1];
rz(-1.7949634) q[1];
sx q[1];
rz(2.3586291) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84945143) q[3];
sx q[3];
rz(-1.9799383) q[3];
sx q[3];
rz(2.9144998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7462848) q[2];
sx q[2];
rz(-1.8751273) q[2];
sx q[2];
rz(0.9300119) q[2];
rz(1.8837455) q[3];
sx q[3];
rz(-1.427429) q[3];
sx q[3];
rz(-0.045684489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3607218) q[0];
sx q[0];
rz(-1.5683132) q[0];
sx q[0];
rz(2.1029396) q[0];
rz(2.5301798) q[1];
sx q[1];
rz(-0.88714209) q[1];
sx q[1];
rz(0.92322737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8002565) q[0];
sx q[0];
rz(-1.0025327) q[0];
sx q[0];
rz(-2.7305014) q[0];
rz(-1.153062) q[2];
sx q[2];
rz(-1.5600403) q[2];
sx q[2];
rz(-0.64537424) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4312909) q[1];
sx q[1];
rz(-1.6552306) q[1];
sx q[1];
rz(-2.2950417) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33780725) q[3];
sx q[3];
rz(-2.4630513) q[3];
sx q[3];
rz(1.9935009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5103147) q[2];
sx q[2];
rz(-0.46773657) q[2];
sx q[2];
rz(-0.54747096) q[2];
rz(0.064420961) q[3];
sx q[3];
rz(-1.8378601) q[3];
sx q[3];
rz(0.87375435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
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
rz(0.11988457) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7728876) q[0];
sx q[0];
rz(-1.794422) q[0];
sx q[0];
rz(0.86754936) q[0];
rz(-pi) q[1];
rz(0.85544678) q[2];
sx q[2];
rz(-2.3513633) q[2];
sx q[2];
rz(-1.2145192) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6189055) q[1];
sx q[1];
rz(-0.22929103) q[1];
sx q[1];
rz(0.35405901) q[1];
x q[2];
rz(-2.0771039) q[3];
sx q[3];
rz(-1.0303921) q[3];
sx q[3];
rz(-3.0960954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.039006058) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53073019) q[0];
sx q[0];
rz(-3.0513638) q[0];
sx q[0];
rz(2.1395785) q[0];
rz(0.37725457) q[1];
sx q[1];
rz(-2.1899624) q[1];
sx q[1];
rz(2.196905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7101038) q[0];
sx q[0];
rz(-2.1444107) q[0];
sx q[0];
rz(-1.3106034) q[0];
x q[1];
rz(2.2425507) q[2];
sx q[2];
rz(-1.0790677) q[2];
sx q[2];
rz(1.1373625) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78996736) q[1];
sx q[1];
rz(-1.891544) q[1];
sx q[1];
rz(-1.411639) q[1];
x q[2];
rz(-1.764749) q[3];
sx q[3];
rz(-0.55046457) q[3];
sx q[3];
rz(-2.6725519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79823309) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(3.0976683) q[2];
rz(-1.515306) q[3];
sx q[3];
rz(-1.1224727) q[3];
sx q[3];
rz(1.6759492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8632904) q[0];
sx q[0];
rz(-0.55178061) q[0];
sx q[0];
rz(2.5569051) q[0];
rz(-2.0901285) q[1];
sx q[1];
rz(-2.3221071) q[1];
sx q[1];
rz(2.6712766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0338351) q[0];
sx q[0];
rz(-1.4246539) q[0];
sx q[0];
rz(0.58476292) q[0];
rz(-pi) q[1];
rz(-2.7696848) q[2];
sx q[2];
rz(-1.1924517) q[2];
sx q[2];
rz(-1.3828948) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.85519771) q[1];
sx q[1];
rz(-0.87955399) q[1];
sx q[1];
rz(1.0950086) q[1];
x q[2];
rz(1.5585654) q[3];
sx q[3];
rz(-1.6837956) q[3];
sx q[3];
rz(1.5025653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26479244) q[2];
sx q[2];
rz(-0.53692836) q[2];
sx q[2];
rz(-0.31141591) q[2];
rz(0.16820678) q[3];
sx q[3];
rz(-1.5752537) q[3];
sx q[3];
rz(2.8581207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6300221) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(2.2054963) q[0];
rz(2.6452737) q[1];
sx q[1];
rz(-0.66718188) q[1];
sx q[1];
rz(-0.47028968) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9966489) q[0];
sx q[0];
rz(-1.2723288) q[0];
sx q[0];
rz(-2.6281283) q[0];
rz(-pi) q[1];
rz(1.3730714) q[2];
sx q[2];
rz(-1.5682967) q[2];
sx q[2];
rz(-3.1001774) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3067842) q[1];
sx q[1];
rz(-2.7925599) q[1];
sx q[1];
rz(1.7743737) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3948729) q[3];
sx q[3];
rz(-1.298549) q[3];
sx q[3];
rz(1.7737089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5101667) q[2];
sx q[2];
rz(-1.3672071) q[2];
sx q[2];
rz(1.6660956) q[2];
rz(2.3462319) q[3];
sx q[3];
rz(-0.16203351) q[3];
sx q[3];
rz(1.2696666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0610166) q[0];
sx q[0];
rz(-2.5503655) q[0];
sx q[0];
rz(3.0812145) q[0];
rz(0.16054842) q[1];
sx q[1];
rz(-1.6146654) q[1];
sx q[1];
rz(2.1626332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6290603) q[0];
sx q[0];
rz(-0.94659014) q[0];
sx q[0];
rz(-1.7140304) q[0];
rz(-0.8884807) q[2];
sx q[2];
rz(-1.7856626) q[2];
sx q[2];
rz(-2.9944098) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97040788) q[1];
sx q[1];
rz(-1.1029585) q[1];
sx q[1];
rz(-2.4572608) q[1];
rz(-pi) q[2];
rz(-1.5185202) q[3];
sx q[3];
rz(-2.5617122) q[3];
sx q[3];
rz(1.1842022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0742566) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(-0.072889797) q[2];
rz(0.59761754) q[3];
sx q[3];
rz(-1.3772734) q[3];
sx q[3];
rz(3.1387175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.446796) q[0];
sx q[0];
rz(-2.1294761) q[0];
sx q[0];
rz(0.52892518) q[0];
rz(-0.18813285) q[1];
sx q[1];
rz(-2.4362322) q[1];
sx q[1];
rz(0.58427748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54160252) q[0];
sx q[0];
rz(-1.8162105) q[0];
sx q[0];
rz(0.26588579) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9599914) q[2];
sx q[2];
rz(-1.9799332) q[2];
sx q[2];
rz(1.7117866) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1288695) q[1];
sx q[1];
rz(-1.979961) q[1];
sx q[1];
rz(-1.6175458) q[1];
x q[2];
rz(1.1951625) q[3];
sx q[3];
rz(-2.4886819) q[3];
sx q[3];
rz(0.87797166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75954413) q[2];
sx q[2];
rz(-2.1376762) q[2];
sx q[2];
rz(0.071852597) q[2];
rz(-2.1182649) q[3];
sx q[3];
rz(-1.5691248) q[3];
sx q[3];
rz(-0.48880997) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95579424) q[0];
sx q[0];
rz(-2.4837942) q[0];
sx q[0];
rz(2.876045) q[0];
rz(2.4664948) q[1];
sx q[1];
rz(-1.5470807) q[1];
sx q[1];
rz(3.0463228) q[1];
rz(-0.13785378) q[2];
sx q[2];
rz(-2.6225435) q[2];
sx q[2];
rz(-0.23487716) q[2];
rz(2.3768589) q[3];
sx q[3];
rz(-0.22976362) q[3];
sx q[3];
rz(0.094658628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
