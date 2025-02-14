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
rz(0.84063831) q[0];
rz(0.12159881) q[1];
sx q[1];
rz(-1.2727979) q[1];
sx q[1];
rz(0.29903856) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39002573) q[0];
sx q[0];
rz(-1.9932858) q[0];
sx q[0];
rz(0.83777352) q[0];
x q[1];
rz(-2.0571124) q[2];
sx q[2];
rz(-1.5019226) q[2];
sx q[2];
rz(-1.8191847) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92235293) q[1];
sx q[1];
rz(-1.427622) q[1];
sx q[1];
rz(0.24210614) q[1];
rz(-1.1160581) q[3];
sx q[3];
rz(-1.699243) q[3];
sx q[3];
rz(2.0609531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.76089871) q[2];
sx q[2];
rz(-1.0502522) q[2];
sx q[2];
rz(-0.32763457) q[2];
rz(1.3753752) q[3];
sx q[3];
rz(-1.7211434) q[3];
sx q[3];
rz(0.49427858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37680092) q[0];
sx q[0];
rz(-1.5400274) q[0];
sx q[0];
rz(-0.85533992) q[0];
rz(0.0080464706) q[1];
sx q[1];
rz(-1.2866373) q[1];
sx q[1];
rz(-2.6904552) q[1];
x q[2];
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
rz(2.9443953) q[2];
sx q[2];
rz(-0.78193808) q[2];
sx q[2];
rz(0.3202084) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0135368) q[1];
sx q[1];
rz(-1.3971796) q[1];
sx q[1];
rz(-0.44349576) q[1];
rz(-pi) q[2];
rz(-1.9782009) q[3];
sx q[3];
rz(-1.1278858) q[3];
sx q[3];
rz(-2.9748084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1796639) q[2];
sx q[2];
rz(-1.1397866) q[2];
sx q[2];
rz(-0.12953225) q[2];
rz(0.1772964) q[3];
sx q[3];
rz(-2.5640021) q[3];
sx q[3];
rz(3.0887443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7496846) q[0];
sx q[0];
rz(-0.41202298) q[0];
sx q[0];
rz(-2.5823197) q[0];
rz(3.0468805) q[1];
sx q[1];
rz(-1.6030703) q[1];
sx q[1];
rz(-2.6643378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6420631) q[0];
sx q[0];
rz(-1.9669878) q[0];
sx q[0];
rz(0.27735932) q[0];
rz(1.4213218) q[2];
sx q[2];
rz(-1.2129307) q[2];
sx q[2];
rz(-2.0400816) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2169358) q[1];
sx q[1];
rz(-2.3291596) q[1];
sx q[1];
rz(-1.2595909) q[1];
rz(-0.84945143) q[3];
sx q[3];
rz(-1.9799383) q[3];
sx q[3];
rz(-0.22709286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39530784) q[2];
sx q[2];
rz(-1.2664653) q[2];
sx q[2];
rz(0.9300119) q[2];
rz(-1.2578472) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(0.045684489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.78087085) q[0];
sx q[0];
rz(-1.5683132) q[0];
sx q[0];
rz(2.1029396) q[0];
rz(0.61141283) q[1];
sx q[1];
rz(-2.2544506) q[1];
sx q[1];
rz(0.92322737) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1425689) q[0];
sx q[0];
rz(-1.9143189) q[0];
sx q[0];
rz(2.1792063) q[0];
rz(3.1298248) q[2];
sx q[2];
rz(-1.9885049) q[2];
sx q[2];
rz(2.2209446) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.065040914) q[1];
sx q[1];
rz(-0.84969798) q[1];
sx q[1];
rz(-3.0290719) q[1];
x q[2];
rz(-2.8037854) q[3];
sx q[3];
rz(-2.4630513) q[3];
sx q[3];
rz(-1.9935009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5103147) q[2];
sx q[2];
rz(-0.46773657) q[2];
sx q[2];
rz(-0.54747096) q[2];
rz(3.0771717) q[3];
sx q[3];
rz(-1.3037325) q[3];
sx q[3];
rz(-2.2678383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9861458) q[0];
sx q[0];
rz(-2.2387945) q[0];
sx q[0];
rz(-1.6814394) q[0];
rz(-1.0001146) q[1];
sx q[1];
rz(-2.2747048) q[1];
sx q[1];
rz(3.0217081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6836707) q[0];
sx q[0];
rz(-2.409465) q[0];
sx q[0];
rz(-1.2326272) q[0];
rz(-0.91953711) q[2];
sx q[2];
rz(-1.0860168) q[2];
sx q[2];
rz(-0.19249053) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5226871) q[1];
sx q[1];
rz(-2.9123016) q[1];
sx q[1];
rz(-2.7875336) q[1];
rz(-pi) q[2];
rz(-2.4618231) q[3];
sx q[3];
rz(-0.7228557) q[3];
sx q[3];
rz(-2.2732609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1025866) q[2];
sx q[2];
rz(-2.234499) q[2];
sx q[2];
rz(0.26068035) q[2];
rz(0.87812224) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53073019) q[0];
sx q[0];
rz(-3.0513638) q[0];
sx q[0];
rz(-2.1395785) q[0];
rz(2.7643381) q[1];
sx q[1];
rz(-2.1899624) q[1];
sx q[1];
rz(-2.196905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7101038) q[0];
sx q[0];
rz(-0.99718192) q[0];
sx q[0];
rz(-1.8309893) q[0];
rz(2.281419) q[2];
sx q[2];
rz(-0.80931907) q[2];
sx q[2];
rz(-3.0391673) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2608883) q[1];
sx q[1];
rz(-0.35683888) q[1];
sx q[1];
rz(0.44512213) q[1];
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
rz(0.79823309) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(3.0976683) q[2];
rz(1.515306) q[3];
sx q[3];
rz(-1.1224727) q[3];
sx q[3];
rz(1.4656434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.8632904) q[0];
sx q[0];
rz(-0.55178061) q[0];
sx q[0];
rz(0.58468753) q[0];
rz(-2.0901285) q[1];
sx q[1];
rz(-0.81948558) q[1];
sx q[1];
rz(-2.6712766) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4620275) q[0];
sx q[0];
rz(-0.60067486) q[0];
sx q[0];
rz(0.26058773) q[0];
rz(2.7696848) q[2];
sx q[2];
rz(-1.9491409) q[2];
sx q[2];
rz(1.7586979) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1086169) q[1];
sx q[1];
rz(-1.931463) q[1];
sx q[1];
rz(0.74957871) q[1];
rz(-1.5585654) q[3];
sx q[3];
rz(-1.6837956) q[3];
sx q[3];
rz(-1.5025653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26479244) q[2];
sx q[2];
rz(-2.6046643) q[2];
sx q[2];
rz(-2.8301767) q[2];
rz(-0.16820678) q[3];
sx q[3];
rz(-1.5663389) q[3];
sx q[3];
rz(-0.28347191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-0.51157057) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(0.93609634) q[0];
rz(-0.49631897) q[1];
sx q[1];
rz(-2.4744108) q[1];
sx q[1];
rz(0.47028968) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054809991) q[0];
sx q[0];
rz(-0.58712372) q[0];
sx q[0];
rz(-2.5820288) q[0];
rz(1.7685212) q[2];
sx q[2];
rz(-1.573296) q[2];
sx q[2];
rz(-3.1001774) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5230323) q[1];
sx q[1];
rz(-1.2292687) q[1];
sx q[1];
rz(-0.073445436) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3948729) q[3];
sx q[3];
rz(-1.8430437) q[3];
sx q[3];
rz(-1.3678838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0805761) q[0];
sx q[0];
rz(-0.59122714) q[0];
sx q[0];
rz(0.06037816) q[0];
rz(-2.9810442) q[1];
sx q[1];
rz(-1.6146654) q[1];
sx q[1];
rz(2.1626332) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2706021) q[0];
sx q[0];
rz(-0.63828642) q[0];
sx q[0];
rz(-0.19564512) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9039745) q[2];
sx q[2];
rz(-0.71014437) q[2];
sx q[2];
rz(-1.4613446) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.95293249) q[1];
sx q[1];
rz(-0.97129909) q[1];
sx q[1];
rz(-2.1486234) q[1];
rz(-pi) q[2];
rz(1.6230725) q[3];
sx q[3];
rz(-2.5617122) q[3];
sx q[3];
rz(1.1842022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0742566) q[2];
sx q[2];
rz(-0.29198519) q[2];
sx q[2];
rz(-0.072889797) q[2];
rz(-2.5439751) q[3];
sx q[3];
rz(-1.7643192) q[3];
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
x q[3];
rz(pi/2) q[3];
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
rz(2.6947967) q[0];
sx q[0];
rz(-1.0121166) q[0];
sx q[0];
rz(-0.52892518) q[0];
rz(-2.9534598) q[1];
sx q[1];
rz(-2.4362322) q[1];
sx q[1];
rz(2.5573152) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0952632) q[0];
sx q[0];
rz(-1.8285311) q[0];
sx q[0];
rz(1.3168174) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7033349) q[2];
sx q[2];
rz(-1.9263679) q[2];
sx q[2];
rz(-0.3027161) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6021298) q[1];
sx q[1];
rz(-1.5279084) q[1];
sx q[1];
rz(2.7320288) q[1];
x q[2];
rz(0.27354555) q[3];
sx q[3];
rz(-2.1714032) q[3];
sx q[3];
rz(0.41714868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.75954413) q[2];
sx q[2];
rz(-1.0039165) q[2];
sx q[2];
rz(0.071852597) q[2];
rz(-2.1182649) q[3];
sx q[3];
rz(-1.5724678) q[3];
sx q[3];
rz(-2.6527827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1857984) q[0];
sx q[0];
rz(-2.4837942) q[0];
sx q[0];
rz(2.876045) q[0];
rz(0.67509782) q[1];
sx q[1];
rz(-1.594512) q[1];
sx q[1];
rz(-0.095269861) q[1];
rz(2.6266392) q[2];
sx q[2];
rz(-1.6390159) q[2];
sx q[2];
rz(-1.9255571) q[2];
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
