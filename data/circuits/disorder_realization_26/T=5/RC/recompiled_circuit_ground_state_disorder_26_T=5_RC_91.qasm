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
rz(1.5344319) q[0];
sx q[0];
rz(-0.91437712) q[0];
sx q[0];
rz(-2.5975511) q[0];
rz(-pi) q[1];
rz(3.063721) q[2];
sx q[2];
rz(-1.0857333) q[2];
sx q[2];
rz(0.21202206) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2192397) q[1];
sx q[1];
rz(-1.7139707) q[1];
sx q[1];
rz(-2.8994865) q[1];
x q[2];
rz(1.1160581) q[3];
sx q[3];
rz(-1.4423496) q[3];
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
rz(-1.4204493) q[3];
sx q[3];
rz(-0.49427858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7647917) q[0];
sx q[0];
rz(-1.5400274) q[0];
sx q[0];
rz(-0.85533992) q[0];
rz(0.0080464706) q[1];
sx q[1];
rz(-1.8549553) q[1];
sx q[1];
rz(-0.45113742) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2349885) q[0];
sx q[0];
rz(-1.7275817) q[0];
sx q[0];
rz(-2.2995728) q[0];
rz(-2.3694384) q[2];
sx q[2];
rz(-1.4322965) q[2];
sx q[2];
rz(-2.0318299) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0135368) q[1];
sx q[1];
rz(-1.744413) q[1];
sx q[1];
rz(2.6980969) q[1];
rz(-pi) q[2];
rz(-1.9782009) q[3];
sx q[3];
rz(-2.0137069) q[3];
sx q[3];
rz(2.9748084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.96192876) q[2];
sx q[2];
rz(-1.1397866) q[2];
sx q[2];
rz(3.0120604) q[2];
rz(2.9642963) q[3];
sx q[3];
rz(-0.57759053) q[3];
sx q[3];
rz(-0.052848335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.6030703) q[1];
sx q[1];
rz(-0.47725484) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6420631) q[0];
sx q[0];
rz(-1.1746049) q[0];
sx q[0];
rz(0.27735932) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7202708) q[2];
sx q[2];
rz(-1.2129307) q[2];
sx q[2];
rz(-1.1015111) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.86384799) q[1];
sx q[1];
rz(-1.3466292) q[1];
sx q[1];
rz(2.3586291) q[1];
rz(-pi) q[2];
rz(-2.2921412) q[3];
sx q[3];
rz(-1.9799383) q[3];
sx q[3];
rz(0.22709286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7462848) q[2];
sx q[2];
rz(-1.8751273) q[2];
sx q[2];
rz(-0.9300119) q[2];
rz(-1.8837455) q[3];
sx q[3];
rz(-1.427429) q[3];
sx q[3];
rz(-3.0959082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3607218) q[0];
sx q[0];
rz(-1.5732795) q[0];
sx q[0];
rz(1.038653) q[0];
rz(0.61141283) q[1];
sx q[1];
rz(-2.2544506) q[1];
sx q[1];
rz(-2.2183653) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99902376) q[0];
sx q[0];
rz(-1.2272738) q[0];
sx q[0];
rz(2.1792063) q[0];
rz(-pi) q[1];
rz(3.1298248) q[2];
sx q[2];
rz(-1.9885049) q[2];
sx q[2];
rz(2.2209446) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.23455305) q[1];
sx q[1];
rz(-0.72826339) q[1];
sx q[1];
rz(1.4437463) q[1];
x q[2];
rz(1.3096894) q[3];
sx q[3];
rz(-2.2045643) q[3];
sx q[3];
rz(1.5719617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5103147) q[2];
sx q[2];
rz(-2.6738561) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.15544686) q[0];
sx q[0];
rz(-0.90279818) q[0];
sx q[0];
rz(-1.6814394) q[0];
rz(2.1414781) q[1];
sx q[1];
rz(-0.8668879) q[1];
sx q[1];
rz(-3.0217081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.125357) q[0];
sx q[0];
rz(-0.88847697) q[0];
sx q[0];
rz(2.8518139) q[0];
x q[1];
rz(-0.91953711) q[2];
sx q[2];
rz(-1.0860168) q[2];
sx q[2];
rz(2.9491021) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5226871) q[1];
sx q[1];
rz(-0.22929103) q[1];
sx q[1];
rz(-0.35405901) q[1];
rz(2.0771039) q[3];
sx q[3];
rz(-1.0303921) q[3];
sx q[3];
rz(-0.04549724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.039006058) q[2];
sx q[2];
rz(-0.90709364) q[2];
sx q[2];
rz(0.26068035) q[2];
rz(-0.87812224) q[3];
sx q[3];
rz(-1.9120646) q[3];
sx q[3];
rz(-2.3584283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6108625) q[0];
sx q[0];
rz(-3.0513638) q[0];
sx q[0];
rz(1.0020142) q[0];
rz(0.37725457) q[1];
sx q[1];
rz(-2.1899624) q[1];
sx q[1];
rz(-0.9446876) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024600895) q[0];
sx q[0];
rz(-0.62378609) q[0];
sx q[0];
rz(2.7626415) q[0];
x q[1];
rz(-2.2425507) q[2];
sx q[2];
rz(-2.0625249) q[2];
sx q[2];
rz(-2.0042302) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3516253) q[1];
sx q[1];
rz(-1.891544) q[1];
sx q[1];
rz(-1.411639) q[1];
rz(-pi) q[2];
rz(2.112859) q[3];
sx q[3];
rz(-1.6717864) q[3];
sx q[3];
rz(-0.93588692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3433596) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(-3.0976683) q[2];
rz(-1.515306) q[3];
sx q[3];
rz(-2.01912) q[3];
sx q[3];
rz(-1.6759492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(0.47031602) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7746587) q[0];
sx q[0];
rz(-0.99307382) q[0];
sx q[0];
rz(1.7455186) q[0];
rz(-pi) q[1];
rz(0.83018556) q[2];
sx q[2];
rz(-2.6174394) q[2];
sx q[2];
rz(2.1955127) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17543381) q[1];
sx q[1];
rz(-0.81635469) q[1];
sx q[1];
rz(0.50558009) q[1];
x q[2];
rz(-0.10736088) q[3];
sx q[3];
rz(-3.0279362) q[3];
sx q[3];
rz(1.74708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.26479244) q[2];
sx q[2];
rz(-0.53692836) q[2];
sx q[2];
rz(-0.31141591) q[2];
rz(2.9733859) q[3];
sx q[3];
rz(-1.5663389) q[3];
sx q[3];
rz(-0.28347191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51157057) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(-2.2054963) q[0];
rz(2.6452737) q[1];
sx q[1];
rz(-2.4744108) q[1];
sx q[1];
rz(0.47028968) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1449438) q[0];
sx q[0];
rz(-1.2723288) q[0];
sx q[0];
rz(0.51346438) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0025493267) q[2];
sx q[2];
rz(-1.373072) q[2];
sx q[2];
rz(1.5288803) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6185604) q[1];
sx q[1];
rz(-1.9123239) q[1];
sx q[1];
rz(0.073445436) q[1];
x q[2];
rz(-1.3948729) q[3];
sx q[3];
rz(-1.298549) q[3];
sx q[3];
rz(1.7737089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63142598) q[2];
sx q[2];
rz(-1.7743856) q[2];
sx q[2];
rz(1.6660956) q[2];
rz(-0.79536074) q[3];
sx q[3];
rz(-0.16203351) q[3];
sx q[3];
rz(1.2696666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0805761) q[0];
sx q[0];
rz(-0.59122714) q[0];
sx q[0];
rz(0.06037816) q[0];
rz(-0.16054842) q[1];
sx q[1];
rz(-1.5269273) q[1];
sx q[1];
rz(-0.9789595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2706021) q[0];
sx q[0];
rz(-0.63828642) q[0];
sx q[0];
rz(-0.19564512) q[0];
rz(2.8674815) q[2];
sx q[2];
rz(-0.90702552) q[2];
sx q[2];
rz(1.2520777) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1886602) q[1];
sx q[1];
rz(-2.1702936) q[1];
sx q[1];
rz(0.99296928) q[1];
rz(2.1500504) q[3];
sx q[3];
rz(-1.5994306) q[3];
sx q[3];
rz(-0.34285173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6947967) q[0];
sx q[0];
rz(-1.0121166) q[0];
sx q[0];
rz(2.6126675) q[0];
rz(2.9534598) q[1];
sx q[1];
rz(-0.70536047) q[1];
sx q[1];
rz(2.5573152) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3007418) q[0];
sx q[0];
rz(-2.7817431) q[0];
sx q[0];
rz(0.76143439) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7033349) q[2];
sx q[2];
rz(-1.9263679) q[2];
sx q[2];
rz(2.8388765) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6021298) q[1];
sx q[1];
rz(-1.5279084) q[1];
sx q[1];
rz(-0.40956386) q[1];
rz(-pi) q[2];
rz(2.8680471) q[3];
sx q[3];
rz(-2.1714032) q[3];
sx q[3];
rz(-0.41714868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75954413) q[2];
sx q[2];
rz(-1.0039165) q[2];
sx q[2];
rz(0.071852597) q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1857984) q[0];
sx q[0];
rz(-0.6577984) q[0];
sx q[0];
rz(-0.26554769) q[0];
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
