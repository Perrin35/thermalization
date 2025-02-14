OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8074789) q[0];
sx q[0];
rz(-0.69598323) q[0];
sx q[0];
rz(-0.88710436) q[0];
rz(2.4298985) q[1];
sx q[1];
rz(-2.0790172) q[1];
sx q[1];
rz(0.33023155) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7831551) q[0];
sx q[0];
rz(-0.81029725) q[0];
sx q[0];
rz(-2.1970219) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55203931) q[2];
sx q[2];
rz(-0.62969724) q[2];
sx q[2];
rz(-2.3384467) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2621999) q[1];
sx q[1];
rz(-2.6385251) q[1];
sx q[1];
rz(-0.21299739) q[1];
rz(-1.7403148) q[3];
sx q[3];
rz(-2.7054686) q[3];
sx q[3];
rz(2.2317302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7547138) q[2];
sx q[2];
rz(-0.42749307) q[2];
sx q[2];
rz(-1.0270366) q[2];
rz(2.2586281) q[3];
sx q[3];
rz(-1.797978) q[3];
sx q[3];
rz(-0.55356717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.7850194) q[0];
sx q[0];
rz(-3.0632601) q[0];
sx q[0];
rz(2.8739492) q[0];
rz(-1.74125) q[1];
sx q[1];
rz(-0.61336556) q[1];
sx q[1];
rz(-0.5563446) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2927383) q[0];
sx q[0];
rz(-1.8195099) q[0];
sx q[0];
rz(-1.066904) q[0];
rz(-pi) q[1];
rz(0.39926932) q[2];
sx q[2];
rz(-0.74490163) q[2];
sx q[2];
rz(0.58303729) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9832335) q[1];
sx q[1];
rz(-2.4438847) q[1];
sx q[1];
rz(-2.1844728) q[1];
x q[2];
rz(2.2675927) q[3];
sx q[3];
rz(-0.79588875) q[3];
sx q[3];
rz(-2.9615642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9050425) q[2];
sx q[2];
rz(-1.4507797) q[2];
sx q[2];
rz(1.2518008) q[2];
rz(1.8768138) q[3];
sx q[3];
rz(-0.71139657) q[3];
sx q[3];
rz(-1.0714162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4816137) q[0];
sx q[0];
rz(-2.5857506) q[0];
sx q[0];
rz(-0.83928338) q[0];
rz(1.6370157) q[1];
sx q[1];
rz(-1.4312276) q[1];
sx q[1];
rz(1.379871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8913075) q[0];
sx q[0];
rz(-1.4472741) q[0];
sx q[0];
rz(3.059292) q[0];
rz(-pi) q[1];
rz(-1.2261058) q[2];
sx q[2];
rz(-2.569103) q[2];
sx q[2];
rz(2.9782611) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8306337) q[1];
sx q[1];
rz(-1.1006025) q[1];
sx q[1];
rz(-2.7814835) q[1];
rz(-pi) q[2];
rz(1.0187418) q[3];
sx q[3];
rz(-1.0699439) q[3];
sx q[3];
rz(-0.93397442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.442349) q[2];
sx q[2];
rz(-1.2171429) q[2];
sx q[2];
rz(2.4231363) q[2];
rz(0.62026223) q[3];
sx q[3];
rz(-0.93554997) q[3];
sx q[3];
rz(-0.81234318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45850596) q[0];
sx q[0];
rz(-1.2309256) q[0];
sx q[0];
rz(-1.7339535) q[0];
rz(2.5900224) q[1];
sx q[1];
rz(-3.0307814) q[1];
sx q[1];
rz(1.3161906) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16770076) q[0];
sx q[0];
rz(-2.1657054) q[0];
sx q[0];
rz(0.3428726) q[0];
x q[1];
rz(1.8502683) q[2];
sx q[2];
rz(-1.5588461) q[2];
sx q[2];
rz(-0.751647) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0312457) q[1];
sx q[1];
rz(-1.7025885) q[1];
sx q[1];
rz(2.736997) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85547282) q[3];
sx q[3];
rz(-1.0294204) q[3];
sx q[3];
rz(-2.4349041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5482285) q[2];
sx q[2];
rz(-0.43539384) q[2];
sx q[2];
rz(1.0080053) q[2];
rz(-0.27901444) q[3];
sx q[3];
rz(-0.81441003) q[3];
sx q[3];
rz(-2.6830955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8654883) q[0];
sx q[0];
rz(-1.3300329) q[0];
sx q[0];
rz(-2.8651067) q[0];
rz(-0.28929389) q[1];
sx q[1];
rz(-1.1951059) q[1];
sx q[1];
rz(2.8791265) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.58231) q[0];
sx q[0];
rz(-1.773343) q[0];
sx q[0];
rz(1.3900533) q[0];
x q[1];
rz(2.1621373) q[2];
sx q[2];
rz(-1.5797414) q[2];
sx q[2];
rz(1.1696377) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20690675) q[1];
sx q[1];
rz(-1.0240203) q[1];
sx q[1];
rz(1.9711167) q[1];
rz(-pi) q[2];
rz(-1.4175884) q[3];
sx q[3];
rz(-0.99780449) q[3];
sx q[3];
rz(0.38148034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7485973) q[2];
sx q[2];
rz(-0.59610072) q[2];
sx q[2];
rz(-1.3470915) q[2];
rz(0.10284452) q[3];
sx q[3];
rz(-1.9072073) q[3];
sx q[3];
rz(1.5170521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0031925072) q[0];
sx q[0];
rz(-1.6380558) q[0];
sx q[0];
rz(0.060977161) q[0];
rz(1.2334476) q[1];
sx q[1];
rz(-1.3950709) q[1];
sx q[1];
rz(-1.6159509) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2642941) q[0];
sx q[0];
rz(-1.0420615) q[0];
sx q[0];
rz(-1.651047) q[0];
rz(1.5360918) q[2];
sx q[2];
rz(-2.3922709) q[2];
sx q[2];
rz(2.9108436) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96887302) q[1];
sx q[1];
rz(-0.27471009) q[1];
sx q[1];
rz(1.9897047) q[1];
rz(2.9750729) q[3];
sx q[3];
rz(-2.2158479) q[3];
sx q[3];
rz(1.695961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47152758) q[2];
sx q[2];
rz(-1.3482886) q[2];
sx q[2];
rz(1.8111551) q[2];
rz(1.1250251) q[3];
sx q[3];
rz(-0.87441134) q[3];
sx q[3];
rz(0.037954656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9576981) q[0];
sx q[0];
rz(-1.6989468) q[0];
sx q[0];
rz(-2.884602) q[0];
rz(0.43486241) q[1];
sx q[1];
rz(-2.3915274) q[1];
sx q[1];
rz(0.90352568) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5241476) q[0];
sx q[0];
rz(-0.092369583) q[0];
sx q[0];
rz(-1.525982) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9952094) q[2];
sx q[2];
rz(-2.3902262) q[2];
sx q[2];
rz(2.9379972) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.74106797) q[1];
sx q[1];
rz(-1.6323106) q[1];
sx q[1];
rz(1.6183302) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1210455) q[3];
sx q[3];
rz(-1.9686832) q[3];
sx q[3];
rz(2.3031421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8187108) q[2];
sx q[2];
rz(-1.633753) q[2];
sx q[2];
rz(1.7203077) q[2];
rz(0.70982248) q[3];
sx q[3];
rz(-2.0028508) q[3];
sx q[3];
rz(0.53409725) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98599559) q[0];
sx q[0];
rz(-0.43102145) q[0];
sx q[0];
rz(-1.0850329) q[0];
rz(2.494508) q[1];
sx q[1];
rz(-1.1754464) q[1];
sx q[1];
rz(-3.0771902) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0376799) q[0];
sx q[0];
rz(-1.2453916) q[0];
sx q[0];
rz(2.4766981) q[0];
rz(-0.66775457) q[2];
sx q[2];
rz(-0.96558324) q[2];
sx q[2];
rz(-0.79082327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6820087) q[1];
sx q[1];
rz(-1.1136076) q[1];
sx q[1];
rz(-2.3383635) q[1];
rz(1.9265979) q[3];
sx q[3];
rz(-1.3094153) q[3];
sx q[3];
rz(1.3047117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.25962466) q[2];
sx q[2];
rz(-0.3825863) q[2];
sx q[2];
rz(3.0982223) q[2];
rz(0.67443332) q[3];
sx q[3];
rz(-1.784227) q[3];
sx q[3];
rz(-1.9969214) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44716537) q[0];
sx q[0];
rz(-0.58037102) q[0];
sx q[0];
rz(0.32661435) q[0];
rz(-1.0940374) q[1];
sx q[1];
rz(-2.2973165) q[1];
sx q[1];
rz(0.48386595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.992386) q[0];
sx q[0];
rz(-1.0629553) q[0];
sx q[0];
rz(-0.81224982) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3442114) q[2];
sx q[2];
rz(-1.1562776) q[2];
sx q[2];
rz(-0.58794566) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6623508) q[1];
sx q[1];
rz(-1.7852946) q[1];
sx q[1];
rz(0.7995601) q[1];
rz(-1.3392162) q[3];
sx q[3];
rz(-0.94144034) q[3];
sx q[3];
rz(-1.6762102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6924374) q[2];
sx q[2];
rz(-2.1026976) q[2];
sx q[2];
rz(0.078744002) q[2];
rz(0.76830831) q[3];
sx q[3];
rz(-1.2945622) q[3];
sx q[3];
rz(2.9505742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.47062) q[0];
sx q[0];
rz(-0.60523954) q[0];
sx q[0];
rz(-3.0249) q[0];
rz(2.281588) q[1];
sx q[1];
rz(-1.8000894) q[1];
sx q[1];
rz(-0.87342993) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54166543) q[0];
sx q[0];
rz(-2.3573236) q[0];
sx q[0];
rz(-1.2265217) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8928746) q[2];
sx q[2];
rz(-1.0828185) q[2];
sx q[2];
rz(-2.8805594) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8598249) q[1];
sx q[1];
rz(-1.072653) q[1];
sx q[1];
rz(-1.7120275) q[1];
x q[2];
rz(0.61140538) q[3];
sx q[3];
rz(-1.6263925) q[3];
sx q[3];
rz(1.9724595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7200835) q[2];
sx q[2];
rz(-2.8953711) q[2];
sx q[2];
rz(-2.1923547) q[2];
rz(1.2468437) q[3];
sx q[3];
rz(-1.6722164) q[3];
sx q[3];
rz(2.5505572) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.128189) q[0];
sx q[0];
rz(-1.0008151) q[0];
sx q[0];
rz(-1.6304954) q[0];
rz(-0.36698256) q[1];
sx q[1];
rz(-1.783168) q[1];
sx q[1];
rz(2.5617243) q[1];
rz(1.0147711) q[2];
sx q[2];
rz(-2.6064739) q[2];
sx q[2];
rz(0.50630416) q[2];
rz(2.4511724) q[3];
sx q[3];
rz(-1.2391117) q[3];
sx q[3];
rz(1.5285872) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
