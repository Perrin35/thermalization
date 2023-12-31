OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(-1.6516049) q[0];
sx q[0];
rz(-2.2111501) q[0];
rz(0.62970495) q[1];
sx q[1];
rz(-2.0071109) q[1];
sx q[1];
rz(-1.1073444) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005172) q[0];
sx q[0];
rz(-2.374875) q[0];
sx q[0];
rz(2.6630152) q[0];
x q[1];
rz(2.2912575) q[2];
sx q[2];
rz(-1.9413661) q[2];
sx q[2];
rz(-0.71760273) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2335637) q[1];
sx q[1];
rz(-1.6828013) q[1];
sx q[1];
rz(1.2909375) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0290306) q[3];
sx q[3];
rz(-0.68800612) q[3];
sx q[3];
rz(2.1306761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.063623108) q[2];
sx q[2];
rz(-2.412553) q[2];
sx q[2];
rz(1.3280274) q[2];
rz(-0.32087457) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(0.13197556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6593453) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(0.88062084) q[0];
rz(1.2940787) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(-0.81726384) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1999232) q[0];
sx q[0];
rz(-2.7669853) q[0];
sx q[0];
rz(1.1767715) q[0];
rz(1.1772637) q[2];
sx q[2];
rz(-0.52429188) q[2];
sx q[2];
rz(2.2146068) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4546928) q[1];
sx q[1];
rz(-2.2001839) q[1];
sx q[1];
rz(1.8886186) q[1];
rz(-pi) q[2];
rz(-2.0440631) q[3];
sx q[3];
rz(-1.4023937) q[3];
sx q[3];
rz(2.4646204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91784224) q[2];
sx q[2];
rz(-0.68085256) q[2];
sx q[2];
rz(2.7775653) q[2];
rz(2.1552127) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(-1.4646437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36689511) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(-2.1752775) q[0];
rz(-2.9486588) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(1.4470709) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7265978) q[0];
sx q[0];
rz(-1.5883755) q[0];
sx q[0];
rz(-0.21501644) q[0];
rz(-pi) q[1];
rz(3.0588805) q[2];
sx q[2];
rz(-1.2809922) q[2];
sx q[2];
rz(0.92083425) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.88168624) q[1];
sx q[1];
rz(-0.73598624) q[1];
sx q[1];
rz(-2.5658539) q[1];
rz(-2.0855911) q[3];
sx q[3];
rz(-2.6714973) q[3];
sx q[3];
rz(-1.6197636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43859279) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(-2.036371) q[2];
rz(-2.3953719) q[3];
sx q[3];
rz(-1.6522224) q[3];
sx q[3];
rz(0.97572774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0063909) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(1.4659457) q[0];
rz(2.8566467) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(0.38898653) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31127351) q[0];
sx q[0];
rz(-1.8000326) q[0];
sx q[0];
rz(-0.29005187) q[0];
rz(-2.94611) q[2];
sx q[2];
rz(-1.0921548) q[2];
sx q[2];
rz(-0.10749707) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4219141) q[1];
sx q[1];
rz(-0.72243566) q[1];
sx q[1];
rz(2.4721485) q[1];
x q[2];
rz(-2.9205434) q[3];
sx q[3];
rz(-0.76023686) q[3];
sx q[3];
rz(1.5787214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(-0.45670613) q[2];
rz(1.6263973) q[3];
sx q[3];
rz(-0.94907343) q[3];
sx q[3];
rz(-2.8592498) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.221955) q[0];
sx q[0];
rz(-2.0018405) q[0];
sx q[0];
rz(2.7815681) q[0];
rz(2.4941764) q[1];
sx q[1];
rz(-1.5292239) q[1];
sx q[1];
rz(-2.6470851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21365989) q[0];
sx q[0];
rz(-1.4466009) q[0];
sx q[0];
rz(0.28161063) q[0];
rz(-pi) q[1];
rz(-1.2890655) q[2];
sx q[2];
rz(-3.0006471) q[2];
sx q[2];
rz(-0.81705392) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9322885) q[1];
sx q[1];
rz(-2.5087025) q[1];
sx q[1];
rz(3.0401405) q[1];
rz(2.6974929) q[3];
sx q[3];
rz(-1.5937623) q[3];
sx q[3];
rz(-0.75957739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4313724) q[2];
sx q[2];
rz(-2.4857095) q[2];
sx q[2];
rz(-2.8430856) q[2];
rz(-2.8295637) q[3];
sx q[3];
rz(-1.3307064) q[3];
sx q[3];
rz(2.503094) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9962149) q[0];
sx q[0];
rz(-2.4185116) q[0];
sx q[0];
rz(0.19590713) q[0];
rz(-0.021082489) q[1];
sx q[1];
rz(-1.3985876) q[1];
sx q[1];
rz(-1.235199) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0876448) q[0];
sx q[0];
rz(-1.2637648) q[0];
sx q[0];
rz(-0.35015492) q[0];
x q[1];
rz(0.70448204) q[2];
sx q[2];
rz(-1.8330169) q[2];
sx q[2];
rz(0.30202497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0732121) q[1];
sx q[1];
rz(-0.74406032) q[1];
sx q[1];
rz(-0.63853227) q[1];
x q[2];
rz(-2.8060693) q[3];
sx q[3];
rz(-0.24970679) q[3];
sx q[3];
rz(1.6186796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6482676) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(-2.7098999) q[2];
rz(1.4124983) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(-0.13599642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7504808) q[0];
sx q[0];
rz(-1.1688122) q[0];
sx q[0];
rz(-3.0294763) q[0];
rz(-2.926459) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(-2.0281866) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8002692) q[0];
sx q[0];
rz(-1.1870664) q[0];
sx q[0];
rz(1.8510438) q[0];
rz(-pi) q[1];
rz(-1.8298803) q[2];
sx q[2];
rz(-1.2652745) q[2];
sx q[2];
rz(-2.6079026) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0943162) q[1];
sx q[1];
rz(-0.89135209) q[1];
sx q[1];
rz(-1.1761155) q[1];
rz(-0.21344276) q[3];
sx q[3];
rz(-1.1323954) q[3];
sx q[3];
rz(1.9813117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1404184) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(-3.0155638) q[2];
rz(-2.0942988) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(-0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66013181) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(1.6814167) q[0];
rz(-1.8966282) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(-1.2089027) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93223244) q[0];
sx q[0];
rz(-2.8673842) q[0];
sx q[0];
rz(-1.2742395) q[0];
x q[1];
rz(2.3102343) q[2];
sx q[2];
rz(-2.0908329) q[2];
sx q[2];
rz(-1.3607197) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.066684494) q[1];
sx q[1];
rz(-1.9004596) q[1];
sx q[1];
rz(0.7633022) q[1];
rz(-1.8672089) q[3];
sx q[3];
rz(-2.489438) q[3];
sx q[3];
rz(1.46331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.37844354) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(-1.8219927) q[2];
rz(0.59213263) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(0.035141703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9086583) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(-1.7804902) q[0];
rz(-1.2110442) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(0.39168721) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2254588) q[0];
sx q[0];
rz(-2.707021) q[0];
sx q[0];
rz(-0.01390121) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6892745) q[2];
sx q[2];
rz(-1.3903793) q[2];
sx q[2];
rz(-0.42186055) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45410941) q[1];
sx q[1];
rz(-0.81067639) q[1];
sx q[1];
rz(-2.0337385) q[1];
x q[2];
rz(-3.016032) q[3];
sx q[3];
rz(-2.6124622) q[3];
sx q[3];
rz(0.76847968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52035511) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(-1.8048145) q[2];
rz(2.3796066) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(-0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-15/(14*pi)) q[0];
sx q[0];
rz(-0.30277345) q[0];
sx q[0];
rz(2.5706932) q[0];
rz(1.4292498) q[1];
sx q[1];
rz(-2.0776904) q[1];
sx q[1];
rz(-2.9796519) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83308342) q[0];
sx q[0];
rz(-2.1662795) q[0];
sx q[0];
rz(-0.34992976) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3822046) q[2];
sx q[2];
rz(-1.5618426) q[2];
sx q[2];
rz(0.21536516) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3028114) q[1];
sx q[1];
rz(-1.4421842) q[1];
sx q[1];
rz(-1.7437115) q[1];
x q[2];
rz(-0.52311388) q[3];
sx q[3];
rz(-1.5418846) q[3];
sx q[3];
rz(-1.6541964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1841715) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(1.7133678) q[2];
rz(-1.9421633) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(1.7709581) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7006871) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(-1.7715001) q[1];
sx q[1];
rz(-2.1961828) q[1];
sx q[1];
rz(-0.97074769) q[1];
rz(-1.7182405) q[2];
sx q[2];
rz(-1.6885919) q[2];
sx q[2];
rz(-0.13292776) q[2];
rz(1.7907501) q[3];
sx q[3];
rz(-1.9913047) q[3];
sx q[3];
rz(-1.2515968) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
