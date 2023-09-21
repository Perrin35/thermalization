OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4176183) q[0];
sx q[0];
rz(-1.4899878) q[0];
sx q[0];
rz(2.2111501) q[0];
rz(0.62970495) q[1];
sx q[1];
rz(-2.0071109) q[1];
sx q[1];
rz(-1.1073444) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21270277) q[0];
sx q[0];
rz(-1.895993) q[0];
sx q[0];
rz(0.70744275) q[0];
rz(-2.1030175) q[2];
sx q[2];
rz(-0.79471171) q[2];
sx q[2];
rz(-2.6796535) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.107923) q[1];
sx q[1];
rz(-0.30089295) q[1];
sx q[1];
rz(1.184102) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0290306) q[3];
sx q[3];
rz(-0.68800612) q[3];
sx q[3];
rz(1.0109166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0779695) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(1.3280274) q[2];
rz(0.32087457) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(-0.13197556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6593453) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(0.88062084) q[0];
rz(1.847514) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(0.81726384) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9416695) q[0];
sx q[0];
rz(-2.7669853) q[0];
sx q[0];
rz(-1.1767715) q[0];
rz(-1.9643289) q[2];
sx q[2];
rz(-0.52429188) q[2];
sx q[2];
rz(2.2146068) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0751794) q[1];
sx q[1];
rz(-1.315409) q[1];
sx q[1];
rz(2.4875719) q[1];
rz(-pi) q[2];
rz(-2.0440631) q[3];
sx q[3];
rz(-1.739199) q[3];
sx q[3];
rz(-2.4646204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2237504) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(2.7775653) q[2];
rz(-2.1552127) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(-1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7746975) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(-2.1752775) q[0];
rz(2.9486588) q[1];
sx q[1];
rz(-1.0886334) q[1];
sx q[1];
rz(-1.6945217) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23611785) q[0];
sx q[0];
rz(-2.9258699) q[0];
sx q[0];
rz(-0.082213684) q[0];
x q[1];
rz(1.3005199) q[2];
sx q[2];
rz(-0.30105653) q[2];
sx q[2];
rz(-1.2031872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5406394) q[1];
sx q[1];
rz(-0.97266957) q[1];
sx q[1];
rz(-2.0289434) q[1];
x q[2];
rz(1.0560016) q[3];
sx q[3];
rz(-0.47009531) q[3];
sx q[3];
rz(1.6197636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7029999) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(2.036371) q[2];
rz(2.3953719) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(-2.1658649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352017) q[0];
sx q[0];
rz(-0.52440301) q[0];
sx q[0];
rz(-1.4659457) q[0];
rz(2.8566467) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(-2.7526061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60880946) q[0];
sx q[0];
rz(-2.7739077) q[0];
sx q[0];
rz(-2.4572548) q[0];
rz(-pi) q[1];
rz(2.94611) q[2];
sx q[2];
rz(-2.0494378) q[2];
sx q[2];
rz(3.0340956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.092408808) q[1];
sx q[1];
rz(-1.0256983) q[1];
sx q[1];
rz(-1.070302) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9205434) q[3];
sx q[3];
rz(-0.76023686) q[3];
sx q[3];
rz(1.5628712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(-2.6848865) q[2];
rz(1.5151954) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91963768) q[0];
sx q[0];
rz(-2.0018405) q[0];
sx q[0];
rz(2.7815681) q[0];
rz(0.64741627) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(0.49450758) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8202782) q[0];
sx q[0];
rz(-1.2914133) q[0];
sx q[0];
rz(1.4415635) q[0];
x q[1];
rz(-1.8525271) q[2];
sx q[2];
rz(-0.14094555) q[2];
sx q[2];
rz(-0.81705392) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8620017) q[1];
sx q[1];
rz(-1.5108567) q[1];
sx q[1];
rz(2.5111591) q[1];
rz(-pi) q[2];
rz(-2.6974929) q[3];
sx q[3];
rz(-1.5937623) q[3];
sx q[3];
rz(0.75957739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71022025) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(-2.8430856) q[2];
rz(-0.31202894) q[3];
sx q[3];
rz(-1.3307064) q[3];
sx q[3];
rz(0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1453778) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(0.19590713) q[0];
rz(3.1205102) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(-1.9063937) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0876448) q[0];
sx q[0];
rz(-1.2637648) q[0];
sx q[0];
rz(-2.7914377) q[0];
x q[1];
rz(-1.2321129) q[2];
sx q[2];
rz(-0.89502305) q[2];
sx q[2];
rz(2.0896926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8630353) q[1];
sx q[1];
rz(-0.99579358) q[1];
sx q[1];
rz(1.0689736) q[1];
rz(-pi) q[2];
rz(-1.6545717) q[3];
sx q[3];
rz(-1.8063074) q[3];
sx q[3];
rz(-4/(1*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6482676) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(-2.7098999) q[2];
rz(-1.4124983) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(-3.0055962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39111185) q[0];
sx q[0];
rz(-1.1688122) q[0];
sx q[0];
rz(0.11211638) q[0];
rz(2.926459) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(2.0281866) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33681413) q[0];
sx q[0];
rz(-1.830173) q[0];
sx q[0];
rz(2.74385) q[0];
rz(-2.4593997) q[2];
sx q[2];
rz(-2.7436514) q[2];
sx q[2];
rz(1.8856018) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4091195) q[1];
sx q[1];
rz(-1.2670244) q[1];
sx q[1];
rz(-0.71883808) q[1];
x q[2];
rz(-1.1464305) q[3];
sx q[3];
rz(-2.6570435) q[3];
sx q[3];
rz(-2.4534006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1404184) q[2];
sx q[2];
rz(-0.73264709) q[2];
sx q[2];
rz(-3.0155638) q[2];
rz(2.0942988) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-1.0943202) q[1];
sx q[1];
rz(1.2089027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9019933) q[0];
sx q[0];
rz(-1.308846) q[0];
sx q[0];
rz(-0.082017935) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4822794) q[2];
sx q[2];
rz(-2.1954143) q[2];
sx q[2];
rz(0.63559947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.804467) q[1];
sx q[1];
rz(-2.2837688) q[1];
sx q[1];
rz(2.0130403) q[1];
rz(-1.2743837) q[3];
sx q[3];
rz(-0.65215462) q[3];
sx q[3];
rz(1.46331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7631491) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(-1.8219927) q[2];
rz(0.59213263) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(-3.106451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2329344) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(1.3611025) q[0];
rz(-1.2110442) q[1];
sx q[1];
rz(-2.2369604) q[1];
sx q[1];
rz(2.7499054) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64205326) q[0];
sx q[0];
rz(-1.5649438) q[0];
sx q[0];
rz(-0.43453479) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45231818) q[2];
sx q[2];
rz(-1.7512133) q[2];
sx q[2];
rz(-2.7197321) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3561331) q[1];
sx q[1];
rz(-1.9003938) q[1];
sx q[1];
rz(2.3258924) q[1];
x q[2];
rz(0.52569315) q[3];
sx q[3];
rz(-1.5075397) q[3];
sx q[3];
rz(0.69378187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52035511) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(-1.8048145) q[2];
rz(-2.3796066) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-15/(14*pi)) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(-2.5706932) q[0];
rz(1.7123429) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(0.16194078) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8853332) q[0];
sx q[0];
rz(-0.67978871) q[0];
sx q[0];
rz(-1.1023561) q[0];
x q[1];
rz(-1.759388) q[2];
sx q[2];
rz(-1.5618426) q[2];
sx q[2];
rz(-2.9262275) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0429749) q[1];
sx q[1];
rz(-0.21511714) q[1];
sx q[1];
rz(0.92623644) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5374244) q[3];
sx q[3];
rz(-2.0936692) q[3];
sx q[3];
rz(3.0748622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.95742115) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(1.7133678) q[2];
rz(-1.9421633) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(1.3706346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4409055) q[0];
sx q[0];
rz(-2.6661243) q[0];
sx q[0];
rz(2.0211924) q[0];
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