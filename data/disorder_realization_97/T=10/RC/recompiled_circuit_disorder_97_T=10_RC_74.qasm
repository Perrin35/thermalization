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
rz(-2.5118877) q[1];
sx q[1];
rz(-1.1344818) q[1];
sx q[1];
rz(1.1073444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21270277) q[0];
sx q[0];
rz(-1.895993) q[0];
sx q[0];
rz(-2.4341499) q[0];
x q[1];
rz(2.6644457) q[2];
sx q[2];
rz(-0.90847441) q[2];
sx q[2];
rz(1.9805816) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0336696) q[1];
sx q[1];
rz(-0.30089295) q[1];
sx q[1];
rz(-1.184102) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9355448) q[3];
sx q[3];
rz(-1.8555292) q[3];
sx q[3];
rz(0.92393827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0779695) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(1.3280274) q[2];
rz(2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(-3.0096171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48224738) q[0];
sx q[0];
rz(-3.0292065) q[0];
sx q[0];
rz(-0.88062084) q[0];
rz(-1.2940787) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(2.3243288) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1399122) q[0];
sx q[0];
rz(-1.4298555) q[0];
sx q[0];
rz(-1.9190448) q[0];
x q[1];
rz(-1.1772637) q[2];
sx q[2];
rz(-0.52429188) q[2];
sx q[2];
rz(-2.2146068) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6868999) q[1];
sx q[1];
rz(-0.94140879) q[1];
sx q[1];
rz(1.252974) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2137787) q[3];
sx q[3];
rz(-0.5001874) q[3];
sx q[3];
rz(-0.57750765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2237504) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(-0.36402738) q[2];
rz(-2.1552127) q[3];
sx q[3];
rz(-1.4168408) q[3];
sx q[3];
rz(1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7746975) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(2.1752775) q[0];
rz(2.9486588) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(1.6945217) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7265978) q[0];
sx q[0];
rz(-1.5532171) q[0];
sx q[0];
rz(-2.9265762) q[0];
x q[1];
rz(-1.3005199) q[2];
sx q[2];
rz(-2.8405361) q[2];
sx q[2];
rz(-1.2031872) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5406394) q[1];
sx q[1];
rz(-2.1689231) q[1];
sx q[1];
rz(-2.0289434) q[1];
rz(-pi) q[2];
rz(0.2451285) q[3];
sx q[3];
rz(-1.9760625) q[3];
sx q[3];
rz(-2.1851636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7029999) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(2.036371) q[2];
rz(0.74622074) q[3];
sx q[3];
rz(-1.6522224) q[3];
sx q[3];
rz(0.97572774) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0063909) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(1.6756469) q[0];
rz(0.28494596) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(2.7526061) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5327832) q[0];
sx q[0];
rz(-2.7739077) q[0];
sx q[0];
rz(-0.68433783) q[0];
rz(-pi) q[1];
rz(-1.0842501) q[2];
sx q[2];
rz(-1.3975189) q[2];
sx q[2];
rz(1.7692406) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0491838) q[1];
sx q[1];
rz(-1.0256983) q[1];
sx q[1];
rz(-1.070302) q[1];
x q[2];
rz(-1.7763406) q[3];
sx q[3];
rz(-2.3081429) q[3];
sx q[3];
rz(1.8635686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(2.6848865) q[2];
rz(-1.6263973) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.221955) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(0.36002457) q[0];
rz(0.64741627) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(0.49450758) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8202782) q[0];
sx q[0];
rz(-1.2914133) q[0];
sx q[0];
rz(-1.4415635) q[0];
x q[1];
rz(1.2890655) q[2];
sx q[2];
rz(-0.14094555) q[2];
sx q[2];
rz(2.3245387) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9322885) q[1];
sx q[1];
rz(-0.63289019) q[1];
sx q[1];
rz(3.0401405) q[1];
rz(-pi) q[2];
rz(-0.0534119) q[3];
sx q[3];
rz(-0.44465372) q[3];
sx q[3];
rz(-0.85944552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71022025) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(-2.8430856) q[2];
rz(-0.31202894) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(-0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1453778) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(2.9456855) q[0];
rz(3.1205102) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(-1.9063937) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93341953) q[0];
sx q[0];
rz(-2.6801077) q[0];
sx q[0];
rz(-2.3955406) q[0];
rz(0.70448204) q[2];
sx q[2];
rz(-1.8330169) q[2];
sx q[2];
rz(-2.8395677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1393309) q[1];
sx q[1];
rz(-1.1552703) q[1];
sx q[1];
rz(-0.63654391) q[1];
x q[2];
rz(2.9052832) q[3];
sx q[3];
rz(-1.4893388) q[3];
sx q[3];
rz(-2.8636275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6482676) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(2.7098999) q[2];
rz(1.7290944) q[3];
sx q[3];
rz(-2.4192211) q[3];
sx q[3];
rz(3.0055962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7504808) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(-3.0294763) q[0];
rz(-2.926459) q[1];
sx q[1];
rz(-1.5810177) q[1];
sx q[1];
rz(2.0281866) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33681413) q[0];
sx q[0];
rz(-1.830173) q[0];
sx q[0];
rz(2.74385) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68219296) q[2];
sx q[2];
rz(-0.39794121) q[2];
sx q[2];
rz(-1.8856018) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0472764) q[1];
sx q[1];
rz(-0.89135209) q[1];
sx q[1];
rz(1.9654771) q[1];
rz(1.1464305) q[3];
sx q[3];
rz(-0.48454912) q[3];
sx q[3];
rz(-2.4534006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.001174288) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(3.0155638) q[2];
rz(-1.0472939) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66013181) q[0];
sx q[0];
rz(-2.3829057) q[0];
sx q[0];
rz(1.460176) q[0];
rz(-1.2449645) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(1.9326899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7891114) q[0];
sx q[0];
rz(-1.6500104) q[0];
sx q[0];
rz(1.3080025) q[0];
x q[1];
rz(2.4822794) q[2];
sx q[2];
rz(-2.1954143) q[2];
sx q[2];
rz(0.63559947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.804467) q[1];
sx q[1];
rz(-0.85782385) q[1];
sx q[1];
rz(2.0130403) q[1];
rz(-2.2015757) q[3];
sx q[3];
rz(-1.3925843) q[3];
sx q[3];
rz(0.34561397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37844354) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(-1.3195999) q[2];
rz(-0.59213263) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(3.106451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2329344) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(-1.3611025) q[0];
rz(-1.2110442) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(-2.7499054) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91613382) q[0];
sx q[0];
rz(-0.43457169) q[0];
sx q[0];
rz(3.1276914) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7462256) q[2];
sx q[2];
rz(-2.6569416) q[2];
sx q[2];
rz(-1.5026827) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.78545953) q[1];
sx q[1];
rz(-1.2411989) q[1];
sx q[1];
rz(2.3258924) q[1];
rz(-pi) q[2];
rz(2.6158995) q[3];
sx q[3];
rz(-1.5075397) q[3];
sx q[3];
rz(-0.69378187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52035511) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(1.8048145) q[2];
rz(2.3796066) q[3];
sx q[3];
rz(-2.8218994) q[3];
sx q[3];
rz(0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.8005463) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(2.5706932) q[0];
rz(1.7123429) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(0.16194078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83308342) q[0];
sx q[0];
rz(-0.97531318) q[0];
sx q[0];
rz(0.34992976) q[0];
rz(-pi) q[1];
rz(-1.5230721) q[2];
sx q[2];
rz(-2.952791) q[2];
sx q[2];
rz(-1.4023086) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.83878126) q[1];
sx q[1];
rz(-1.6994085) q[1];
sx q[1];
rz(1.3978811) q[1];
rz(-pi) q[2];
x q[2];
rz(0.057823618) q[3];
sx q[3];
rz(-2.6177546) q[3];
sx q[3];
rz(3.0081089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.95742115) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(1.7133678) q[2];
rz(1.1994294) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(1.3706346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7006871) q[0];
sx q[0];
rz(-2.6661243) q[0];
sx q[0];
rz(2.0211924) q[0];
rz(1.7715001) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(-1.7182405) q[2];
sx q[2];
rz(-1.6885919) q[2];
sx q[2];
rz(-0.13292776) q[2];
rz(0.45392848) q[3];
sx q[3];
rz(-0.47149999) q[3];
sx q[3];
rz(-0.75054689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
