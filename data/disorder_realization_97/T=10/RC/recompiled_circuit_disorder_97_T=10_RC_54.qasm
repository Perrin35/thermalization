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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5168092) q[0];
sx q[0];
rz(-0.9073572) q[0];
sx q[0];
rz(-1.1532564) q[0];
rz(-0.47714699) q[2];
sx q[2];
rz(-2.2331182) q[2];
sx q[2];
rz(-1.9805816) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5109374) q[1];
sx q[1];
rz(-1.8488548) q[1];
sx q[1];
rz(3.0250938) q[1];
rz(-pi) q[2];
rz(-2.2060478) q[3];
sx q[3];
rz(-1.8555292) q[3];
sx q[3];
rz(-0.92393827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.063623108) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(-1.3280274) q[2];
rz(-2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(3.0096171) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6593453) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(2.2609718) q[0];
rz(-1.847514) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(-0.81726384) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1999232) q[0];
sx q[0];
rz(-0.37460735) q[0];
sx q[0];
rz(1.1767715) q[0];
rz(-pi) q[1];
rz(-0.21821071) q[2];
sx q[2];
rz(-2.0514224) q[2];
sx q[2];
rz(2.6618119) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6868999) q[1];
sx q[1];
rz(-2.2001839) q[1];
sx q[1];
rz(1.8886186) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18873429) q[3];
sx q[3];
rz(-1.1047603) q[3];
sx q[3];
rz(-2.1621494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91784224) q[2];
sx q[2];
rz(-0.68085256) q[2];
sx q[2];
rz(-2.7775653) q[2];
rz(-0.98637995) q[3];
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
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7746975) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(-0.96631518) q[0];
rz(-2.9486588) q[1];
sx q[1];
rz(-1.0886334) q[1];
sx q[1];
rz(-1.4470709) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7265978) q[0];
sx q[0];
rz(-1.5532171) q[0];
sx q[0];
rz(-2.9265762) q[0];
rz(-pi) q[1];
rz(1.3005199) q[2];
sx q[2];
rz(-2.8405361) q[2];
sx q[2];
rz(-1.9384055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24070534) q[1];
sx q[1];
rz(-1.9449688) q[1];
sx q[1];
rz(0.64970533) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8964642) q[3];
sx q[3];
rz(-1.9760625) q[3];
sx q[3];
rz(0.95642904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43859279) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(-1.1052216) q[2];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352017) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(1.6756469) q[0];
rz(2.8566467) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(-2.7526061) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8143512) q[0];
sx q[0];
rz(-1.8530493) q[0];
sx q[0];
rz(1.3319356) q[0];
rz(-2.94611) q[2];
sx q[2];
rz(-1.0921548) q[2];
sx q[2];
rz(3.0340956) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71967857) q[1];
sx q[1];
rz(-0.72243566) q[1];
sx q[1];
rz(-0.6694442) q[1];
x q[2];
rz(0.22104927) q[3];
sx q[3];
rz(-2.3813558) q[3];
sx q[3];
rz(-1.5787214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4185562) q[2];
sx q[2];
rz(-1.1636461) q[2];
sx q[2];
rz(-0.45670613) q[2];
rz(-1.6263973) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(-2.8592498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.221955) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(-0.36002457) q[0];
rz(0.64741627) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(-2.6470851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3798824) q[0];
sx q[0];
rz(-2.8344791) q[0];
sx q[0];
rz(2.7193927) q[0];
rz(1.2890655) q[2];
sx q[2];
rz(-0.14094555) q[2];
sx q[2];
rz(-0.81705392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33489409) q[1];
sx q[1];
rz(-2.1999199) q[1];
sx q[1];
rz(-1.4966399) q[1];
x q[2];
rz(-3.0881808) q[3];
sx q[3];
rz(-2.6969389) q[3];
sx q[3];
rz(-0.85944552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4313724) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(-0.29850706) q[2];
rz(-0.31202894) q[3];
sx q[3];
rz(-1.3307064) q[3];
sx q[3];
rz(-2.503094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1453778) q[0];
sx q[0];
rz(-2.4185116) q[0];
sx q[0];
rz(0.19590713) q[0];
rz(-3.1205102) q[1];
sx q[1];
rz(-1.3985876) q[1];
sx q[1];
rz(-1.9063937) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7346749) q[0];
sx q[0];
rz(-1.2376681) q[0];
sx q[0];
rz(1.245265) q[0];
rz(-pi) q[1];
rz(0.70448204) q[2];
sx q[2];
rz(-1.3085758) q[2];
sx q[2];
rz(2.8395677) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0022618) q[1];
sx q[1];
rz(-1.9863223) q[1];
sx q[1];
rz(-0.63654391) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9052832) q[3];
sx q[3];
rz(-1.4893388) q[3];
sx q[3];
rz(2.8636275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49332508) q[2];
sx q[2];
rz(-0.92612925) q[2];
sx q[2];
rz(-2.7098999) q[2];
rz(1.4124983) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(3.0055962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39111185) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(3.0294763) q[0];
rz(-2.926459) q[1];
sx q[1];
rz(-1.5810177) q[1];
sx q[1];
rz(-1.1134061) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68590251) q[0];
sx q[0];
rz(-0.4710353) q[0];
sx q[0];
rz(2.5409565) q[0];
x q[1];
rz(0.68219296) q[2];
sx q[2];
rz(-0.39794121) q[2];
sx q[2];
rz(1.2559909) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4091195) q[1];
sx q[1];
rz(-1.8745683) q[1];
sx q[1];
rz(-2.4227546) q[1];
rz(-1.1464305) q[3];
sx q[3];
rz(-2.6570435) q[3];
sx q[3];
rz(0.68819203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.001174288) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(-0.12602885) q[2];
rz(2.0942988) q[3];
sx q[3];
rz(-1.8208561) q[3];
sx q[3];
rz(-0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.4814608) q[0];
sx q[0];
rz(-2.3829057) q[0];
sx q[0];
rz(-1.460176) q[0];
rz(-1.2449645) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(-1.2089027) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2395994) q[0];
sx q[0];
rz(-1.308846) q[0];
sx q[0];
rz(-3.0595747) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4822794) q[2];
sx q[2];
rz(-0.94617832) q[2];
sx q[2];
rz(-0.63559947) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.066684494) q[1];
sx q[1];
rz(-1.2411331) q[1];
sx q[1];
rz(0.7633022) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94001694) q[3];
sx q[3];
rz(-1.7490083) q[3];
sx q[3];
rz(-0.34561397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37844354) q[2];
sx q[2];
rz(-1.903879) q[2];
sx q[2];
rz(-1.3195999) q[2];
rz(0.59213263) q[3];
sx q[3];
rz(-1.7254646) q[3];
sx q[3];
rz(3.106451) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2329344) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(1.7804902) q[0];
rz(1.2110442) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(-0.39168721) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93145934) q[0];
sx q[0];
rz(-1.1362695) q[0];
sx q[0];
rz(-1.5772485) q[0];
rz(0.39536706) q[2];
sx q[2];
rz(-2.6569416) q[2];
sx q[2];
rz(-1.5026827) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6874832) q[1];
sx q[1];
rz(-2.3309163) q[1];
sx q[1];
rz(-2.0337385) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.016032) q[3];
sx q[3];
rz(-2.6124622) q[3];
sx q[3];
rz(-2.373113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6212375) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(1.8048145) q[2];
rz(2.3796066) q[3];
sx q[3];
rz(-2.8218994) q[3];
sx q[3];
rz(-2.3412162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(-0.57089943) q[0];
rz(1.7123429) q[1];
sx q[1];
rz(-2.0776904) q[1];
sx q[1];
rz(-0.16194078) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8853332) q[0];
sx q[0];
rz(-0.67978871) q[0];
sx q[0];
rz(-1.1023561) q[0];
x q[1];
rz(1.5230721) q[2];
sx q[2];
rz(-2.952791) q[2];
sx q[2];
rz(-1.7392841) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.098617741) q[1];
sx q[1];
rz(-2.9264755) q[1];
sx q[1];
rz(-0.92623644) q[1];
rz(-pi) q[2];
rz(3.083769) q[3];
sx q[3];
rz(-2.6177546) q[3];
sx q[3];
rz(-3.0081089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95742115) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(-1.4282248) q[2];
rz(1.9421633) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(1.3706346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-1.4233521) q[2];
sx q[2];
rz(-1.4530008) q[2];
sx q[2];
rz(3.0086649) q[2];
rz(-1.3508425) q[3];
sx q[3];
rz(-1.9913047) q[3];
sx q[3];
rz(-1.2515968) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
