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
rz(-2.5118877) q[1];
sx q[1];
rz(-1.1344818) q[1];
sx q[1];
rz(-2.0342483) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6247834) q[0];
sx q[0];
rz(-0.9073572) q[0];
sx q[0];
rz(-1.9883363) q[0];
rz(-pi) q[1];
rz(0.8503352) q[2];
sx q[2];
rz(-1.2002266) q[2];
sx q[2];
rz(2.4239899) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.107923) q[1];
sx q[1];
rz(-2.8406997) q[1];
sx q[1];
rz(1.184102) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9355448) q[3];
sx q[3];
rz(-1.2860635) q[3];
sx q[3];
rz(2.2176544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.063623108) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(1.3280274) q[2];
rz(0.32087457) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(3.0096171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48224738) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(-0.88062084) q[0];
rz(-1.847514) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(2.3243288) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5215223) q[0];
sx q[0];
rz(-1.9154473) q[0];
sx q[0];
rz(-0.14981139) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9233819) q[2];
sx q[2];
rz(-2.0514224) q[2];
sx q[2];
rz(2.6618119) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4546928) q[1];
sx q[1];
rz(-2.2001839) q[1];
sx q[1];
rz(-1.8886186) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9278139) q[3];
sx q[3];
rz(-0.5001874) q[3];
sx q[3];
rz(-2.564085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.91784224) q[2];
sx q[2];
rz(-0.68085256) q[2];
sx q[2];
rz(0.36402738) q[2];
rz(-0.98637995) q[3];
sx q[3];
rz(-1.4168408) q[3];
sx q[3];
rz(1.4646437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7746975) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(-0.96631518) q[0];
rz(2.9486588) q[1];
sx q[1];
rz(-1.0886334) q[1];
sx q[1];
rz(-1.6945217) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9896302) q[0];
sx q[0];
rz(-1.3558136) q[0];
sx q[0];
rz(1.5887898) q[0];
rz(-pi) q[1];
rz(-3.0588805) q[2];
sx q[2];
rz(-1.8606004) q[2];
sx q[2];
rz(0.92083425) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88168624) q[1];
sx q[1];
rz(-2.4056064) q[1];
sx q[1];
rz(2.5658539) q[1];
rz(-1.0560016) q[3];
sx q[3];
rz(-2.6714973) q[3];
sx q[3];
rz(-1.521829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43859279) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(1.1052216) q[2];
rz(2.3953719) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(0.97572774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0063909) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(1.4659457) q[0];
rz(-0.28494596) q[1];
sx q[1];
rz(-1.0703215) q[1];
sx q[1];
rz(2.7526061) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31127351) q[0];
sx q[0];
rz(-1.3415601) q[0];
sx q[0];
rz(2.8515408) q[0];
x q[1];
rz(-2.0573425) q[2];
sx q[2];
rz(-1.3975189) q[2];
sx q[2];
rz(1.372352) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3868689) q[1];
sx q[1];
rz(-1.1479953) q[1];
sx q[1];
rz(2.5368284) q[1];
rz(-pi) q[2];
rz(0.74794482) q[3];
sx q[3];
rz(-1.4191295) q[3];
sx q[3];
rz(-0.15350728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(0.45670613) q[2];
rz(1.6263973) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(-0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91963768) q[0];
sx q[0];
rz(-2.0018405) q[0];
sx q[0];
rz(-2.7815681) q[0];
rz(-2.4941764) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(-2.6470851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7617103) q[0];
sx q[0];
rz(-0.30711353) q[0];
sx q[0];
rz(2.7193927) q[0];
x q[1];
rz(-1.8525271) q[2];
sx q[2];
rz(-0.14094555) q[2];
sx q[2];
rz(-0.81705392) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9322885) q[1];
sx q[1];
rz(-2.5087025) q[1];
sx q[1];
rz(0.10145213) q[1];
rz(-0.0534119) q[3];
sx q[3];
rz(-2.6969389) q[3];
sx q[3];
rz(-2.2821471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.71022025) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(-0.29850706) q[2];
rz(-2.8295637) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9962149) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(0.19590713) q[0];
rz(3.1205102) q[1];
sx q[1];
rz(-1.3985876) q[1];
sx q[1];
rz(-1.235199) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053947833) q[0];
sx q[0];
rz(-1.2637648) q[0];
sx q[0];
rz(2.7914377) q[0];
rz(-pi) q[1];
rz(0.39288315) q[2];
sx q[2];
rz(-0.74379951) q[2];
sx q[2];
rz(-1.5768029) q[2];
rz(-pi) q[3];
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
rz(2.5050487) q[1];
rz(-pi) q[2];
rz(2.8060693) q[3];
sx q[3];
rz(-0.24970679) q[3];
sx q[3];
rz(1.5229131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49332508) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(2.7098999) q[2];
rz(1.7290944) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(0.13599642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39111185) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(-0.11211638) q[0];
rz(-2.926459) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(1.1134061) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8002692) q[0];
sx q[0];
rz(-1.1870664) q[0];
sx q[0];
rz(1.2905489) q[0];
x q[1];
rz(1.3117123) q[2];
sx q[2];
rz(-1.8763181) q[2];
sx q[2];
rz(-0.53369001) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
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
rz(0.68819203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1404184) q[2];
sx q[2];
rz(-0.73264709) q[2];
sx q[2];
rz(-3.0155638) q[2];
rz(-1.0472939) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66013181) q[0];
sx q[0];
rz(-2.3829057) q[0];
sx q[0];
rz(-1.460176) q[0];
rz(1.2449645) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(-1.9326899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35248127) q[0];
sx q[0];
rz(-1.6500104) q[0];
sx q[0];
rz(-1.8335901) q[0];
rz(-2.4822794) q[2];
sx q[2];
rz(-2.1954143) q[2];
sx q[2];
rz(-0.63559947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.804467) q[1];
sx q[1];
rz(-2.2837688) q[1];
sx q[1];
rz(2.0130403) q[1];
rz(-pi) q[2];
rz(0.21945159) q[3];
sx q[3];
rz(-0.95153522) q[3];
sx q[3];
rz(-1.0964364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7631491) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(-1.8219927) q[2];
rz(-0.59213263) q[3];
sx q[3];
rz(-1.7254646) q[3];
sx q[3];
rz(-3.106451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2329344) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(1.3611025) q[0];
rz(-1.9305485) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(-0.39168721) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93145934) q[0];
sx q[0];
rz(-1.1362695) q[0];
sx q[0];
rz(1.5772485) q[0];
rz(-pi) q[1];
rz(-0.39536706) q[2];
sx q[2];
rz(-0.48465109) q[2];
sx q[2];
rz(1.6389099) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6874832) q[1];
sx q[1];
rz(-2.3309163) q[1];
sx q[1];
rz(1.1078542) q[1];
rz(-2.6158995) q[3];
sx q[3];
rz(-1.634053) q[3];
sx q[3];
rz(-0.69378187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6212375) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(1.3367782) q[2];
rz(-0.76198602) q[3];
sx q[3];
rz(-2.8218994) q[3];
sx q[3];
rz(-2.3412162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8005463) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(0.57089943) q[0];
rz(1.7123429) q[1];
sx q[1];
rz(-2.0776904) q[1];
sx q[1];
rz(-0.16194078) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2019661) q[0];
sx q[0];
rz(-1.8585748) q[0];
sx q[0];
rz(-0.94596011) q[0];
rz(-pi) q[1];
rz(-1.3822046) q[2];
sx q[2];
rz(-1.57975) q[2];
sx q[2];
rz(-2.9262275) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83878126) q[1];
sx q[1];
rz(-1.4421842) q[1];
sx q[1];
rz(1.3978811) q[1];
rz(-3.083769) q[3];
sx q[3];
rz(-2.6177546) q[3];
sx q[3];
rz(-0.1334838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1841715) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(1.4282248) q[2];
rz(1.1994294) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(1.7709581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7006871) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(-1.3700925) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(-3.0225171) q[2];
sx q[2];
rz(-1.7172114) q[2];
sx q[2];
rz(-1.7211771) q[2];
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