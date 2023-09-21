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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9288899) q[0];
sx q[0];
rz(-1.2455997) q[0];
sx q[0];
rz(-0.70744275) q[0];
rz(-pi) q[1];
rz(-1.0385752) q[2];
sx q[2];
rz(-0.79471171) q[2];
sx q[2];
rz(2.6796535) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5109374) q[1];
sx q[1];
rz(-1.8488548) q[1];
sx q[1];
rz(-3.0250938) q[1];
x q[2];
rz(0.9355448) q[3];
sx q[3];
rz(-1.2860635) q[3];
sx q[3];
rz(0.92393827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0779695) q[2];
sx q[2];
rz(-2.412553) q[2];
sx q[2];
rz(1.8135653) q[2];
rz(2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(0.13197556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48224738) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(-2.2609718) q[0];
rz(-1.847514) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(0.81726384) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0016804455) q[0];
sx q[0];
rz(-1.4298555) q[0];
sx q[0];
rz(1.2225479) q[0];
x q[1];
rz(1.0802644) q[2];
sx q[2];
rz(-1.3776407) q[2];
sx q[2];
rz(-2.1527388) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4546928) q[1];
sx q[1];
rz(-0.94140879) q[1];
sx q[1];
rz(1.252974) q[1];
rz(1.2137787) q[3];
sx q[3];
rz(-0.5001874) q[3];
sx q[3];
rz(0.57750765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91784224) q[2];
sx q[2];
rz(-0.68085256) q[2];
sx q[2];
rz(2.7775653) q[2];
rz(-0.98637995) q[3];
sx q[3];
rz(-1.4168408) q[3];
sx q[3];
rz(1.4646437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36689511) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(2.1752775) q[0];
rz(0.19293383) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(1.4470709) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23611785) q[0];
sx q[0];
rz(-2.9258699) q[0];
sx q[0];
rz(0.082213684) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3005199) q[2];
sx q[2];
rz(-2.8405361) q[2];
sx q[2];
rz(1.9384055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88168624) q[1];
sx q[1];
rz(-0.73598624) q[1];
sx q[1];
rz(-2.5658539) q[1];
x q[2];
rz(-1.1544187) q[3];
sx q[3];
rz(-1.3458816) q[3];
sx q[3];
rz(2.625536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7029999) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(2.036371) q[2];
rz(-2.3953719) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(-0.97572774) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8303191) q[0];
sx q[0];
rz(-1.8000326) q[0];
sx q[0];
rz(0.29005187) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0573425) q[2];
sx q[2];
rz(-1.3975189) q[2];
sx q[2];
rz(-1.372352) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3868689) q[1];
sx q[1];
rz(-1.9935973) q[1];
sx q[1];
rz(-2.5368284) q[1];
x q[2];
rz(-2.3936478) q[3];
sx q[3];
rz(-1.7224632) q[3];
sx q[3];
rz(-2.9880854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(-0.45670613) q[2];
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
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.221955) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(0.36002457) q[0];
rz(2.4941764) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(2.6470851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3213145) q[0];
sx q[0];
rz(-1.2914133) q[0];
sx q[0];
rz(1.4415635) q[0];
x q[1];
rz(-1.4353384) q[2];
sx q[2];
rz(-1.5317305) q[2];
sx q[2];
rz(2.1087697) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2795909) q[1];
sx q[1];
rz(-1.5108567) q[1];
sx q[1];
rz(-0.63043352) q[1];
rz(-2.6974929) q[3];
sx q[3];
rz(-1.5478304) q[3];
sx q[3];
rz(-0.75957739) q[3];
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
rz(-2.8295637) q[3];
sx q[3];
rz(-1.3307064) q[3];
sx q[3];
rz(-0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.9962149) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(2.9456855) q[0];
rz(-0.021082489) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(-1.9063937) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93341953) q[0];
sx q[0];
rz(-2.6801077) q[0];
sx q[0];
rz(0.74605201) q[0];
rz(-pi) q[1];
rz(1.9094798) q[2];
sx q[2];
rz(-0.89502305) q[2];
sx q[2];
rz(2.0896926) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27855733) q[1];
sx q[1];
rz(-2.1457991) q[1];
sx q[1];
rz(-2.0726191) q[1];
x q[2];
rz(-2.8060693) q[3];
sx q[3];
rz(-2.8918859) q[3];
sx q[3];
rz(-1.6186796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49332508) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(-0.43169272) q[2];
rz(-1.4124983) q[3];
sx q[3];
rz(-2.4192211) q[3];
sx q[3];
rz(-0.13599642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39111185) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(3.0294763) q[0];
rz(-0.21513367) q[1];
sx q[1];
rz(-1.5810177) q[1];
sx q[1];
rz(1.1134061) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8047785) q[0];
sx q[0];
rz(-1.830173) q[0];
sx q[0];
rz(-0.39774261) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31539519) q[2];
sx q[2];
rz(-1.3239667) q[2];
sx q[2];
rz(0.95755267) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73247319) q[1];
sx q[1];
rz(-1.2670244) q[1];
sx q[1];
rz(2.4227546) q[1];
rz(2.0180842) q[3];
sx q[3];
rz(-1.3778068) q[3];
sx q[3];
rz(-2.6393294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1404184) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(3.0155638) q[2];
rz(-2.0942988) q[3];
sx q[3];
rz(-1.8208561) q[3];
sx q[3];
rz(-2.4333911) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66013181) q[0];
sx q[0];
rz(-2.3829057) q[0];
sx q[0];
rz(-1.6814167) q[0];
rz(1.2449645) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(-1.9326899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9019933) q[0];
sx q[0];
rz(-1.8327466) q[0];
sx q[0];
rz(3.0595747) q[0];
rz(-pi) q[1];
rz(-0.86645855) q[2];
sx q[2];
rz(-2.2668215) q[2];
sx q[2];
rz(2.8528086) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9637451) q[1];
sx q[1];
rz(-2.3235465) q[1];
sx q[1];
rz(0.45958105) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9221411) q[3];
sx q[3];
rz(-2.1900574) q[3];
sx q[3];
rz(2.0451562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7631491) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(1.3195999) q[2];
rz(2.54946) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(-0.035141703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9086583) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(-1.3611025) q[0];
rz(-1.9305485) q[1];
sx q[1];
rz(-2.2369604) q[1];
sx q[1];
rz(0.39168721) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2101333) q[0];
sx q[0];
rz(-1.1362695) q[0];
sx q[0];
rz(1.5772485) q[0];
x q[1];
rz(2.7462256) q[2];
sx q[2];
rz(-2.6569416) q[2];
sx q[2];
rz(-1.6389099) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0605676) q[1];
sx q[1];
rz(-2.276366) q[1];
sx q[1];
rz(-0.4391567) q[1];
rz(-pi) q[2];
rz(0.12556062) q[3];
sx q[3];
rz(-2.6124622) q[3];
sx q[3];
rz(-2.373113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52035511) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(-1.8048145) q[2];
rz(-0.76198602) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(2.3412162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(15/(14*pi)) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(-0.57089943) q[0];
rz(-1.7123429) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(2.9796519) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3085092) q[0];
sx q[0];
rz(-2.1662795) q[0];
sx q[0];
rz(-2.7916629) q[0];
rz(-pi) q[1];
rz(1.5230721) q[2];
sx q[2];
rz(-2.952791) q[2];
sx q[2];
rz(1.4023086) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.83878126) q[1];
sx q[1];
rz(-1.4421842) q[1];
sx q[1];
rz(1.7437115) q[1];
rz(-1.5374244) q[3];
sx q[3];
rz(-2.0936692) q[3];
sx q[3];
rz(0.066730412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95742115) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(-1.7133678) q[2];
rz(1.9421633) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(1.7709581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7006871) q[0];
sx q[0];
rz(-2.6661243) q[0];
sx q[0];
rz(2.0211924) q[0];
rz(-1.3700925) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(0.8926819) q[2];
sx q[2];
rz(-2.953139) q[2];
sx q[2];
rz(2.1072731) q[2];
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