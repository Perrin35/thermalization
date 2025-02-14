OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27788568) q[0];
sx q[0];
rz(-2.7122893) q[0];
sx q[0];
rz(2.8306146) q[0];
rz(-1.3485981) q[1];
sx q[1];
rz(-2.0452979) q[1];
sx q[1];
rz(-0.57408339) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.047217) q[0];
sx q[0];
rz(-1.2021016) q[0];
sx q[0];
rz(-1.1218907) q[0];
rz(-0.98784222) q[2];
sx q[2];
rz(-1.541809) q[2];
sx q[2];
rz(-2.0984954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52770268) q[1];
sx q[1];
rz(-0.37606323) q[1];
sx q[1];
rz(-1.8164053) q[1];
rz(3.0776377) q[3];
sx q[3];
rz(-1.9465145) q[3];
sx q[3];
rz(-2.2967754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.2744039) q[2];
sx q[2];
rz(-1.2588661) q[2];
sx q[2];
rz(1.3192419) q[2];
rz(3.0684209) q[3];
sx q[3];
rz(-1.3061482) q[3];
sx q[3];
rz(-2.3225972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0935593) q[0];
sx q[0];
rz(-1.1807384) q[0];
sx q[0];
rz(-2.4617885) q[0];
rz(-2.3381084) q[1];
sx q[1];
rz(-1.7796703) q[1];
sx q[1];
rz(2.5018073) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073994324) q[0];
sx q[0];
rz(-0.72866458) q[0];
sx q[0];
rz(-1.0858028) q[0];
rz(-pi) q[1];
rz(1.670094) q[2];
sx q[2];
rz(-0.84337762) q[2];
sx q[2];
rz(2.119144) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8315994) q[1];
sx q[1];
rz(-0.91504708) q[1];
sx q[1];
rz(0.80226957) q[1];
x q[2];
rz(1.3469668) q[3];
sx q[3];
rz(-1.2696243) q[3];
sx q[3];
rz(-0.45066842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.60376072) q[2];
sx q[2];
rz(-0.49850264) q[2];
sx q[2];
rz(0.70408386) q[2];
rz(-3.0294561) q[3];
sx q[3];
rz(-1.4487368) q[3];
sx q[3];
rz(1.0224379) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10244399) q[0];
sx q[0];
rz(-2.0344489) q[0];
sx q[0];
rz(2.802134) q[0];
rz(-2.0231694) q[1];
sx q[1];
rz(-2.2148841) q[1];
sx q[1];
rz(0.035645398) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46860002) q[0];
sx q[0];
rz(-0.22449271) q[0];
sx q[0];
rz(0.39546449) q[0];
x q[1];
rz(2.8025134) q[2];
sx q[2];
rz(-1.2847752) q[2];
sx q[2];
rz(-1.7947444) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6468387) q[1];
sx q[1];
rz(-0.6907845) q[1];
sx q[1];
rz(0.34091807) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6142637) q[3];
sx q[3];
rz(-1.302056) q[3];
sx q[3];
rz(0.72878557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3493335) q[2];
sx q[2];
rz(-2.4534241) q[2];
sx q[2];
rz(2.2815857) q[2];
rz(-2.3490014) q[3];
sx q[3];
rz(-0.20321295) q[3];
sx q[3];
rz(2.121076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.129313) q[0];
sx q[0];
rz(-1.7565933) q[0];
sx q[0];
rz(2.5087575) q[0];
rz(1.9889779) q[1];
sx q[1];
rz(-0.8747789) q[1];
sx q[1];
rz(-1.4702183) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2180207) q[0];
sx q[0];
rz(-1.583719) q[0];
sx q[0];
rz(-0.92148975) q[0];
rz(-pi) q[1];
rz(-1.7196413) q[2];
sx q[2];
rz(-2.334377) q[2];
sx q[2];
rz(2.5850353) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6286575) q[1];
sx q[1];
rz(-2.3875065) q[1];
sx q[1];
rz(-1.9505461) q[1];
x q[2];
rz(2.4864462) q[3];
sx q[3];
rz(-1.0628502) q[3];
sx q[3];
rz(-0.20482132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70964885) q[2];
sx q[2];
rz(-0.96082965) q[2];
sx q[2];
rz(0.43080899) q[2];
rz(2.4591947) q[3];
sx q[3];
rz(-1.6513446) q[3];
sx q[3];
rz(0.94990134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34669852) q[0];
sx q[0];
rz(-1.8759202) q[0];
sx q[0];
rz(1.9516113) q[0];
rz(2.8311912) q[1];
sx q[1];
rz(-2.0605395) q[1];
sx q[1];
rz(0.50500542) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22689817) q[0];
sx q[0];
rz(-2.3932308) q[0];
sx q[0];
rz(-0.90004653) q[0];
x q[1];
rz(-1.2799576) q[2];
sx q[2];
rz(-2.7105936) q[2];
sx q[2];
rz(-2.9484826) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0481927) q[1];
sx q[1];
rz(-1.6580771) q[1];
sx q[1];
rz(-2.643226) q[1];
x q[2];
rz(1.6564063) q[3];
sx q[3];
rz(-2.2946649) q[3];
sx q[3];
rz(1.5392661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.40671047) q[2];
sx q[2];
rz(-0.70047417) q[2];
sx q[2];
rz(0.90275466) q[2];
rz(2.086575) q[3];
sx q[3];
rz(-2.2746634) q[3];
sx q[3];
rz(0.46190754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2257253) q[0];
sx q[0];
rz(-2.4287455) q[0];
sx q[0];
rz(1.073904) q[0];
rz(-1.3794927) q[1];
sx q[1];
rz(-0.72934377) q[1];
sx q[1];
rz(0.097537907) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3341971) q[0];
sx q[0];
rz(-1.0369891) q[0];
sx q[0];
rz(-1.3601204) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9729207) q[2];
sx q[2];
rz(-1.9430117) q[2];
sx q[2];
rz(1.8255359) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8733858) q[1];
sx q[1];
rz(-1.5097575) q[1];
sx q[1];
rz(-0.048442099) q[1];
rz(-pi) q[2];
rz(-0.82322466) q[3];
sx q[3];
rz(-0.71079094) q[3];
sx q[3];
rz(-1.2430134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.11233106) q[2];
sx q[2];
rz(-1.4533726) q[2];
sx q[2];
rz(-0.41109273) q[2];
rz(-1.4136275) q[3];
sx q[3];
rz(-2.3634383) q[3];
sx q[3];
rz(1.5467862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3635062) q[0];
sx q[0];
rz(-0.030487617) q[0];
sx q[0];
rz(1.7331069) q[0];
rz(2.1259437) q[1];
sx q[1];
rz(-1.3280832) q[1];
sx q[1];
rz(-1.4782864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41970872) q[0];
sx q[0];
rz(-2.1745178) q[0];
sx q[0];
rz(1.2092071) q[0];
x q[1];
rz(-2.1514966) q[2];
sx q[2];
rz(-1.0063356) q[2];
sx q[2];
rz(1.9866458) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.46131733) q[1];
sx q[1];
rz(-1.6467386) q[1];
sx q[1];
rz(-2.0148177) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8084548) q[3];
sx q[3];
rz(-2.6293652) q[3];
sx q[3];
rz(2.81306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0279072) q[2];
sx q[2];
rz(-2.5137641) q[2];
sx q[2];
rz(-2.9713463) q[2];
rz(1.8611192) q[3];
sx q[3];
rz(-1.4743285) q[3];
sx q[3];
rz(1.1580275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4729507) q[0];
sx q[0];
rz(-2.1773715) q[0];
sx q[0];
rz(0.52841312) q[0];
rz(-0.93327418) q[1];
sx q[1];
rz(-0.92527881) q[1];
sx q[1];
rz(0.55567137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66131573) q[0];
sx q[0];
rz(-2.4877036) q[0];
sx q[0];
rz(1.8956899) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3494928) q[2];
sx q[2];
rz(-1.6413771) q[2];
sx q[2];
rz(-2.5472699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5720776) q[1];
sx q[1];
rz(-1.5387427) q[1];
sx q[1];
rz(-0.7116913) q[1];
rz(1.2597841) q[3];
sx q[3];
rz(-2.2385983) q[3];
sx q[3];
rz(0.0086812191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3282503) q[2];
sx q[2];
rz(-2.3395061) q[2];
sx q[2];
rz(-2.8343406) q[2];
rz(0.79636374) q[3];
sx q[3];
rz(-2.4107404) q[3];
sx q[3];
rz(0.097631924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5873213) q[0];
sx q[0];
rz(-1.048943) q[0];
sx q[0];
rz(1.8121423) q[0];
rz(-2.9216271) q[1];
sx q[1];
rz(-2.797762) q[1];
sx q[1];
rz(1.8230009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7466965) q[0];
sx q[0];
rz(-1.3550955) q[0];
sx q[0];
rz(-1.8442276) q[0];
x q[1];
rz(2.9302915) q[2];
sx q[2];
rz(-2.6323279) q[2];
sx q[2];
rz(-0.76274055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4107549) q[1];
sx q[1];
rz(-1.9111834) q[1];
sx q[1];
rz(-1.0761258) q[1];
x q[2];
rz(-2.6133282) q[3];
sx q[3];
rz(-0.79341989) q[3];
sx q[3];
rz(1.8509282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69507504) q[2];
sx q[2];
rz(-1.4277642) q[2];
sx q[2];
rz(-1.2850777) q[2];
rz(0.88589969) q[3];
sx q[3];
rz(-1.748184) q[3];
sx q[3];
rz(-1.6133962) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1000243) q[0];
sx q[0];
rz(-2.4055241) q[0];
sx q[0];
rz(-2.8247483) q[0];
rz(-0.44081229) q[1];
sx q[1];
rz(-0.79439729) q[1];
sx q[1];
rz(-2.7712834) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.646831) q[0];
sx q[0];
rz(-1.6017822) q[0];
sx q[0];
rz(1.9520743) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4914527) q[2];
sx q[2];
rz(-2.1863089) q[2];
sx q[2];
rz(0.56734771) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84853828) q[1];
sx q[1];
rz(-1.8369743) q[1];
sx q[1];
rz(-1.8817503) q[1];
rz(-pi) q[2];
rz(-1.6746117) q[3];
sx q[3];
rz(-1.5201638) q[3];
sx q[3];
rz(-1.4736444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73654282) q[2];
sx q[2];
rz(-0.49095792) q[2];
sx q[2];
rz(1.6323818) q[2];
rz(-0.80471188) q[3];
sx q[3];
rz(-1.5038306) q[3];
sx q[3];
rz(-3.0338083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1152773) q[0];
sx q[0];
rz(-1.4125217) q[0];
sx q[0];
rz(-1.9052292) q[0];
rz(-0.23030494) q[1];
sx q[1];
rz(-2.8254012) q[1];
sx q[1];
rz(0.91513035) q[1];
rz(2.1565746) q[2];
sx q[2];
rz(-2.14902) q[2];
sx q[2];
rz(2.4075748) q[2];
rz(-2.0143916) q[3];
sx q[3];
rz(-0.66772912) q[3];
sx q[3];
rz(-2.276788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
