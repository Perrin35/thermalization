OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8925722) q[0];
sx q[0];
rz(-2.1828716) q[0];
sx q[0];
rz(0.22554654) q[0];
rz(-1.072999) q[1];
sx q[1];
rz(3.5993242) q[1];
sx q[1];
rz(9.2158894) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6207558) q[0];
sx q[0];
rz(-1.136047) q[0];
sx q[0];
rz(1.8692383) q[0];
rz(0.3772119) q[2];
sx q[2];
rz(-2.4754731) q[2];
sx q[2];
rz(-1.3782901) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.75336136) q[1];
sx q[1];
rz(-1.8750138) q[1];
sx q[1];
rz(2.6016584) q[1];
x q[2];
rz(2.8257583) q[3];
sx q[3];
rz(-1.2441934) q[3];
sx q[3];
rz(2.0706319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0245584) q[2];
sx q[2];
rz(-0.25091761) q[2];
sx q[2];
rz(-0.92570242) q[2];
rz(0.77445817) q[3];
sx q[3];
rz(-1.9175994) q[3];
sx q[3];
rz(1.2831877) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16933146) q[0];
sx q[0];
rz(-0.65003482) q[0];
sx q[0];
rz(1.6768804) q[0];
rz(-0.88893923) q[1];
sx q[1];
rz(-0.89193946) q[1];
sx q[1];
rz(-1.3533786) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1036677) q[0];
sx q[0];
rz(-2.1620181) q[0];
sx q[0];
rz(1.9210438) q[0];
rz(2.775506) q[2];
sx q[2];
rz(-1.5502068) q[2];
sx q[2];
rz(2.1799157) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1802952) q[1];
sx q[1];
rz(-0.63034454) q[1];
sx q[1];
rz(-0.35758361) q[1];
rz(-pi) q[2];
rz(0.7781182) q[3];
sx q[3];
rz(-1.3740842) q[3];
sx q[3];
rz(2.9784378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9951524) q[2];
sx q[2];
rz(-1.012994) q[2];
sx q[2];
rz(-0.97274441) q[2];
rz(-1.3257596) q[3];
sx q[3];
rz(-1.9917859) q[3];
sx q[3];
rz(-0.61746922) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72508088) q[0];
sx q[0];
rz(-1.4840115) q[0];
sx q[0];
rz(-0.02956477) q[0];
rz(-1.0211396) q[1];
sx q[1];
rz(-0.28426668) q[1];
sx q[1];
rz(-2.8254642) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3409948) q[0];
sx q[0];
rz(-1.5506239) q[0];
sx q[0];
rz(1.5490269) q[0];
rz(-pi) q[1];
rz(-3.0477858) q[2];
sx q[2];
rz(-1.7711763) q[2];
sx q[2];
rz(2.0721638) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7252032) q[1];
sx q[1];
rz(-1.3969743) q[1];
sx q[1];
rz(1.9085788) q[1];
rz(1.49316) q[3];
sx q[3];
rz(-2.2901553) q[3];
sx q[3];
rz(0.84859728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0315087) q[2];
sx q[2];
rz(-1.8434593) q[2];
sx q[2];
rz(-2.617344) q[2];
rz(-2.5414069) q[3];
sx q[3];
rz(-0.46349183) q[3];
sx q[3];
rz(-2.0704827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4635187) q[0];
sx q[0];
rz(-1.2200032) q[0];
sx q[0];
rz(2.779261) q[0];
rz(-0.70829567) q[1];
sx q[1];
rz(-0.34168044) q[1];
sx q[1];
rz(1.1452311) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0247097) q[0];
sx q[0];
rz(-1.6779461) q[0];
sx q[0];
rz(0.080579455) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0037654) q[2];
sx q[2];
rz(-0.7793372) q[2];
sx q[2];
rz(-0.61160751) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.64001361) q[1];
sx q[1];
rz(-2.3928436) q[1];
sx q[1];
rz(-1.1925936) q[1];
rz(-1.5208552) q[3];
sx q[3];
rz(-1.4068713) q[3];
sx q[3];
rz(-2.0633179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6142673) q[2];
sx q[2];
rz(-2.3136487) q[2];
sx q[2];
rz(-3.1329727) q[2];
rz(2.4873867) q[3];
sx q[3];
rz(-0.71627408) q[3];
sx q[3];
rz(0.27230984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24097405) q[0];
sx q[0];
rz(-2.8130377) q[0];
sx q[0];
rz(-2.3531438) q[0];
rz(-2.12766) q[1];
sx q[1];
rz(-1.8221816) q[1];
sx q[1];
rz(3.0533155) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4409613) q[0];
sx q[0];
rz(-0.92871237) q[0];
sx q[0];
rz(-2.2034646) q[0];
rz(-2.5779326) q[2];
sx q[2];
rz(-1.6573326) q[2];
sx q[2];
rz(1.5147839) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8535117) q[1];
sx q[1];
rz(-1.8520903) q[1];
sx q[1];
rz(-1.049925) q[1];
x q[2];
rz(2.0658947) q[3];
sx q[3];
rz(-1.5468239) q[3];
sx q[3];
rz(-2.3119462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91904116) q[2];
sx q[2];
rz(-1.7307948) q[2];
sx q[2];
rz(-2.6338573) q[2];
rz(-0.12568411) q[3];
sx q[3];
rz(-1.1832184) q[3];
sx q[3];
rz(-2.3465033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3313726) q[0];
sx q[0];
rz(-2.3851244) q[0];
sx q[0];
rz(0.08547011) q[0];
rz(-2.5147009) q[1];
sx q[1];
rz(-1.9848928) q[1];
sx q[1];
rz(-1.289182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4816718) q[0];
sx q[0];
rz(-0.22511765) q[0];
sx q[0];
rz(1.0876924) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7250546) q[2];
sx q[2];
rz(-2.3084894) q[2];
sx q[2];
rz(-1.8869893) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5252555) q[1];
sx q[1];
rz(-0.29381815) q[1];
sx q[1];
rz(2.1928685) q[1];
x q[2];
rz(-0.014628476) q[3];
sx q[3];
rz(-1.3335776) q[3];
sx q[3];
rz(2.5166409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.36279303) q[2];
sx q[2];
rz(-2.1976566) q[2];
sx q[2];
rz(-1.8865406) q[2];
rz(3.0826027) q[3];
sx q[3];
rz(-1.9374282) q[3];
sx q[3];
rz(-2.1571933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4208218) q[0];
sx q[0];
rz(-1.2728007) q[0];
sx q[0];
rz(0.76501784) q[0];
rz(1.7401241) q[1];
sx q[1];
rz(-2.2547289) q[1];
sx q[1];
rz(2.9041362) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7516605) q[0];
sx q[0];
rz(-2.513305) q[0];
sx q[0];
rz(1.6310255) q[0];
x q[1];
rz(1.6638905) q[2];
sx q[2];
rz(-2.4590465) q[2];
sx q[2];
rz(1.0034221) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2903257) q[1];
sx q[1];
rz(-1.8614166) q[1];
sx q[1];
rz(1.5253011) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1116099) q[3];
sx q[3];
rz(-0.87453547) q[3];
sx q[3];
rz(-3.1122623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6852297) q[2];
sx q[2];
rz(-2.3380029) q[2];
sx q[2];
rz(2.2881499) q[2];
rz(-1.8026836) q[3];
sx q[3];
rz(-1.5693376) q[3];
sx q[3];
rz(0.16551031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6859739) q[0];
sx q[0];
rz(-1.9299262) q[0];
sx q[0];
rz(-2.9562505) q[0];
rz(2.7475884) q[1];
sx q[1];
rz(-1.7454255) q[1];
sx q[1];
rz(2.0448763) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029848969) q[0];
sx q[0];
rz(-1.9557448) q[0];
sx q[0];
rz(-2.0864331) q[0];
rz(0.11285891) q[2];
sx q[2];
rz(-0.90870171) q[2];
sx q[2];
rz(-0.87967949) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.19272596) q[1];
sx q[1];
rz(-2.1565686) q[1];
sx q[1];
rz(-1.1652258) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6043209) q[3];
sx q[3];
rz(-2.4747765) q[3];
sx q[3];
rz(1.3992573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77157053) q[2];
sx q[2];
rz(-3.0574419) q[2];
sx q[2];
rz(0.30878511) q[2];
rz(1.8874946) q[3];
sx q[3];
rz(-1.2718688) q[3];
sx q[3];
rz(-2.0103644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9517188) q[0];
sx q[0];
rz(-1.8872486) q[0];
sx q[0];
rz(0.21981123) q[0];
rz(-1.9728569) q[1];
sx q[1];
rz(-1.7857779) q[1];
sx q[1];
rz(-1.8692325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6913805) q[0];
sx q[0];
rz(-0.60465971) q[0];
sx q[0];
rz(1.613841) q[0];
x q[1];
rz(0.51080334) q[2];
sx q[2];
rz(-1.496721) q[2];
sx q[2];
rz(-1.6965716) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5516806) q[1];
sx q[1];
rz(-1.9724692) q[1];
sx q[1];
rz(-1.3856616) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79325647) q[3];
sx q[3];
rz(-2.4783387) q[3];
sx q[3];
rz(2.548717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9610338) q[2];
sx q[2];
rz(-0.84332931) q[2];
sx q[2];
rz(-2.6832704) q[2];
rz(-0.35442963) q[3];
sx q[3];
rz(-1.7731881) q[3];
sx q[3];
rz(-1.1752769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.71427041) q[0];
sx q[0];
rz(-1.9388119) q[0];
sx q[0];
rz(-0.74139968) q[0];
rz(1.4330014) q[1];
sx q[1];
rz(-2.7129136) q[1];
sx q[1];
rz(3.0523849) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2994381) q[0];
sx q[0];
rz(-2.0426867) q[0];
sx q[0];
rz(0.20345511) q[0];
x q[1];
rz(2.600895) q[2];
sx q[2];
rz(-0.6650228) q[2];
sx q[2];
rz(-2.2407209) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20232378) q[1];
sx q[1];
rz(-1.2319733) q[1];
sx q[1];
rz(-3.0975049) q[1];
rz(-pi) q[2];
rz(-0.44762917) q[3];
sx q[3];
rz(-2.3530686) q[3];
sx q[3];
rz(0.41390362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1409113) q[2];
sx q[2];
rz(-3.0912283) q[2];
sx q[2];
rz(-2.4939406) q[2];
rz(-0.40376136) q[3];
sx q[3];
rz(-1.253456) q[3];
sx q[3];
rz(-0.13127479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972926) q[0];
sx q[0];
rz(-2.1680752) q[0];
sx q[0];
rz(-2.0808921) q[0];
rz(0.79509673) q[1];
sx q[1];
rz(-0.92082321) q[1];
sx q[1];
rz(-0.77650741) q[1];
rz(-2.2262103) q[2];
sx q[2];
rz(-1.6554828) q[2];
sx q[2];
rz(0.2632904) q[2];
rz(0.078362502) q[3];
sx q[3];
rz(-1.6607124) q[3];
sx q[3];
rz(2.511047) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
