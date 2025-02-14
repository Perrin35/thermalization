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
rz(-0.072862536) q[0];
sx q[0];
rz(-2.1003567) q[0];
sx q[0];
rz(0.6676724) q[0];
rz(-0.11198894) q[1];
sx q[1];
rz(-1.656124) q[1];
sx q[1];
rz(0.30593425) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4914305) q[0];
sx q[0];
rz(-1.9823649) q[0];
sx q[0];
rz(1.4812019) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61688136) q[2];
sx q[2];
rz(-1.1921669) q[2];
sx q[2];
rz(-2.4509492) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.48253912) q[1];
sx q[1];
rz(-1.6371563) q[1];
sx q[1];
rz(0.37239667) q[1];
x q[2];
rz(3.0612321) q[3];
sx q[3];
rz(-0.70801641) q[3];
sx q[3];
rz(-2.3637259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2970807) q[2];
sx q[2];
rz(-1.7252555) q[2];
sx q[2];
rz(-0.97999209) q[2];
rz(2.8196107) q[3];
sx q[3];
rz(-1.201509) q[3];
sx q[3];
rz(0.15596786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62671536) q[0];
sx q[0];
rz(-1.5276271) q[0];
sx q[0];
rz(-0.8901341) q[0];
rz(2.0251677) q[1];
sx q[1];
rz(-2.2869488) q[1];
sx q[1];
rz(1.9568303) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6044004) q[0];
sx q[0];
rz(-1.0774776) q[0];
sx q[0];
rz(1.7866578) q[0];
rz(-pi) q[1];
rz(-2.7689887) q[2];
sx q[2];
rz(-1.0516343) q[2];
sx q[2];
rz(-2.2302996) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0240619) q[1];
sx q[1];
rz(-1.6819997) q[1];
sx q[1];
rz(0.11604138) q[1];
rz(-1.5873648) q[3];
sx q[3];
rz(-1.9958766) q[3];
sx q[3];
rz(-2.1147902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4778135) q[2];
sx q[2];
rz(-0.16872831) q[2];
sx q[2];
rz(-1.4230049) q[2];
rz(-2.6723828) q[3];
sx q[3];
rz(-1.6459295) q[3];
sx q[3];
rz(-2.7928228) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94217268) q[0];
sx q[0];
rz(-1.9439789) q[0];
sx q[0];
rz(-1.0823826) q[0];
rz(-0.17467817) q[1];
sx q[1];
rz(-1.7593316) q[1];
sx q[1];
rz(-0.47294012) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.32047) q[0];
sx q[0];
rz(-1.1468915) q[0];
sx q[0];
rz(-0.17571077) q[0];
rz(-1.9835112) q[2];
sx q[2];
rz(-0.7447401) q[2];
sx q[2];
rz(-1.2117529) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.14607695) q[1];
sx q[1];
rz(-0.56997609) q[1];
sx q[1];
rz(1.5289299) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93227902) q[3];
sx q[3];
rz(-2.4014056) q[3];
sx q[3];
rz(0.96961752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.48367721) q[2];
sx q[2];
rz(-1.4853073) q[2];
sx q[2];
rz(-1.8243054) q[2];
rz(-0.45051908) q[3];
sx q[3];
rz(-2.2299288) q[3];
sx q[3];
rz(1.3594422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6555742) q[0];
sx q[0];
rz(-2.0756606) q[0];
sx q[0];
rz(2.4701212) q[0];
rz(2.2010522) q[1];
sx q[1];
rz(-0.42010072) q[1];
sx q[1];
rz(-2.1817575) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.335054) q[0];
sx q[0];
rz(-1.949652) q[0];
sx q[0];
rz(-0.23083487) q[0];
rz(-pi) q[1];
rz(-1.3449538) q[2];
sx q[2];
rz(-0.91691559) q[2];
sx q[2];
rz(-2.8381062) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.48573182) q[1];
sx q[1];
rz(-1.8645687) q[1];
sx q[1];
rz(1.345977) q[1];
rz(0.22573443) q[3];
sx q[3];
rz(-2.6896296) q[3];
sx q[3];
rz(2.0737157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.81022108) q[2];
sx q[2];
rz(-2.0600617) q[2];
sx q[2];
rz(-0.01037154) q[2];
rz(1.2822019) q[3];
sx q[3];
rz(-2.5344262) q[3];
sx q[3];
rz(0.31266323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4717167) q[0];
sx q[0];
rz(-0.80729055) q[0];
sx q[0];
rz(-0.5759936) q[0];
rz(-2.5895789) q[1];
sx q[1];
rz(-1.8355398) q[1];
sx q[1];
rz(0.920151) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3475889) q[0];
sx q[0];
rz(-2.2146642) q[0];
sx q[0];
rz(-2.3809303) q[0];
rz(0.29187496) q[2];
sx q[2];
rz(-1.5426965) q[2];
sx q[2];
rz(2.4296938) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.53729) q[1];
sx q[1];
rz(-1.3370945) q[1];
sx q[1];
rz(-1.2140732) q[1];
rz(-pi) q[2];
rz(2.4534485) q[3];
sx q[3];
rz(-0.16040186) q[3];
sx q[3];
rz(-0.80996603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0612023) q[2];
sx q[2];
rz(-1.3866501) q[2];
sx q[2];
rz(-0.45026711) q[2];
rz(3.0007512) q[3];
sx q[3];
rz(-1.9167269) q[3];
sx q[3];
rz(-2.5950477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30215728) q[0];
sx q[0];
rz(-1.7385229) q[0];
sx q[0];
rz(2.9490525) q[0];
rz(2.4977066) q[1];
sx q[1];
rz(-1.5637249) q[1];
sx q[1];
rz(0.083316915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9238038) q[0];
sx q[0];
rz(-1.4155651) q[0];
sx q[0];
rz(-0.96331994) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4519789) q[2];
sx q[2];
rz(-1.6518704) q[2];
sx q[2];
rz(0.18866779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.5789448) q[1];
sx q[1];
rz(-1.399843) q[1];
sx q[1];
rz(-1.2557545) q[1];
x q[2];
rz(1.313435) q[3];
sx q[3];
rz(-1.93522) q[3];
sx q[3];
rz(2.6217883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4087499) q[2];
sx q[2];
rz(-2.7538444) q[2];
sx q[2];
rz(-0.55633083) q[2];
rz(-0.88931102) q[3];
sx q[3];
rz(-1.5143062) q[3];
sx q[3];
rz(-2.6591532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0675875) q[0];
sx q[0];
rz(-2.6370625) q[0];
sx q[0];
rz(2.0489847) q[0];
rz(-0.69674528) q[1];
sx q[1];
rz(-0.69843355) q[1];
sx q[1];
rz(-2.4901857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6987885) q[0];
sx q[0];
rz(-2.4518928) q[0];
sx q[0];
rz(-2.882363) q[0];
x q[1];
rz(-0.90187855) q[2];
sx q[2];
rz(-1.9142391) q[2];
sx q[2];
rz(2.4879188) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40245956) q[1];
sx q[1];
rz(-0.62109438) q[1];
sx q[1];
rz(-2.0517292) q[1];
x q[2];
rz(-0.015373165) q[3];
sx q[3];
rz(-1.5173994) q[3];
sx q[3];
rz(-1.3863409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6394028) q[2];
sx q[2];
rz(-0.71523372) q[2];
sx q[2];
rz(-2.5779842) q[2];
rz(0.11858502) q[3];
sx q[3];
rz(-1.010681) q[3];
sx q[3];
rz(-2.3732869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3350155) q[0];
sx q[0];
rz(-1.7105569) q[0];
sx q[0];
rz(0.12120506) q[0];
rz(-2.99446) q[1];
sx q[1];
rz(-1.9857261) q[1];
sx q[1];
rz(-2.3187231) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3980902) q[0];
sx q[0];
rz(-1.8202298) q[0];
sx q[0];
rz(1.477308) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6001094) q[2];
sx q[2];
rz(-2.1357059) q[2];
sx q[2];
rz(0.82290111) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.71119324) q[1];
sx q[1];
rz(-2.3533513) q[1];
sx q[1];
rz(1.3543717) q[1];
rz(-pi) q[2];
rz(-1.9328281) q[3];
sx q[3];
rz(-2.9144807) q[3];
sx q[3];
rz(2.1696573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11468981) q[2];
sx q[2];
rz(-2.5451886) q[2];
sx q[2];
rz(-0.58843311) q[2];
rz(2.9450997) q[3];
sx q[3];
rz(-1.7076219) q[3];
sx q[3];
rz(2.4832723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999417) q[0];
sx q[0];
rz(-0.64084941) q[0];
sx q[0];
rz(0.95630056) q[0];
rz(0.20005964) q[1];
sx q[1];
rz(-1.450489) q[1];
sx q[1];
rz(0.53093451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0414687) q[0];
sx q[0];
rz(-1.6371797) q[0];
sx q[0];
rz(-2.0602577) q[0];
x q[1];
rz(-0.63636698) q[2];
sx q[2];
rz(-1.8312935) q[2];
sx q[2];
rz(-0.72450984) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8968205) q[1];
sx q[1];
rz(-1.7339118) q[1];
sx q[1];
rz(-0.54508026) q[1];
rz(-pi) q[2];
rz(0.45147272) q[3];
sx q[3];
rz(-1.6992308) q[3];
sx q[3];
rz(2.4344545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7134573) q[2];
sx q[2];
rz(-1.0604481) q[2];
sx q[2];
rz(1.1327845) q[2];
rz(-2.1652341) q[3];
sx q[3];
rz(-0.68621695) q[3];
sx q[3];
rz(-1.8648225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79162947) q[0];
sx q[0];
rz(-3.0296453) q[0];
sx q[0];
rz(-1.5891225) q[0];
rz(-0.38996977) q[1];
sx q[1];
rz(-1.5733122) q[1];
sx q[1];
rz(-0.76036298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2137879) q[0];
sx q[0];
rz(-2.669691) q[0];
sx q[0];
rz(-0.15165059) q[0];
x q[1];
rz(-1.5993714) q[2];
sx q[2];
rz(-1.1090389) q[2];
sx q[2];
rz(-0.28959488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.37959209) q[1];
sx q[1];
rz(-1.6093644) q[1];
sx q[1];
rz(-1.9412575) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.032234) q[3];
sx q[3];
rz(-1.3038774) q[3];
sx q[3];
rz(2.0570786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4752263) q[2];
sx q[2];
rz(-2.1658289) q[2];
sx q[2];
rz(0.59398061) q[2];
rz(-2.5178759) q[3];
sx q[3];
rz(-1.2902322) q[3];
sx q[3];
rz(-0.82998961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8754616) q[0];
sx q[0];
rz(-0.13127357) q[0];
sx q[0];
rz(1.776478) q[0];
rz(0.020137067) q[1];
sx q[1];
rz(-0.29300856) q[1];
sx q[1];
rz(-1.551052) q[1];
rz(0.0082848194) q[2];
sx q[2];
rz(-1.1859672) q[2];
sx q[2];
rz(-2.7993508) q[2];
rz(0.094811335) q[3];
sx q[3];
rz(-2.9651793) q[3];
sx q[3];
rz(1.7549979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
