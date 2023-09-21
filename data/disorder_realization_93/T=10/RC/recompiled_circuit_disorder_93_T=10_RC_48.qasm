OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(1.8811037) q[0];
rz(2.1029544) q[1];
sx q[1];
rz(-1.3488052) q[1];
sx q[1];
rz(-2.2178712) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38102725) q[0];
sx q[0];
rz(-2.0187223) q[0];
sx q[0];
rz(2.7678124) q[0];
rz(-pi) q[1];
x q[1];
rz(2.997424) q[2];
sx q[2];
rz(-1.8509794) q[2];
sx q[2];
rz(-0.35137128) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1957789) q[1];
sx q[1];
rz(-1.073277) q[1];
sx q[1];
rz(2.0583908) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63763036) q[3];
sx q[3];
rz(-2.5507567) q[3];
sx q[3];
rz(-1.6319815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2279921) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(2.4285765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0682003) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(-1.9447928) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-0.49566832) q[1];
sx q[1];
rz(-1.4555567) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8287795) q[0];
sx q[0];
rz(-1.7377186) q[0];
sx q[0];
rz(-2.1524327) q[0];
x q[1];
rz(2.7669719) q[2];
sx q[2];
rz(-1.6472367) q[2];
sx q[2];
rz(2.3945216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5160408) q[1];
sx q[1];
rz(-1.5572773) q[1];
sx q[1];
rz(-1.0479755) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.060828408) q[3];
sx q[3];
rz(-0.41234499) q[3];
sx q[3];
rz(2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7130647) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(1.3519752) q[2];
rz(2.9591566) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(-0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882554) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(2.341111) q[0];
rz(3.1128186) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(-1.9690537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4370255) q[0];
sx q[0];
rz(-2.4798658) q[0];
sx q[0];
rz(-2.6252803) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2433979) q[2];
sx q[2];
rz(-0.84257579) q[2];
sx q[2];
rz(2.3936405) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6781569) q[1];
sx q[1];
rz(-1.0532182) q[1];
sx q[1];
rz(-0.28097681) q[1];
rz(-pi) q[2];
rz(1.3200687) q[3];
sx q[3];
rz(-1.6276976) q[3];
sx q[3];
rz(-2.1123321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0671493) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(2.2303936) q[2];
rz(2.1905812) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38055414) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(3.0134841) q[0];
rz(-3.065486) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(-0.52350837) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9273705) q[0];
sx q[0];
rz(-0.71338755) q[0];
sx q[0];
rz(0.58332304) q[0];
x q[1];
rz(0.33004327) q[2];
sx q[2];
rz(-1.4903307) q[2];
sx q[2];
rz(2.6835494) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7586786) q[1];
sx q[1];
rz(-0.74712336) q[1];
sx q[1];
rz(1.9488571) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.773049) q[3];
sx q[3];
rz(-0.60086717) q[3];
sx q[3];
rz(2.3893389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52544242) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(0.564044) q[2];
rz(2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600835) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(1.6500641) q[0];
rz(-0.87961698) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(-0.99194828) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0569699) q[0];
sx q[0];
rz(-1.7449433) q[0];
sx q[0];
rz(1.6790381) q[0];
rz(-pi) q[1];
rz(1.5993824) q[2];
sx q[2];
rz(-2.8386142) q[2];
sx q[2];
rz(-2.6945393) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.93047749) q[1];
sx q[1];
rz(-2.8301297) q[1];
sx q[1];
rz(-0.13179563) q[1];
rz(-0.81327849) q[3];
sx q[3];
rz(-1.6653898) q[3];
sx q[3];
rz(-2.7300342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0118959) q[2];
sx q[2];
rz(-0.36964881) q[2];
sx q[2];
rz(0.27080718) q[2];
rz(-2.9233542) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4145684) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(-1.3487934) q[0];
rz(-2.7596966) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(-1.4250925) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3544918) q[0];
sx q[0];
rz(-2.0658501) q[0];
sx q[0];
rz(-1.3160734) q[0];
x q[1];
rz(3.0503057) q[2];
sx q[2];
rz(-1.316615) q[2];
sx q[2];
rz(-2.9448178) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35364321) q[1];
sx q[1];
rz(-2.0197581) q[1];
sx q[1];
rz(3.06762) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85961996) q[3];
sx q[3];
rz(-2.4058127) q[3];
sx q[3];
rz(1.7914634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(0.18051906) q[3];
sx q[3];
rz(-1.6231977) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6329704) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(0.82180506) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9200631) q[0];
sx q[0];
rz(-2.1863345) q[0];
sx q[0];
rz(-0.91381844) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1173382) q[2];
sx q[2];
rz(-0.80596906) q[2];
sx q[2];
rz(-2.1794127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2591178) q[1];
sx q[1];
rz(-1.7977409) q[1];
sx q[1];
rz(0.74101733) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58737289) q[3];
sx q[3];
rz(-1.2578739) q[3];
sx q[3];
rz(1.8250993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8043148) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(0.24469963) q[2];
rz(3.0120567) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(2.3186671) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(-1.8364505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37874052) q[0];
sx q[0];
rz(-1.2438602) q[0];
sx q[0];
rz(-0.61600323) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70456409) q[2];
sx q[2];
rz(-1.8578055) q[2];
sx q[2];
rz(-1.9412083) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95272428) q[1];
sx q[1];
rz(-2.3844068) q[1];
sx q[1];
rz(0.56307478) q[1];
rz(-1.6373789) q[3];
sx q[3];
rz(-2.801602) q[3];
sx q[3];
rz(-0.41223994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(1.6513599) q[2];
rz(-2.0643318) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(3.0100477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7548783) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(2.819678) q[0];
rz(-1.6053258) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(2.4386491) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87899938) q[0];
sx q[0];
rz(-1.2125373) q[0];
sx q[0];
rz(-2.7260289) q[0];
rz(-pi) q[1];
rz(-1.0614971) q[2];
sx q[2];
rz(-1.5361538) q[2];
sx q[2];
rz(0.73355567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.34300464) q[1];
sx q[1];
rz(-1.7354021) q[1];
sx q[1];
rz(2.1144457) q[1];
rz(-2.8076257) q[3];
sx q[3];
rz(-2.1595862) q[3];
sx q[3];
rz(-2.2759618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(1.8019603) q[2];
rz(0.30570269) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(-1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16383485) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(-1.2257858) q[0];
rz(-2.2380791) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7077431) q[0];
sx q[0];
rz(-1.187547) q[0];
sx q[0];
rz(-1.7953403) q[0];
x q[1];
rz(2.2655728) q[2];
sx q[2];
rz(-0.48601905) q[2];
sx q[2];
rz(1.8234058) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.66099) q[1];
sx q[1];
rz(-0.27462474) q[1];
sx q[1];
rz(1.8250699) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7077984) q[3];
sx q[3];
rz(-1.9926096) q[3];
sx q[3];
rz(-1.0292366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.998385) q[2];
rz(-3.0269567) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(-2.519683) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(2.2767699) q[2];
sx q[2];
rz(-2.13158) q[2];
sx q[2];
rz(-2.0958015) q[2];
rz(-3.0631089) q[3];
sx q[3];
rz(-2.2208636) q[3];
sx q[3];
rz(0.29490864) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
