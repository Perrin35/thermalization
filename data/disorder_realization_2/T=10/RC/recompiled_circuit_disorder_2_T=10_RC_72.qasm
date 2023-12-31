OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6931273) q[0];
sx q[0];
rz(5.7603523) q[0];
sx q[0];
rz(8.8011959) q[0];
rz(0.29016718) q[1];
sx q[1];
rz(-2.4224412) q[1];
sx q[1];
rz(-0.50049385) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96007632) q[0];
sx q[0];
rz(-1.5746563) q[0];
sx q[0];
rz(0.5998248) q[0];
x q[1];
rz(3.0481553) q[2];
sx q[2];
rz(-1.4215901) q[2];
sx q[2];
rz(-2.0051533) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11195586) q[1];
sx q[1];
rz(-1.8957378) q[1];
sx q[1];
rz(-0.8128266) q[1];
x q[2];
rz(2.0088828) q[3];
sx q[3];
rz(-1.4287018) q[3];
sx q[3];
rz(2.4349468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3502675) q[2];
sx q[2];
rz(-1.2056377) q[2];
sx q[2];
rz(1.9187437) q[2];
rz(1.6932999) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(2.1712415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(0.37110776) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(1.0789385) q[0];
rz(1.7547912) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(-2.4761377) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6404214) q[0];
sx q[0];
rz(-1.2873532) q[0];
sx q[0];
rz(-0.47407504) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3999248) q[2];
sx q[2];
rz(-1.8881646) q[2];
sx q[2];
rz(-0.7764118) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.14109719) q[1];
sx q[1];
rz(-0.54447237) q[1];
sx q[1];
rz(0.46846868) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86136787) q[3];
sx q[3];
rz(-1.5966166) q[3];
sx q[3];
rz(2.4840419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60454303) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(1.6710619) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(0.69141928) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55494088) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(-2.3828322) q[0];
rz(1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(1.741515) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5132719) q[0];
sx q[0];
rz(-1.3043881) q[0];
sx q[0];
rz(-0.4011641) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7064352) q[2];
sx q[2];
rz(-1.7559768) q[2];
sx q[2];
rz(-0.26091012) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.53576614) q[1];
sx q[1];
rz(-1.8927791) q[1];
sx q[1];
rz(0.95226007) q[1];
rz(-pi) q[2];
rz(0.16617822) q[3];
sx q[3];
rz(-1.2697392) q[3];
sx q[3];
rz(0.64134669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65163461) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(0.65762323) q[2];
rz(1.970132) q[3];
sx q[3];
rz(-1.4843342) q[3];
sx q[3];
rz(1.4191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3477429) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(1.6963652) q[0];
rz(1.4472648) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(2.7935374) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30663438) q[0];
sx q[0];
rz(-1.277703) q[0];
sx q[0];
rz(2.4819863) q[0];
rz(-pi) q[1];
rz(-2.9910827) q[2];
sx q[2];
rz(-1.9410656) q[2];
sx q[2];
rz(0.18266695) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5665633) q[1];
sx q[1];
rz(-1.4815147) q[1];
sx q[1];
rz(2.1818698) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7054889) q[3];
sx q[3];
rz(-1.4105984) q[3];
sx q[3];
rz(-0.45405162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41670123) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(0.42281881) q[2];
rz(-2.4041798) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(3.056934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87930644) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(2.6111531) q[0];
rz(0.92492217) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(1.8431429) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4209375) q[0];
sx q[0];
rz(-1.4268095) q[0];
sx q[0];
rz(-1.7341341) q[0];
x q[1];
rz(1.8680044) q[2];
sx q[2];
rz(-0.98515918) q[2];
sx q[2];
rz(-2.1538018) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4358208) q[1];
sx q[1];
rz(-0.37322361) q[1];
sx q[1];
rz(-0.2758287) q[1];
rz(0.97814822) q[3];
sx q[3];
rz(-1.1453298) q[3];
sx q[3];
rz(0.076171906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9412781) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(0.47362622) q[2];
rz(-0.099362699) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(-2.30106) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6376003) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(2.1437058) q[0];
rz(-0.87431327) q[1];
sx q[1];
rz(-2.120178) q[1];
sx q[1];
rz(-0.46674892) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1186819) q[0];
sx q[0];
rz(-0.70436275) q[0];
sx q[0];
rz(0.33606152) q[0];
rz(2.7141063) q[2];
sx q[2];
rz(-1.3689405) q[2];
sx q[2];
rz(1.0496548) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0007243) q[1];
sx q[1];
rz(-1.5233526) q[1];
sx q[1];
rz(-0.15111698) q[1];
x q[2];
rz(-0.096188992) q[3];
sx q[3];
rz(-1.4078762) q[3];
sx q[3];
rz(-0.067226203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2816887) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(-1.5768645) q[2];
rz(0.62670296) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(0.68157649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8108114) q[0];
sx q[0];
rz(-0.89650506) q[0];
sx q[0];
rz(-2.7217641) q[0];
rz(0.22142521) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(1.5931169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67328582) q[0];
sx q[0];
rz(-1.562403) q[0];
sx q[0];
rz(3.054137) q[0];
x q[1];
rz(2.9485547) q[2];
sx q[2];
rz(-1.2912573) q[2];
sx q[2];
rz(2.2045731) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35510264) q[1];
sx q[1];
rz(-0.10663248) q[1];
sx q[1];
rz(2.998888) q[1];
rz(1.1057304) q[3];
sx q[3];
rz(-1.2601868) q[3];
sx q[3];
rz(0.066699337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.55591136) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(-1.7377724) q[2];
rz(2.3445271) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(-1.2169303) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7109011) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(0.28924334) q[0];
rz(0.62492433) q[1];
sx q[1];
rz(-1.1661252) q[1];
sx q[1];
rz(-1.3141059) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02567357) q[0];
sx q[0];
rz(-0.81633767) q[0];
sx q[0];
rz(1.4195819) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35347519) q[2];
sx q[2];
rz(-2.3620053) q[2];
sx q[2];
rz(-1.0970864) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96771679) q[1];
sx q[1];
rz(-1.6248676) q[1];
sx q[1];
rz(-1.6915583) q[1];
rz(-0.81359158) q[3];
sx q[3];
rz(-1.3435257) q[3];
sx q[3];
rz(-3.1384625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73359314) q[2];
sx q[2];
rz(-1.5780129) q[2];
sx q[2];
rz(1.4286263) q[2];
rz(-2.1777878) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(-2.111964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6190417) q[0];
sx q[0];
rz(-1.6864809) q[0];
sx q[0];
rz(1.8956986) q[0];
rz(0.11101162) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(2.5949809) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1221736) q[0];
sx q[0];
rz(-1.1328567) q[0];
sx q[0];
rz(2.7287448) q[0];
rz(-pi) q[1];
rz(-1.1999646) q[2];
sx q[2];
rz(-1.9413345) q[2];
sx q[2];
rz(1.3818936) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4709028) q[1];
sx q[1];
rz(-2.7109475) q[1];
sx q[1];
rz(-1.9914658) q[1];
x q[2];
rz(-0.9162174) q[3];
sx q[3];
rz(-1.9700053) q[3];
sx q[3];
rz(-2.6938714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1116011) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(1.4477504) q[2];
rz(-1.1374121) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85957134) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(0.71722537) q[0];
rz(-1.2099129) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(0.70770121) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4004138) q[0];
sx q[0];
rz(-2.3614863) q[0];
sx q[0];
rz(-2.0692503) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7949445) q[2];
sx q[2];
rz(-1.4321616) q[2];
sx q[2];
rz(-2.5773406) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.087079436) q[1];
sx q[1];
rz(-2.2594249) q[1];
sx q[1];
rz(-2.7282532) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0503065) q[3];
sx q[3];
rz(-2.4647053) q[3];
sx q[3];
rz(-0.53938473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5132961) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(0.90325242) q[2];
rz(-1.603027) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(0.85047754) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.306504) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(2.3256336) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(-1.1076526) q[2];
sx q[2];
rz(-1.1504428) q[2];
sx q[2];
rz(3.0124315) q[2];
rz(2.49414) q[3];
sx q[3];
rz(-0.068131937) q[3];
sx q[3];
rz(1.6905231) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
